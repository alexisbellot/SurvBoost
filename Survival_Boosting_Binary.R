Survival_Boosting_Binary = function(data, tree_depth = 10, n_rounds = 50, verbose = TRUE,phi=0.2,power=1,
                        control = NULL, horizon=NULL,subsample=FALSE,frac=0.8){
  
  """
  Fit prediction model to time-to-event data with binary error function
  Inputs:
  - data: R data frame with time to event information coded -Survival- and event indicator coded -Status- 
  - tree_depth: the maximum depth of each tree in the ensemble
  - n_rounds: number of boosting rounds
  - control: rpart.control object to further specify the construction of each tree
  - horizon: the time horizon to be optimized, if NULL all time horizons considered
  - subsample: whether to subsample the data before fitting trees
  - frac: subsampling fraction

  Outputs:
  - fitted model with all estimated parameters
  """
  
  ## packages needed
  library(pec)
  library(party)
  library(survival)
  library(rpart)
  
  # check data types
  if(!all(data$Survival > 0))
    stop("Survival time must be > 0")
  
  # check for presence of rpart control
  if(is.null(control)){
    # very important  to adjust these settings carefully to avoid errors in tree construction
    # ie avoid 1 instance per leaf
    control=rpart.control(minsplit = 20, minbucket = round(20/3), cp = 0.01, 
                          maxcompete = 0, maxsurrogate = 1, usesurrogate = 2, xval = 10,
                          surrogatestyle = 0, maxdepth = tree_depth)
  } 
  
  # initialize output placeholders
  trees                = ws = es = list()
  betas                = phis = adj_es = n_root_nodes =numeric()
  preds                = matrix(nrow=nrow(data),ncol=n_rounds)
  importance           = matrix(0,nrow=n_rounds,ncol=ncol(data)-2)
  colnames(importance) = colnames(subset(data,select=-c(Survival,Status)))
  
  # ipcw for prediction error estimation
  times   = seq(from=min(data$Survival), to=max(data$Survival), length.out = 10)
  formula = as.formula(paste("Surv(Survival,Status)", paste(colnames(subset(data,select=-c(Survival,Status))), collapse=" + "), sep=" ~ "))
  if(!is.null(horizon)){
    brier_weights = ipcw(formula,data=data,method="cox",times=horizon,what="IPCW.times")$IPCW.times
  }
  else {brier_weights = ipcw(formula,data=data,method="cox",times=times,what="IPCW.times")$IPCW.times
  status  = matrix(as.numeric(sapply(as.list(data$Survival),function(x)x > times)),nrow=length(data$Survival))}#n*10 mat
  
  # mean error as a function of phi
  ev_phi = function(a){mean(as.numeric(brier_weights*(as.numeric(data$Survival>horizon)-pred)^2  > a))}
  
  # initialize weights
  w = ifelse(rep(subsample,nrow(data)),rep(1/nrow(data), nrow(data)),rep(1, nrow(data)))
  
  
  ## BOOSTING ITERATION --------------------------------------------------------------------------------------
  ##
  for(i in seq(n_rounds)){
    
    # fit survival tree
    if(subsample){
      sub_sample    = sample(1:nrow(data),round(frac*nrow(data)),replace=F, prob = w) 
      tree          = pecRpart(Surv(Survival,Status)~.,data=data[sub_sample,], control = control)# fitted model
    }
    else{tree = pecRpart(Surv(Survival,Status)~.,data=data, weights = w, control = control)} # fitted model
    
    
    # prediction and performance evauation
    if(!is.null(horizon)){
      pred = predictSurvProb(tree, newdata=data,times=horizon) # predictions insample
      if(i==1 && is.null(phi)){
        s   = seq(from=0.17,to=0.22,by=0.01)
        phi = s[which.min(abs(c(ev_phi(s[1]),ev_phi(s[2]),ev_phi(s[3]),ev_phi(s[4]),ev_phi(s[5]),ev_phi(s[6]))-0.4))]}
      e   = as.numeric(brier_weights*(as.numeric(data$Survival>horizon)-pred)^2  > phi)# error is w_i*(delta_i(t)-S_i(t))^2
    }
    else{
      pred    = predictSurvProb(tree, newdata=data,times=times) # predictions insample
      e       = as.numeric(rowSums(brier_weights*(status-pred)^2)/length(data$Survival) > phi)} # error integrated over time
    
    
    # If tree perfectly gets data, boosting terminates
    if(sum(e,na.rm=T) < 2){
      # handle the case where first base classifier fits data perfectly
      if(i == 1){
        trees[[i]]  = tree
        betas[[i]]  = 1
        terms       = tree$terms
        break
      }
      break
    }
    
    e[is.na(e)] = 1
    # adjusted error based on weights
    adj_e = sum(e*w) 
    beta  = adj_e
    # updated weights
    w     = w*(e + (1-e)*beta) 
    w     = pmin(5/nrow(data),pmax(w/sum(w),0.2/nrow(data)))
    
    # store each tree and predictions
    trees[[i]]      = tree; n_root_nodes[i] = length(tree$levels); betas[i] = beta; phis[i] = phi
    importance[i,]  = replace(importance[i,],names(tree$rpart$variable.importance),as.numeric(tree$rpart$variable.importance))
    adj_es[i]       = adj_e; ws[[i]] = w;es[[i]] = e; preds[,i] = pred
    
    if(verbose & (i %% 10 == 0))
      cat("Iteration: ", i, "\n")
  }
  
  out        = list(betas = unlist(betas), tree_depth = tree_depth, phis=phis,variable_importance=importance,
                    adj_es=adj_es,ws=ws,es=es,preds=preds, trees = trees,n_root_nodes = n_root_nodes)
  class(out) = "adaboost"
  out
}