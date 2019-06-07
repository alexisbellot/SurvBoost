Survival_Boosting_Continuous = function(data, tree_depth = 10, n_rounds = 50, verbose = TRUE,
                             control = NULL, horizon=median(data$Survival),subsample=TRUE,frac=0.8){
  """
  Fit prediction model to time-to-event data with continuous error function
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
    
  } else if(control$maxdepth != tree_depth){
    warning(paste('tree_depth set to: ', control$maxdepth))
  }
  
  # initialize output placeholders
  trees         = ws = es = list()
  betas         = adj_es = n_terminal_nodes =numeric()
  preds         = matrix(nrow=nrow(data),ncol=n_rounds)
  importance    = matrix(0,nrow=n_rounds,ncol=ncol(data)-2)
  colnames(importance) = colnames(subset(data,select=-c(Survival,Status)))
  
  linMap = function(x, from, to){(x - min(x)) / max(x - min(x)) * (to - from) + from}
  
  w = ifelse(rep(subsample,nrow(data)),rep(1/nrow(data), nrow(data)),rep(1, nrow(data)))
  #w             = rep(1, nrow(data)) # weights
  
  times   = seq(from=min(data$Survival), to=max(data$Survival), length.out = 10)
  formula = as.formula(paste("Surv(Survival,Status)", paste(colnames(subset(data,select=-c(Survival,Status))), collapse=" + "), sep=" ~ "))
  
  # ipcw
  if(!is.null(horizon)){
    brier_weights = ipcw(formula,data=data,method="cox",times=horizon,what="IPCW.times")$IPCW.times
  }
  else {brier_weights = ipcw(formula,data=data,method="cox",times=times,what="IPCW.times")$IPCW.times
  status  = matrix(as.numeric(sapply(as.list(data$Survival),function(x)x > times)),nrow=length(data$Survival))}#n*10 mat
  
  # Boosting iterations
  for(i in seq(n_rounds)){
    
    # subsample manually fraction of the data of data with prob. prop to e
    if(subsample){
      sub_sample    = sample(1:nrow(data),round(frac*nrow(data)),replace=F, prob = w) 
      tree          = pecRpart(Surv(Survival,Status)~.,data=data[sub_sample,], control = control)# fitted model
    }
    else{tree = pecRpart(Surv(Survival,Status)~.,data=data, weights = w, control = control)} # fitted model
    
    # evaluation at a defined time horizon (time dependent brier)
    if(!is.null(horizon)){
      pred = predictSurvProb(tree, newdata=data,times=horizon) # predictions insample
      e    = brier_weights*(as.numeric(data$Survival>horizon)-pred)^2  # error is w_i*(delta_i(t)-S_i(t))^2
      e[is.na(e)] = mean(e,na.rm=T)
    }
    # evaluation integrated over time (time dependent brier)
    else{pred    = predictSurvProb(tree, newdata=data,times=times) # predictions insample
    e       = rowSums(brier_weights*(status-pred)^2)/length(data$Survival)}
    
    # If tree perfectly gets data, boosting terminates
    if(abs(mean(e)) < 1e-8){
      # handle the case where first base classifier fits data perfectly
      if(i == 1){
        trees[[i]]  = tree
        alphas[[i]] = 1
        terms       = tree$terms
        break
      }
      break
    }
    
    adj_e = sum(e*w/sum(w)) # adjusted error based on weights
    if(adj_e>0.5){ break }
    beta  = adj_e/(1-adj_e) # limit to 0.5
    w     = w*as.numeric(beta^(1-e)) # updated weights
    w     = as.numeric(round(linMap(w,from=1,to=10),digits=0)) # limit weights for regularisation
    w     = w/sum(w)
    #w     = pmin(5/nrow(data),pmax(w/sum(w),0.2/nrow(data)))

    
    # store for predictions
    trees[[i]]      = tree; n_terminal_nodes[i] = length(tree$levels); betas[i] = beta
    importance[i,]  = replace(importance[i,],names(tree$rpart$variable.importance),as.numeric(tree$rpart$variable.importance))
    adj_es[i]       = adj_e; ws[[i]] = w;es[[i]] = e; preds[,i] = pred
    
    if(verbose & (i %% 10 == 0))
      cat("Iteration: ", i, "\n")
  }
  out        = list(betas = unlist(betas), trees = trees, tree_depth = tree_depth, variable_importance=importance,
                    n_terminal_nodes=n_terminal_nodes,adj_es=adj_es,ws=ws,preds=preds)
  class(out) = "adaboost"
  out
}
