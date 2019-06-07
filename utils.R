
## utilities for prediction and performance computation

library(pec)


## Prediction
## -----------------------------------------------------------------------------------------------------------------------------

predict.adaboost = function(object, X, times, n_tree = NULL, ...){
  # handle args
  if(is.null(n_tree)){
    tree_seq = seq_along(object$betas)
  } else{
    if(n_tree > length(object$betas))
      stop('n_tree must be less than the number of trees used in fit')
    tree_seq = seq(1, n_tree)
  }
  
  # evaluate score function on sample
  f = 0
  for(i in tree_seq){
    tree = object$trees[[i]]
    pred = predictSurvProb(tree, newdata=data.frame(X),times=times)
    ## create an index of missing values 
    index <- which(is.na(pred)|is.nan(pred), arr.ind = TRUE) 
    ## calculate the row means and "duplicate" them to assign to appropriate cells 
    pred[index] <- colMeans(pred, na.rm = TRUE)[index[, "col"]] 
    f    = f + log(1/object$betas[i])*pred
  }
  # return weighted average
  return(f/sum(log(1/object$betas[tree_seq])))
}



## Performance computation
## -----------------------------------------------------------------------------------------------------------------------------


performance = function(data,k = 2,n_rounds=250,subsample=TRUE,frac=0.8) {
  
  n     = nrow(data)
  folds = split(sample(n), seq_len(k))
  
  xval.fold = function(fold) {
    
    unique.times  = mean(data$Survival)
    ## Fit model
    ada           = Survival_Boosting(data=data[-fold,],tree_depth=3, subsample=subsample,frac=frac,
                                      horizon=mean(data$Survival), n_rounds=n_rounds)
    ## Survival predictions
    predictions   = predict(ada, data[fold,],times=unique.times)
    ## C indeces
    cindex.ada    = cindex(object=as.matrix(predictions),formula=Surv(Survival,Status)~1, data=data[fold,],eval.times=unique.times)
    brier.ada     = pec(object=as.matrix(predictions),formula=Surv(Survival,Status)~1,data=data[fold,],start=unique.times[1],maxtime=unique.times[length(unique.times)],
                        exact=FALSE, times=unique.times,reference=FALSE,splitMethod = "none")
    
    return(data.frame('cindex.ada'=cindex.ada$AppCindex,'brier.ada'=brier.ada$AppErr,'unique.times'=unique.times))
  }#XVAL.FOLD
  
  results = lapply(folds, xval.fold)
  return(do.call("rbind", results))
  
}#XVAL

perf.est = function(data,k=2,n_rounds=250,subsample=TRUE,frac=0.8){
  xval.out = performance(data,k,n_rounds,subsample,frac)
  return(data.frame(ada.cindex=tapply(xval.out[,1],as.factor(xval.out[,3]),mean),
                    sd.ada.cindex=tapply(xval.out[,1],as.factor(xval.out[,3]),sd),
                    ada.brier=tapply(xval.out[,2],as.factor(xval.out[,3]),mean),
                    sd.ada.brier=tapply(xval.out[,2],as.factor(xval.out[,3]),sd)))
}#PERFORMANCE ESTIMATE
