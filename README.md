# Boosted Trees for Risk Prognosis

This is an R implementation of the paper ["Boosted Trees for Risk Prognosis"](http://proceedings.mlr.press/v85/bellot18a.html). 

In this project we present a new approach to ensemble learning for risk prognosis in heterogeneous medical populations. Our aim is to improve overall prognosis by focusing on under-represented patient subgroups with an atypical disease presentation; with current prognostic tools, these subgroups are being consistently mis-estimated. Our method proceeds sequentially by learning nonparametric survival estimators which iteratively learn to improve predictions of previously misdiagnosed patients - a process called boosting. This results in fully nonparametric survival estimates, that is, constrained neither by assumptions regarding the baseline hazard nor assumptions regarding the underlying covariate interactions - and thus differentiating our approach from existing boosting methods for survival analysis. 

*Please cite the above paper if this resource is used in any publication*

## Requirements

* R version 3.5
* Packages: "pec","rpart","party","survival"

## First steps
To get started, check demo.R which will guide you through an application of our algorithm. 

If you have questions or comments about anything regarding this work, please do not hesitate to contact [Alexis](https://alexisbellot.github.io/Website/)
