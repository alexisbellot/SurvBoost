"""
Example performance estimation on the publicly available GBSG2 cancer data 
"""

## Required packages

library('survival')
library('pec')
library('party')
library('rpart')


## Load model and utility function for performance computation

source('/utils.R')
source('/Survival_Boosting_Continuous.R')

## Load and preprocess data

data(GBSG2)
data = GBSG2
colnames(data)[10]<- 'Status'
data$horTh= as.numeric(data$horTh) -1
data$menostat = as.numeric(data$menostat) -1
data$tgrade1 = as.numeric(data$tgrade=="I")
data$tgrade2 = as.numeric(data$tgrade=="II")
#data$tgrade3 = as.numeric(data$tgrade=="III")
data$time = data$time/7
data = subset(data, select=-c(tgrade))
colnames(data)[1] = "Treatment"
colnames(data)[8] = "Survival"


## Run algorithm and compute performance

perf = perf.est(data,k=3,n_rounds=250,subsample=0.8)

