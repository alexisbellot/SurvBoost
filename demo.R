"""
Example performance estimation on the publicly available GBSG2 cancer data 
"""

## Required packages

library('survival')
library('pec')
library('party')
library('rpart')
library('ggplot2')


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

## Run algorithm and plot predictions
ada = Survival_Boosting_Continuous(data=data,horizon=median(data$Survival))
predictions = predict.adaboost(ada,data[1:10,],times = seq(1,100,length.out = 100))

predictions_frame = data.frame(times = seq(1,100,length.out = 100), patient_1=predictions[1,1:100],
                               patient_2=predictions[2,1:100])

ggplot(predictions_frame, aes(times))+ xlim(-1,100) + ylim(0.75,1)+
  geom_line(aes(y=patient_1,col="1"),size=1)+
  geom_line(aes(y=patient_2, col="2"),size=1)+
  scale_colour_manual(name="",values=c("1"="orange", "2"="blue"))+
  theme(legend.position="",legend.text=element_text(size=10),
        panel.grid.major = element_line(colour = "gray90"),
        panel.grid.minor = element_line(colour = "gray90"),
        panel.background = element_rect(fill='white'))+labs(x = "Time (Weeks)", y="Survival probability")


## Run algorithm and compute performance

perf = perf.est(data,k=3,n_rounds=250,subsample=0.8)

