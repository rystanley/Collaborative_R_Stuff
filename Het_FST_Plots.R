### This plot will make a combined frequency distribution plot of FST and HZ of Green Crab Data -------

#load libraries
library(ggplot2)

#load data
dat=read.csv("Locus_Fst.csv",sep=",")

dat_Het=dat[,c("Locus","Het")];colnames(dat_Het)[2]="Value";dat_Het$metric="Hetero"
dat_FST=dat[,c("Locus","FST")];colnames(dat_FST)[2]="Value";dat_FST$metric="FST"

#set the break levels for the histogram
b1=seq(0,0.7,by=0.02)

#Plot the density of observations
pDens=ggplot(dat_FST,aes(x=Value))+
  geom_histogram(aes(y=..density..),breaks=b1,fill="white",colour="black")+
  theme_bw()+geom_density(data=dat_Het,lwd=0.75)+
  labs(x="Observed value",y="Density of observations")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank());pDens

##Plot the frequency of observations
#Note that this needs some hard coding to work. See in the breaks function (b1) for the split 
#levels (by=) and the number of observations. This split level and number of observations must
#be entered in the geom_density() layer to match the data. 

nrow(dat) #number of observations

#Plot the frequency of observations
pFreq=ggplot(dat_FST,aes(x=Value))+
  geom_histogram(aes(y=..count..),fill="white",breaks=b1,colour="black")+
  geom_density(data=dat_Het,aes(y=..density..*9137*0.02),lwd=0.75)+
  labs(x="Observed value",y="Frequency of observations")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank());pFreq

## save the plots 
ggsave("Het_FST_Density.png",pDens)
ggsave("Het_FST_Frequency.png",pFreq)
