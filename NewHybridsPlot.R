## This file will create a structure like plot for "New Hybrids" ouput

##load libraries -----------
library(ggplot2)
library(dplyr)
library(reshape2)
library(grid)

## load data -----------
nhdata <- read.csv("Outlier_NewHybsresults.csv")

#Change column names
names(nhdata)[2] <- "Sweeps"

#Melt the data for plotting
nDat <- melt(nhdata[,],id=c("Pop","Sweeps"))

#Here you need to fill in the rest (i.e. levels or order you want your plot to be rastered). If not specified ggplot will
#use alphabetical.

#nDat$Pop=factor(nDat$Pop,levels=c("TKT","BDB","CBI", .. you fill out the rest)
  
#Plot the data
p1=ggplot(nDat,aes(x=factor(Sweeps),y=value,fill=variable,group=Pop)) + 
    geom_bar(position = "stack",stat = "identity",width=1)+
    theme_bw()+scale_y_continuous(expand = c(0,0),labels = percent)+
    facet_grid(~Pop,scales="free")+
    theme(axis.text.x = element_blank(),
          axis.ticks=element_blank(),
          panel.margin = unit(0.1, "lines"),
          legend.position="bottom")+
    labs(x="Population",y="",fill="");p1

#Save plot  
ggsave("NewHybrids_Output.png",
       width=8.5,height=7.5,p1)
  

