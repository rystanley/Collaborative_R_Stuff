#Libraries
library(ggplot2)

## load data
outliers <- read.table("c:/Users/test/OneDrive/PostDoc/DFO/Nick/GenePop/CrabFinal_outliers.txt",sep="\t",header=T)
outliers <- as.vector(as.numeric(outliers[,1]))
ciData <- read.csv("c:/Users/test/Desktop/Nick Arlequin/CI_Intervals.csv",sep=",",header=T)
colnames(ciData)[c(1,length(ciData))] <-c("Obs_Het_BP", "FST")
ciData$colour="black"
ciData <- subset(ciData,Obs_Het_BP<0.55 & Obs_Het_BP>=0.1)


aData <- read.table("C:/Users/test/Desktop/Nick Arlequin/fdist2_ObsOut.txt",sep="\t",header=T)
colnames(aData) <- c("Locus","Obs_Het_BP","FST","pval","quantile")


aData$colour <- "black"
aData[outliers,"colour"] <- "red"

    
    ggplot(aData,aes(x=Obs_Het_BP,y=FST,colour=colour,size=colour))+geom_point()+
      scale_colour_manual(values=c("Black","red"))+scale_size_manual(values=c(1,3))+
      theme_bw()+geom_line(data = ciData,aes(guide=FALSE),lty=2)+
      theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      labs(x="Observed Heterozygosity",y=expression("F"["st"]))
    

  