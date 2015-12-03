## load libraries --------  
library(VennDiagram)

## load data -------
OutData <- read.csv("Outliers.csv",header=T,sep=",")

#Parse data --------
BayeScan <- as.vector(OutData[which(OutData$Method=="BayeScan"),"Loci"])
Arelquin <- as.vector(OutData[which(OutData$Method=="Arl"),"Loci"])
ArelquinIsland <- as.vector(OutData[which(OutData$Method=="ArlHeir"),"Loci"])

# number of overlapping predictions among all methods
length(Reduce(intersect,list(BayeScan,Arelquin,ArelquinIsland)))

##Create venn diagram ---------
venn.diagram(list(BayeScan,ArelquinIsland,Arelquin),
             "Outlier_Venn.tiff",
             fill=c("white","grey90","lightblue"),
             category=c("BS","AH","ANH"),
             rotation.degree=0.45)

