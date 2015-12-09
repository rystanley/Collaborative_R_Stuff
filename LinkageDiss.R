#Code to filter New Hybrids file for linkage disequilibrium

#load libraries
library(dplyr)

#load data
FilterData <- read.csv("Top500_LD.csv",sep=",",row.names=1)
FilterData=as.matrix(FilterData) # convert to matrix
diag(FilterData) <- 0 # replace the values of 1 with 0 in the diagonal
FilterData=as.data.frame(FilterData) # convert back to dataframe

Temp1=apply(FilterData,2,max,na.rm=1) # get the maximum r2 for each column (SNP)
Names <- gsub("X","",colnames(t(Temp1))) # get the SNP names (because they are read in as numeric R adds an X to them)

Temp2=data.frame(SNPS=Names,Linkage=as.vector(Temp1)) #create an output dataframe

Out=filter(Temp2,Linkage<0.25) #filter the output dataframe for only those who have max r2 less than 0.25

write.table(Out,file="c:/Users/test/OneDrive/PostDoc/DFO/Nick/New Hybrids/LinkageDis.csv",sep=",",row.names=F) #save file. 
