### code for data prep - Structure <- bgc format ---------

    #Following the example code by:
    #https://popgencode.wordpress.com/2015/04/17/preparing-data-file-bgc-genotype-count/

##load libraries -------------
library(adegenet)

#set working directory (makes for cleaner code) ------------------
setwd("c:/Users/RyanStanley/OneDrive/PostDoc/DFO/Nick/BGC/")

#Load in geneteic data
CrabNorth <- import2genind('CrabOutliers_PureNorth.str', onerowperind=F, n.ind=132, n.loc=117, row.marknames=1, col.lab=1, col.pop=2, ask=F)
CrabSouth <- import2genind('CrabOutliers_PureSouth.str', onerowperind=F, n.ind=66, n.loc=117, row.marknames=1, col.lab=1, col.pop=2, ask=F)
CrabMixed <- import2genind('Outliers_AdmixedIn.str', onerowperind=F, n.ind=43, n.loc=117, row.marknames=1, col.lab=1, col.pop=2, ask=F)

#convert to adegenet format
CrabNorth <- genind2genpop(CrabNorth)
CrabSouth <- genind2genpop(CrabSouth)
CrabMixed <- genind2genpop(CrabMixed) ### I think you have a monomophic allele or something. Not an even number of alleles. 

#Grab the allele frequency format data * the code on the web is missing something
CN_Data=t(CrabNorth@tab)
CrabNorth_allele <- data.frame(Allele=rownames(CN_Data))
CrabNorth_allele <- cbind(CrabNorth_allele,CN_Data)
rownames(CrabNorth_allele) <- NULL

CS_Data=t(CrabSouth@tab)
CrabSouth_allele <- data.frame(Allele=rownames(CS_Data))
CrabSouth_allele <- cbind(CrabSouth_allele,CS_Data)
rownames(CrabSouth_allele) <- NULL

CM_Data=t(CrabMixed@tab)
CrabMixed_allele <- data.frame(Allele=rownames(CM_Data))
CrabMixed_allele <- cbind(CrabMixed_allele,CM_Data)
rownames(CrabMixed_allele) <- NULL ### ** issue here. There seems to be a missing allele for one of the loci

#remove unneeded data from this step
rm(CN_Data,CS_Data,CM_Data)

#Get locus names without the allele appended value
LocusNorth <- as.character(CrabNorth_allele$Allele)
LocusNorth <- lapply(strsplit(LocusNorth, '.', fixed=TRUE), '[[', 1)
LocusNorth <- LocusNorth[rep(c(TRUE,FALSE),nrow(CrabNorth_allele)/2)]
LocusNorth <- unlist(LocusNorth)

LocusSouth <- as.character(CrabSouth_allele$Allele)
LocusSouth <- lapply(strsplit(LocusSouth, '.', fixed=TRUE), '[[', 1)
LocusSouth <- LocusSouth[rep(c(TRUE,FALSE),nrow(CrabSouth_allele)/2)]
LocusSouth <- unlist(LocusSouth)

LocusMixed <- as.character(CrabMixed_allele$Allele)
LocusMixed <- lapply(strsplit(LocusMixed, '.', fixed=TRUE), '[[', 1)
LocusMixed <- LocusMixed[rep(c(TRUE,FALSE),nrow(CrabMixed_allele)/2)]
LocusMixed <- unlist(LocusMixed)

#Get the data for each allele in each population (every other row as denoted by even and odd row integers)
North_Odd <- CrabNorth_allele[seq(1,length(CrabNorth_allele$Allele),2),]
North_Even <- CrabNorth_allele[seq(2,length(CrabNorth_allele$Allele),2),]

South_Odd <- CrabSouth_allele[seq(1,length(CrabSouth_allele$Allele),2),]
South_Even <- CrabSouth_allele[seq(2,length(CrabSouth_allele$Allele),2),]

Mixed_Odd <- CrabMixed_allele[seq(1,length(CrabMixed_allele$Allele),2),]
Mixed_Even <- CrabMixed_allele[seq(2,length(CrabMixed_allele$Allele),2),]

#Create dataframes ## here is where I am confused. Is there supposed to be a new text file for each population?

#For population labeled #1 in the crab data
NorthData <- data.frame(Loci=LocusNorth,filler=rep('#',length(LocusNorth)),PopOdd=North_Odd$`1`,PopEven=North_Even$`1`)

#for each pop
head(CrabNorth_allele)
write.table(NorthData,"test.txt", row.names=F, col.names=F, quote=F, sep='\t')

#North Populations
for(i in colnames(CrabNorth_allele)[2:length(CrabNorth_allele)])
  {
  write.table(data.frame(Loci=LocusNorth,filler=rep('#',length(LocusNorth)),PopOdd=North_Odd[,i],PopEven=North_Even[,i]),
              paste("BGC Files/North_pop",i,"_allct_bgc.txt",sep=""), row.names=F, col.names=F, quote=F, sep='\t')
}

#South Populations
for(i in colnames(CrabSouth_allele)[2:length(CrabSouth_allele)])
{
  write.table(data.frame(Loci=LocusSouth,filler=rep('#',length(LocusSouth)),PopOdd=South_Odd[,i],PopEven=South_Even[,i]),
              paste("BGC Files/South_pop",i,"_allct_bgc.txt",sep=""), row.names=F, col.names=F, quote=F, sep='\t')
}

#Mixed Populations * failes due to monomorphic allele issue
# for(i in colnames(CrabMixed_allele)[2:length(CrabMixed_allele)])
# {
#   write.table(data.frame(Loci=LocusMixed,filler=rep('#',length(LocusMixed)),PopOdd=Mixed_Odd[,i],PopEven=Mixed_Even[,i]),
#               paste("Mixed_pop",i,"_allct_bgc.txt",sep=""), row.names=F, col.names=F, quote=F, sep='\t')
# }


