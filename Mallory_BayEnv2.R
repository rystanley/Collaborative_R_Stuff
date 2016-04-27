#set wd
setwd("~/School/MASTER_Masters/Genetic_Oceano/BayEnv_associations_results/Data/FINAL DATA/")
setwd("~/Desktop/Mallory_BayEnv_Aug2015")

####Env Data ####
#load env data (not transposed)
STAND_master_oceano_data_NoCorrTest_FINAL_FixedNames <- read.csv("~/School/MASTER_Masters/Oceanographic_Data/FINAL DATA/NoCorrelationTest/STAND_master_oceano_data_NoCorrTest_FINAL_NoNA_FixedNames.csv")
#rownames(STAND_master_oceano_data_NoCorrTest_FINAL_FixedNames) <-STAND_master_oceano_data_NoCorrTest_FINAL_FixedNames$X
#STAND_master_oceano_data_NoCorrTest_FINAL_FixedNames <- STAND_master_oceano_data_NoCorrTest_FINAL_FixedNames[,-1]

##transpose oceano data
BayEnv_oceano_data_stand_allvars <- as.data.frame(t(STAND_master_oceano_data_NoCorrTest_FINAL_FixedNames))
write.csv(BayEnv_oceano_data_stand_allvars, "BayEnv_oceano_data_stand_allvars.csv")

####edit oceano data in excel to create files of different env variables



#####Gen Data ####


raw_allele_data <- read.csv("finaldata_HWE_BayEnv_inprog_TESTING.csv", header=T,na.strings=c(" ","","NA"))
rownames(raw_allele_data) <- raw_allele_data$X
raw_allele_data <- raw_allele_data[,-1]

trans_raw_allele_data <- as.data.frame(t(raw_allele_data))

#import list of loci
locus_list <- read.csv("~/School/MASTER_Masters/Genetic_Data/June_2014_finals/All_loci.csv", header=F)
locus_list <- as.vector(locus_list)

#make a list of each locus listed twice
locus_list2 <- as.vector(rep(locus_list$V1, each=2))

#make list of individuals
ind_list <- read.csv("~/School/MASTER_Masters/Genetic_Data/June_2014_finals/All_individuals.csv", header=F)
ind_list <- as.vector(ind_list)

#make list of individuals each listed twice
ind_list2 <- as.vector(rep(ind_list$V1, each=2))

#add individual column in trans_raw_allele_data
trans_raw_allele_data$Ind <- ind_list2

#reorder
trans_raw_allele_data <- trans_raw_allele_data[,c(7164,1:7163)]

#issues coming up here because some loci have 2 alleles and a blank, some are missing that blank
# will try adding a blank row to the bottom of test3 and see if that fixes it

 # create a one-row matrix the same length as data

 temprow <- matrix(c(rep.int(NA,length(data))),nrow=1,ncol=length(trans_raw_allele_data))

# # make it a data.frame and give cols the same names as data

 newrow <- data.frame(temprow)
 colnames(newrow) <- colnames(trans_raw_allele_data)
 # rbind the empty row to data

trans_raw_allele_data <- rbind(trans_raw_allele_data,newrow)

###data frame is ready - trans_raw_allele_data
#transpose the dataframe and export to excel for use in final counts
re_trans_raw_allele_data <- as.data.frame(t(trans_raw_allele_data))

write.csv(re_trans_raw_allele_data,"re_trans_raw_allele_data.csv")


#make a nwe data frame with 4 columns - Locus, Allele1, Allele2, Extra

Locus_allele <- as.data.frame(locus_list)

Locus_allele$Allele1 <- NA
Locus_allele$Allele2 <- NA
Locus_allele$Extra <- NA

rownames(Locus_allele) <- locus_list
Locus_allele <- Locus_allele[,-1]

##add unique alleles per locus to each row in the data frame
for (i in 2:7164) {
    Locus_allele[(i-1),] <- unique(trans_raw_allele_data[,i])

}

###it worked!!  now to figure out how to get ride of the NAs and only have 2 cols of alleles

Locus_allele$TESTING <- paste(Locus_allele$Allele1, Locus_allele$Allele2, Locus_allele$Extra, sep=":")

#export this to Excel to quickly get rid of NAs using text-to-column
write.csv(Locus_allele, "Locus_allele.csv")
#fixed up in excel

#reimport
Locus_allele_final <- read.csv("~/School/MASTER_Masters/Genetic_Data/June_2014_finals/All_Locus_allele.csv")
colnames(Locus_allele_final) <- c("Locus", "Allele1", "Allele2")


###I have a list of Locus and each unique allele

#make a data frame with each locus listed twice (rows), and 12 populations (columns)
allele_data_BayEnv_prog <- as.data.frame(locus_list2)
allele_data_BayEnv_prog[,2:13] <- NA
colnames(allele_data_BayEnv_prog) <- c("Locus","SUN","LTB","MGD","NTS","PSB","BOF","SSM","GMI","SSB","GMO","GEO","MDA")

#export again for more excel fun!
write.csv(allele_data_BayEnv_prog, "allele_data_BayEnv_prog.csv")

###finished in Excel using countif()

####reimport
BayEnv_allele_data_final_HWE <- read.csv("~/School/MASTER_Masters/Genetic_Oceano/BayEnv2Outliers/FINAL DATA/BayEnv_allele_data_final_HWE.csv")


####FINAL Gen and Env DATA####
#files are: BayEnv_allele_data_final_HWE.csv and BayEnv_oceano_data_stand_allvars.csv

###final R data Frames BayEnv_allele_data_final_HWE, BayEnv_oceano_data_stand_allvars


#####create subset for matrix estimation####
#make list of all loci with numbers
locus_list_with_numbers <- locus_list
locus_list_with_numbers$Number <- 1:7163
write.csv(locus_list_with_numbers, "locus_list_with_numbers.csv")

###MATRIX A###
#subset of 100 loci- just do numbers, the use those numbers to match with locus_list
subset_locus_numbers_MATRIX_A <- sample(1:7163, 100, replace=F)
subset_locus_numbers_MATRIX_A <- sort(subset_locus_numbers_MATRIX_A)
subset_locus_df_MATRIX_A <- as.data.frame(subset_locus_numbers_MATRIX_A)
subset_locus_df_MATRIX_A$Locus <- NA

#this is really complicated in R, export to use index match in excel - match loci to subset numbers
write.csv(subset_locus_df_MATRIX_A, "subset_locus_df_MATRIX_A.csv")

#reimport subset_locus_list
subset_locus_df_MATRIX_A <- read.csv("subset_locus_df_MATRIX_A.csv", header=T)

###make a vector of subset locus names
subset_locus_names_MATRIX_A <- as.vector(subset_locus_df_MATRIX_A$Locus)

####subset final data using that list
BayEnv_subset_allele_data_final_HWE_MATRIX_A <- subset(BayEnv_allele_data_final_HWE, BayEnv_allele_data_final_HWE$Locus %in% subset_locus_names_MATRIX_A)

BayEnv_allele_data_minus_subset_final_HWE_MATRIX_A <- subset(BayEnv_allele_data_final_HWE, !BayEnv_allele_data_final_HWE$Locus %in% subset_locus_names_MATRIX_A)

###export!
write.csv(BayEnv_subset_allele_data_final_HWE_MATRIX_A, "BayEnv_subset_allele_data_final_HWE_MATRIX_A.csv")
write.csv(BayEnv_allele_data_minus_subset_final_HWE_MATRIX_A, "BayEnv_allele_data_minus_subset_final_HWE_MATRIX_A.csv")


###MATRIX B###
#subset of 100 loci- just do numbers, the use those numbers to match with locus_list
subset_locus_numbers_MATRIX_B <- sample(1:7163, 100, replace=F)
subset_locus_numbers_MATRIX_B <- sort(subset_locus_numbers_MATRIX_B)
subset_locus_df_MATRIX_B <- as.data.frame(subset_locus_numbers_MATRIX_B)
subset_locus_df_MATRIX_B$Locus <- NA

#this is really complicated in R, export to use index match in excel - match loci to subset numbers
write.csv(subset_locus_df_MATRIX_B, "subset_locus_df_MATRIX_B.csv")

#reimport subset_locus_list
subset_locus_df_MATRIX_B <- read.csv("subset_locus_df_MATRIX_B.csv", header=T)

###make a vector of subset locus names
subset_locus_names_MATRIX_B <- as.vector(subset_locus_df_MATRIX_B$Locus)

####subset final data using that list
BayEnv_subset_allele_data_final_HWE_MATRIX_B <- subset(BayEnv_allele_data_final_HWE, BayEnv_allele_data_final_HWE$Locus %in% subset_locus_names_MATRIX_B)

BayEnv_allele_data_minus_subset_final_HWE_MATRIX_B <- subset(BayEnv_allele_data_final_HWE, !BayEnv_allele_data_final_HWE$Locus %in% subset_locus_names_MATRIX_B)

###export!
write.csv(BayEnv_subset_allele_data_final_HWE_MATRIX_B, "BayEnv_subset_allele_data_final_HWE_MATRIX_B.csv")
write.csv(BayEnv_allele_data_minus_subset_final_HWE_MATRIX_B, "BayEnv_allele_data_minus_subset_final_HWE_MATRIX_B.csv")


###MATRIX C###
#subset of 100 loci- just do numbers, the use those numbers to match with locus_list
subset_locus_numbers_MATRIX_C <- sample(1:7163, 100, replace=F)
subset_locus_numbers_MATRIX_C <- sort(subset_locus_numbers_MATRIX_C)
subset_locus_df_MATRIX_C <- as.data.frame(subset_locus_numbers_MATRIX_C)
subset_locus_df_MATRIX_C$Locus <- NA

#this is really complicated in R, export to use index match in excel - match loci to subset numbers
write.csv(subset_locus_df_MATRIX_C, "subset_locus_df_MATRIX_C.csv")

#reimport subset_locus_list
subset_locus_df_MATRIX_C <- read.csv("subset_locus_df_MATRIX_C.csv", header=T)

###make a vector of subset locus names
subset_locus_names_MATRIX_C <- as.vector(subset_locus_df_MATRIX_C$Locus)

####subset final data using that list
BayEnv_subset_allele_data_final_HWE_MATRIX_C <- subset(BayEnv_allele_data_final_HWE, BayEnv_allele_data_final_HWE$Locus %in% subset_locus_names_MATRIX_C)

BayEnv_allele_data_minus_subset_final_HWE_MATRIX_C <- subset(BayEnv_allele_data_final_HWE, !BayEnv_allele_data_final_HWE$Locus %in% subset_locus_names_MATRIX_C)

###export!
write.csv(BayEnv_subset_allele_data_final_HWE_MATRIX_C, "BayEnv_subset_allele_data_final_HWE_MATRIX_C.csv")
write.csv(BayEnv_allele_data_minus_subset_final_HWE_MATRIX_C, "BayEnv_allele_data_minus_subset_final_HWE_MATRIX_C.csv")

###MATRIX D###
#subset of 100 loci- just do numbers, the use those numbers to match with locus_list
subset_locus_numbers_MATRIX_D <- sample(1:7163, 100, replace=F)
subset_locus_numbers_MATRIX_D <- sort(subset_locus_numbers_MATRIX_D)
subset_locus_df_MATRIX_D <- as.data.frame(subset_locus_numbers_MATRIX_D)
subset_locus_df_MATRIX_D$Locus <- NA

#this is really complicated in R, export to use index match in excel - match loci to subset numbers
write.csv(subset_locus_df_MATRIX_D, "subset_locus_df_MATRIX_D.csv")

#reimport subset_locus_list
subset_locus_df_MATRIX_D <- read.csv("subset_locus_df_MATRIX_D.csv", header=T)

###make a vector of subset locus names
subset_locus_names_MATRIX_D <- as.vector(subset_locus_df_MATRIX_D$Locus)

####subset final data using that list
BayEnv_subset_allele_data_final_HWE_MATRIX_D <- subset(BayEnv_allele_data_final_HWE, BayEnv_allele_data_final_HWE$Locus %in% subset_locus_names_MATRIX_D)

BayEnv_allele_data_minus_subset_final_HWE_MATRIX_D <- subset(BayEnv_allele_data_final_HWE, !BayEnv_allele_data_final_HWE$Locus %in% subset_locus_names_MATRIX_D)

###export!
write.csv(BayEnv_subset_allele_data_final_HWE_MATRIX_D, "BayEnv_subset_allele_data_final_HWE_MATRIX_D.csv")
write.csv(BayEnv_allele_data_minus_subset_final_HWE_MATRIX_D, "BayEnv_allele_data_minus_subset_final_HWE_MATRIX_D.csv")

###MATRIX E###
#subset of 100 loci- just do numbers, the use those numbers to match with locus_list
subset_locus_numbers_MATRIX_E <- sample(1:7163, 100, replace=F)
subset_locus_numbers_MATRIX_E <- sort(subset_locus_numbers_MATRIX_E)
subset_locus_df_MATRIX_E <- as.data.frame(subset_locus_numbers_MATRIX_E)
subset_locus_df_MATRIX_E$Locus <- NA

#this is really complicated in R, export to use index match in excel - match loci to subset numbers
write.csv(subset_locus_df_MATRIX_E, "subset_locus_df_MATRIX_E.csv")

#reimport subset_locus_list
subset_locus_df_MATRIX_E <- read.csv("subset_locus_df_MATRIX_E.csv", header=T)

###make a vector of subset locus names
subset_locus_names_MATRIX_E <- as.vector(subset_locus_df_MATRIX_E$Locus)

####subset final data using that list
BayEnv_subset_allele_data_final_HWE_MATRIX_E <- subset(BayEnv_allele_data_final_HWE, BayEnv_allele_data_final_HWE$Locus %in% subset_locus_names_MATRIX_E)

BayEnv_allele_data_minus_subset_final_HWE_MATRIX_E <- subset(BayEnv_allele_data_final_HWE, !BayEnv_allele_data_final_HWE$Locus %in% subset_locus_names_MATRIX_E)

###export!
write.csv(BayEnv_subset_allele_data_final_HWE_MATRIX_E, "BayEnv_subset_allele_data_final_HWE_MATRIX_E.csv")
write.csv(BayEnv_allele_data_minus_subset_final_HWE_MATRIX_E, "BayEnv_allele_data_minus_subset_final_HWE_MATRIX_E.csv")


####10% neutral matrix####

#subset locus list so that it's only neutral loci - do in excel

#create list of 700 neutral loci to use
#subset of 100 loci- just do numbers, the use those numbers to match with locus_map_numbers_neutral.csv
subset_locus_numbers_NeutralMatrix <- sample(1:7051, 700, replace=F)
subset_locus_numbers_NeutralMatrix <- sort(subset_locus_numbers_NeutralMatrix)
subset_locus_numbers_NeutralMatrix <- as.data.frame(subset_locus_numbers_NeutralMatrix)
subset_locus_numbers_NeutralMatrix$Locus <- NA

#this is really complicated in R, export to use index match in excel - match loci to subset numbers
write.csv(subset_locus_numbers_NeutralMatrix, "subset_locus_numbers_NeutralMatrix.csv")
#reimport subset_locus_list
subset_locus_numbers_NeutralMatrix <- read.csv("subset_locus_numbers_NeutralMatrix.csv", header=T)

###make a vector of subset locus names
subset_locus_names_NeutralMatrix <- as.vector(subset_locus_numbers_NeutralMatrix$Locus)

####subset final data using that list
#Just the 700 neutral loci  
BayEnv_subset_allele_data_final_HWE_NeutralMatrix <- subset(BayEnv_allele_data_final_HWE, BayEnv_allele_data_final_HWE$Locus %in% subset_locus_names_NeutralMatrix)
#all other loci (including outliers)
BayEnv_allele_data_minus_subset_final_HWE_NeutralMatrix <- subset(BayEnv_allele_data_final_HWE, !BayEnv_allele_data_final_HWE$Locus %in% subset_locus_names_NeutralMatrix)

###export!
write.csv(BayEnv_subset_allele_data_final_HWE_NeutralMatrix, "BayEnv_subset_allele_data_final_HWE_NeutralMatrix.csv")
write.csv(BayEnv_allele_data_minus_subset_final_HWE_NeutralMatrix, "BayEnv_allele_data_minus_subset_final_HWE_NeutralMatrix.csv")


#####MATRIX COMPARISONS#####
library("gplots", lib.loc="~/R/win-library/3.2")
pops <- c("SUN","LTB","MGD","NTS","PSB","BOF","SSM","GMI","SSB","GMO","GEO","MDA")

###use BayEnv to create covariance matrices, then use final matrix to create NeutralMatrix_final file
###create a few different covar matrices for the neutral data, then compare to ensure it's estimated well?
###see codelines_Mallory_toDo_Halifax_Nov14-21.txt


#load covariance matrices into R

matrixA_final <- as.matrix(read.csv("matrixA_final.csv", header=F))
matrixB_final <- as.matrix(read.csv("matrixB_final.csv", header=F))
matrixC_final <- as.matrix(read.csv("matrixC_final.csv", header=F))
matrixD_final <- as.matrix(read.csv("matrixD_final.csv", header=F))
matrixE_final <- as.matrix(read.csv("matrixE_final.csv", header=F))
neutralMatrix_final <- as.matrix(read.csv("neutralMatrix_final.csv", header=F))

#convert covariance matrices into correlation matrices using cov2cor (see BayEnv2 manual)
matrixA_final_correlation <- cov2cor(matrixA_final)
rownames(matrixA_final_correlation) <- pops
colnames(matrixA_final_correlation) <- pops

matrixB_final_correlation <- cov2cor(matrixB_final)
rownames(matrixB_final_correlation) <- pops
colnames(matrixB_final_correlation) <- pops

matrixC_final_correlation <- cov2cor(matrixC_final)
rownames(matrixC_final_correlation) <- pops
colnames(matrixC_final_correlation) <- pops

matrixD_final_correlation <- cov2cor(matrixD_final)
rownames(matrixD_final_correlation) <- pops
colnames(matrixD_final_correlation) <- pops

matrixE_final_correlation <- cov2cor(matrixE_final)
rownames(matrixE_final_correlation) <- pops
colnames(matrixE_final_correlation) <- pops

neutralMatrix_final_correlation <- cov2cor(neutralMatrix_final)
rownames(neutralMatrix_final_correlation) <- pops
colnames(neutralMatrix_final_correlation) <- pops

#load pop pairwise Fst
pairwise_Fst <- read.csv("C:/Users/Mallory/Documents/School/MASTER_Masters/Genetic_Data/June_2014_finals/Fst_Gst_pop_pair/Pairwise_Fst_arlequin_allloci.csv")
rownames(pairwise_Fst) <- pairwise_Fst$X
pairwise_Fst <- pairwise_Fst[,-1]
pairwise_Fst <- as.matrix(pairwise_Fst)

neutral_pairwise_Fst <-  read.csv("C:/Users/Mallory/Documents/School/MASTER_Masters/Genetic_Data/June_2014_finals/Fst_Gst_pop_pair/Pairwise_Fst_arlequin_neutralLoci.csv")
rownames(neutral_pairwise_Fst) <- neutral_pairwise_Fst$X
neutral_pairwise_Fst <- neutral_pairwise_Fst[,-1]
neutral_pairwise_Fst <- as.matrix(neutral_pairwise_Fst)

#heatmaps for the matrices


heatmap.2(pairwise_Fst, Rowv = FALSE, Colv="RowV", trace="none", main="Fst")
heatmap.2(neutral_pairwise_Fst, Rowv = FALSE, Colv="RowV", trace="none", main="Neutral Fst")
heatmap.2(matrixA_final_correlation, Rowv = FALSE, Colv="RowV", trace="none", main="A")
heatmap.2(matrixB_final_correlation, Rowv = FALSE, Colv="RowV", trace="none", main="B")
heatmap.2(matrixC_final_correlation, Rowv = FALSE, Colv="RowV", trace="none", main="C")
heatmap.2(matrixD_final_correlation, Rowv = FALSE, Colv="RowV", trace="none", main="D")
heatmap.2(matrixE_final_correlation, Rowv = FALSE, Colv="RowV", trace="none", main="E")
heatmap.2(neutralMatrix_final_correlation, Rowv = FALSE, Colv="RowV", trace="none", main="Neutral")


heatmap.2(matrixA_final, Rowv = FALSE, Colv="RowV", trace="none", main="A cov")
heatmap.2(matrixB_final, Rowv = FALSE, Colv="RowV", trace="none", main="B cov")
heatmap.2(matrixC_final, Rowv = FALSE, Colv="RowV", trace="none", main="C cov")
heatmap.2(matrixD_final, Rowv = FALSE, Colv="RowV", trace="none", main="D cov")
heatmap.2(matrixE_final, Rowv = FALSE, Colv="RowV", trace="none", main="E cov")
heatmap.2(neutralMatrix_final, Rowv = FALSE, Colv="RowV", trace="none", main="Neutral cov")

####PROCESSING####
#BayEnv2 run through R command line (see codelines_Mallory.txt)
#BF results pulled from .bf files into excel file
#matched with correct loci, etc.

#set wd
setwd("~/School/MASTER_Masters/Genetic_Oceano/BayEnv_associations_results/Final_RESULTS")
setwd("~/School/MASTER_Masters/Genetic_Oceano/BayEnv_associations_results")

#load data for allele freq clustering
allele_freq_bins <- read.csv("BayEnv_freq_bins.csv")
bins <- c("A","B","C","D","E")
BinList <- list()
BinList$A <- as.vector(subset(allele_freq_bins$Locus, allele_freq_bins$Bin == "A"))
BinList$B <- as.vector(subset(allele_freq_bins$Locus, allele_freq_bins$Bin == "B"))
BinList$C <- as.vector(subset(allele_freq_bins$Locus, allele_freq_bins$Bin == "C"))
BinList$D <- as.vector(subset(allele_freq_bins$Locus, allele_freq_bins$Bin == "D"))
BinList$E <- as.vector(subset(allele_freq_bins$Locus, allele_freq_bins$Bin == "E"))



####HeatMapPCs####
#load data
HeatMapPCs_MatrixA <- read.csv("BayEnv2Results_HeatMapPCs_MatrixA.csv")
HeatMapPCs_MatrixB <- read.csv("BayEnv2Results_HeatMapPCs_MatrixB.csv")
HeatMapPCs_MatrixC <- read.csv("BayEnv2Results_HeatMapPCs_MatrixC.csv")
HeatMapPCs_MatrixD <- read.csv("BayEnv2Results_HeatMapPCs_MatrixD.csv")
HeatMapPCs_MatrixE <- read.csv("BayEnv2Results_HeatMapPCs_MatrixE.csv")
HeatMapsPCs_NeutralMatrix <- read.csv("BayEnv2Results_HeatMapPCs_NeutralMatrix.csv")
#
# #steps:
# For each matrix:
# Make range of BF values for each Env variable
# Determine the top 5% of that range of values
# Pull out the name of the loci that have BFs above that 5% cutoff

#pull out the names of the env vars
PCvars <- colnames(HeatMapPCs_MatrixA[,c(3,5,7,9,11,13,15)])

###Matrix A###
HeatMapPCs_MatrixA_Results <- list()
for (j in bins){
  HeatMapPCs_MatrixA_Results[[j]] <-list()
  for(i in PCvars){
  locusperbin <- subset(HeatMapPCs_MatrixA, HeatMapPCs_MatrixA$Locus %in% BinList[[j]])
   max <- max(locusperbin[[i]], na.rm=TRUE)
   min <- min(locusperbin[[i]], na.rm=TRUE)
   fivepercentcutoff <- max-(0.05*(max-min))
   HeatMapPCs_MatrixA_Results[[j]][[i]] <- list()
   HeatMapPCs_MatrixA_Results[[j]][[i]] <- as.vector(subset(locusperbin$Locus, 
                                                                 locusperbin[[i]] > fivepercentcutoff))
   remove(max,min,fivepercentcutoff,locusperbin)
  }
  remove(i)
}
remove(j)




#This makes a list called HeatMapPCs_MatrixA_Results, 
#and each item in the list is the env var with a single "outlier" locus listed
#double check in excel that the correct number of loci has been selected

#MatrixB
HeatMapPCs_MatrixB_Results <- list()
for (j in bins){
  HeatMapPCs_MatrixB_Results[[j]] <-list()
  for(i in PCvars){
    locusperbin <- subset(HeatMapPCs_MatrixB, HeatMapPCs_MatrixB$Locus %in% BinList[[j]])
    max <- max(locusperbin[[i]], na.rm=TRUE)
    min <- min(locusperbin[[i]], na.rm=TRUE)
    fivepercentcutoff <- max-(0.05*(max-min))
    HeatMapPCs_MatrixB_Results[[j]][[i]] <- list()
    HeatMapPCs_MatrixB_Results[[j]][[i]] <- as.vector(subset(locusperbin$Locus, 
                                                             locusperbin[[i]] > fivepercentcutoff))
    remove(max,min,fivepercentcutoff,locusperbin)
  }
  remove(i)
}
remove(j)

#MatrixC
HeatMapPCs_MatrixC_Results <- list()
for (j in bins){
  HeatMapPCs_MatrixC_Results[[j]] <-list()
  for(i in PCvars){
    locusperbin <- subset(HeatMapPCs_MatrixC, HeatMapPCs_MatrixC$Locus %in% BinList[[j]])
    max <- max(locusperbin[[i]], na.rm=TRUE)
    min <- min(locusperbin[[i]], na.rm=TRUE)
    fivepercentcutoff <- max-(0.05*(max-min))
    HeatMapPCs_MatrixC_Results[[j]][[i]] <- list()
    HeatMapPCs_MatrixC_Results[[j]][[i]] <- as.vector(subset(locusperbin$Locus, 
                                                             locusperbin[[i]] > fivepercentcutoff))
    remove(max,min,fivepercentcutoff,locusperbin)
  }
  remove(i)
}
remove(j)

#MatrixD
HeatMapPCs_MatrixD_Results <- list()
for (j in bins){
  HeatMapPCs_MatrixD_Results[[j]] <-list()
  for(i in PCvars){
    locusperbin <- subset(HeatMapPCs_MatrixD, HeatMapPCs_MatrixD$Locus %in% BinList[[j]])
    max <- max(locusperbin[[i]], na.rm=TRUE)
    min <- min(locusperbin[[i]], na.rm=TRUE)
    fivepercentcutoff <- max-(0.05*(max-min))
    HeatMapPCs_MatrixD_Results[[j]][[i]] <- list()
    HeatMapPCs_MatrixD_Results[[j]][[i]] <- as.vector(subset(locusperbin$Locus, 
                                                             locusperbin[[i]] > fivepercentcutoff))
    remove(max,min,fivepercentcutoff,locusperbin)
  }
  remove(i)
}
remove(j)

#MatrixE
HeatMapPCs_MatrixE_Results <- list()
for (j in bins){
  HeatMapPCs_MatrixE_Results[[j]] <-list()
  for(i in PCvars){
    locusperbin <- subset(HeatMapPCs_MatrixE, HeatMapPCs_MatrixE$Locus %in% BinList[[j]])
    max <- max(locusperbin[[i]], na.rm=TRUE)
    min <- min(locusperbin[[i]], na.rm=TRUE)
    fivepercentcutoff <- max-(0.05*(max-min))
    HeatMapPCs_MatrixE_Results[[j]][[i]] <- list()
    HeatMapPCs_MatrixE_Results[[j]][[i]] <- as.vector(subset(locusperbin$Locus, 
                                                             locusperbin[[i]] > fivepercentcutoff))
    remove(max,min,fivepercentcutoff,locusperbin)
  }
  remove(i)
}
remove(j)


##NeutralMatrix##
HeatMapPCs_NeutralMatrix_Results <- list()
for (j in bins){
  HeatMapPCs_NeutralMatrix_Results[[j]] <-list()
  for(i in PCvars){
    locusperbin <- subset(HeatMapsPCs_NeutralMatrix, HeatMapsPCs_NeutralMatrix$Locus %in% BinList[[j]])
    max <- max(locusperbin[[i]], na.rm=TRUE)
    min <- min(locusperbin[[i]], na.rm=TRUE)
    fivepercentcutoff <- max-(0.05*(max-min))
    HeatMapPCs_NeutralMatrix_Results[[j]][[i]] <- list()
    HeatMapPCs_NeutralMatrix_Results[[j]][[i]] <- as.vector(subset(locusperbin$Locus, 
                                                             locusperbin[[i]] > fivepercentcutoff))
    remove(max,min,fivepercentcutoff,locusperbin)
  }
  remove(i)
}
remove(j)

remove(PCvars)

####TempSalCluster####
#load data
TempSalCluster_MatrixA <- read.csv("BayEnv2Results_TempSalCluster_MatrixA.csv")
TempSalCluster_MatrixB <- read.csv("BayEnv2Results_TempSalCluster_MatrixB.csv")
TempSalCluster_MatrixC <- read.csv("BayEnv2Results_TempSalCluster_MatrixC.csv")
TempSalCluster_MatrixD <- read.csv("BayEnv2Results_TempSalCluster_MatrixD.csv")
TempSalCluster_MatrixE <- read.csv("BayEnv2Results_TempSalCluster_MatrixE.csv")
#

#pull out the names of the env vars
TempSalClustervars <- colnames(TempSalCluster_MatrixA[,c(3:17)])

###Matrix A###
TempSalCluster_MatrixA_Results <- list()
for (j in bins){
  TempSalCluster_MatrixA_Results[[j]] <-list()
  for(i in TempSalClustervars){
    locusperbin <- subset(TempSalCluster_MatrixA, TempSalCluster_MatrixA$Locus %in% BinList[[j]])
    max <- max(locusperbin[[i]], na.rm=TRUE)
    min <- min(locusperbin[[i]], na.rm=TRUE)
    fivepercentcutoff <- max-(0.05*(max-min))
    TempSalCluster_MatrixA_Results[[j]][[i]] <- list()
    TempSalCluster_MatrixA_Results[[j]][[i]] <- as.vector(subset(locusperbin$Locus, 
                                                             locusperbin[[i]] > fivepercentcutoff))
    remove(max,min,fivepercentcutoff,locusperbin)
  }
  remove(i)
}
remove(j)

#This makes a list called TempSalCluster_MatrixA_Results, 
#and each item in the list is the env var with a single "outlier" locus listed
#double check in excel that the correct number of loci has been selected

#MatrixB
TempSalCluster_MatrixB_Results <- list()
for (j in bins){
  TempSalCluster_MatrixB_Results[[j]] <-list()
  for(i in TempSalClustervars){
    locusperbin <- subset(TempSalCluster_MatrixB, TempSalCluster_MatrixB$Locus %in% BinList[[j]])
    max <- max(locusperbin[[i]], na.rm=TRUE)
    min <- min(locusperbin[[i]], na.rm=TRUE)
    fivepercentcutoff <- max-(0.05*(max-min))
    TempSalCluster_MatrixB_Results[[j]][[i]] <- list()
    TempSalCluster_MatrixB_Results[[j]][[i]] <- as.vector(subset(locusperbin$Locus, 
                                                             locusperbin[[i]] > fivepercentcutoff))
    remove(max,min,fivepercentcutoff,locusperbin)
  }
  remove(i)
}
remove(j)

#MatrixC
TempSalCluster_MatrixC_Results <- list()
for (j in bins){
  TempSalCluster_MatrixC_Results[[j]] <-list()
  for(i in TempSalClustervars){
    locusperbin <- subset(TempSalCluster_MatrixC, TempSalCluster_MatrixC$Locus %in% BinList[[j]])
    max <- max(locusperbin[[i]], na.rm=TRUE)
    min <- min(locusperbin[[i]], na.rm=TRUE)
    fivepercentcutoff <- max-(0.05*(max-min))
    TempSalCluster_MatrixC_Results[[j]][[i]] <- list()
    TempSalCluster_MatrixC_Results[[j]][[i]] <- as.vector(subset(locusperbin$Locus, 
                                                             locusperbin[[i]] > fivepercentcutoff))
    remove(max,min,fivepercentcutoff,locusperbin)
  }
  remove(i)
}
remove(j)

#MatrixD
TempSalCluster_MatrixD_Results <- list()
for (j in bins){
  TempSalCluster_MatrixD_Results[[j]] <-list()
  for(i in TempSalClustervars){
    locusperbin <- subset(TempSalCluster_MatrixD, TempSalCluster_MatrixD$Locus %in% BinList[[j]])
    max <- max(locusperbin[[i]], na.rm=TRUE)
    min <- min(locusperbin[[i]], na.rm=TRUE)
    fivepercentcutoff <- max-(0.05*(max-min))
    TempSalCluster_MatrixD_Results[[j]][[i]] <- list()
    TempSalCluster_MatrixD_Results[[j]][[i]] <- as.vector(subset(locusperbin$Locus, 
                                                             locusperbin[[i]] > fivepercentcutoff))
    remove(max,min,fivepercentcutoff,locusperbin)
  }
  remove(i)
}
remove(j)

#MatrixE
TempSalCluster_MatrixE_Results <- list()
for (j in bins){
  TempSalCluster_MatrixE_Results[[j]] <-list()
  for(i in TempSalClustervars){
    locusperbin <- subset(TempSalCluster_MatrixE, TempSalCluster_MatrixE$Locus %in% BinList[[j]])
    max <- max(locusperbin[[i]], na.rm=TRUE)
    min <- min(locusperbin[[i]], na.rm=TRUE)
    fivepercentcutoff <- max-(0.05*(max-min))
    TempSalCluster_MatrixE_Results[[j]][[i]] <- list()
    TempSalCluster_MatrixE_Results[[j]][[i]] <- as.vector(subset(locusperbin$Locus, 
                                                             locusperbin[[i]] > fivepercentcutoff))
    remove(max,min,fivepercentcutoff,locusperbin)
  }
  remove(i)
}
remove(j)

remove(TempSalClustervars)

####ChlASalTemp_nomonths####
#load data
nomonths_MatrixA <- read.csv("BayEnv2Results_nomonths_MatrixA.csv")
nomonths_MatrixB <- read.csv("BayEnv2Results_nomonths_MatrixB.csv")
nomonths_MatrixC <- read.csv("BayEnv2Results_nomonths_MatrixC.csv")
nomonths_MatrixD <- read.csv("BayEnv2Results_nomonths_MatrixD.csv")
nomonths_MatrixE <- read.csv("BayEnv2Results_nomonths_MatrixE.csv")
ChlASalTemp_nomonths_NeutralMatrix <- read.csv("BayEnv2Results_ChlASalTemp_nomonths_NeutralMatrix.csv")
#

#pull out the names of the env vars
ncol(nomonths_MatrixA)
nomonthsvars <- colnames(nomonths_MatrixA[,c(3:38)])

###Matrix A###
nomonths_MatrixA_Results <- list()
for (j in bins){
  nomonths_MatrixA_Results[[j]] <-list()
  for(i in nomonthsvars){
    locusperbin <- subset(nomonths_MatrixA, nomonths_MatrixA$Locus %in% BinList[[j]])
    max <- max(locusperbin[[i]], na.rm=TRUE)
    min <- min(locusperbin[[i]], na.rm=TRUE)
    fivepercentcutoff <- max-(0.05*(max-min))
    nomonths_MatrixA_Results[[j]][[i]] <- list()
    nomonths_MatrixA_Results[[j]][[i]] <- as.vector(subset(locusperbin$Locus, 
                                                                 locusperbin[[i]] > fivepercentcutoff))
    remove(max,min,fivepercentcutoff,locusperbin)
  }
  remove(i)
}
remove(j)

#This makes a list called nomonths_MatrixA_Results, 
#and each item in the list is the env var with a single "outlier" locus listed
#double check in excel that the correct number of loci has been selected

#MatrixB
nomonths_MatrixB_Results <- list()
for (j in bins){
  nomonths_MatrixB_Results[[j]] <-list()
  for(i in nomonthsvars){
    locusperbin <- subset(nomonths_MatrixB, nomonths_MatrixB$Locus %in% BinList[[j]])
    max <- max(locusperbin[[i]], na.rm=TRUE)
    min <- min(locusperbin[[i]], na.rm=TRUE)
    fivepercentcutoff <- max-(0.05*(max-min))
    nomonths_MatrixB_Results[[j]][[i]] <- list()
    nomonths_MatrixB_Results[[j]][[i]] <- as.vector(subset(locusperbin$Locus, 
                                                                 locusperbin[[i]] > fivepercentcutoff))
    remove(max,min,fivepercentcutoff,locusperbin)
  }
  remove(i)
}
remove(j)

#MatrixC
nomonths_MatrixC_Results <- list()
for (j in bins){
  nomonths_MatrixC_Results[[j]] <-list()
  for(i in nomonthsvars){
    locusperbin <- subset(nomonths_MatrixC, nomonths_MatrixC$Locus %in% BinList[[j]])
    max <- max(locusperbin[[i]], na.rm=TRUE)
    min <- min(locusperbin[[i]], na.rm=TRUE)
    fivepercentcutoff <- max-(0.05*(max-min))
    nomonths_MatrixC_Results[[j]][[i]] <- list()
    nomonths_MatrixC_Results[[j]][[i]] <- as.vector(subset(locusperbin$Locus, 
                                                                 locusperbin[[i]] > fivepercentcutoff))
    remove(max,min,fivepercentcutoff,locusperbin)
  }
  remove(i)
}
remove(j)

#MatrixD
nomonths_MatrixD_Results <- list()
for (j in bins){
  nomonths_MatrixD_Results[[j]] <-list()
  for(i in nomonthsvars){
    locusperbin <- subset(nomonths_MatrixD, nomonths_MatrixD$Locus %in% BinList[[j]])
    max <- max(locusperbin[[i]], na.rm=TRUE)
    min <- min(locusperbin[[i]], na.rm=TRUE)
    fivepercentcutoff <- max-(0.05*(max-min))
    nomonths_MatrixD_Results[[j]][[i]] <- list()
    nomonths_MatrixD_Results[[j]][[i]] <- as.vector(subset(locusperbin$Locus, 
                                                                 locusperbin[[i]] > fivepercentcutoff))
    remove(max,min,fivepercentcutoff,locusperbin)
  }
  remove(i)
}
remove(j)

#MatrixE
nomonths_MatrixE_Results <- list()
for (j in bins){
  nomonths_MatrixE_Results[[j]] <-list()
  for(i in nomonthsvars){
    locusperbin <- subset(nomonths_MatrixE, nomonths_MatrixE$Locus %in% BinList[[j]])
    max <- max(locusperbin[[i]], na.rm=TRUE)
    min <- min(locusperbin[[i]], na.rm=TRUE)
    fivepercentcutoff <- max-(0.05*(max-min))
    nomonths_MatrixE_Results[[j]][[i]] <- list()
    nomonths_MatrixE_Results[[j]][[i]] <- as.vector(subset(locusperbin$Locus, 
                                                                 locusperbin[[i]] > fivepercentcutoff))
    remove(max,min,fivepercentcutoff,locusperbin)
  }
  remove(i)
}
remove(j)

###NeutralMatrix###
ChlASalTemp_nomonths_NeutralMatrix_Results <- list()
for (j in bins){
  ChlASalTemp_nomonths_NeutralMatrix_Results[[j]] <-list()
  for(i in nomonthsvars){
    locusperbin <- subset(ChlASalTemp_nomonths_NeutralMatrix, 
                          ChlASalTemp_nomonths_NeutralMatrix$Locus %in% BinList[[j]])
    max <- max(locusperbin[[i]], na.rm=TRUE)
    min <- min(locusperbin[[i]], na.rm=TRUE)
    fivepercentcutoff <- max-(0.05*(max-min))
    ChlASalTemp_nomonths_NeutralMatrix_Results[[j]][[i]] <- list()
    ChlASalTemp_nomonths_NeutralMatrix_Results[[j]][[i]] <- as.vector(subset(locusperbin$Locus, 
                                                           locusperbin[[i]] > fivepercentcutoff))
    remove(max,min,fivepercentcutoff,locusperbin)
  }
  remove(i)
}
remove(j)

remove(nomonthsvars)


#####AllEnv NoMonths#####
BayEnv2Results_AllEnv_nomonths_NeutralMatrix.csv

AllEnvNoMonths_NeutralMatrix <- read.csv("BayEnv2Results_AllEnv_nomonths_NeutralMatrix.csv")

# #steps:
# For each matrix:
# Make range of BF values for each Env variable
# Determine the top 5% of that range of values
# Pull out the name of the loci that have BFs above that 5% cutoff

#pull out the names of the env vars
AllEnvNoMonthvars <- colnames(AllEnvNoMonths_NeutralMatrix[,c(3:92)])

###Neutral Matrix###
AllEnvNoMonths_NeutralMatrix_Results <- list()
for (j in bins){
  AllEnvNoMonths_NeutralMatrix_Results[[j]] <-list()
  for(i in AllEnvNoMonthvars){
    locusperbin <- subset(AllEnvNoMonths_NeutralMatrix, AllEnvNoMonths_NeutralMatrix$Locus %in% BinList[[j]])
    max <- max(locusperbin[[i]], na.rm=TRUE)
    min <- min(locusperbin[[i]], na.rm=TRUE)
    fivepercentcutoff <- max-(0.05*(max-min))
    AllEnvNoMonths_NeutralMatrix_Results[[j]][[i]] <- list()
    AllEnvNoMonths_NeutralMatrix_Results[[j]][[i]] <- as.vector(subset(locusperbin$Locus, 
                                                             locusperbin[[i]] > fivepercentcutoff))
    remove(max,min,fivepercentcutoff,locusperbin)
  }
  remove(i)
}
remove(j)

remove(AllEnvNoMonthvars)

View(AllEnvNoMonths_NeutralMatrix_Results)


###Combine results from same vars into dataframe to export####
View(HeatMapPCs_MatrixE_Results)
HeatMapPCs_allresults_bins <- rbind(HeatMapPCs_MatrixA_Results, HeatMapPCs_MatrixB_Results, HeatMapPCs_MatrixC_Results,
              HeatMapPCs_MatrixD_Results,HeatMapPCs_MatrixE_Results)

TempSalCluster_allresults_bins <- rbind(TempSalCluster_MatrixA_Results, TempSalCluster_MatrixB_Results, TempSalCluster_MatrixC_Results,
                               TempSalCluster_MatrixD_Results,TempSalCluster_MatrixE_Results)


nomonths_allresults_bins <- rbind(nomonths_MatrixA_Results, nomonths_MatrixB_Results, nomonths_MatrixC_Results,
                                   nomonths_MatrixD_Results,nomonths_MatrixE_Results)

ChlASalTemp_nomonths_NeutralMatrix_Results_bins <- rbind(ChlASalTemp_nomonths_NeutralMatrix_Results$A, ChlASalTemp_nomonths_NeutralMatrix_Results$B,
              ChlASalTemp_nomonths_NeutralMatrix_Results$C, ChlASalTemp_nomonths_NeutralMatrix_Results$D, 
              ChlASalTemp_nomonths_NeutralMatrix_Results$E)

HeatMapPCs_NeutralMatrix_Results_bins <- rbind(HeatMapPCs_NeutralMatrix_Results$A,HeatMapPCs_NeutralMatrix_Results$B,
                                               HeatMapPCs_NeutralMatrix_Results$C,HeatMapPCs_NeutralMatrix_Results$D,
                                               HeatMapPCs_NeutralMatrix_Results$E)

AllEnvNoMonths_NeutralMatrix_Results_bins <- rbind(AllEnvNoMonths_NeutralMatrix_Results$A,AllEnvNoMonths_NeutralMatrix_Results$B,
                                               AllEnvNoMonths_NeutralMatrix_Results$C,AllEnvNoMonths_NeutralMatrix_Results$D,
                                               AllEnvNoMonths_NeutralMatrix_Results$E)

#####export!####
write.csv(HeatMapPCs_allresults_bins, "HeatMapPCs_allresults_bins.csv")
write.csv(TempSalCluster_allresults_bins, "TempSalCluster_allresults_bins.csv")
write.csv(nomonths_allresults_bins, "nnomonths_allresults_bins.csv")


write.csv(ChlASalTemp_nomonths_NeutralMatrix_Results_bins, "ChlASalTemp_nomonths_NeutralMatrix_Results_bins.csv")
write.csv(HeatMapPCs_NeutralMatrix_Results_bins, "HeatMapPCs_NeutralMatrix_Results_bins.csv")
write.csv(AllEnvNoMonths_NeutralMatrix_Results_bins, "AllEnvNoMonths_NeutralMatrix_Results_bins.csv")

#####Save final#####
save.image("~/School/MASTER_Masters/R/Genetic_Environmental/BayEnv2_workspace.RData")
save.image("~/Desktop/Mallory_BayEnv_Aug2015/BayEnv2_workspace.RData")


