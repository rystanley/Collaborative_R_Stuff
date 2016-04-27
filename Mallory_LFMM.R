library("LEA", lib.loc="~/R/win-library/3.2")
library("LEA", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.1")

source("http://bioconductor.org/biocLite.R")
biocLite("LEA")
require("LEA")

load("C:/Users/Mallory/Documents/School/MASTER_Masters/R/Genetic_Environmental/LFMM_workspace.RData")


#####NOTES####
#K= number of clusters in the data, we already know we've got 2!

#####AllLociAllEnvNoMonths####
#setwd
setwd("~/School/MASTER_Masters/Genetic_Oceano/Outliers/LFMM_associations_results/AllLociAllEnv_NoMonths_Windows")

#linux
setwd("~/Desktop/AllLociAllEnv_NoMonths_Mallory/")
#gen file - All_loci_HWE.lfmm
#env file - STAND_oceano_data_NoCorrTest_noGeo_FINAL_NoNA_LFMM.env

allloci_allenvNoMonths.lfmm  <- lfmm(input.file="All_loci_HWE.lfmm", 
                             environment.file="STAND_oceano_data_NoCorrTest_noGeo_FINAL_NoNA_NoMonths_LFMM.env", 
                             K = 1:4, rep = 3, burnin=5000, iterations=15000, project = "new", CPU=30)



allloci_allenvNoMonths_windows.lfmm  <- lfmm(input.file="All_loci_HWE.lfmm", 
                                     environment.file="STAND_oceano_data_NoCorrTest_noGeo_FINAL_NoNA_NoMonths_LFMM.env", 
                                     K = 2, rep = 3, burnin=5000, iterations=15000, project = "new", CPU=1)



#####AllLociChlASalTemp####
#setwd
setwd("~/School/MASTER_Masters/Genetic_Oceano/LFMM_associations_results/AllLociChlASalTemp/")

#linux
setwd("~/Desktop/Mallory/LFMM_associations_Aug11_2015_1245/AllLociChlASalTemp/")
#gen file - All_loci_HWE.lfmm
#env file - STAND_oceano_data_ChlATempSal_nomonths_LFMM.env

allloci_ChlASalTemp.lfmm  <- lfmm(input.file="All_loci_HWE.lfmm", 
                                  environment.file="STAND_oceano_data_ChlATempSal_nomonths_LFMM.env", 
                                  K = 1:6, rep = 3, burnin=5000, iterations=15000, project = "new",
                                  CPU=20)

####AllLociHeatMapTempSalCluster####
#setwd
setwd("~/School/MASTER_Masters/Genetic_Oceano/LFMM_associations_results/AllLociHeatMapTempSalCluster/")

#linux
setwd("~/Desktop/Mallory/LFMM_associations_Aug11_2015_1245/AllLociHeatMapTempSalCluster/")
#gen file - All_loci_HWE.lfmm
#env file - STAND_oceano_data_ChlASalTemp_HeatMap_TempSalCluster_LFMM.env

allloci_HeatMapTempSalCluster.lfmm  <- lfmm(input.file="All_loci_HWE.lfmm", 
                                            environment.file="STAND_oceano_data_ChlASalTemp_HeatMap_TempSalCluster_LFMM.env", 
                                            K = 1:6, rep = 3, burnin=5000, iterations=15000, project = "new", CPU=20)

#####AllLociPCs####
#setwd
setwd("~/School/MASTER_Masters/Genetic_Oceano/LFMM_associations_results/AllLociPCs/")

#linux
#setwd("~/Desktop/Mallory/LFMM_associations_Aug11_2015_1245/AllLociPCs/")

#gen file - All_loci_HWE.lfmm
#env file - STAND_oceano_data_ChlASalTemp_NoMonth_HeatMapPCs_LFMM.env

allloci_PCs.lfmm  <- lfmm(input.file="All_loci_HWE.lfmm", 
                          environment.file="STAND_oceano_data_ChlASalTemp_NoMonth_HeatMapPCs_LFMM.env", 
                          K = 1:6, rep = 3, burnin=5000, iterations=15000, project = "new")





####Processing results#####

#####AllLociAllEnv --> didn't run properly, can't process anything because it didn't create the project file####

setwd("~/School/MASTER_Masters/Genetic_Oceano/LFMM_associations_results/AllLociAllEnv")

AllLociAllEnv_project = load.lfmmProject("All_loci_HWE_STAND_oceano_data_ChlATempSal_nomonths_LFMM.lfmmProject")


###working to see if I can get the d forloop nested within the K forloop
for (j in 1:6){
  assign(paste0("test.k", j), list())
  for (i in 1:7){
    name <- paste0("d",i)
    kname <- paste0(j)
    temp_score <- z.scores(AllLociPCs_project, K=j, d=i)
    paste0("test.k", j)[[name]] <- temp_score
    remove(temp_score, name,kname)
  }
  remove(i)    
}
remove(j)


remove(test.k1,test.k2,test.k4,test.k3,test.k5,test.k6)

#####AllLociAllEnvNoMonths#####
setwd("~/School/MASTER_Masters/Genetic_Oceano/Outliers/LFMM_associations_results/AllLociAllEnv_NoMonths_Windows/")

AllLociAllEnvNoMonths_project = load.lfmmProject("All_loci_HWE_STAND_oceano_data_NoCorrTest_noGeo_FINAL_NoNA_NoMonths_LFMM.lfmmProject")

###z-scores###
#make a list of 7163 (#loci) z.scores for each env var (d) for a specific K

###k=2##
AllLociAllEnvNoMonths_project.zs.k2 <- list()
for (i in 1:90){
  name <- paste0("d",i)
  temp_score <- z.scores(AllLociAllEnvNoMonths_project, K=2, d=i)
  AllLociAllEnvNoMonths_project.zs.k2[[name]] <- temp_score
  remove(temp_score, name)
}
remove(i)
##combine z-scores from the 3 runs##
AllLociAllEnvNoMonths_project.combinedzs.k2 <- list()
AllLociAllEnvNoMonths_project.combinedzs.k2$d1 <- apply(AllLociAllEnvNoMonths_project.zs.k2$d1, MARGIN = 1, median)
#...see AllLociAllEnvNoMonths_CodeLines.xlsx



###then determine correct # of K using GIF
#calculate genomic inflation factor from z-scores
##here, I want one data frame, 90 rows, 1 columns, K=2 as columns, d=1:90 as rows
AllLociAllEnvNoMonths_project.lambda <- data.frame(matrix(nrow=90,ncol=1))
colnames(AllLociAllEnvNoMonths_project.lambda) <- c("K2")
AllLociAllEnvNoMonths_project.lambda[1,1] <- median(AllLociAllEnvNoMonths_project.combinedzs.k1$d1^2)/.456
AllLociAllEnvNoMonths_project.lambda[2,1] <- median(AllLociAllEnvNoMonths_project.combinedzs.k1$d3^2)/.456
#...see AllLociAllEnvNoMonths_CodeLines.xlsx

View(AllLociAllEnvNoMonths_project.lambda)

#average GIF over all env vars?
AllLociAllEnvNoMonths_project.averagelambda <- data.frame(matrix(nrow=1,ncol=1))
AllLociAllEnvNoMonths_project.averagelambda[1,1] <- mean(AllLociAllEnvNoMonths_project.lambda[1:90,1])

###K=2 is the best###

#calculate adjusted p-values
AllLociAllEnvNoMonths_project.adjpvals.k2 <- list()
AllLociAllEnvNoMonths_project.adjpvals.k2$d1 <- pchisq((AllLociAllEnvNoMonths_project.combinedzs.k2$d1^2)/AllLociAllEnvNoMonths_project.lambda[1,1],df=1, lower=FALSE)
#...see AllLociAllEnvNoMonths_CodeLines.xlsx

#list of candidate loci, with correction for multiple testing (returned as AllLociAllEnvNoMonths.d1.candidates)

#for d=1:90
L = 7163
#return a list of candidates with an expected FDR of alpha
w = which(sort(AllLociAllEnvNoMonths_project.adjpvals.k2$d2) < 0.05 * (1:L) / L)
AllLociAllEnvNoMonths.d2.candidates = order(AllLociAllEnvNoMonths_project.adjpvals.k2$d2) [w]
AllLociAllEnvNoMonths.d2.candidates


#####AllLociChlASalTemp####
setwd("~/School/MASTER_Masters/Genetic_Oceano/LFMM_associations_results/AllLociChlASalTemp/")

AllLociChlASalTemp_project = load.lfmmProject("All_loci_HWE_STAND_oceano_data_ChlATempSal_nomonths_LFMM.lfmmProject")

###z-scores###
#make a list of 7163 (#loci) z.scores for each env var (d) for a specific K

###k=1###
AllLociChlASalTemp_project.zs.k1 <- list()
for (i in 1:36){
  name <- paste0("d",i)
  temp_score <- z.scores(AllLociChlASalTemp_project, K=1, d=i)
  AllLociChlASalTemp_project.zs.k1[[name]] <- temp_score
  remove(temp_score, name)
}
remove(i)
##combine z-scores from the 3 runs##
AllLociChlASalTemp_project.combinedzs.k1 <- list()
AllLociChlASalTemp_project.combinedzs.k1$d1 <- apply(AllLociChlASalTemp_project.zs.k1$d1, MARGIN = 1, median)
#...see AllLociChlASalTemp_CodeLines.xlsx

###k=2##
AllLociChlASalTemp_project.zs.k2 <- list()
for (i in 1:36){
  name <- paste0("d",i)
  temp_score <- z.scores(AllLociChlASalTemp_project, K=2, d=i)
  AllLociChlASalTemp_project.zs.k2[[name]] <- temp_score
  remove(temp_score, name)
}
remove(i)
##combine z-scores from the 3 runs##
AllLociChlASalTemp_project.combinedzs.k2 <- list()
AllLociChlASalTemp_project.combinedzs.k2$d1 <- apply(AllLociChlASalTemp_project.zs.k2$d1, MARGIN = 1, median)
#...see AllLociChlASalTemp_CodeLines.xlsx

###k3##
AllLociChlASalTemp_project.zs.k3 <- list()
for (i in 1:36){
  name <- paste0("d",i)
  temp_score <- z.scores(AllLociChlASalTemp_project, K=3, d=i)
  AllLociChlASalTemp_project.zs.k3[[name]] <- temp_score
  remove(temp_score, name)
}
remove(i)
##combine z-scores from the 3 runs##
AllLociChlASalTemp_project.combinedzs.k3 <- list()
AllLociChlASalTemp_project.combinedzs.k3$d1 <- apply(AllLociChlASalTemp_project.zs.k3$d1, MARGIN = 1, median)
#...see AllLociChlASalTemp_CodeLines.xlsx


###k4###
AllLociChlASalTemp_project.zs.k4 <- list()
for (i in 1:36){
  name <- paste0("d",i)
  temp_score <- z.scores(AllLociChlASalTemp_project, K=4, d=i)
  AllLociChlASalTemp_project.zs.k4[[name]] <- temp_score
  remove(temp_score, name)
}
remove(i)
##combine z-scores from the 3 runs##
AllLociChlASalTemp_project.combinedzs.k4 <- list()
AllLociChlASalTemp_project.combinedzs.k4$d1 <- apply(AllLociChlASalTemp_project.zs.k4$d1, MARGIN = 1, median)
#...see AllLociChlASalTemp_CodeLines.xlsx

###k5###
AllLociChlASalTemp_project.zs.k5 <- list()
for (i in 1:36){
  name <- paste0("d",i)
  temp_score <- z.scores(AllLociChlASalTemp_project, K=5, d=i)
  AllLociChlASalTemp_project.zs.k5[[name]] <- temp_score
  remove(temp_score, name)
}
remove(i)
##combine z-scores from the 3 runs##
AllLociChlASalTemp_project.combinedzs.k5 <- list()
AllLociChlASalTemp_project.combinedzs.k5$d1 <- apply(AllLociChlASalTemp_project.zs.k5$d1, MARGIN = 1, median)
#...see AllLociChlASalTemp_CodeLines.xlsx

###k6###
AllLociChlASalTemp_project.zs.k6 <- list()
for (i in 1:36){
  name <- paste0("d",i)
  temp_score <- z.scores(AllLociChlASalTemp_project, K=6, d=i)
  AllLociChlASalTemp_project.zs.k6[[name]] <- temp_score
  remove(temp_score, name)
}
remove(i)
##combine z-scores from the 3 runs##
AllLociChlASalTemp_project.combinedzs.k6 <- list()
AllLociChlASalTemp_project.combinedzs.k6$d1 <- apply(AllLociChlASalTemp_project.zs.k6$d1, MARGIN = 1, median)
#...see AllLociChlASalTemp_CodeLines.xlsx


###then determine correct # of K using GIF
#calculate genomic inflation factor from z-scores
##here, I want one data frame, 36 rows, 6 columns, K=1:6 as columns, d=1:36 as rows
AllLociChlASalTemp_project.lambda <- data.frame(matrix(nrow=36,ncol=6))
colnames(AllLociChlASalTemp_project.lambda) <- c("K1","K2","K3","K4","K5","K6")
AllLociChlASalTemp_project.lambda[1,1] <- median(AllLociChlASalTemp_project.combinedzs.k1$d1^2)/.456
AllLociChlASalTemp_project.lambda[2,1] <- median(AllLociChlASalTemp_project.combinedzs.k1$d3^2)/.456
#...see AllLociChlASalTemp_CodeLines.xlsx

View(AllLociChlASalTemp_project.lambda)

#average GIF over all env vars?
AllLociChlASalTemp_project.averagelambda <- data.frame(matrix(nrow=1,ncol=6))
AllLociChlASalTemp_project.averagelambda[1,1] <- mean(AllLociChlASalTemp_project.lambda[1:36,1])
AllLociChlASalTemp_project.averagelambda[1,2] <- mean(AllLociChlASalTemp_project.lambda[1:36,2])
AllLociChlASalTemp_project.averagelambda[1,3] <- mean(AllLociChlASalTemp_project.lambda[1:36,3])
AllLociChlASalTemp_project.averagelambda[1,4] <- mean(AllLociChlASalTemp_project.lambda[1:36,4])
AllLociChlASalTemp_project.averagelambda[1,5] <- mean(AllLociChlASalTemp_project.lambda[1:36,5])
AllLociChlASalTemp_project.averagelambda[1,6] <- mean(AllLociChlASalTemp_project.lambda[1:36,6])

###K=2 is the best###

#calculate adjusted p-values
AllLociChlASalTemp_project.adjpvals.k2 <- list()
AllLociChlASalTemp_project.adjpvals.k2$d1 <- pchisq((AllLociChlASalTemp_project.combinedzs.k2$d1^2)/AllLociChlASalTemp_project.lambda[1,2],df=1, lower=FALSE)
#...see AllLociChlASalTemp_CodeLines.xlsx

#list of candidate loci, with correction for multiple testing (returned as AllLociChlASalTemp.d1.candidates)

L = 7163
#return a list of candidates with an expected FDR of alpha
w = which(sort(AllLociChlASalTemp_project.adjpvals.k2$d36) < 0.05 * (1:L) / L)
AllLociChlASalTemp.d36.candidates = order(AllLociChlASalTemp_project.adjpvals.k2$d36) [w]
AllLociChlASalTemp.d36.candidates
  

####AllLociHeatMapTempSalCluster####
setwd("~/School/MASTER_Masters/Genetic_Oceano/LFMM_associations_results/AllLociHeatMapTempSalCluster/")

#load project
AllLociHeatMapTempSalCluster_project = load.lfmmProject("All_loci_HWE_STAND_oceano_data_ChlASalTemp_HeatMap_TempSalCluster_LFMM.lfmmProject")

###z-scores###
#make a list of 7163 (#loci) z.scores for each env var (d) for a specific K

###k=1###
AllLociHeatMapTempSalCluster_project.zs.k1 <- list()
for (i in 1:15){
  name <- paste0("d",i)
  temp_score <- z.scores(AllLociHeatMapTempSalCluster_project, K=1, d=i)
  AllLociHeatMapTempSalCluster_project.zs.k1[[name]] <- temp_score
  remove(temp_score, name)
}
remove(i)
##combine z-scores from the 3 runs##
AllLociHeatMapTempSalCluster_project.combinedzs.k1 <- list()
AllLociHeatMapTempSalCluster_project.combinedzs.k1$d1 <- apply(AllLociHeatMapTempSalCluster_project.zs.k1$d1, MARGIN = 1, median)
#...see AllLociHeatMapTempSalCluster_CodeLines.xlsx

###k=2##
AllLociHeatMapTempSalCluster_project.zs.k2 <- list()
for (i in 1:15){
  name <- paste0("d",i)
  temp_score <- z.scores(AllLociHeatMapTempSalCluster_project, K=2, d=i)
  AllLociHeatMapTempSalCluster_project.zs.k2[[name]] <- temp_score
  remove(temp_score, name)
}
remove(i)
##combine z-scores from the 3 runs##
AllLociHeatMapTempSalCluster_project.combinedzs.k2 <- list()
AllLociHeatMapTempSalCluster_project.combinedzs.k2$d1 <- apply(AllLociHeatMapTempSalCluster_project.zs.k2$d1, MARGIN = 1, median)
#...see AllLociHeatMapTempSalCluster_CodeLines.xlsx

###k3##
AllLociHeatMapTempSalCluster_project.zs.k3 <- list()
for (i in 1:15){
  name <- paste0("d",i)
  temp_score <- z.scores(AllLociHeatMapTempSalCluster_project, K=3, d=i)
  AllLociHeatMapTempSalCluster_project.zs.k3[[name]] <- temp_score
  remove(temp_score, name)
}
remove(i)
##combine z-scores from the 3 runs##
AllLociHeatMapTempSalCluster_project.combinedzs.k3 <- list()
AllLociHeatMapTempSalCluster_project.combinedzs.k3$d1 <- apply(AllLociHeatMapTempSalCluster_project.zs.k3$d1, MARGIN = 1, median)
#...see AllLociHeatMapTempSalCluster_CodeLines.xlsx


###k4###
AllLociHeatMapTempSalCluster_project.zs.k4 <- list()
for (i in 1:15){
  name <- paste0("d",i)
  temp_score <- z.scores(AllLociHeatMapTempSalCluster_project, K=4, d=i)
  AllLociHeatMapTempSalCluster_project.zs.k4[[name]] <- temp_score
  remove(temp_score, name)
}
remove(i)
##combine z-scores from the 3 runs##
AllLociHeatMapTempSalCluster_project.combinedzs.k4 <- list()
AllLociHeatMapTempSalCluster_project.combinedzs.k4$d1 <- apply(AllLociHeatMapTempSalCluster_project.zs.k4$d1, MARGIN = 1, median)
#...see AllLociHeatMapTempSalCluster_CodeLines.xlsx

###k5###
AllLociHeatMapTempSalCluster_project.zs.k5 <- list()
for (i in 1:15){
  name <- paste0("d",i)
  temp_score <- z.scores(AllLociHeatMapTempSalCluster_project, K=5, d=i)
  AllLociHeatMapTempSalCluster_project.zs.k5[[name]] <- temp_score
  remove(temp_score, name)
}
remove(i)
##combine z-scores from the 3 runs##
AllLociHeatMapTempSalCluster_project.combinedzs.k5 <- list()
AllLociHeatMapTempSalCluster_project.combinedzs.k5$d1 <- apply(AllLociHeatMapTempSalCluster_project.zs.k5$d1, MARGIN = 1, median)
#...see AllLociHeatMapTempSalCluster_CodeLines.xlsx

###k6###
AllLociHeatMapTempSalCluster_project.zs.k6 <- list()
for (i in 1:15){
  name <- paste0("d",i)
  temp_score <- z.scores(AllLociHeatMapTempSalCluster_project, K=6, d=i)
  AllLociHeatMapTempSalCluster_project.zs.k6[[name]] <- temp_score
  remove(temp_score, name)
}
remove(i)
##combine z-scores from the 3 runs##
AllLociHeatMapTempSalCluster_project.combinedzs.k6 <- list()
AllLociHeatMapTempSalCluster_project.combinedzs.k6$d1 <- apply(AllLociHeatMapTempSalCluster_project.zs.k6$d1, MARGIN = 1, median)
#...see AllLociHeatMapTempSalCluster_CodeLines.xlsx

###then determine correct # of K using GIF
#calculate genomic inflation factor from z-scores
##here, I want one data frame, 15 rows, 6 columns, K=1:6 as columns, d=1:15 as rows
AllLociHeatMapTempSalCluster_project.lambda <- data.frame(matrix(nrow=15,ncol=6))
colnames(AllLociHeatMapTempSalCluster_project.lambda) <- c("K1","K2","K3","K4","K5","K6")
AllLociHeatMapTempSalCluster_project.lambda[1,1] <- median(AllLociHeatMapTempSalCluster_project.combinedzs.k1$d1^2)/.456
AllLociHeatMapTempSalCluster_project.lambda[2,1] <- median(AllLociHeatMapTempSalCluster_project.combinedzs.k1$d2^2)/.456
#...see AllLociHeatMapTempSalCluster_CodeLines.xlsx

View(AllLociHeatMapTempSalCluster_project.lambda)

#average GIF over all env vars?
AllLociHeatMapTempSalCluster_project.averagelambda <- data.frame(matrix(nrow=1,ncol=6))
AllLociHeatMapTempSalCluster_project.averagelambda[1,1] <- mean(AllLociHeatMapTempSalCluster_project.lambda[1:15,1])
AllLociHeatMapTempSalCluster_project.averagelambda[1,2] <- mean(AllLociHeatMapTempSalCluster_project.lambda[1:15,2])
AllLociHeatMapTempSalCluster_project.averagelambda[1,3] <- mean(AllLociHeatMapTempSalCluster_project.lambda[1:15,3])
AllLociHeatMapTempSalCluster_project.averagelambda[1,4] <- mean(AllLociHeatMapTempSalCluster_project.lambda[1:15,4])
AllLociHeatMapTempSalCluster_project.averagelambda[1,5] <- mean(AllLociHeatMapTempSalCluster_project.lambda[1:15,5])
AllLociHeatMapTempSalCluster_project.averagelambda[1,6] <- mean(AllLociHeatMapTempSalCluster_project.lambda[1:15,6])

###K=2 is the best###


#calculate adjusted p-values
AllLociHeatMapTempSalCluster_project.adjpvals.k2 <- list()
AllLociHeatMapTempSalCluster_project.adjpvals.k2$d1 <- pchisq((AllLociHeatMapTempSalCluster_project.combinedzs.k2$d1^2)/AllLociHeatMapTempSalCluster_project.lambda[1,2],df=1, lower=FALSE)
#...see AllLociHeatMapTempSalCluster_CodeLines.xlsx

#list of candidate loci, with correction for multiple testing (returned as AllLociHeatMapTempSalCluster.d1.candidates)

L = 7163
#return a list of candidates with an expected FDR of alpha
w = which(sort(AllLociHeatMapTempSalCluster_project.adjpvals.k2$d15) < 0.05 * (1:L) / L)
AllLociHeatMapTempSalCluster.d15.candidates = order(AllLociHeatMapTempSalCluster_project.adjpvals.k2$d15) [w]
AllLociHeatMapTempSalCluster.d15.candidates

#####AllLociPCs####
setwd("~/School/MASTER_Masters/Genetic_Oceano/LFMM_associations_results/AllLociPCs/")

#load project
AllLociPCs_project = load.lfmmProject("All_loci_HWE_STAND_oceano_data_ChlASalTemp_NoMonth_HeatMapPCs_LFMM.lfmmProject")


###z-scores###
#make a list of 7163 (#loci) z.scores for each env var (d) for a specific K

###k=1###
AllLociPCs_project.zs.k1 <- list()
for (i in 1:7){
  name <- paste0("d",i)
  temp_score <- z.scores(AllLociPCs_project, K=1, d=i)
  AllLociPCs_project.zs.k1[[name]] <- temp_score
  remove(temp_score, name)
}
remove(i)
##combine z-scores from the 3 runs##
AllLociPCs_project.combinedzs.k1 <- list()
AllLociPCs_project.combinedzs.k1$d1 <- apply(AllLociPCs_project.zs.k1$d1, MARGIN = 1, median)
AllLociPCs_project.combinedzs.k1$d2 <- apply(AllLociPCs_project.zs.k1$d2, MARGIN = 1, median)
AllLociPCs_project.combinedzs.k1$d3 <- apply(AllLociPCs_project.zs.k1$d3, MARGIN = 1, median)
AllLociPCs_project.combinedzs.k1$d4 <- apply(AllLociPCs_project.zs.k1$d4, MARGIN = 1, median)
AllLociPCs_project.combinedzs.k1$d5 <- apply(AllLociPCs_project.zs.k1$d5, MARGIN = 1, median)
AllLociPCs_project.combinedzs.k1$d6 <- apply(AllLociPCs_project.zs.k1$d6, MARGIN = 1, median)
AllLociPCs_project.combinedzs.k1$d7 <- apply(AllLociPCs_project.zs.k1$d7, MARGIN = 1, median)
#...see AllLociPCs_CodeLines.xlsx

###k=2##
AllLociPCs_project.zs.k2 <- list()
for (i in 1:7){
  name <- paste0("d",i)
  temp_score <- z.scores(AllLociPCs_project, K=2, d=i)
  AllLociPCs_project.zs.k2[[name]] <- temp_score
  remove(temp_score, name)
}
remove(i)
##combine z-scores from the 3 runs##
AllLociPCs_project.combinedzs.k2 <- list()
AllLociPCs_project.combinedzs.k2$d1 <- apply(AllLociPCs_project.zs.k2$d1, MARGIN = 1, median)
AllLociPCs_project.combinedzs.k2$d2 <- apply(AllLociPCs_project.zs.k2$d2, MARGIN = 1, median)
AllLociPCs_project.combinedzs.k2$d3 <- apply(AllLociPCs_project.zs.k2$d3, MARGIN = 1, median)
AllLociPCs_project.combinedzs.k2$d4 <- apply(AllLociPCs_project.zs.k2$d4, MARGIN = 1, median)
AllLociPCs_project.combinedzs.k2$d5 <- apply(AllLociPCs_project.zs.k2$d5, MARGIN = 1, median)
AllLociPCs_project.combinedzs.k2$d6 <- apply(AllLociPCs_project.zs.k2$d6, MARGIN = 1, median)
AllLociPCs_project.combinedzs.k2$d7 <- apply(AllLociPCs_project.zs.k2$d7, MARGIN = 1, median)

#...see AllLociPCs_CodeLines.xlsx

###k3##
AllLociPCs_project.zs.k3 <- list()
for (i in 1:7){
  name <- paste0("d",i)
  temp_score <- z.scores(AllLociPCs_project, K=3, d=i)
  AllLociPCs_project.zs.k3[[name]] <- temp_score
  remove(temp_score, name)
}
remove(i)
##combine z-scores from the 3 runs##
AllLociPCs_project.combinedzs.k3 <- list()
AllLociPCs_project.combinedzs.k3$d1 <- apply(AllLociPCs_project.zs.k3$d1, MARGIN = 1, median)
AllLociPCs_project.combinedzs.k3$d2 <- apply(AllLociPCs_project.zs.k3$d2, MARGIN = 1, median)
AllLociPCs_project.combinedzs.k3$d3 <- apply(AllLociPCs_project.zs.k3$d3, MARGIN = 1, median)
AllLociPCs_project.combinedzs.k3$d4 <- apply(AllLociPCs_project.zs.k3$d4, MARGIN = 1, median)
AllLociPCs_project.combinedzs.k3$d5 <- apply(AllLociPCs_project.zs.k3$d5, MARGIN = 1, median)
AllLociPCs_project.combinedzs.k3$d6 <- apply(AllLociPCs_project.zs.k3$d6, MARGIN = 1, median)
AllLociPCs_project.combinedzs.k3$d7 <- apply(AllLociPCs_project.zs.k3$d7, MARGIN = 1, median)

#...see AllLociPCs_CodeLines.xlsx


###k4###
AllLociPCs_project.zs.k4 <- list()
for (i in 1:7){
  name <- paste0("d",i)
  temp_score <- z.scores(AllLociPCs_project, K=4, d=i)
  AllLociPCs_project.zs.k4[[name]] <- temp_score
  remove(temp_score, name)
}
remove(i)
##combine z-scores from the 3 runs##
AllLociPCs_project.combinedzs.k4 <- list()
AllLociPCs_project.combinedzs.k4$d1 <- apply(AllLociPCs_project.zs.k4$d1, MARGIN = 1, median)
AllLociPCs_project.combinedzs.k4$d2 <- apply(AllLociPCs_project.zs.k4$d2, MARGIN = 1, median)
AllLociPCs_project.combinedzs.k4$d3 <- apply(AllLociPCs_project.zs.k4$d3, MARGIN = 1, median)
AllLociPCs_project.combinedzs.k4$d4 <- apply(AllLociPCs_project.zs.k4$d4, MARGIN = 1, median)
AllLociPCs_project.combinedzs.k4$d5 <- apply(AllLociPCs_project.zs.k4$d5, MARGIN = 1, median)
AllLociPCs_project.combinedzs.k4$d6 <- apply(AllLociPCs_project.zs.k4$d6, MARGIN = 1, median)
AllLociPCs_project.combinedzs.k4$d7 <- apply(AllLociPCs_project.zs.k4$d7, MARGIN = 1, median)

#...see AllLociPCs_CodeLines.xlsx

###k5###
AllLociPCs_project.zs.k5 <- list()
for (i in 1:7){
  name <- paste0("d",i)
  temp_score <- z.scores(AllLociPCs_project, K=5, d=i)
  AllLociPCs_project.zs.k5[[name]] <- temp_score
  remove(temp_score, name)
}
remove(i)
##combine z-scores from the 3 runs##
AllLociPCs_project.combinedzs.k5 <- list()
AllLociPCs_project.combinedzs.k5$d1 <- apply(AllLociPCs_project.zs.k5$d1, MARGIN = 1, median)
AllLociPCs_project.combinedzs.k5$d2 <- apply(AllLociPCs_project.zs.k5$d2, MARGIN = 1, median)
AllLociPCs_project.combinedzs.k5$d3 <- apply(AllLociPCs_project.zs.k5$d3, MARGIN = 1, median)
AllLociPCs_project.combinedzs.k5$d4 <- apply(AllLociPCs_project.zs.k5$d4, MARGIN = 1, median)
AllLociPCs_project.combinedzs.k5$d5 <- apply(AllLociPCs_project.zs.k5$d5, MARGIN = 1, median)
AllLociPCs_project.combinedzs.k5$d6 <- apply(AllLociPCs_project.zs.k5$d6, MARGIN = 1, median)
AllLociPCs_project.combinedzs.k5$d7 <- apply(AllLociPCs_project.zs.k5$d7, MARGIN = 1, median)

#...see AllLociPCs_CodeLines.xlsx

###k6###
AllLociPCs_project.zs.k6 <- list()
for (i in 1:7){
  name <- paste0("d",i)
  temp_score <- z.scores(AllLociPCs_project, K=6, d=i)
  AllLociPCs_project.zs.k6[[name]] <- temp_score
  remove(temp_score, name)
}
remove(i)
##combine z-scores from the 3 runs##
AllLociPCs_project.combinedzs.k6 <- list()
AllLociPCs_project.combinedzs.k6$d1 <- apply(AllLociPCs_project.zs.k6$d1, MARGIN = 1, median)
AllLociPCs_project.combinedzs.k6$d2 <- apply(AllLociPCs_project.zs.k6$d2, MARGIN = 1, median)
AllLociPCs_project.combinedzs.k6$d3 <- apply(AllLociPCs_project.zs.k6$d3, MARGIN = 1, median)
AllLociPCs_project.combinedzs.k6$d4 <- apply(AllLociPCs_project.zs.k6$d4, MARGIN = 1, median)
AllLociPCs_project.combinedzs.k6$d5 <- apply(AllLociPCs_project.zs.k6$d5, MARGIN = 1, median)
AllLociPCs_project.combinedzs.k6$d6 <- apply(AllLociPCs_project.zs.k6$d6, MARGIN = 1, median)
AllLociPCs_project.combinedzs.k6$d7 <- apply(AllLociPCs_project.zs.k6$d7, MARGIN = 1, median)

#...see AllLociPCs_CodeLines.xlsx

###then determine correct # of K using GIF
#calculate genomic inflation factor from z-scores
##here, I want one data frame, 7 rows, 6 columns, K=1:6 as columns, d=1:7 as rows
AllLociPCs_project.lambda <- data.frame(matrix(nrow=7,ncol=6))
colnames(AllLociPCs_project.lambda) <- c("K1","K2","K3","K4","K5","K6")
AllLociPCs_project.lambda[1,1] <- median(AllLociPCs_project.combinedzs.k1$d1^2)/.456
AllLociPCs_project.lambda[2,1] <- median(AllLociPCs_project.combinedzs.k1$d2^2)/.456
#...see AllLociPCs_CodeLines.xlsx

View(AllLociPCs_project.lambda)

#average GIF over all env vars?
AllLociPCs_project.averagelambda <- data.frame(matrix(nrow=1,ncol=6))
AllLociPCs_project.averagelambda[1,1] <- mean(AllLociPCs_project.lambda[1:7,1])
AllLociPCs_project.averagelambda[1,2] <- mean(AllLociPCs_project.lambda[1:7,2])
AllLociPCs_project.averagelambda[1,3] <- mean(AllLociPCs_project.lambda[1:7,3])
AllLociPCs_project.averagelambda[1,4] <- mean(AllLociPCs_project.lambda[1:7,4])
AllLociPCs_project.averagelambda[1,5] <- mean(AllLociPCs_project.lambda[1:7,5])
AllLociPCs_project.averagelambda[1,6] <- mean(AllLociPCs_project.lambda[1:7,6])

View(AllLociPCs_project.averagelambda)

###K=2 is the best###


#calculate adjusted p-values
AllLociPCs_project.adjpvals.k2 <- list()
AllLociPCs_project.adjpvals.k2$d1 <- pchisq((AllLociPCs_project.combinedzs.k2$d1^2)/AllLociPCs_project.lambda[1,2],df=1, lower=FALSE)
AllLociPCs_project.adjpvals.k2$d2 <- pchisq((AllLociPCs_project.combinedzs.k2$d2^2)/AllLociPCs_project.lambda[2,2],df=1, lower=FALSE)
AllLociPCs_project.adjpvals.k2$d3 <- pchisq((AllLociPCs_project.combinedzs.k2$d3^2)/AllLociPCs_project.lambda[3,2],df=1, lower=FALSE)
AllLociPCs_project.adjpvals.k2$d4 <- pchisq((AllLociPCs_project.combinedzs.k2$d4^2)/AllLociPCs_project.lambda[4,2],df=1, lower=FALSE)
AllLociPCs_project.adjpvals.k2$d5 <- pchisq((AllLociPCs_project.combinedzs.k2$d5^2)/AllLociPCs_project.lambda[5,2],df=1, lower=FALSE)
AllLociPCs_project.adjpvals.k2$d6 <- pchisq((AllLociPCs_project.combinedzs.k2$d6^2)/AllLociPCs_project.lambda[6,2],df=1, lower=FALSE)
AllLociPCs_project.adjpvals.k2$d7 <- pchisq((AllLociPCs_project.combinedzs.k2$d7^2)/AllLociPCs_project.lambda[7,2],df=1, lower=FALSE)

#...see AllLociPCs_CodeLines.xlsx

#list of candidate loci, with correction for multiple testing (returned as AllLociPCs.d1.candidates)

L = 7163
#return a list of candidates with an expected FDR of alpha
w = which(sort(AllLociPCs_project.adjpvals.k2$d15) < 0.05 * (1:L) / L)
AllLociPCs.d15.candidates = order(AllLociPCs_project.adjpvals.k2$d15) [w]
AllLociPCs.d15.candidates


#####FINAL SAVE####
save.image("~/School/MASTER_Masters/R/Genetic_Environmental/LFMM_workspace.RData")
save.image("~/Desktop/AllLociAllEnv_Mallory/LFMM_workspace.RData")
