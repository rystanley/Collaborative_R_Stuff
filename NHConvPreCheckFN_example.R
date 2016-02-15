## Check for convergence FN
PreDir <- "/Users/brendanwringe/Desktop/DFO Aquaculture Interaction/South West Rivers Analysis/Frequency Based Sim/West/TestCheck/"
  PreDir <- "/Users/brendanwringe/Desktop/DFO Aquaculture Interaction/South West Rivers Analysis/Frequency Based Sim/West/West Correct Wild Fst Top Loci/West - ReRun Feb 9 2016/"

  
preCheckFN <- function(PreDir){

  ### NewHybrids has a sweet habit of failing to converge from time to time. The issue seems to be that sometimes the seeds given to it cause the MCMC chains to begin so far
    ## out in the boonies not even the Littlest Hobo could get them home. 
  ## When convergence fails in this manner, it is very easy to diagnose - NH says, "to hell with it, everyone's an F2 and I'm out" 
  ## Knowing this, can check to see if there are waaaaay too many F2s to make sense, and then flag those runs as likely needing to be re-run.(https://www.youtube.com/watch?v=lgGKSjiw0HQ)
  
  #Check to make sure the packages required are there and if not install them
  packages <- c("dplyr", "tidyr", "stringr","ggplot2","reshape2","grid","scales")
  if (length(setdiff(packages, rownames(installed.packages()))) > 0) { 
    install.packages(setdiff(packages, rownames(installed.packages())))
  } 
  
  #load each library
  require(dplyr)
  require(ggplot2)
  require(tidyr)
  require(stringr)
  require(reshape2)
  require(grid)
  require(scales)

  
  ## 
  tbCheck <- list.files(PreDir)

  possibleProbs_amChecking <- NULL
  for(i in tbCheck){
    
    LiamChecking <- list.files(paste0(PreDir, i))
    pzCheckingFind <- LiamChecking[grep("PofZ", LiamChecking)]
    amChecking <- read.table(paste0(PreDir, i, "/", pzCheckingFind), header = T)
    
    #possibleProbs_amChecking <- c(possibleProbs_amChecking, pzCheckingFind)
    
    if(sum(amChecking[,6])>sum(sum(amChecking[,3]),sum(amChecking[,4]))){
      possibleProbs_amChecking <- c(possibleProbs_amChecking, pzCheckingFind)
      }
    
   # if(sum(amChecking$X0.250.0.250.0.250.0.250)>sum(sum(amChecking$X1.000.0.000.0.000.0.000),sum(amChecking$X0.000.0.000.0.000.1.000))){
    #  possibleProps_amChecking <- c(possibleProps_amChecking, pzCheckingFind)
     # }
    
  }
  
  if(length(possibleProbs_amChecking) < 1){
    print("Looks good bud, giv'er")
  }
  if(length(possibleProbs_amChecking) >= 1){
    print(paste("Bud, you might want to have a closer look at", possibleProbs_amChecking, sep = " "))
  }
  
  
  
}