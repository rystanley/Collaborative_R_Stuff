calc.prop.correctFN <- function(NHResults, print.results = TRUE, all.hyb = FALSE){

num.sim <- nrow(NHResults)/(length(NHResults)-1)

## there is the potential for NH to assign different groups to Pop1 and Pop2 between runs since functionally there is no differrence
  ## which population it designates one and 2 
      ## the only issue is that to calculate the mean proportion that have been assigned correctly, need to have the individuals given
      ## be in the same population - also have to then correc the BCs

## if the first population given to NH has been assinged to be Pop1 by NH, then the probabilities that the individuals of that population
  ## are in Pop1 will be greater than the probability they are in fact in Pop2 
    ## so if this is true, the first Population given is Pop1, so, calculate based on that - also have to adjust so that BC1 is associatd
    ### with Pop1
if(sum(NHResults$Pure1[1:200]) > sum(NHResults$Pure2[1:200])){
  
   ## calculate the number of indiviudals in the known populations whose probability of being in the correct population exceeds
        ## given by NH exceeds the given level of stringency / the total number of individuals in that population
  
  ## proportion correct at p = 0.5
  prop.corr.farm.p.5 <- length(which(NHResults$Pure1[1:200] > 0.5))/num.sim
  prop.corr.wild.p.5 <- length(which(NHResults$Pure2[201:400] > 0.5))/num.sim
  prop.corr.F1.p.5 <- length(which(NHResults$F1[401:600] > 0.5))/num.sim
  prop.corr.F2.p.5 <- length(which(NHResults$F2[601:800] > 0.5))/num.sim
  prop.corr.farmBC.p.5 <- length(which(NHResults$BC1[801:1000] > 0.5))/num.sim
  prop.corr.wildBC.p.5 <- length(which(NHResults$BC2[1001:1200] > 0.5))/num.sim
  
  ## proportion correct at p=0.75
  
  prop.corr.farm.p.75 <- length(which(NHResults$Pure1[1:200] > 0.75))/num.sim
  prop.corr.wild.p.75 <- length(which(NHResults$Pure2[201:400] > 0.75))/num.sim
  prop.corr.F1.p.75 <- length(which(NHResults$F1[401:600] > 0.75))/num.sim
  prop.corr.F2.p.75 <- length(which(NHResults$F2[601:800] > 0.75))/num.sim
  prop.corr.farmBC.p.75 <- length(which(NHResults$BC1[801:1000] > 0.75))/num.sim
  prop.corr.wildBC.p.75 <- length(which(NHResults$BC2[1001:1200] > 0.75))/num.sim
  
  ## proportion correct at p=0.9
  
  prop.corr.farm.p.9 <- length(which(NHResults$Pure1[1:200] > 0.9))/num.sim
  prop.corr.wild.p.9 <- length(which(NHResults$Pure2[201:400] > 0.9))/num.sim
  prop.corr.F1.p.9 <- length(which(NHResults$F1[401:600] > 0.9))/num.sim
  prop.corr.F2.p.9 <- length(which(NHResults$F2[601:800] > 0.9))/num.sim
  prop.corr.farmBC.p.9 <- length(which(NHResults$BC1[801:1000] > 0.9))/num.sim
  prop.corr.wildBC.p.9 <- length(which(NHResults$BC2[1001:1200] > 0.9))/num.sim
  
  
  ## should the proportion of all indiviuals known to be hyybrids which have been identified as hybrids be calculated at each
    ## level of stringency
  if(all.hyb = TRUE){
    hyb.det.ALL.p.5 <- length(which(NHResults$pHyb[401:1200] > 0.5))/(num.sim*4)
    hyb.det.ALL.p.75 <- length(which(NHResults$pHyb[401:1200] > 0.75))/(num.sim*4)
    hyb.det.ALL.p.9 <- length(which(NHResults$pHyb[401:1200] > 0.9))/(num.sim*4)
  }
  
}

### if the First population given to NH has been denoted Pop2 by NH
if(sum(NHResults$Pure1[1:200]) < sum(NHResults$Pure2[1:200])){
 
  ## proportion correct at p=0.5
  
   prop.corr.farm.p.5 <- length(which(NHResults$Pure1[201:400] > 0.5))/num.sim
  prop.corr.wild.p.5 <- length(which(NHResults$Pure2[1:200] > 0.5))/num.sim
  prop.corr.F1.p.5 <- length(which(NHResults$F1[401:600] > 0.5))/num.sim
  prop.corr.F2.p.5 <- length(which(NHResults$F2[601:800] > 0.5))/num.sim
  prop.corr.farmBC.p.5 <- length(which(NHResults$BC1[1001:1200] > 0.5))/num.sim
  prop.corr.wildBC.p.5 <- length(which(NHResults$BC2[801:1000] > 0.5))/num.sim
  
  ## proportion correct at p=0.75
  
  prop.corr.farm.p.75 <- length(which(NHResults$Pure1[201:400] > 0.75))/num.sim
  prop.corr.wild.p.75 <- length(which(NHResults$Pure2[1:200] > 0.75))/num.sim
  prop.corr.F1.p.75 <- length(which(NHResults$F1[401:600] > 0.75))/num.sim
  prop.corr.F2.p.75 <- length(which(NHResults$F2[601:800] > 0.75))/num.sim
  prop.corr.farmBC.p.75 <- length(which(NHResults$BC1[1001:1200] > 0.75))/num.sim
  prop.corr.wildBC.p.75 <- length(which(NHResults$BC2[801:1000] > 0.75))/num.sim
  
   ## proportion correct at p=0.9
  
  prop.corr.farm.p.9 <- length(which(NHResults$Pure1[201:400] > 0.9))/num.sim
  prop.corr.wild.p.9 <- length(which(NHResults$Pure2[1:200] > 0.9))/num.sim
  prop.corr.F1.p.9 <- length(which(NHResults$F1[401:600] > 0.9))/num.sim
  prop.corr.F2.p.9 <- length(which(NHResults$F2[601:800] > 0.9))/num.sim
  prop.corr.farmBC.p.9 <- length(which(NHResults$BC1[1001:1200] > 0.9))/num.sim
  prop.corr.wildBC.p.9 <- length(which(NHResults$BC2[801:1000] > 0.9))/num.sim
  
  if(all.hyb = TRUE){
    hyb.det.ALL.p.5 <- length(which(NHResults$pHyb[401:1200] > 0.5))/(num.sim*4)
    hyb.det.ALL.p.75 <- length(which(NHResults$pHyb[401:1200] > 0.75))/(num.sim*4)
    hyb.det.ALL.p.9 <- length(which(NHResults$pHyb[401:1200] > 0.9))/(num.sim*4)
  }
  
}

if(print.results = TRUE){
print(paste("Farm Prop Correct 0.5", prop.corr.farm.p.5))
print(paste("Wild Prop Correct 0.5", prop.corr.wild.p.5))
print(paste("F1 Prop Correct 0.5", prop.corr.F1.p.5))
print(paste("F2 Prop Correct 0.5", prop.corr.F2.p.5))
print(paste("Farm BC Prop Correct 0.5", prop.corr.farmBC.p.5))
print(paste("Wild BC Prop Correct 0.5", prop.corr.wildBC.p.5))
if(all.hyb = TRUE){
  print(paste("Hybrid Detected Regardless of Class 0.5", hyb.det.ALL.p.5))
}
print("")
print(paste("Farm Prop Correct 0.75", prop.corr.farm.p.75))
print(paste("Wild Prop Correct 0.75", prop.corr.wild.p.75))
print(paste("F1 Prop Correct 0.75", prop.corr.F1.p.75))
print(paste("F2 Prop Correct 0.75", prop.corr.F2.p.75))
print(paste("Farm BC Prop Correct 0.75", prop.corr.farmBC.p.75))
print(paste("Wild BC Prop Correct 0.75", prop.corr.wildBC.p.75))
if(all.hyb = TRUE){
  print(paste("Hybrid Detected Regardless of Class 0.75", hyb.det.ALL.p.75))
}
print("")
print(paste("Farm Prop Correct 0.9", prop.corr.farm.p.9))
print(paste("Wild Prop Correct 0.9", prop.corr.wild.p.9))
print(paste("F1 Prop Correct 0.9", prop.corr.F1.p.9))
print(paste("F2 Prop Correct 0.9", prop.corr.F2.p.9))
print(paste("Farm BC Prop Correct 0.9", prop.corr.farmBC.p.9))
print(paste("Wild BC Prop Correct 0.9", prop.corr.wildBC.p.9))
if(all.hyb = TRUE){
  print(paste("Hybrid Detected Regardless of Class 0.9", hyb.det.ALL.p.9))
}

} ## end of whether or not to print the results if test


## if the proportion of all hybs assigned correctly has been calculated, output

if(all.hyb = TRUE){
p.5 <- rbind(prop.corr.farm.p.5, prop.corr.wild.p.5, prop.corr.F1.p.5, prop.corr.F2.p.5, prop.corr.farmBC.p.5,
  prop.corr.wildBC.p.5,  hyb.det.ALL.p.5)
p.75 <- rbind(prop.corr.farm.p.75, prop.corr.wild.p.75, prop.corr.F1.p.75, prop.corr.F2.p.75, prop.corr.farmBC.p.75,
  prop.corr.wildBC.p.75,  hyb.det.ALL.p.75)
p.9 <- rbind(prop.corr.farm.p.9, prop.corr.wild.p.9, prop.corr.F1.p.9, prop.corr.F2.p.9, prop.corr.farmBC.p.9,
  prop.corr.wildBC.p.9,  hyb.det.ALL.p.9)

corr.names <- c("Farm","Wild", "F1", "F2", "FarmBC", "WildBC", "All.Hyb")
prop.corr.output <- cbind(corr.names, p.5, p.75, p.9)
colnames(prop.corr.output)[2:4] <- c("p0.5", "p0.75", "p0.9")
}

## if the proportion of all individuals known to be hybs have not been calculated -- output
if(all.hyb = FALSE){
p.5 <- rbind(prop.corr.farm.p.5, prop.corr.wild.p.5, prop.corr.F1.p.5, prop.corr.F2.p.5, prop.corr.farmBC.p.5,
  prop.corr.wildBC.p.5)
p.75 <- rbind(prop.corr.farm.p.75, prop.corr.wild.p.75, prop.corr.F1.p.75, prop.corr.F2.p.75, prop.corr.farmBC.p.75,
  prop.corr.wildBC.p.75)
p.9 <- rbind(prop.corr.farm.p.9, prop.corr.wild.p.9, prop.corr.F1.p.9, prop.corr.F2.p.9, prop.corr.farmBC.p.9,
  prop.corr.wildBC.p.9)

corr.names <- c("Farm","Wild", "F1", "F2", "FarmBC", "WildBC")
prop.corr.output <- cbind(corr.names, p.5, p.75, p.9)
colnames(prop.corr.output)[2:4] <- c("p0.5", "p0.75", "p0.9")
}


name.assign <- "Res.prop.corr"
print(name.assign)
assign(x = name.assign, value = prop.corr.output, envir = globalenv())



}