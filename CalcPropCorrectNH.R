calc.prop.correctFN <- function(NHResults){

num.sim <- nrow(NHResults)/(length(NHResults)-1)


if(sum(NHResults$Pure1[1:200]) > sum(NHResults$Pure2[1:200])){
  prop.corr.farm <- length(which(NHResults$Pure1[1:200] > 0.5))/num.sim
  prop.corr.wild <- length(which(NHResults$Pure2[201:400] > 0.5))/num.sim
  prop.corr.F1 <- length(which(NHResults$F1[401:600] > 0.5))/num.sim
  prop.corr.F2 <- length(which(NHResults$F2[601:800] > 0.5))/num.sim
  prop.corr.farmBC <- length(which(NHResults$BC1[801:1000] > 0.5))/num.sim
  prop.corr.wildBC <- length(which(NHResults$BC2[1001:1200] > 0.5))/num.sim
}


if(sum(NHResults$Pure1[1:200]) < sum(NHResults$Pure2[1:200])){
  prop.corr.farm <- length(which(NHResults$Pure1[201:400] > 0.5))/num.sim
  prop.corr.wild <- length(which(NHResults$Pure2[1:200] > 0.5))/num.sim
  prop.corr.F1 <- length(which(NHResults$F1[401:600] > 0.5))/num.sim
  prop.corr.F2 <- length(which(NHResults$F2[601:800] > 0.5))/num.sim
  prop.corr.farmBC <- length(which(NHResults$BC1[1001:1200] > 0.5))/num.sim
  prop.corr.wildBC <- length(which(NHResults$BC2[801:1000] > 0.5))/num.sim
}


print(paste("Farm Prop Correct", prop.corr.farm))
print(paste("Wild Prop Correct", prop.corr.wild))
print(paste("F1 Prop Correct", prop.corr.F1))
print(paste("F2 Prop Correct", prop.corr.F2))
print(paste("Farm BC Prop Correct", prop.corr.farmBC))
print(paste("Wild BC Prop Correct", prop.corr.wildBC))


}