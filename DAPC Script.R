#Discriminant Analysis of Principal Components (DAPC)
library(adegenet)

dat<-read.structure("FullStructure.str") #Interactive Step

grp1<-find.clusters(dat) #Interactive Step, choose all PCs and best K value (lowest BIC)

#show the number of individuals assigned to each population
table(pop(dat),grp1$grp)
grp1 #shows some useful stats

#DAPC TIME
dapc1<-dapc(dat, grp1$grp) #the grp1$grp gives the DAPC the prior of your clusters found in the previous step
#interactive steps

#summary stats
dapc1 
summary(dapc1)
#show the assignment of each invidiual to a group
dapc1$posterior

#Plot your DAPC, modify as you like
myCol <- c("red","blue")

png(filename = "DAPC_Scatterplot2.png", 
    width = 2400, height = 2000, res = 300, bg="transparent")
scatter(dapc1, posi.da="bottomright", bg="white",pch=17:22, cstar=0, col=myCol, scree.pca=FALSE,legend=TRUE,txt.leg=c("North","South"),solid=.4)
dev.off()
#OR
#I like this one
scatter(dapc1, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, col=myCol, solid=.4,
        cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:2))
#Now plot which invidiuals are in which cluster, with a subset

assignplot(dapc1, subset=1:50)

#Finally, plot a STRUCTURE like plot, with individuals on the x-axis and membership prob. on the Y
compoplot(dapc1, posi="bottomright",
         txt.leg=paste("Cluster", 1:5), lab="",
         ncol=1, xlab="individuals", col=funky(6))

###Save final workspace
save.image("C:/Users/Nick/Desktop/Postdoc/Green Crab/RAD/DAPC/All/DAPC_Full.RData")
