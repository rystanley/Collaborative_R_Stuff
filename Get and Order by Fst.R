### Editing the output of GenoDive so that it is workable in R

Genotest <- read.csv("~/Desktop/DFO Aquaculture Interaction/Batch_Merge/Genotest.csv")

levels(Genotest$Source.of.Variation)
which(Genotest$Source.of.Variation == "")
Genotest = Genotest[-which(Genotest$Source.of.Variation == ""),]


which(Genotest$Locus != "")
hold.locus.name = Genotest$Locus[which(Genotest$Locus != "")]

Salmon.SC.Fst = Genotest[which(Genotest$F.stat == "F_st"),]
Salmon.SC.Fst$Locus = hold.locus.name

Salmon.SC.Fst = droplevels(Salmon.SC.Fst)

summary(Salmon.SC.Fst$F.value)
hist(Salmon.SC.Fst$F.value)

SC_and_Aqua <- read.csv("~/Desktop/DFO Aquaculture Interaction/Batch_Merge/Salmon - South Coast + Aquaculture Pops.csv")

SC_and_Aqua = SC_and_Aqua[-which(SC_and_Aqua$Source.of.Variation == ""),]
sc.hold.locus = SC_and_Aqua$Locus[which(SC_and_Aqua$Locus != "")]
SC_and_Aqua = SC_and_Aqua[which(SC_and_Aqua$F.stat == "F_st"),]
SC_and_Aqua$Locus = sc.hold.locus

summary(SC_and_Aqua$F.value)
hist(SC_and_Aqua$F.value)

plot(as.factor(1:nrow(SC_and_Aqua)), SC_and_Aqua[order(-SC_and_Aqua$F.value),]$F.value)


View( SC_and_Aqua[order(-SC_and_Aqua$F.value),])
View( SC_and_Aqua[order(-SC_and_Aqua$F.value),])[1:100]


loci.to.get = SC_and_Aqua[order(-SC_and_Aqua$F.value),]$Locus[1:100]
loci.to.get = droplevels(loci.to.get)
write.csv(loci.to.get, "top100loci.csv")