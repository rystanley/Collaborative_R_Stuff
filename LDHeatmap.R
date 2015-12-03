
#load packages
library("mixOmics")
library("stats")
library("gplots")
library("pheatmap")
library("RColorBrewer")


#import loci
#not necessary for heatmap but retained from Mallory's code
loci_list <- read.csv("C:/Users/Mallory/Documents/School/MASTER_Masters/Genetic_Data/June_2014_finals/VanWyngaarden_finaldata_HWE_FinalData_LocusList.csv")

outlier_loci <- as.vector(loci_list[,1])



#plot color pallettes
my_palette <- colorRampPalette(c("royalblue4", "palegoldenrod", "brown"))(n = 11)


my_palette_2 <- colorRampPalette(rev(brewer.pal(99, "Spectral")), space="Lab")

#outlier
#setwd
setwd("C:/Users/Nick/Desktop/Postdoc/Green Crab/Plink Files/")


#load data
outlier_LD_r2_matrix <- read.csv("Outlier_LD.csv", row.names=1)
loci <- row.names(outlier_LD_r2_matrix)
colnames(outlier_LD_r2_matrix) <- loci

outlier_LD_r2_matrix <- as.matrix(outlier_LD_r2_matrix)

#colors default is "heat.colors"
png(filename="LDHeatmap3.png",height=2400,width=2400,res=400,bg="white")
heatmap.2(outlier_LD_r2_matrix,
          Rowv=FALSE, #rows should be reordered as required
          Colv = "Rowv", #columns should be treated as rows
          dendrogram="none", #no trees
          scale="none",
          breaks=100, #number of break points used to bin into colours
          col=my_palette_2,
          trace="none", #whether lines should be drawn between cols, rows,
          margins=c(1,1),#margins for column names and row names
          labRow= " ",
          labCol= " ",
          #cexCol=0.4, #column label size
          #cexRow=0.4,#row label size
          #srtCol = 90,#col label angle, degrees from horizontal
          #srtRow = 90,
          key=TRUE,
          keysize = 1,
          key.title="",
          density.info="none",
          lmat = rbind(c(3,1),c(4,2)), #order of elements (1 is heatmap, 4 is key) #makes the plot a 2x2 grid, each "cell" has 1 element - row tree, key, col tree, and heatmap
          lhei = c(20,5), #row height for plot elements
          lwid = c(8,30)  #column width for plot elements
)

dev.off()

#####final save####
save.image("C:/Users/Nick/Desktop/Postdoc/Green Crab/Plink Files/LDHeatmap.RData")




