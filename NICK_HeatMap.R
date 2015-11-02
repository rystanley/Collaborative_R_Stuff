## This code is to make a heatmap + cluster dendrogram of greencrab data.  

### Install package Heatplus
  library(gdata)
  library(gplots)
  library(ggplot2)
  library(reshape)
  library(grid)
  library(ggdendro)
  library(gridExtra)

### Input data ----------
    nData <- read.table("outlier_fst.txt",sep="\t",header=T) #read in data
    rownames(nData) <- nData[,1] #make the row names = to the first column
    fst <- nData[,-1] # remove the first column, now reduntant
    rm(nData) # delete the data which is read in
    
    upperTriangle(fst) <- lowerTriangle(fst) #create a 'square' data matrix instead of triangle dissimilarity
    fst <- data.matrix(fst) # convert to matrix format

## Data check ---------
    names(fst)==rownames(fst) # check ot make sure all the names are correct (should all be TRUE)

# Make a data frame from a matrix
    fst2 <- melt(fst)
    names(fst2) <- c("popA","popB","value")
    
    fst3=subset(fst2, subset=(fst2$value!=0.000)) #remove paired to one location fst values
    summary(fst3)


# Set grouping variables based on fst percentiles ------------

    #cut into values based on 20th percentiles
    fst2$SLICE  <-  "0"
    
    fst2[which(fst2$value>0 & fst2$value<quantile(fst2$value,0.20)),"SLICE"] <- 
      paste("[0.001-",round(quantile(fst2$value,0.2),3),"]",sep="")
    
    fst2[which(fst2$value>quantile(fst2$value,0.20) & fst2$value<quantile(fst2$value,0.40)),"SLICE"] <- 
      paste("[",round(quantile(fst2$value,0.20),3),"-",
                round(quantile(fst2$value,0.4),3),"]",sep="")
    
    fst2[which(fst2$value>quantile(fst2$value,0.40) & fst2$value<quantile(fst2$value,0.60)),"SLICE"] <- 
      paste("[",round(quantile(fst2$value,0.40),3),"-",
                round(quantile(fst2$value,0.6),3),"]",sep="")
    
    fst2[which(fst2$value>quantile(fst2$value,0.60) & fst2$value<quantile(fst2$value,0.80)),"SLICE"] <- 
      paste("[",round(quantile(fst2$value,0.60),3),"-",
                round(quantile(fst2$value,0.8),3),"]",sep="")

    fst2[which(fst2$value>quantile(fst2$value,0.80)),"SLICE"] <- 
      paste("[",round(quantile(fst2$value,0.80),3),"-",
                round(max(fst2$value),3),"]",sep="")
    
    #Set factor levels which ggplot uses to order the input
    fst2$SLICE  <-  factor(fst2$SLICE, levels=c(names(table(fst2$SLICE))[c(6,1:5)]))

## Create dendrogram -------------
    hc <- hclust(dist(fst))
    
    p1  <-  ggdendrogram(hc, theme_dendro=F, color= "black", size=6)+
      theme(axis.text.y=element_blank(), 
            axis.text.x= element_text(size=12,face="bold",colour="black"),
            axis.ticks=element_blank(),
            axis.title.y=element_blank(),
            axis.title.x = element_blank(),
            panel.background=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    
## Reorder populations as it appears on the dendrogram ----------
    fst2$popA=factor(fst2$popA,levels=(hc$labels[hc$order]),ordered=T)
    fst2$popB=factor(fst2$popB,levels=rev((hc$labels[hc$order])),ordered=T)
    
## Make the heatmap ---------------

  #Colour pallet (white for zero then smooth scale between yellow and dark red)
      FillCols=c("#FFFFFF",colorpanel(5,low="yellow",high="darkred"))
      ColLevs=c("",levels(fst2$SLICE)[2:6]) # do not put anything for zero (noted in figure caption)
  
      hm=qplot(fst2$popA, fst2$popB, fill=fst2$SLICE, data=fst2,geom="tile",colour="black")+
      scale_fill_manual(values = FillCols, labels=ColLevs)+
      theme(panel.border = element_rect(colour="black", fill=NA, size=1),
            panel.background = element_rect(fill = "white"),
            axis.title=element_text(size=12,colour="black",face="bold"),
            legend.title=element_text(size=12,colour="black",face="bold"),
            legend.text=element_text(colour="black",size=12,face="bold"),
            axis.text.x=element_text(size=12,colour="black",face="bold", angle=90),
            axis.title.y=element_text(size=12,colour="black",face="bold"),
            axis.text.y=element_text(size=12,colour="black",face="bold"),
            legend.position="bottom")+
      guides(fill=guide_legend(title=expression(paste(F[st], " value",sep=""))), keyheight=10,colour=FALSE,face="bold")+
      scale_x_discrete(breaks=NULL)+ 
        theme(axis.title.x = element_blank(), axis.title.y = element_blank());hm

### convert plots into grob objects -------------
    gp1<-ggplotGrob(p1)
    gp2<-ggplotGrob(hm)

## make grobs have same width dimensions (so they overlap properly) -----------
    maxWidth <-  grid::unit.pmax(gp1$widths[1:5], gp2$widths[1:5])
    gp1$widths[1:5] <- as.list(maxWidth)
    gp2$widths[1:5] <- as.list(maxWidth)

## Save plot ----------------
    #* Note you can change the specifics of this to a png, jpeg, bmp, tiff, pdf, etc. and play with settings
    
    png("GreenCrab_Dendro_HeatMap.png",width=600,height=480)
    
        grid.arrange(gp1, gp2, heights=c(1/5,4/5), ncol=1) 
    
    dev.off()
