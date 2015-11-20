##Mallory say relax

## Scallop analysis ###

# Libraries -----------------
library(ggplot2)
library(dplyr)
library(boa)
library(scales)
library(VennDiagram)
library(RCurl)
#library(gstudio)

# Functions -----------------------

#function to grab R code of Github using raw urls
SourceGitFunc <- function(url)
{
  
  require(RCurl)
  script <- getURL(url, ssl.verifypeer = FALSE)
  eval(parse(text = script),envir=.GlobalEnv)
}

#Source the function which will create the plot
SourceGitFunc("https://raw.githubusercontent.com/rystanley/RAD_R_Functions/master/AlleleFreq.R")


#Funciton for recoding entries in a vector
recoderFunc <- function(data, oldvalue, newvalue) {
  if (is.factor(data))     data     <- as.character(data)
  if (is.factor(oldvalue)) oldvalue <- as.character(oldvalue)
  if (is.factor(newvalue)) newvalue <- as.character(newvalue)
  newvec <- data
  for (i in unique(oldvalue)) newvec[data == i] <- newvalue[oldvalue == i]
  newvec
}
# Load data ------------------

GenePopData <- read.table("GenePopFixed.txt",
                          header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)


Outliers <- read.table("AllPopsOutliers.txt",
                       header = TRUE, sep = "\t")

Outliers <- as.character(Outliers$Loci)

Neutral <- read.table("AllPopsNeutral.txt",
                      header=TRUE, sep="\t")

Neutral  <- as.character(Neutral$Loci)


## Read in the population name updates ---------------

NameUpdates <- read.csv("PopNameFix.csv",sep=",")
colnames(NameUpdates) <- c("old","new")

# Get the order of populations --------------
PopOrder <- rev(c("bon","ltb","mgd","nts","psb","bof","ssm","gmi","ssb","gmo","geo","mda"))


### Create the allele frequency plot --------------

    #Get the data instead of plot
    Scallop_Heat_Outlier=AlleleFreqHeatMap(GenePopData,subs = Outliers,keep = TRUE,POP="CHAR",refPop="bon",
                                           OrderPops=PopOrder,optimizer = TRUE,standardize = TRUE,plot=FALSE)
    
    # Fix the population names assigned based on the sample IDs
    for (i in 1:nrow(NameUpdates))
      {
      Scallop_Heat_Outlier[,"Pop"] <- recoderFunc(Scallop_Heat_Outlier[,"Pop"],NameUpdates[i,"old"],NameUpdates[i,"new"])
    }
    
    #Create the plot
    p1 <- ggplot(Scallop_Heat_Outlier,aes(x=SNP,y=Pop,fill=FreqStand))+
      geom_tile()+ scale_y_discrete(expand = c(0,0))+scale_x_discrete(expand = c(0,0))+
      theme_bw()+theme(legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1))+
      scale_fill_gradient(low="blue",high=muted("red"))+
      labs(y="Population",x="SNP",fill="Standardized allele frequency")
    
    Scallop_Heat_Outlier <- p1
    Scallop_Heat_Outlier# show the plot output

    ggsave("Scallop_Outlier_AlleleFreqHeatMap.png",Scallop_Heat_Outlier,
           width=12,height=8)

## Create the allele frequency plot
Scallop_Heat_Neutral=AlleleFreqHeatMap(GenePopData,subs = Neutral,keep = TRUE,POP="CHAR",refPop="bon",
                                       OrderPops=PopOrder,optimizer = TRUE,standardize = TRUE,plot=FALSE)

# Fix the population names assigned based on the sample IDs
for (i in 1:nrow(NameUpdates))
{
  Scallop_Heat_Neutral[,"Pop"] <- recoderFunc(Scallop_Heat_Neutral[,"Pop"],NameUpdates[i,"old"],NameUpdates[i,"new"])
}

#Create the plot
p2 <- ggplot(Scallop_Heat_Neutral,aes(x=SNP,y=Pop,fill=FreqStand))+
  geom_tile()+ scale_y_discrete(expand = c(0,0))+scale_x_discrete(expand = c(0,0))+
  theme_bw()+theme(legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_gradient(low="blue",high=muted("red"))+
  labs(y="Population",x="SNP",fill="Standardized allele frequency")

Scallop_Heat_Neutral <- p2
Scallop_Heat_Neutral# show the plot output

ggsave("Scallop_Neutral_AlleleFreqHeatMap.png",Scallop_Heat_Neutral,
       width=12,height=8)



### Use the clinal Neutral markers -------------
Cline_Neutral <- read.csv("Cline_Neutral_Loci.csv")
Cline_Neutral <- as.character(Cline_Neutral$Neutral)
Cline_Neutral <- gsub("[a-zA-Z]+","",Cline_Neutral) # remove the allele letters

Scallop_Heat_NeutralCline=AlleleFreqHeatMap(GenePopData,subs = Cline_Neutral,keep = TRUE,POP="CHAR",refPop="bon",
                                            OrderPops=PopOrder,optimizer = TRUE,standardize = TRUE,plot=FALSE)

# Fix the population names assigned based on the sample IDs
for (i in 1:nrow(NameUpdates))
{
  Scallop_Heat_NeutralCline[,"Pop"] <- recoderFunc(Scallop_Heat_NeutralCline[,"Pop"],NameUpdates[i,"old"],NameUpdates[i,"new"])
}

#Create the plot
p3 <- ggplot(Scallop_Heat_NeutralCline,aes(x=SNP,y=Pop,fill=FreqStand))+
  geom_tile()+ scale_y_discrete(expand = c(0,0))+scale_x_discrete(expand = c(0,0))+
  theme_bw()+theme(legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_gradient(low="blue",high=muted("red"))+
  labs(y="Population",x="SNP",fill="Standardized allele frequency")

Scallop_Heat_NeutralCline <- p3
Scallop_Heat_NeutralCline# show the plot output

ggsave("Scallop_Neutral_Cline_AlleleFreqHeatMap.png",
       Scallop_Heat_NeutralCline,width=12,height=8)

### Use the clinal outlier markers -------------
Cline_Outlier <- read.csv("Cline_Outlier_Loci.csv")
Cline_Outlier <- as.character(Cline_Outlier$Outlier)
Cline_Outlier <- gsub("[a-zA-Z]+","",Cline_Outlier) # remove the allele letters

Scallop_Heat_OutlierCline=AlleleFreqHeatMap(GenePopData,subs = Cline_Outlier,keep = TRUE,POP="CHAR",refPop="bon",
                                            OrderPops=PopOrder,optimizer = TRUE,standardize = TRUE,plot=FALSE)

# Fix the population names assigned based on the sample IDs
for (i in 1:nrow(NameUpdates))
{
  Scallop_Heat_OutlierCline[,"Pop"] <- recoderFunc(Scallop_Heat_OutlierCline[,"Pop"],NameUpdates[i,"old"],NameUpdates[i,"new"])
}

#Create the plot
p4 <- ggplot(Scallop_Heat_OutlierCline,aes(x=SNP,y=Pop,fill=FreqStand))+
  geom_tile()+ scale_y_discrete(expand = c(0,0))+scale_x_discrete(expand = c(0,0))+
  theme_bw()+theme(legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_gradient(low="blue",high=muted("red"))+
  labs(y="Population",x="SNP",fill="Standardized allele frequency")

Scallop_Heat_OutlierCline <- p4
Scallop_Heat_OutlierCline# show the plot output


Scallop_Heat_OutlierCline # show the plot output

ggsave("Scallop_Outlier_Cline_AlleleFreqHeatMap.png",
       Scallop_Heat_OutlierCline,width=12,height=8)
