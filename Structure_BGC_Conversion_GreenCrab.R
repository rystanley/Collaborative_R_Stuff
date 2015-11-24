### code for data prep - Structure <- bgc format ---------

    #Following the example code by:
    #https://popgencode.wordpress.com/2015/04/17/preparing-data-file-bgc-genotype-count/

##load libraries -------------
        library(adegenet)
        library(dplyr)

#set working directory (makes for cleaner code) ------------------
# e.g. (setwd("c:/Users/RyanStanley/OneDrive/PostDoc/DFO/Nick/BGC/"))
        setwd("~BGC/")

#Load in geneteic data
        CrabNorth <- import2genind('CrabOutliers_PureNorth.str', onerowperind=F, n.ind=132, n.loc=117, row.marknames=1, col.lab=1, col.pop=2, ask=F)
        CrabSouth <- import2genind('CrabOutliers_PureSouth.str', onerowperind=F, n.ind=66, n.loc=117, row.marknames=1, col.lab=1, col.pop=2, ask=F)

#convert to adegenet format
        CrabNorth <- genind2genpop(CrabNorth)
        CrabSouth <- genind2genpop(CrabSouth)

#Grab the allele frequency format data * the code on the web is missing something
        CN_Data=t(CrabNorth@tab)
        CrabNorth_allele <- data.frame(Allele=rownames(CN_Data))
        CrabNorth_allele <- cbind(CrabNorth_allele,CN_Data)
        rownames(CrabNorth_allele) <- NULL
        
        CS_Data=t(CrabSouth@tab)
        CrabSouth_allele <- data.frame(Allele=rownames(CS_Data))
        CrabSouth_allele <- cbind(CrabSouth_allele,CS_Data)
        rownames(CrabSouth_allele) <- NULL

#remove unneeded data from this step
        rm(CN_Data,CS_Data)

#Create Population dataframes 

      #Function to wrangle data in BGC format according to table 1 in BGC manual
          BGCFormat=function(x,Pop)
          {
            #extract the Loci names
            LociNames <- lapply(strsplit(as.character(x$Allele), '.', fixed=TRUE), '[[', 1)
            x$LociNames <- unlist(LociNames)
    
            temp <- x[,c("Allele","LociNames",Pop)] # subset out just the data to a specific population
            
            tempData <- data.frame(Loci=temp[seq(1,nrow(temp),2),"LociNames"],
                                   Allele=paste(temp[seq(1,nrow(temp),2),Pop],
                                                             temp[seq(2,nrow(temp),2),Pop],sep=" "))
            tempData$Allele=as.character(tempData$Allele) # convert the joined populatoin based allele count for a given locus to character
            
            #format for BGC output accordig to example table 1 in the manual
            Output <- NULL
            for (i in 1:nrow(tempData))
            {
            Output <- c(Output,paste("locus ",tempData[i,"Loci"],sep=""),tempData[i,"Allele"])
            }
            
            return(Output)
          }
        
        
  #Create population files for the north and south groupings 
      
        for (i in colnames(CrabNorth_allele[2:length(CrabNorth_allele)]))
        {
          BGCData <- BGCFormat(CrabNorth_allele,i)
          write.table(BGCData,paste("BGC Files/North_pop",i,"_bgc.txt",sep=""), 
                      row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
        }
        
        for (i in colnames(CrabSouth_allele[2:length(CrabSouth_allele)]))
        {
          BGCData <- BGCFormat(CrabSouth_allele,i)
          write.table(BGCData,paste("BGC Files/South_pop",i,"_bgc.txt",sep=""), 
                      row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
        }
     
        rm(BGCData) # clean workspace

#Create the admixed file format --------------
        
        # function to assign a binary 0 1 or 2 scoring where 2 is the larger allele (2 or 4) and 1 is the smaller allele (1 or 3)
        aCount <- function(x){
          if(length(unique(as.numeric(x[c("allele1","allele2")])))==1 &
             sum(as.numeric(x[c("allele1","allele2")])!= -9)>1 &
             unique(as.numeric(x["allele1"])==as.numeric(x["alleleMin"]))){return(c(1,1))} # same alleles (small)
          
          if(length(unique(as.numeric(x[c("allele1","allele2")])))==1 &
             sum(as.numeric(x[c("allele1","allele2")])!= -9)>1 &
             unique(as.numeric(x["allele1"])==as.numeric(x["alleleMax"]))){return(c(2,2))} # same alleles (small)
          
          if(as.numeric(x[,"allele1"])>as.numeric(x[,"allele2"])){return(c(2,0))}
          if(as.numeric(x[,"allele1"])<as.numeric(x[,"allele2"])){return(c(0,2))}
          if(unique(as.numeric(x[c("allele1","allele2")]))== -9){return(c(0,0))}
        }
      
      #Load mixed data
        CrabMixed <- import2genind('Outliers_AdmixedIn.str', onerowperind=F, n.ind=43, n.loc=117, row.marknames=1, col.lab=1, col.pop=2, ask=F)
        CrabMixed <- genind2genpop(CrabMixed) ### I think you have a monomophic allele or something. Not an even number of alleles. 
        head(CrabMixed) # should be ab 1 x 230 matrix 
        
      #format to the tabular output
        CM_Data=t(CrabMixed@tab)
        CrabMixed_allele <- data.frame(Allele=rownames(CM_Data))
        CrabMixed_allele <- cbind(CrabMixed_allele,CM_Data)
        rownames(CrabMixed_allele) <- NULL ### ** issue here. There seems to be a missing allele for one of the loci
        
      #Grab the names of loci in the admixed data
        LocusMixed1 <- as.character(CrabMixed_allele$Allele)
        LocusMixed1 <- lapply(strsplit(LocusMixed1, '.', fixed=TRUE), '[[', 1)
        
      #Flag loci which have only one allele
        table(unlist(LocusMixed1))[which(as.vector(table(unlist(LocusMixed1)))==1)] # here we see the issue with locus 3657 5488 6798 7068  946 all only have one allele
        names(table(unlist(LocusMixed1))[which(as.vector(table(unlist(LocusMixed1)))==1)]) #names of the monomorphic alleles
        
      #names of the loci with only one allele (need to be removed)
        MonoAlleles <- names(table(unlist(LocusMixed1))[which(as.vector(table(unlist(LocusMixed1)))==1)])
        LocusMixed2 <- unlist(LocusMixed1)
        
        LocusMixed3 <- LocusMixed2[!(LocusMixed2 %in% MonoAlleles)]
        table(unlist(LocusMixed3))[which(as.vector(table(unlist(LocusMixed3)))==1)] #check to see if it worked (Should be nothing)
        
        #remove repeats
        LocusMixed <- LocusMixed3[rep(c(TRUE,FALSE),(nrow(CrabMixed_allele)-length(MonoAlleles))/2)]
        
        # fixed up structure file where I added a header for the crab ID and the population
        MixedStruct1 <- read.table("Outliers_AdmixedIn.txt",header=T)
        
        Loci <- gsub("X","",as.character(colnames(MixedStruct1)[3:length(MixedStruct1)])) # remove the "X" which R adds to the numeric column names
        colnames(MixedStruct1)[3:length(MixedStruct1)] <- Loci #adjust column names
        
        MixedStruct=MixedStruct1[,c(colnames(MixedStruct1)[1:2],LocusMixed)] #subset out the loci which only have one allele
        

     #Data wrangle and construct the format for admixture populatons as per the instructions and table 2 in the BGC manual
        MixedData <- NULL
        for(col in names(MixedStruct)[3:length(MixedStruct)])
          {
          temp1 <- MixedStruct[,c("ID","Pop",col)]
        
          Locushold <- paste("locus ",col,sep="") # Start the locus data
            
            for (i in unique(MixedStruct$Pop)) # each populatoin
              {
              temp2 <- filter(temp1,Pop==i) # subset for the population
              
              #Reformat the data for one row for each individaul (ID, Pop, Allele1, Allele2)
              temp3 <- data.frame(ID=temp2[seq(1,nrow(temp2),2),"ID"],
                                  Pop=temp2[seq(1,nrow(temp2),2),"Pop"],
                                  allele1=temp2[seq(1,nrow(temp2),2),col],
                                  allele2=temp2[seq(2,nrow(temp2),2),col],
                                  alleleMax=max(temp2[,col],na.rm=T),
                                  alleleMin=min(temp2[which(temp2[,col]>(-1)),col],na.rm=T))
              
              temp4 <- as.data.frame(temp3%>%group_by(ID)%>%do(col1=aCount(.)[1],col2=aCount(.)[2])%>%ungroup())
              
              temp5 <- paste(temp4[,"col1"],temp4[,"col2"],sep=" ")
              
              Locushold <- c(Locushold,paste("pop ",i,sep=""),temp5)
            } #end of population loop
          
            MixedData <- c(MixedData,Locushold) # add each successive locus 
          } #end of locus loop
  
      #Save output for BGC
      write.table(MixedData,"BGC Files/MixedData_Crab.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
      
      
      
      
## To create the files as per the tutorial online ----------------
#       
#       #Get locus names without the allele appended value
#       LocusNorth <- as.character(CrabNorth_allele$Allele)
#       LocusNorth <- lapply(strsplit(LocusNorth, '.', fixed=TRUE), '[[', 1)
#       LocusNorth <- LocusNorth[rep(c(TRUE,FALSE),nrow(CrabNorth_allele)/2)]
#       LocusNorth <- unlist(LocusNorth)
#       
#       LocusSouth <- as.character(CrabSouth_allele$Allele)
#       LocusSouth <- lapply(strsplit(LocusSouth, '.', fixed=TRUE), '[[', 1)
#       LocusSouth <- LocusSouth[rep(c(TRUE,FALSE),nrow(CrabSouth_allele)/2)]
#       LocusSouth <- unlist(LocusSouth)
#       
#       #Get the data for each allele in each population (every other row as denoted by even and odd row integers)
#       North_Odd <- CrabNorth_allele[seq(1,length(CrabNorth_allele$Allele),2),]
#       North_Even <- CrabNorth_allele[seq(2,length(CrabNorth_allele$Allele),2),]
#       
#       South_Odd <- CrabSouth_allele[seq(1,length(CrabSouth_allele$Allele),2),]
#       South_Even <- CrabSouth_allele[seq(2,length(CrabSouth_allele$Allele),2),]
      
#       NorthData <- data.frame(Loci=LocusNorth,filler=rep('#',length(LocusNorth)),PopOdd=North_Odd$`1`,PopEven=North_Even$`1`)
#       
#       #for each pop
#       head(CrabNorth_allele)
#       write.table(NorthData,"test.txt", row.names=F, col.names=F, quote=F, sep='\t')
#       
#       #North Populations
#       for(i in colnames(CrabNorth_allele)[2:length(CrabNorth_allele)])
#       {
#         write.table(data.frame(Loci=LocusNorth,filler=rep('#',length(LocusNorth)),PopOdd=North_Odd[,i],PopEven=North_Even[,i]),
#                     paste("BGC Files/North_pop",i,"_allct_bgc.txt",sep=""), row.names=F, col.names=F, quote=F, sep='\t')
#       }
#       
#       #South Populations
#       for(i in colnames(CrabSouth_allele)[2:length(CrabSouth_allele)])
#       {
#         write.table(data.frame(Loci=LocusSouth,filler=rep('#',length(LocusSouth)),PopOdd=South_Odd[,i],PopEven=South_Even[,i]),
#                     paste("BGC Files/South_pop",i,"_allct_bgc.txt",sep=""), row.names=F, col.names=F, quote=F, sep='\t')
#       }
