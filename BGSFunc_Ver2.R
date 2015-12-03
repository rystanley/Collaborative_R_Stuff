BGCFunc <- function(p1,p2,admixed,dir="",name=""){

  #Function details -------------
  
  #This function will convert a structure files from two parental locations (p1 and p2) and an
  #admixed population into BGC format for Bayesian Genomic Cline analysis
  #this function will by default remove any monomorphic alleles as per the requirements of BGC
  #A list of monomorphic alleles detected will be created for each population (p1,p2 and admixed)
  
  #p1 p2 & admixed <- paths to conventional structure files (first column is the sample ID and second is the population) 
  #                   Note that only one population should be specified in the second column of the structure files for parental data
  
  # dir <- is a directory where the output txt files will be saved. Note that this should be labeled specifically
  #       to the BGC analysis which is being undertaken. This must be specified
  
  # name <- is a flag to be added to the file names to make them more specific. Default is NULL so if not specified
  #         the output txt files saved are genetic (BGC_p1.txt, BGC_p2.txt, BGC_admixture.txt and BGC_monomorphic.txt)
  
  # Example usage
  #setwd("Documents/Files/BGC_RAW")
  #BGCFunc(p1 = 'Outliers_NorthPops_.str',p2 = 'Outliers_SouthPops_.str',
          #admixed = 'Outliers_AdmixedIn.str',dir="Documents/Files/BGC_RAW/processed",name="NorthSouthAnalysis")
  
  #returns the files: 'NorthSouthAnalysis_Admixed_BGC.txt' 'NorthSouthAnalysis_MonomorphicAlleleList.txt'
  #                   'NorthSouthAnalysis_Parental1_BGC' 'NorthSouthAnalysis_Parental2_BGC'
  
  # in the "Documents/Files/BGC_RAW/processed" folder
  
  ## Custom data wrangling functions needed in this analysis ---------------
  
  #this function will read in a structure file as a dataframe.
  StructureRead <- function(dir){
    
    # This two step process helps with the issues around there being missing
    #column names for the first two files of structure. 
    snpData <- read.table(dir,skip=1) #skip the first column because it doesn't align
    SNPnames <- as.character(as.matrix(read.table(dir, nrow=1)))
    SNPnames <- c("ID","Pop",SNPnames)
    colnames(snpData) <- SNPnames
    return(snpData)
  }
  
  #Function to assign a binary 0 1 or 2 scoring where 2 is the larger allele (2 or 4) and 1 is the smaller allele (1 or 3)
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
 

   
##load required libraries -------------
  
  #Check to make sure the packages required are there
  packages <- "dplyr"
  if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    install.packages(setdiff(packages, rownames(installed.packages())))  
  } 
  
#Load library
        require(dplyr)

#load the genetic data
        P1_raw <- StructureRead(p1)
        P2_raw <- StructureRead(p2)
        P3_raw <- StructureRead(admixed)
        
        # names of the snps
        snpnames <- names(P1_raw[3:length(P1_raw)]) 
        
        #map to be used for missing alleles (if any present across loci)
        Allele_Map <- data.frame(SNP=snpnames,
                                 Allele1=rep(999,length(snpnames)),
                                 Allele2=rep(999,length(snpnames))) # 999 is a dummy placeholder
        
        for(i in 1:length(snpnames)){
            #unique alleles for a given snp (locus)
            alleleVals <- as.data.frame(table(as.character(c(P1_raw[,snpnames[i]],P2_raw[,snpnames[i]],P3_raw[,snpnames[i]]))))
            
            # if there is missing data (-9) delete it as a possibe allele
            if(length(which(alleleVals[,1]==(-9)))>0){
              alleleVals <- alleleVals[-which(alleleVals[,1]==(-9)),]
              } 
            
            Allele_Map[i,"Allele1"]=as.character(alleleVals[1,1])
            Allele_Map[i,"Allele2"]=as.character(alleleVals[2,1])
        }
        
        #NULL vectors
        P1_BGC <- NULL
        P2_BGC <- NULL
        Admixed_BGC <- NULL
        
        for(i in snpnames){
          # grab vector of alleles and delete replace missing values (-9) with NA
          P1_alleles <- P1_raw[,i];P1_alleles[which(P1_alleles==-9)]=NA
          P2_alleles <- P2_raw[,i];P2_alleles[which(P2_alleles==-9)]=NA 
          P3_alleles <- P3_raw[,i];P3_alleles[which(P3_alleles==-9)]=NA
          
          #If the population only has one allele for a given locus then a zero and the allele have be be added
          if(length(table(P1_alleles))==1){
            hold <- as.data.frame(table(P1_alleles))
            hold[,1] <- as.character(hold[,1]) 
            hold <- rbind(hold,c(setdiff(as.numeric(Allele_Map[which(Allele_Map$SNP==i),c("Allele1","Allele2")]),hold[1,1]),0)) #add in the extra value
            hold <- hold[order(hold[,1]),] #sort the right order from a conventional table output
            P1_alleles <- hold[,2]
            rm(hold)
          } else {P1_alleles <- as.character(as.data.frame(table(P1_alleles))[,2])}
          
          if(length(table(P2_alleles))==1){
            hold <- as.data.frame(table(P2_alleles))
            hold[,1] <- as.character(hold[,1]) 
            hold <- rbind(hold,c(setdiff(as.numeric(Allele_Map[which(Allele_Map$SNP==i),c("Allele1","Allele2")]),hold[1,1]),0)) #add in the extra value
            hold <- hold[order(hold[,1]),] #sort the right order from a conventional table output
            P2_alleles <- hold[,2]
            rm(hold)
          } else {P2_alleles <- as.character(as.data.frame(table(P2_alleles))[,2])}
        
         
          #for a given locus get the format for BGC
          P1_temp <- c(paste("locus_",i,sep=""),paste(P1_alleles[1],P1_alleles[2],sep=" "))
          P2_temp <- c(paste("locus_",i,sep=""),paste(P2_alleles[1],P2_alleles[2],sep=" "))
          
          #Combine output sequentially for each locus
          P1_BGC <- c(P1_BGC,P1_temp)
          P2_BGC <- c(P2_BGC,P2_temp)
        }

      
#Convert the admixed data to BGC format --------------
        
        MixedStruct=P3_raw

       #Data wrangle and construct the format for admixture populatons as per the instructions and table 2 in the BGC manual
          MixedData <- NULL
          for(col in names(MixedStruct)[3:length(MixedStruct)])
            {
            
            temp1 <- MixedStruct[,c("ID","Pop",col)]
            Locushold <- paste("locus_",col,sep="") # Start the locus data
              
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
                
                Locushold <- c(Locushold,paste("pop_",i,sep=""),temp5)
                
              } #end of population loop
            
              MixedData <- c(MixedData,Locushold) # add each successive locus 
            } #end of locus loop
  
  
##Save output for BGC formated for the parental and mixed populations ------------
      
      write.table(x = P1_BGC,file=paste(dir,"/",name,"_Parental1_BGC.txt",sep=""),
                  sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
      
      write.table(x = P2_BGC,file=paste(dir,"/",name,"_Parental2_BGC.txt",sep=""),
                  sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)  
      
      write.table(x = MixedData,file=paste(dir,"/",name,"_Admixed_BGC.txt",sep=""),
                  sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
      
} #end of function