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
  
  # this function will identify monomorphic alleles and return a list of loci headings
  MonoFunc <- function(x){
    
    # x is the genepop allele counts per population created using the 'genind2genepop' function from adegenet
    
    Allele_Data=t(x@tab)
    s1 <- data.frame(Allele=rownames(Allele_Data))
    x_allele <- cbind(s1,Allele_Data)
    rownames(x_allele) <- NULL #remove the rownames which are redundant
    
    #Grab the names of loci in the admixed data
    x1 <- as.character(x_allele$Allele)
    x1 <- lapply(strsplit(x1, '.', fixed=TRUE), '[[', 1)
    
    #names of the loci with only one allele (need to be removed)
    MonoNames <- names(table(unlist(x1))[which(as.vector(table(unlist(x1)))==1)])
    
    if(length(MonoNames)>0){return(MonoNames)}
    if(length(MonoNames)==0){return(NA)}
  }

  #Function to wrangle data in BGC format according to table 1 in BGC manual
  BGCFormat=function(x){
    #extract the Loci names
    LociNames <- lapply(strsplit(as.character(x$Allele), '.', fixed=TRUE), '[[', 1)
    x$LociNames <- unlist(LociNames)
    
    #temp <- x[,c("Allele","LociNames",Pop)] # subset out just the data to a specific population
    
    tempData <- data.frame(Loci=x[seq(1,nrow(x),2),"LociNames"],
                           Allele=paste(x[seq(1,nrow(x),2),2],
                                        x[seq(2,nrow(x),2),2],sep=" "))
    tempData$Allele=as.character(tempData$Allele) # convert the joined populatoin based allele count for a given locus to character
    
    #format for BGC output accordig to example table 1 in the manual
    Output <- NULL
    for (i in 1:nrow(tempData))
    {
      Output <- c(Output,paste("locus ",tempData[i,"Loci"],sep=""),tempData[i,"Allele"])
    }
    
    return(Output)
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
  packages <- c("dplyr", "adegenet")
  if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    install.packages(setdiff(packages, rownames(installed.packages())))  
  } 
  
        require(adegenet)
        require(dplyr)

#load the genetic data
        P1_raw <- StructureRead(p1)
        P2_raw <- StructureRead(p2)
        P3_raw <- StructureRead(admixed)
  
#Import genind objects and use the file specifics from the raw genetic data (e.g. # individuals/loci)
        P1data <- import2genind(file = p1, onerowperind=F, n.ind=length(unique(P1_raw$ID)), 
                            n.loc=(length(P1_raw)-2), row.marknames=1, col.lab=1, col.pop=2, ask=F)
        
        P2data <- import2genind(file= p2, onerowperind=F, n.ind=length(unique(P2_raw$ID)), 
                            n.loc=(length(P2_raw)-2), row.marknames=1, col.lab=1, col.pop=2, ask=F)
        
        ADMdata <- import2genind(file=admixed, onerowperind=F, n.ind=length(unique(P3_raw$ID)), 
                             n.loc=(length(P3_raw)-2), row.marknames=1, col.lab=1, col.pop=2, ask=F)
        
#convert genotypes data (genind) into alleles counts per population (genpop).
        P1 <- genind2genpop(P1data)
        P2 <- genind2genpop(P2data)
        P3 <- genind2genpop(ADMdata)

#Grab the allele frequency format data * the code on the web is missing something
        P1_Data=t(P1@tab) # parental data 1
        P1_allele <- data.frame(Allele=rownames(P1_Data))
        P1_allele <- cbind(P1_allele,P1_Data)
        rownames(P1_allele) <- NULL
        
        P2_Data=t(P1@tab) # parental data 2
        P2_allele <- data.frame(Allele=rownames(P2_Data))
        P2_allele <- cbind(P2_allele,P2_Data)
        rownames(P2_allele) <- NULL
        
        P3_Data=t(P1@tab) #admixed populations
        P3_allele <- data.frame(Allele=rownames(P3_Data))
        P3_allele <- cbind(P3_allele,P3_Data)
        rownames(P3_allele) <- NULL

#remove unneeded data from this step
        rm(P1_Data,P2_Data,P3_Data)
        
# identify monomorphic alleles 
        
      MonoAlleles_P1 <- MonoFunc(P1)
      MonoAlleles_P2 <- MonoFunc(P2)
      MonoAlleles_P3 <- MonoFunc(P3)
      
      #Data frame of monomorphic alleles which will be saved in output
      Monodf <- data.frame(Pops=c(rep("P1",length(MonoAlleles_P1)),
                                  rep("P2",length(MonoAlleles_P2)),
                                  rep("admixed",length(MonoAlleles_P3))),
                           snps=c(MonoAlleles_P1,MonoAlleles_P2,MonoAlleles_P3))
      
      write.table(x = Monodf,file=paste(dir,"/",name,"_MonomorphicAlleleList.txt",sep=""),
                  sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
      
      #List of monomorophic alleles to be removed from each dataset
      MonoAlleles <- c(MonoAlleles_P1,MonoAlleles_P2,MonoAlleles_P3)
      MonoAlleles <- MonoAlleles[!is.na(MonoAlleles)] #remove NAs
      MonoAlleles <- unique(MonoAlleles) #remove repeates
      
      #Clean up the dataoutputs to remove the alleles which are monomorphic
      
      #Parental populations will have the allele counts per pop files (genind2genpop) trimmed
      P1_snps <- lapply(strsplit(as.character(P1_allele$Allele), '.', fixed=TRUE), '[[', 1);
      P1_snps <- unlist(P1_snps)
      P1_allele_sub <- P1_allele[!(P1_snps %in% MonoAlleles),]
      
      P2_snps <- lapply(strsplit(as.character(P2_allele$Allele), '.', fixed=TRUE), '[[', 1);
      P2_snps <- unlist(P2_snps)
      P2_allele_sub <- P2_allele[!(P2_snps %in% MonoAlleles),]
      
      #add mixed population will have the raw output from Structure trimmed
      MixedStruct <- P3_raw[,!(colnames(P3_raw)%in%MonoAlleles)]
      
#Convert the Population Data to BGC format ---------------
      
      P1_BGC <- BGCFormat(P1_allele_sub)
      P2_BGC <- BGCFormat(P2_allele_sub)
      
#Convert the admixed data to BGC format --------------

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
  
  
##Save output for BGC formated for the parental and mixed populations ------------
      
      write.table(x = P1_BGC,file=paste(dir,"/",name,"_Parental1_BGC.txt",sep=""),
                  sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
      
      write.table(x = P2_BGC,file=paste(dir,"/",name,"_Parental2_BGC.txt",sep=""),
                  sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)  
      
      write.table(x = MixedData,file=paste(dir,"/",name,"_Admixed_BGC.txt",sep=""),
                  sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
      
} #end of function