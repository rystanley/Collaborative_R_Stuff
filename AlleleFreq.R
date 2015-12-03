AlleleFreqHeatMap <- function(GenePop,subs=NULL,keep=TRUE,POP="CHAR",refPop,OrderPops=NULL,standardize=TRUE,optimizer=TRUE,plot=TRUE){
  
  #Function details -------------
  
  #This function will convert a GenePop file into a tabular format for the calculation of 
  #allele frequencies of any outliers
  
  ## GenePop <-  the genepop file with loci either in 4 digit (e.g. 1010, 1020) or 
  ## 6 digit (e.g. 110120, 130140) format
  
  ## subs    <-  the loci names of interest or a vector which corresponds the the order of which
  ## they appear in the genepop file. 
  ## These can be either the order by which they occur or the exact name of the loci
  ## (i.e. subs <- c(1,2,3,4) would return the first 4 loci 
  ##       subs <- c("190-56","145_21",456_12") would return loci with defined names
  
  ## keep    <-  logical vector which defines whether you want to remove the loci or keep them. 
  ## the default is to keep them (keep=TRUE) assuming you are removing neutral markers 
  ## and only keeping the subs ("Outliers")
  
  ## POP     <-  is a vector which defines how the populations will be defined
  ## two options: "NUM" for numeric or "CHAR" for character (default)
  ## for the default "CHAR" option the output woudl be the non-numeric component of each
  ## sample id (e.g. "BON00012" would result in a strata called "BON")
  
  ## refPop <- is the reference population to which you choose the major allele frequency which
  ## will be reastered for the remaining populations. 
  ## (e.g. the furthest south or north if looking at a latitudinal gradient)
  
  ## OrderPops <- is a dataframe which contains the latitute and longitude of each sample location. 
  ## if no order is specified the funtion will use the default alphabetical order
  
  ## standardize <- logical vector (default: TRUE) which specifies whether (TRUE) or not (FALSE)
  ## to standardize the allele frequencies to 0-1 within each SNP among populations
  
  ## optimizer <- logical vector (default: TRUE) which specifies whether (TRUE) or not (FALSE)
  ## to try to optimize the differene between your gradient. 
  ## Difference optimizer. The major allele in the ref population might now have the largest
  ## allele frequency compared to the rest of the regions. This creates a situation where the 
  ## gradient is inverted at times leading to a difficult interpretation. This can be resolved
  ## using the reference population to invert situations where the population on one end of the 
  ## spectrum is low and inverts it to high. Note this only occurs in conjuction with standardize
  
  ## plot <- is a logical vector (default: TRUE), 
  ## TRUE: return the ggplot plot object 
  ## FALSE: return the HeatMapData
  
  
  #Libraries ----------
  #Check to make sure the packages required are there
  packages <- c("dplyr", "tidyr", "stringr","ggplot2","scales")
  if (length(setdiff(packages, rownames(installed.packages()))) > 0) { 
    install.packages(setdiff(packages, rownames(installed.packages())))  
  } 

  #load each library
    require(dplyr)
    require(tidyr)
    require(stringr)
    require(ggplot2)
    require(scales)
  
  #Functions -----------
  
    #Function to obtain allele information from a SNP vector. Can return Major and Minor allele
    #frequencies as well as the proportion of data which is missing data. 
  

      AlleleFreq <- function(x,Freq)
      {
        # Freq vector of 2 ("Major" or "Minor" or "missing") allele frequencie or percent missing
        
        if(nchar(as.character(x[1]))==4) # four character locus - 2 digit allele code
          {
          x1 <- as.numeric(substring(x,1,2)) # get first allele
          x2 <- as.numeric(substring(x,3,4)) # get second allele
        }
        
        if(nchar(as.character(x[1]))==6) # six character locus - 3 digit allele code
        {
          x1 <- as.numeric(substring(x,1,3)) # get first allele
          x2 <- as.numeric(substring(x,4,6)) # get second allele
        }
   
        x3 <- c(x1,x2) #combine together in one string
        if(length(unique(x3))>2){x4 <- x3[-which(x3==0)]} # delete the zeros (missing values)
        if(length(unique(x3))<3){x4 <- x3}
        # warning if data is coded wrong (too many alleles or 0 isn't the reference to )
        if(length(unique(x4))>2){print(paste("More than two alleles present. ",
                                             length(unique(x4))," unique alleles in vector:",
                                             paste(unique(x4),collapse=","),sep=""))}
        
        x5 <- as.data.frame(table(x4)/sum(length(x4))) # calculate frequency of each
        colnames(x5)=c("allele","Freq")
        
        #For missing values
        x6 <- as.data.frame(table(x3)/sum(length(x3))) # calculate frequency of each
        
        # if it is 50 - 50 return the first allele
        if(x5[1,2]==0.5){x5[1,2]=x5[1,2]+0.01} # add 0.01% to break the tie
        
        if(Freq=="Major"){return(x5[which(x5$Freq==max(x5$Freq)),])}
        if(Freq=="Minor"){return(x5[which(x5$Freq==min(x5$Freq)),])}
        if(Freq=="Missing"){return(x6[which(x6[,1]==0),])}
        
      }
    
  #Function to return allele frequency information for a vector of SNPs specific to the specifed allele
  
      AlleleFreqSubset <- function(x,allele=NULL,dig=6)
      {
        
        # SNP is the column name of the SNP of interest
        # allele is the specific allele of interest
        
        allele=as.character(allele)
        
        if(nchar(as.character(x[1]))==4) # four character locus - 2 digit allele code
        {
          x1 <- as.numeric(substring(x,1,2)) # get first allele
          x2 <- as.numeric(substring(x,3,4)) # get second allele
        }
        
        if(nchar(as.character(x[1]))==6) # six character locus - 3 digit allele code
        {
          x1 <- as.numeric(substring(x,1,3)) # get first allele
          x2 <- as.numeric(substring(x,4,6)) # get second allele
        }
        
        x3 <- c(x1,x2) #combine together in one string
        if(length(unique(x3))>2){x4 <- x3[-which(x3==0)]} # delete the zeros (missing values)
        if(length(unique(x3))<3){x4 <- x3}
        # warning if data is coded wrong (too many alleles or 0 isn't the reference to )
        if(length(unique(x4))>2){print(paste("More than two alleles present. ",
                                             length(unique(x4))," unique alleles in vector:",
                                             paste(unique(x4),collapse=","),sep=""))}
        
        x5 <- as.data.frame(table(x4)/sum(length(x4))) # calculate frequency of each
        colnames(x5) <- c("allele","Freq")
        
        # if the input allele is 0 just use the larges allele in terms of it's 
        # numeric code. This will ensure all populations are standardized. 
        
        if(allele==0){allele <- x5[which(x5$allele==max(as.numeric(as.character(x5$allele)))),"allele"]} 
        
        x5$allele <- as.character(x5$allele) # set as character for matching
        
        # no allele data at locus: returns NA
        if(sum(unique(x3))==0){return(NA)} 
        
        # only alternate alelle present: returns 0
        if(!(allele %in% unique(x3)) & sum(unique(x3))!=0){return(0)} 
        
        #allele present at locus: returns allele frequency
        if(allele %in% unique(x3) & sum(unique(x3))!=0){return(x5[which(x5$allele==allele),"Freq"])} 
          
      }
      
      #function to rescale a vector range to span 0 and 1
      scale01 <- function(x){(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))}

##Process the data -------------
      
      ## Stacks version information
          stacks.version <- GenePop[1,] # this could be blank or any other source. First row is ignored by GenePop
      
      #Remove first label of the stacks version
          GenePop <- as.vector(GenePop)
          GenePop <- GenePop[-1,]
      
      #Add an index column to Genepop and format as a dataframe
          GenePop <- data.frame(data=GenePop,ind=1:length(GenePop))
          GenePop$data <- as.character(GenePop$data)
      
      #ID the rows which flag the Populations
          Pops  <-  which(GenePop$data == "Pop" | GenePop$data =="pop" | GenePop$data == "POP")
          npops  <-  1:length(Pops)
      
      ## Seperate the data into the column headers and the rest
          ColumnData <- GenePop[1:(Pops[1]-1),"data"]
          snpData <- GenePop[Pops[1]:NROW(GenePop),]
      
      #Get a datafile with just the snp data no pops
          tempPops <- which(snpData$data=="Pop"| snpData$data =="pop" | snpData$data == "POP") ## Changed because we allowed
      ## alternate spelling on line 48, so had to change this so it would identify properly and not make an empty DF
          snpData <- snpData[-tempPops,]
      
      #Seperate the snpdata
      #First we pull out the population data which follows
      #"TEXT ,  "
          temp <- separate(snpData,data,into=c("Pops","snps"),sep=",")
          temp$snps <- substring(temp$snps,3) # delete the extra spaces at the beginning
          temp2 <- as.data.frame(do.call(rbind, strsplit(temp$snps," "))) #split characters by spaces
        
          #Contingency to see if R read in the top line as the "stacks version"
          if (length(temp2)!=length(ColumnData)){colnames(temp2) <- c(stacks.version,ColumnData)}
          if (length(temp2)==length(ColumnData)){colnames(temp2) <- ColumnData}
          if (length(temp2)!=length(ColumnData)){stacks.version="No stacks version specified"}
          
      ## Get the Alpha names from the 
          NamePops=temp[,1] # Names of each
          NameExtract=str_extract(NamePops, "[a-zA-Z]+" ) # extract the text from the individuals names to denote population
      
      ## Now add the population tags using npops (number of populations and Pops for the inter differences)
          tPops <- c(Pops,NROW(GenePop))
            PopIDs <- NULL
                for (i in 2:length(tPops)){
                  hold <- tPops[i]-tPops[i-1]-1
                  if(i==length(tPops)){hold=hold+1}
                  pophold <- rep(npops[i-1],hold)
                  PopIDs <- c(PopIDs,pophold)
                }
          
            ## add the population names
          if(POP=="NUM"){temp2$Pop <- PopIDs}
          if(POP=="CHAR"){temp2$Pop <- NameExtract}

#Subset the loci of interest ------------
   
    if(keep)
      {
        if(length(subs)>0) #is there anything to subset
        {
          
          if(is.numeric(subs)) # are they column entry (numeric)
          {
            temp2 <- temp2[,c(subs,length(temp2))]
          }
          
          if(!is.numeric(subs)) # are they column names (character)
          {
            temp2 <- temp2[,c(subs,"Pop")]
          }
        }
      }
    
    if(!keep)
    {
      if(length(subs)>0) #is there anything to subset
      {
        
        if(is.numeric(subs)) # are they column entry (numeric)
        {
          temp2 <- temp2[,-subs]
        }
        
        if(!is.numeric(subs)) # are they column names (character)
        {
          temp2 <- temp2[,-which(names(temp2) %in% subs)]
        }
      }
    }
    
## Get reference Allele for each SNP --------------
    
    ## Subset out the data for the refernce population
      refData <- subset(temp2,Pop==refPop)
      
    ## for each SNP find the major allele for the reference population
      AlleleReference=NULL
      for(i in 1:(length(refData)-1)){
        snpvec <- refData[,i]
        tempOut <- AlleleFreq(snpvec,Freq = "Major")[1]
        target_allele <- data.frame(SNP=colnames(refData)[i],allele=tempOut)
        AlleleReference <- rbind(AlleleReference,target_allele)
      }
   
## Calculate the SNP and allele population specific frequencies according popRef Major allele Freq ---------
      
      PopList <- unique(temp2$Pop)# List of all populations
      
      #create dataframe which will be populated with the allele frequencies *999 is a placeholder
      HeatMapData <- data.frame(Pop=rep(PopList,each=length(temp2)-1),
                                SNP=rep(AlleleReference$SNP,length(PopList)),
                                Freq=rep(999,length(PopList)*length(AlleleReference$SNP)))
      
      #Populate the dummy dataframe with the frequencies
      for (p in PopList) #Subset for each population
      {
        for (s in AlleleReference$SNP) #Calculate the allele frequency for each SNP or interest
        {
            a <- AlleleReference[which(AlleleReference$SNP==s),"allele"] #allele of interest based on refPop
            subData=temp2[which(temp2$Pop==p),s] #pull out the SNP data for a given population
            HeatMapData[which(HeatMapData$Pop==p & HeatMapData$SNP==s),
                        "Freq"] <- AlleleFreqSubset(subData,a)
        }
      }
    
    #standardize the allele frequency
      if(standardize) #0-1 scaled (TRUE)
      {
        HeatMapData <- as.data.frame(HeatMapData%>%group_by(SNP)%>%mutate(FreqStand=scale01(Freq))%>%ungroup())
      }
      
   
      
      #Order populations by latitude if specified
      if(length(OrderPops)!=0){HeatMapData$Pop=factor(HeatMapData$Pop,levels=OrderPops)}
      
      #Difference optimizer. The major allele in the ref population might now have the largest
      #allele frequency compared to the rest of the regions. This creates a situation where the 
      #gradient is inverted at times leading to a difficult interpretation. This can be resolved
      #using the reference population 
      
      #if the optimizer is goint to work we will have to play around with some general assumptions.
      #Here we assume a missing allele in the reference population (Freq == NA) will represent areas
      #low frequencies 
      HeatMapData$opt_FreqStand <- HeatMapData$FreqStand
      HeatMapData[which(is.na(HeatMapData$opt_FreqStand)),"opt_FreqStand"] <- 0
      
      if(standardize & optimizer)
        {
          #Identify SNP which need to be inverted
           invSNPS <- HeatMapData%>%filter(Pop==refPop,opt_FreqStand<0.5)%>%select_(.,"SNP")
           invSNPS <- as.vector(invSNPS[,1]) #convert to a vector of characters
           
           #invert the standardized allele frequency
           for(i in invSNPS)
             {
             HeatMapData[which(HeatMapData$SNP==i),"FreqStand"]=
             1-HeatMapData[which(HeatMapData$SNP==i),"FreqStand"]
             }
        }
      
#Create the heatmap --------------
      
      if(standardize){p1=ggplot(HeatMapData,aes(x=SNP,y=Pop,fill=FreqStand))} #0-1 scaled (TRUE)
      if(!standardize){p1=ggplot(HeatMapData,aes(x=Pop,y=SNP,fill=Freq))} # 0-1 scaled (FALSE)
      
      p1=p1+geom_tile()+
        scale_y_discrete(expand = c(0,0))+scale_x_discrete(expand = c(0,0))+
        theme_bw()+
        theme(legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1))+
        scale_fill_gradient(low="blue",high=muted("red"))
      
      if(standardize){p1=p1+labs(y="Population",x="SNP",fill="Standardized allele frequency")}
      if(!standardize){p1=p1+labs(y="Population",x="SNP",fill="Allele frequency")}

      
#Return function output -----------
      if(plot){return(p1)}
      if(!plot){return(HeatMapData)}

} #End function