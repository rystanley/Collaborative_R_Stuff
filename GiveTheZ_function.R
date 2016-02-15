NH_Zcore <- function(GetstheZdir, multiapplyZvec=NULL, applyuniqueZvec=NULL){

  ## this function will insert a column of known genotype categories <- Zvec into a NewHybrids formatted dataset
      ## for information on the known genotype categories, refer to the NewHybrids documentation (e.g.. Pure1 = z0, Pure2 = z1, F1 = z2 etc. [numbering follows that of NewHybrids])
    ## fucntion will match a vector of categories (Zvec) to individuals, then output a NH format file, with a name indicating that the category vector
      ## the user must provide the Zvec(s) to the function in the form of a two column .csv file(s) <- NOTE Zvec file must have column names, and the first column must be named "Individual" and
        ## contain the names of individuals to which categories must be given that match those in the first column of the NewHybrids dataset
    
  ## This is done in a batch fashion to every NewHybrids format file within a working folder - 
      #### IF a single file of Zvecs is to be applied to each file, USE multiapplyZvec 
          ###  multiapplyZvec must be specified as a file path + file name (e.g. "~/HoldenUniversity/CoolEthansComputer/HairDolls/Zvexxx.csv")
  
      #### IF each NewHybrids file is to be given a UNIQUE vector of Zvecs use applyuniqueZvec ---- NOTE <- to apply UNIQUE Zvecs, the Zvec files must all be placed in a single folder
          ### separate from the NewHybrids files AND they MUST follow the file name convention "NHFileName_Zvec.csv" where NHFileName is the same name as the file to which the Zvec is
          ## to be applied - consequently, the number of NewHybrids files = number of Zvec files
              ### applyuniqueZvec must be spedifed as a file path to the folder in which the Zvec files live (e.g. "~/HarrisonUniversity/DEANGordonPritchard/LettersOfRecommendation/")
  
  ## The vector of genotype categories NEED NOT be the same length as the NewHybrids file - thus, time can be saved 
  
  
## Check required packages are installed - if not, download
packages <- c("tidyr", "stringr")
  if (length(setdiff(packages, rownames(installed.packages()))) > 0) { 
    install.packages(setdiff(packages, rownames(installed.packages())))
      } ## end of check package loop
## initiate the required packages
library(tidyr)
require(stringr)

## make sure that both multiapplyZvec and applyuniqueZvec haven't both been called

if((length(multiapplyZvec) >0 ) && length(applyuniqueZvec) > 0){
  stop("One or the other bud. You can have multiapplyZvec, or applyuniqueZvec. You aren't a two-year-old at their birthday, so don't get greedy")
}


## files to be converted should all be in a single folder 
## get a list of the files to be converted - then initiate a loop to change each in turn
  NHdataGet <- list.files(GetstheZdir)
      ## assumes you have given files in standard NewHybrids format
  
  ## check if applyuniqueZvec has been called - and if so, error check the shit out of it
  if(length(applyuniqueZvec > 0)){
    ZvecList <- list.files(applyuniqueZvec)
    ## check that an indivual Zvec has been specified for each NHfile
    if(length(NHdataGet) != length(ZvecList)){
      stop("Come on bud, number of Zvecs must equal number of NewHybrids Files") ## if it doesn't let 'em know
    }
    ## check that the files are name right/same so can effectively match the Zvec to the NewHybrids file
        ## need to remove some stuff from teh file names so can check
    NHdataGetCheck <- mapply(FUN = gsub, x = NHdataGet, pattern = ".txt", replacement = "") ## remove from all names in list
    ZvecListCheck <- mapply(FUN = gsub, x = ZvecList, pattern = "_Zvec.csv", replacement = "") ## remove from all names in list

    ZvecMatchLength <- length(sapply(ZvecListCheck, function(x) any(sapply(NHdataGetCheck, str_detect, string = x)))) ## get the lenght of the number of Zvec files that match names in NH files
    if(length(NHdataGet) != ZvecMatchLength){
      stop("Bud, have a look at the files your're trying to give me over here. Either you haven't named 'em right, or you can't count.") # if doesn't match, let 'em know
    }
    
    
  }
  
  for(i in NHdataGet){
    ## get the ith file in the folder
    NHdata <- read.table(paste0(GetstheZdir, i), header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)
    
    ## the first 5 rows contain information about the # of indviduals, loci, etc. need to remove, to work on the genotype data, but 
    #Save the first 5 rows to be added back in at the end.
    addbackin<-NHdata[1:5,]
    
    ## read the ith file back in, called somethign different, without the first 5 rows so can work on genotype data
    NHdata2 <- read.table(paste0(GetstheZdir, i), header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE, skip = 5)
    
    ## in order for separate to work properly, it must be given the nubmer of columns to divide - this is the number of loci included in NH file
        ## number of loci is given in the first 5 rows, and can be deterimined using a regex
    ## get number of loci using regex
    SNPrep <- str_extract(string = addbackin, pattern = "NumLoci [:digit:]{2,5}") ## this does not cleanly take the number, so some cleaning is required
      SNPrep <- as.numeric(gsub(x = SNPrep, pattern = "NumLoci ", replacement = "")) 
      SNPrep <- SNPrep[which(is.na(SNPrep)==FALSE)]

  ## give a column name to allow separate to funciton    
  names(NHdata2) <- "Data2"
  ## separate and give a new name - separates the first column which is the indiviudal IDs from the genotype/loci data - the known genotype category assignments (Zs) must be 
    ## inserted between the indiviual IDs and their corresponding genotypes
  Nhdata3<-separate(data = NHdata2,col=Data2,into = c("Individual",rep("SNP",SNPrep)),sep = " ")


  ## import the csv of z codes and indivudal IDs ; Then merge these with the corresponding genotype data

  if(length(multiapplyZvec) > 0 ){
    Zscorevector <- read.csv(multiapplyZvec)
  }
  if(length(applyuniqueZvec) > 0 ){
    ith_NHdata <- gsub(x = i, pattern = ".txt", replacement = "")
    ZvecListCheck <- mapply(FUN = gsub, x = ZvecList, pattern = "_Zvec.csv", replacement = "") ## remove from all names in list
    
    Zscorevector <- read.csv(paste0(applyuniqueZvec, ZvecList[which(str_detect(ZvecListCheck, pattern = ith_NHdata))]))
      
  }
  
#Zscorevector<-read.csv("/Users/brendanwringe/Desktop/DFO Aquaculture Interaction/South West Rivers Analysis/Frequency Based Sim/West/Zvector.csv")
Zscorevector$Individual <- as.character(Zscorevector$Individual)
#Merge
NHfinal<- merge(y=Nhdata3, x=Zscorevector, by="Individual", all=TRUE)
NHfinal <- NHfinal[order(as.numeric(NHfinal$Individual)),]

#head(NHfinal)

#Now make this an actual file readable by NewHybrids
Loci <- do.call(paste,c(NHfinal[,], sep=" "))
Loci<-gsub("NA"," ",x = Loci)
Loci2<-c(addbackin,Loci)


  NHdataZed <- gsub(x = i, pattern = ".txt", replacement = "_Zed.txt")
  
  
write.table(Loci2, file = paste0(GetstheZdir, NHdataZed),col.names=FALSE,row.names=FALSE,quote=FALSE)


  } ## End of loop
  
} ## End of function


