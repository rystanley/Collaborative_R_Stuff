require(stringr)
require(plyr)
require(tidyr)

rm(list = ls() )

GenePopData <- read.table("GRR.IB.to.SIM.txt", header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)
GenePop <- GenePopData
out.name <- "GRR_FreqSim_N502"
pop.groups <- c("FRM", "WLD")
sample.size <- 50




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
    ColumnData <- GenePop[1:(Pops[1]-1),"data"]  ### SNP Names
    snpData <- GenePop[Pops[1]:NROW(GenePop),]  ### Genotypes

#Get a datafile with just the snp data no pops
    tempPops <- which(snpData$data=="Pop"| snpData$data =="pop" | snpData$data == "POP")
    snpData <- snpData[-tempPops,]

#Seperate the snpdata
#First we pull out the population data which follows "TEXT ,  "
    temp <- separate(snpData,data,into=c("Pops","snps"),sep=",")
    temp$snps <- substring(temp$snps,3) # delete the extra spaces at the beginning
    temp2 <- data.frame(do.call(rbind, str_extract_all(temp$snps, "[0-9]{3}"))) 
   
## Going to have to break the two alleles of the SNPS apart - this will thus double the number of columns
    ## SO <- will want to have SNP_A and SNP_A2 
ColumnData2 <- ColumnData ## Duplicatet the SNP names
ColumnData2 <- paste(ColumnData2, "2", sep = ".")    ## add .2 to each duplicated name
     
## can't just append the duplicated names to the end of the original names - have to intersperse them
places = rep(1:length(ColumnData)*2) ### creates a list of even numbers 2X as long as the number of columns
  ##i.e. the lenght of the original plus the duplicated names 
  ## - this will also mark the position to insert the duplicates
ColumnData.Dup = rep(NA, times = length(ColumnData)*2) ### make an object to feed names into
for(i in 1:length(ColumnData)){ ### for loop to add original, then duplicated name
  a = places[i]-1 ## Original names go first, and are in the odd positions
  b = places[i] ### Duplicated names go second and are in the even posoitons
  Col.name.orig = ColumnData[i] ## Get the name in the ith position
  Col.name.plus2 = ColumnData2[i] ## get the name in the ith positon
  ColumnData.Dup[a] = Col.name.orig ## add the original name to the new vector
  ColumnData.Dup[b] = Col.name.plus2 ## add the duplicate name to the new vector
}


    #Contingency to see if R read in the top line as the "stacks version" -- modified to deal with the duplicated SNP names
    if (length(temp2)!=length(ColumnData.Dup)){colnames(temp2) <- c(stacks.version, paste(stacks.version, "2", sep = "."),ColumnData.Dup)}
    if (length(temp2)==length(ColumnData.Dup)){colnames(temp2) <- ColumnData.Dup}
    #if (length(temp2)/2!=length(ColumnData)){stacks.version="No stacks version specified"}
    
## Get the Alpha names 
    NamePops=temp[,1] # Names of each
  
   if(length(pop.groups) == 0){ ### If unique grouping IDs ≠ number of "Pop" user must give vector of groupings 
                                ### equal to number of "Pop" or else the function will fail
    NameExtract=str_extract(NamePops, "[A-z]{3}" ) ### if looking at higher order grouping (i.e. pops in 
        # regions) can have more unique coding than "Pop" - will want to remove original names so can
        ## keep track of which unique groupings cross. i.e. Cross by "Pop", but remember ID of parents
   }
    
   
  # extract the text from the individuals names to denote population
  ## Now add the population tags using npops (number of populations and Pops for the inter differences)
    tPops <- c(Pops,NROW(GenePop))
      PopIDs <- NULL
          for (i in 2:length(tPops)){
            hold <- tPops[i]-tPops[i-1]-1
            if(i==length(tPops)){hold=hold+1}
            pophold <- rep(npops[i-1],hold)
            PopIDs <- c(PopIDs,pophold)
          }
    
    temp2$Pop <- PopIDs;
    
    
     if(length(pop.groups)!=0){
     hold.names=str_extract(NamePops, "[A-z]{3}" ) ## This may need to be improved in published version
        for(i in 1:length(unique(PopIDs))){
          u.ID.no <- unique(PopIDs)[i]
          to <- min(which(PopIDs==u.ID.no))
          from <- max(which(PopIDs==u.ID.no))
      hold.names[to:from] = paste(pop.groups[i], hold.names[to:from], sep=".")
    }
    NameExtract <- hold.names
     }
    
     ## get the nubmer of indivudals within each "Pop" grouping
    PopLengths <- table(temp2$Pop)
    
     ## Need to be able to tell what row each individual is in, and what population it is
    ind.vector = c(1:nrow(temp)) ### make a vector that is the number of individuals
    ind.matrix = data.frame(temp2$Pop, ind.vector) ## add populatuions to that
    
temp.split <- split(x = temp2, f = temp2$Pop)
    
pop.recall <- NULL
    for(i in 1:length(temp.split)){
      popn <- paste(pop.groups[i], "pop", sep = "_")
      temp.split.hold = temp.split[[i]]
      temp.split.hold = temp.split.hold[-which(names(temp.split.hold) == "Pop")]
      assign(x = popn, value = temp.split.hold, envir = globalenv())
      pop.recall <- c(pop.recall, popn)
      
    }
i = 2
k=1
mat.name.recall <- NULL
for(i in 1:length(pop.recall)){
    temp.mat <- data.frame(matrix(vector(), 2, length(temp2)/2))
    pop.get <- get(pop.recall[i])
    
      temp.mat.hold <- NULL
      for(k in 1:nrow(pop.get)){
        ind.hold <- pop.get[k,]
        temp.mat[1,] <- t(t(ind.hold[c(T,F)]))
        temp.mat[2,] <-  t(t(ind.hold[c(F,T)]))
        
        temp.mat.hold <- rbind(temp.mat.hold, temp.mat)
    
          }
    
      mat.out.name <- paste(pop.recall[i],"matrix", sep = "_")
      assign(x = mat.out.name, value = temp.mat.hold, envir = globalenv())
      mat.name.recall <- c(mat.name.recall, mat.out.name)
      
}


#if(sim.pure == "Yes"){
k=1
### MAKE PURE CROSS
pure.name.recall <- NULL
for(k in 1:length(pop.groups)){

pop1 <- get(mat.name.recall[k])
pop2 <- get(mat.name.recall[k])
off.interspersed.out <- NULL
for(i in 1:sample.size){

hold.off.pop1 <- apply(pop1, FUN = sample, 2, 1)
hold.off.pop2 <- apply(pop2, FUN = sample, 2, 1)
hold.off.interspersed <- data.frame(c(rbind(hold.off.pop1, hold.off.pop2)))

off.interspersed.out <- rbind(off.interspersed.out, t(hold.off.interspersed))


}

pure.name <- paste("Pure", pop.groups[k], sep = "_")
pure.name.recall <- c(pure.name.recall, pure.name)

assign(x = pure.name, value = off.interspersed.out, envir = globalenv())

}

#}




inv.pure.name.recall <- NULL
for(i in 1:length(pop.recall)){
    temp.mat <- data.frame(matrix(vector(), 2, length(temp2)/2))
    pop.get <- get(pure.name.recall[i])
    
      temp.mat.hold <- NULL
      for(k in 1:nrow(pop.get)){
        ind.hold <- pop.get[k,]
        temp.mat[1,] <- t(t(ind.hold[c(T,F)]))
        temp.mat[2,] <-  t(t(ind.hold[c(F,T)]))
        temp.mat.hold <- rbind(temp.mat.hold, temp.mat)
    
    }
    
      inv.pure.out.name <- paste(pure.name.recall[i],"inv", sep = "_")
      assign(x = inv.pure.out.name, value = temp.mat.hold, envir = globalenv())
      inv.pure.name.recall <- c(inv.pure.name.recall, inv.pure.out.name)
      
}




### MAKE F1 CROSS

pop1 <- get(inv.pure.name.recall[1])
pop2 <- get(inv.pure.name.recall[2])
F1.out <- NULL
for(i in 1:sample.size){

hold.off.pop1 <- apply(pop1, FUN = sample, 2, 1)
hold.off.pop2 <- apply(pop2, FUN = sample, 2, 1)
hold.off.interspersed <- data.frame(c(rbind(hold.off.pop1, hold.off.pop2)))

F1.out <- rbind(F1.out, t(hold.off.interspersed))


}


temp.mat <- data.frame(matrix(vector(), 2, length(temp2)/2))
    pop.get <- F1.out
    
      inv.F1 <- NULL
      for(k in 1:nrow(pop.get)){
        ind.hold <- pop.get[k,]
        temp.mat[1,] <- t(t(ind.hold[c(T,F)]))
        temp.mat[2,] <-  t(t(ind.hold[c(F,T)]))
        inv.F1 <- rbind(inv.F1, temp.mat)
    
    }


### MAKE F2 CROSS

pop1 <- inv.F1
pop2 <- inv.F1
F2.out <- NULL
for(i in 1:sample.size){

hold.off.pop1 <- apply(pop1, FUN = sample, 2, 1)
hold.off.pop2 <- apply(pop2, FUN = sample, 2, 1)
hold.off.interspersed <- data.frame(c(rbind(hold.off.pop1, hold.off.pop2)))

F2.out <- rbind(F2.out, t(hold.off.interspersed))


}


### MAKE Back CROSS


BC.name.recall <- NULL
for(k in 1:length(pop.groups)){

pop1 <- get(inv.pure.name.recall[k])
pop2 <- inv.F1
off.interspersed.out <- NULL
for(i in 1:sample.size){

hold.off.pop1 <- apply(pop1, FUN = sample, 2, 1)
hold.off.pop2 <- apply(pop2, FUN = sample, 2, 1)
hold.off.interspersed <- data.frame(c(rbind(hold.off.pop1, hold.off.pop2)))

off.interspersed.out <- rbind(off.interspersed.out, t(hold.off.interspersed))


}

BC.name <- paste("BC", pop.groups[k], sep = "_")
BC.name.recall <- c(BC.name.recall, BC.name)

assign(x = BC.name, value = off.interspersed.out, envir = globalenv())

}



pure.name.recall
for(i in 1:length(pure.name.recall)){
  off.name <- paste(pure.name.recall[i], c(1:sample.size), sep="_")
  hold.dat <- get(pure.name.recall[i])
  hold.dat <- data.frame(off.name, hold.dat)
  ColumnData.Dup.insert = c("ID", ColumnData.Dup)
  colnames(hold.dat) = ColumnData.Dup.insert
  assign(x = pure.name.recall[i], value = hold.dat, envir = globalenv())
}


f1.off.name <- paste("F1", c(1:sample.size), sep = "_")
F1.out <- data.frame(f1.off.name, F1.out)
colnames(F1.out) <- c("ID", ColumnData.Dup)

f2.off.name <- paste("F2", c(1:sample.size), sep = "_")
F2.out <-  data.frame(f2.off.name, F2.out)
colnames(F2.out) <- c("ID", ColumnData.Dup)

BC.name.recall
for(i in 1:length(BC.name.recall)){
  off.name <- paste(BC.name.recall[i], c(1:sample.size), sep="_")
  hold.dat <- get(BC.name.recall[i])
  hold.dat <- data.frame(off.name, hold.dat)
  ColumnData.Dup.insert = c("ID", ColumnData.Dup)
  colnames(hold.dat) = ColumnData.Dup.insert
  assign(x = BC.name.recall[i], value = hold.dat, envir = globalenv())
}


pure.name.recall
BC.name.recall
b=1
for(b in 1:length(pure.name.recall)){
  
  fam.to.bind.name <- pure.name.recall[b]
  fam.to.bind <- get(fam.to.bind.name)
  indiv.hold <- fam.to.bind[,1]
  loci.bind <- which(str_detect(string = colnames(fam.to.bind), pattern = "\\.2")==TRUE)
    col.out <- NULL
    
    for(s in 1:length(loci.bind)){
      place.1 <- (loci.bind[s]-1)
      place.2 <- loci.bind[s]
      hold.col <- paste0(fam.to.bind[,place.1], fam.to.bind[,place.2])
      col.out <- cbind(col.out, hold.col)

        }
  
  fam.reord <- cbind(indiv.hold,col.out)
  colnames(fam.reord) <- c(colnames(fam.to.bind[1]), colnames(fam.to.bind[c((loci.bind-1))]))
  assign(x = fam.to.bind.name, fam.reord, envir = globalenv())
  

for(b in 1:length(pure.name.recall)){
  
  fam.to.remove.untyped.name <- pure.name.recall[b]
  
  fam.to.remove.untyped <- get(fam.to.remove.untyped.name)
  fam.to.remove.untyped[which(str_detect(string = fam.to.remove.untyped, pattern = "000")==TRUE)] = "000000"
  assign(x = fam.to.remove.untyped.name, value = fam.to.remove.untyped, envir = globalenv())
}
  
    }

## F1
  
  fam.to.bind.name <- "F1.out"
  
  fam.to.bind <- get(fam.to.bind.name)
  indiv.hold <- fam.to.bind[,1]
  loci.bind <- which(str_detect(string = colnames(fam.to.bind), pattern = "\\.2")==TRUE)
    col.out <- NULL
    
    for(s in 1:length(loci.bind)){
      place.1 <- (loci.bind[s]-1)
      place.2 <- loci.bind[s]
      hold.col <- paste0(fam.to.bind[,place.1], fam.to.bind[,place.2])
      col.out <- cbind(col.out, hold.col)

        }
  
  fam.reord <- cbind(indiv.hold,col.out)
  colnames(fam.reord) <- c(colnames(fam.to.bind[1]), colnames(fam.to.bind[c((loci.bind-1))]))
  assign(x = fam.to.bind.name, fam.reord, envir = globalenv())

  
  fam.to.remove.untyped.name <- "F1.out"
  
  fam.to.remove.untyped <- get(fam.to.remove.untyped.name)
  fam.to.remove.untyped[which(str_detect(string = fam.to.remove.untyped, pattern = "000")==TRUE)] = "000000"
  assign(x = fam.to.remove.untyped.name, value = fam.to.remove.untyped, envir = globalenv())


  
  fam.to.bind.name <- "F2.out"
  
  fam.to.bind <- get(fam.to.bind.name)
  indiv.hold <- fam.to.bind[,1]
  loci.bind <- which(str_detect(string = colnames(fam.to.bind), pattern = "\\.2")==TRUE)
    col.out <- NULL
    
    for(s in 1:length(loci.bind)){
      place.1 <- (loci.bind[s]-1)
      place.2 <- loci.bind[s]
      hold.col <- paste0(fam.to.bind[,place.1], fam.to.bind[,place.2])
      col.out <- cbind(col.out, hold.col)

        }
  
  fam.reord <- cbind(indiv.hold,col.out)
  colnames(fam.reord) <- c(colnames(fam.to.bind[1]), colnames(fam.to.bind[c((loci.bind-1))]))
  assign(x = fam.to.bind.name, fam.reord, envir = globalenv())

  
  fam.to.remove.untyped.name <- "F2.out"
  
  fam.to.remove.untyped <- get(fam.to.remove.untyped.name)
  fam.to.remove.untyped[which(str_detect(string = fam.to.remove.untyped, pattern = "000")==TRUE)] = "000000"
  assign(x = fam.to.remove.untyped.name, value = fam.to.remove.untyped, envir = globalenv())

  b=1
for(b in 1:length(BC.name.recall)){
  
  fam.to.bind.name <- BC.name.recall[b]
  
  fam.to.bind <- get(fam.to.bind.name)
  indiv.hold <- fam.to.bind[,1]
  loci.bind <- which(str_detect(string = colnames(fam.to.bind), pattern = "\\.2")==TRUE)
    col.out <- NULL
    
    for(s in 1:length(loci.bind)){
      place.1 <- (loci.bind[s]-1)
      place.2 <- loci.bind[s]
      hold.col <- paste0(fam.to.bind[,place.1], fam.to.bind[,place.2])
      col.out <- cbind(col.out, hold.col)

        }
  
  fam.reord <- cbind(indiv.hold,col.out)
  colnames(fam.reord) <- c(colnames(fam.to.bind[1]), colnames(fam.to.bind[c((loci.bind-1))]))
  assign(x = fam.to.bind.name, fam.reord, envir = globalenv())
  
}

for(b in 1:length(BC.name.recall)){
  
  fam.to.remove.untyped.name <- BC.name.recall[b]
  
  fam.to.remove.untyped <- get(fam.to.remove.untyped.name)
  fam.to.remove.untyped[which(str_detect(string = fam.to.remove.untyped, pattern = "000")==TRUE)] = "000000"
  assign(x = fam.to.remove.untyped.name, value = fam.to.remove.untyped, envir = globalenv())
}

#Now recompile the GenePop format


## Function for inserting rows
  insert.vals <- function(Vec,breaks,newVal){
    break.space <- 1:(length(breaks))
    breaks <- breaks+break.space-1 #To space out the insertion points.
    newvec <- rep(NA,length(Vec)+length(breaks)) #Preallocate memory by creating final dataframe.
    for(i in 1:length(breaks)){newvec[breaks[i]]=newVal} #Insert added rows into new dataframe>
    x <- 1:length(newvec)
    x <- x[-(breaks)] #Finding the rows of the new dataframe that will receive old rows
    for(i in 1:length(Vec)){newvec[x[i]]=Vec[i]} 
    return(newvec)}
    

  
pop.names <- c(pure.name.recall, "F1.out", "F2.out", BC.name.recall)  

for(i in 1:length(pop.names)){
  temp.hold <- get(pop.names[i])
  temp.hold[,1] <- paste0(pop.names[i], "_", c(1:nrow(temp.hold)), " ,  ")
  assign(x = pop.names[i], value = temp.hold, envir = globalenv())
  
}

sim.out <- NULL
for(i in 1:length(pop.names)){
  hold.pop <- get(pop.names[i])
  sim.out <- rbind(sim.out, hold.pop)
}

df <- as.data.frame(sim.out[,]) #conver the data to a dataframe 
PopNames <- as.character(df$ID) #remove the populaton ID data
sim.out2 <- df[,2:length(df)] #create second vector which is just the loci data
LociNames <- colnames(sim.out2) #grab the loci names
spacer <- rep(" ,  ",length(PopNames)) # create a spacer vector (this could aslo just use paste but this is easier for bug checking)
Loci.sim <- do.call(paste, c(sim.out2[,], sep = " ")) #paste all the loci together with a " " as a sperator
out.df <- as.data.frame(cbind(PopNames,spacer,Loci.sim)) #bind the names, predefined spacer array and single loci strings togther
df2 <- do.call(paste,c(out.df[,],sep="")) #paste all the columsn into a vector of length(pop.names)*sample.size long

# find the positions by which "Pop" should be added in
PopLengths <- rep(sample.size, length=length(pop.names))
PopPosition <- NULL
PopPosition <- c(PopLengths[1]+1,rep(NA,length(pop.names)-2))
for (i in 2:(length(PopPosition))){
  PopPosition[i] <- PopLengths[i]+PopPosition[i-1]
}

#Insert the value of "Pop" which partitions the data among populations #only if more than one population
Loci <- insert.vals(Vec=df2,breaks=PopPosition,newVal="Pop")

#combine all stacks version, loci labels, first pop entry, and loci data together in one long vector
Loci.out <- c(stacks.version,LociNames,"Pop",Loci)

#save the output as text file. 
write.table(Loci.out,"testing123.txt",col.names=FALSE,row.names=FALSE,quote=FALSE)
    

