Hybridpower_comparison <-function(dir,filetag="",Thresholds=c(0.5,0.6,0.7,0.8,0.9),addThresh=FALSE,samplesize=200){
  
## this function will estimate the power of assignment success and associated variability from New Hybrids for 6 possible classes
#  Pure1, Pure2, F1, F2, BC1, BC2. The code is based on different simulations of individuals and repeat runs through New Hybrids (NH)
# for each simulation. The output is a series of diagnostic plots and data summaries. 
  
## dir - path directory which holds the output from different runs through New Hybrids (e.g. 3 simulations with 3 replicate runs each through NH)
#        note that this directory should only hold the output folders.

## filetag - this is a name tag which will be added to the plots

## Thresholds - this is a vector of thesholds (default c(0.6,0.7,0.8,0.9)) which will to compare among different thesholds to compare among different numbers of SNPs

## addThresh - a logical vector (default: FALSE) which specifies whether the threshold values should be added to the summary plot.     
  
## samplesize - is the number of fish per NH class (default = 200). The assumption is there is equal numbers for each class. If not there must be a 0 or NA added hold the data structure


  #Check to make sure the packages required are there and if not install them
  packages <- c("dplyr", "tidyr", "stringr","ggplot2","reshape2","grid","scales")
  if (length(setdiff(packages, rownames(installed.packages()))) > 0) { 
    install.packages(setdiff(packages, rownames(installed.packages())))
  } 
  
  #load each library
  require(dplyr)
  require(ggplot2)
  require(tidyr)
  require(stringr)
  require(reshape2)
  require(grid)
  require(scales)


      #set directory for which holds the New Hybrids output folders
      filedir <- dir

      lfiles <- setdiff(list.files(dir),"Figures") #ignores Figures folder in case this is run more than once
     
      if(length(which(list.files(dir)=="Figures"))==0){dir.create(paste0(dir,"Figures"))} # if there isn't a 'Figures' folder for output create one
      if(length(which(list.files(paste0(dir,"Figures"))=="pdf"))==0){dir.create(paste0(dir,"Figures/pdf"))} #create a folder for pdfs
      if(length(which(list.files(paste0(dir,"Figures"))=="jpg"))==0){dir.create(paste0(dir,"Figures/jpg"))} #create a folder for jpgs 
      
  
    # Collate the output from New Hybrids together ('p of z' files)
        output <- NULL
        for (i in lfiles)
        {
              tempfiles <- list.files(paste0(filedir,i))
              pzfile <- tempfiles[grep("PofZ",tempfiles)]
              tempfile <- read.table(paste0(filedir,i,"/",pzfile),head=T)
              
              LociandAlleles <- tempfiles[grep("LociAndAlleles", tempfiles)]
              LandAfile <- readChar(paste0(filedir, i, "/", LociandAlleles), file.info(paste0(filedir, i, "/", LociandAlleles))$size)
              numLociExt <- str_extract(string = LandAfile, pattern = paste0("from ", "[:digit:]{1,5}", " loci")) 
              numLociWorking <- gsub(x = numLociExt, pattern = "from ", replacement = "")
              numLociWorking <- as.numeric(gsub(x = numLociWorking, pattern = " loci", replacement = ""))
              
              #identify the simulation and repeat info
                S_ident <- gsub("_","",str_extract(pzfile,paste0("_S","[:digit:]{1}","R","[:digit:]{1}","_")))
                tempfile$sim <- substring(S_ident,1,2)
                tempfile$rep <- substring(S_ident,3,4)
                tempfile$nLoci <- numLociWorking
                
                tempfile <- tempfile[,-grep("IndivName",colnames(tempfile))] #delete IndivName
                
              #rename the columns
                colnames(tempfile) <- c("Indv","Pure1","Pure2","F1","F2","BC1","BC2","sim","rep","nLoci")
                tempfile=tempfile[,c("Indv","sim","rep","nLoci","Pure1","Pure2","F1","F2","BC1","BC2")]# reorder
                
              #common order
                if(sum(tempfile[1:samplesize,"Pure1"],na.rm=T)<sum(tempfile[1:samplesize,"Pure2"],na.rm=T)){
                  pure1 <- tempfile$Pure2;pure2 <- tempfile$Pure1
                  bc1 <- tempfile$BC2;bc2 <- tempfile$BC1
                  
                  tempfile$Pure1 <- pure1;tempfile$Pure2 <- pure2
                  tempfile$BC1 <- bc1;tempfile$BC2 <- bc2
                }
              
              #flag bad runs
              if(sum(tempfile$F2)>sum(sum(tempfile$Pure1),sum(tempfile$Pure2))){tempfile[,5:length(tempfile)]=NA}
              
              output <- rbind(output,tempfile)
          
          }#end of for loop


    ## average and SD the  replicate runs of each simulation in New Hybrids
      sim_data <- as.data.frame(output%>%group_by(nLoci,sim,Indv)%>%summarise(Pure1_sd=sd(Pure1),Pure1=mean(Pure1),
                                                                       Pure2_sd=sd(Pure2),Pure2=mean(Pure2),
                                                                       F1_sd=sd(F1),F1=mean(F1),
                                                                       F2_sd=sd(F2),F2=mean(F2),
                                                                       BC1_sd=sd(BC1),BC1=mean(BC1),
                                                                       BC2_sd=sd(BC2),BC2=mean(BC2))%>%ungroup())
    
    #pull out just the means of the replicates
      sim_means <- sim_data[,-grep("_sd",colnames(sim_data))]
    
    ## assign the classes to the data 
      nIndv <- nrow(sim_means)/6/length(unique(sim_means$sim))/length(unique(sim_means$nLoci)) #number of simulated individuals (assumes the same number for each class)
      sim_means$class=rep(rep(c("Pure1","Pure2","F1","F2","BC1","BC2"),each=nIndv),times=3*length(unique(sim_means$nLoci)))


#Compare the simulations using boxplots
    boxdata <- NULL
    for (i in unique(sim_means$nLoci))
    {
      for (j in unique(sim_means$class))
      {
      temp <- filter(sim_means,class==j,nLoci==i)
      tout <- data.frame(sim=temp$sim,Indv=temp$Indv,class=j,nLoci=i,value=temp[,j])
      boxdata <- rbind(boxdata,tout)
      }
    }
    
    boxdata$nLoci=factor(boxdata$nLoci)

    # Create pot
    p1=ggplot(boxdata,aes(x=nLoci,y=value,fill=sim))+geom_boxplot(alpha=0.8,outlier.size = 0)+theme_bw()+facet_wrap(~class,nrow=3,scales="free_y")+
      labs(y="Class probability",x="Number of loci")+scale_fill_manual(values=c("white","white","white"))+
      theme(strip.background = element_rect(fill="white"),legend.position="none")
   
    #save plot
    if(filetag!=""){ggsave(paste0(dir,"Figures/pdf/",filetag,"_AssignmentSuccess~simulation-nSNPs.pdf"),p1,height = 8,width = 10)}else
    {ggsave(paste0(dir,"Figures/pdf/AssignmentSuccess~simulation-nSNPs.pdf"),p1,height = 8,width = 10)}
    
    if(filetag!=""){ggsave(paste0(dir,"Figures/jpg/",filetag,"_AssignmentSuccess~simulation-nSNPs.jpg"),p1,height = 8,width = 10)}else
    {ggsave(paste0(dir,"Figures/jpg/AssignmentSuccess~simulation-nSNPs.jpg"),p1,height = 8,width = 10)}
    
    #Combined loci
    sim_means2 <- sim_means
    sim_means2$hybrid <- rowSums(sim_means2[,c("F1","F2","BC1","BC2")])
    sim_means2[which(sim_means2$class =="Pure1"),"hybrid"]=sim_means2[which(sim_means2$class =="Pure1"),"Pure1"] #add values of the Pure
    sim_means2[which(sim_means2$class =="Pure2"),"hybrid"]=sim_means2[which(sim_means2$class =="Pure2"),"Pure2"] #add values of the Pure
    
    sim_means2$hclass <- "Hybrid"
    sim_means2[which(sim_means$class=="Pure1"),"hclass"] <- "Pure1"
    sim_means2[which(sim_means$class=="Pure2"),"hclass"] <- "Pure2"
    
    sim_means2$hclass <- factor(sim_means2$hclass,levels=c("Pure1","Pure2","Hybrid"))
    
    h1 <- ggplot(sim_means2,aes(x=factor(nLoci),y=hybrid,fill=sim))+geom_boxplot(alpha=0.8,outlier.size = 0)+theme_bw()+facet_wrap(~hclass,nrow=3,scales="free_y")+
      labs(y="Class probability",x="Number of loci")+scale_fill_manual(values=c("white","white","white"))+
      theme(strip.background = element_rect(fill="white"),legend.position="none")+scale_y_continuous(limits=c(0.6,1))
   
    #Save plot
    if(filetag!=""){ggsave(paste0(dir,"Figures/pdf/",filetag,"_AssignmentSuccess~simulation-nSNPs_Hybrid.pdf"),h1,height = 8,width = 8)}else
    {ggsave(paste0(dir,"Figures/pdf/AssignmentSuccess~simulation-nSNPs_Hybrid.pdf"),h1,height = 8,width = 8)}
    
    if(filetag!=""){ggsave(paste0(dir,"Figures/jpg/",filetag,"_AssignmentSuccess~simulation-nSNPs_Hybrid.jpg"),h1,height = 8,width = 8)}else
    {ggsave(paste0(dir,"Figures/jpg/AssignmentSuccess~simulation-nSNPs_Hybrid.jpg"),h1,height = 8,width = 8)}
    

## Look at assignment success as a function of threshold probability
    num.sim <- length(which(sim_means$sim=="S1"))/6/length(unique(sim_means$nLoci))
    
    ProbOutput <- NULL
    for (s in unique(sim_means$nLoci)){
        lsub <- filter(sim_means,nLoci == s)
        for(i in unique(sim_means$sim)){
          tempsub <- filter(lsub,sim==i)
            for(q in 50:99/100){ # probability of 50 - 99%
             
              farm.p <- length(which(tempsub$Pure1[1:samplesize] > q))/num.sim
              wild.p <- length(which(tempsub$Pure2[(samplesize+1):(2*samplesize)] > q))/num.sim
              F1.p <- length(which(tempsub$F1[((2*samplesize)+1):(3*samplesize)] > q))/num.sim
              F2.p <- length(which(tempsub$F2[((3*samplesize)+1):(4*samplesize)] > q))/num.sim
              farmBC.p <- length(which(tempsub$BC1[((4*samplesize)+1):(5*samplesize)] > q))/num.sim
              wildBC.p <- length(which(tempsub$BC2[((5*samplesize)+1):(6*samplesize)] > q))/num.sim
              tempout <- data.frame(nLoci=s,sim=i,level=q,prob=c(farm.p, wild.p,F1.p,F2.p,farmBC.p,wildBC.p),
                                    class=c("Pure1","Pure2","F1","F2","BC1","BC2")) 
              ProbOutput <- rbind(ProbOutput,tempout)
              
            } # end q loop
          } # end i loop
        } # end s loop
    
    #combined hybrid probabilities
    ProbOutput2 <- NULL
    for (s in unique(sim_means$nLoci)){
      lsub <- filter(sim_means,nLoci == s)
      for(i in unique(sim_means$sim)){
        tempsub <- filter(lsub,sim==i)
        tempsub$phyb <- rowSums(tempsub[,c("F1","F2","BC1","BC2")])
        
        for(q in 50:99/100){ # probability of 50 - 99%
          farm.p <- length(which(tempsub$Pure1[1:samplesize] > q))/num.sim
          wild.p <- length(which(tempsub$Pure2[(samplesize+1):(2*samplesize)] > q))/num.sim
          Hybrid <- length(which(tempsub$phyb[((2*samplesize)+1):(6*samplesize)]>q))/(num.sim*4)
          tempout <- data.frame(nLoci=s,sim=i,level=q,prob=c(farm.p, wild.p,Hybrid),
                                class=c("Pure1","Pure2","Hybrid")) 
          ProbOutput2 <- rbind(ProbOutput2,tempout)
        } # end q loop
      } # end i loop
    } # end s loop

# get the mean and standard error for the estimates of assignment succes based on NH probabilty among simulations    
      FinalData <- data.frame(ProbOutput%>%group_by(nLoci,level,class)%>%summarise(mprob = mean(prob,na.rm=T), 
                                                                            sdprob = sd(prob,na.rm=T))%>%ungroup())
      FinalData$class <- factor(FinalData$class, levels=c("Pure1","Pure2","F1","F2","BC1","BC2")) # NH class
      
    # set plotting levels
      FinalData$group <- "Pure"
      FinalData[which(FinalData$class %in% c("BC1","BC2")),"group"] <- "Back-cross"
      FinalData[which(FinalData$class %in% c("F1","F2")),"group"] <- "Generational hybrids"
      
      FinalData$group <-  factor(FinalData$group,levels=c("Pure","Generational hybrids","Back-cross"))
      FinalData$class <- factor(FinalData$class,levels=c("Pure1","Pure2","F1","F2","BC1","BC2"))
    
    #plot the class groupings
      p3 <- ggplot(FinalData,aes(x=level,y=mprob,col=class))+geom_line(lwd=1.25)+theme_bw()+facet_grid(group~nLoci,scales="free_y")+
        theme(legend.position="bottom",strip.background = element_rect(fill="white",colour = "black"))+
        scale_color_brewer(palette = "Dark2")+
        labs(x="Probability threshold",y="Assignment success",col="Classification")
      
      #Save plot
      if(filetag!=""){ggsave(paste0(dir,"Figures/pdf/",filetag,"_AssinmentSuccess~level-class.pdf"),p3,height = 10,width = 8)} else 
      {ggsave(paste0(dir,"Figures/pdf/AssinmentSuccess~level-class.pdf"),p3,height = 10,width = 8)}
      
      if(filetag!=""){ggsave(paste0(dir,"Figures/jpg/",filetag,"_AssinmentSuccess~level-class.jpg"),p3,height = 10,width = 8)} else 
      {ggsave(paste0(dir,"Figures/jpg/AssinmentSuccess~level-class.jpg"),p3,height = 10,width = 8)}
      
      #ComboHybrids
      FinalData2 <- data.frame(ProbOutput2%>%group_by(nLoci,level,class)%>%summarise(mprob = mean(prob,na.rm=T), 
                                                                               sdprob = sd(prob,na.rm=T))%>%ungroup())
      FinalData2$class <- factor(FinalData2$class, levels=c("Pure1","Pure2","Hybrid")) # set plotting levels
      
      h3 <- ggplot(FinalData2)+
        geom_line(aes(x=level,y=mprob,col=class),lwd=1.25)+
        geom_line(aes(x=level,y=mprob+sdprob,col=class),lty=2)+
        geom_line(aes(x=level,y=mprob-sdprob,col=class),lty=2)+
        theme_bw()+
        facet_grid(~nLoci)+
        theme(legend.position="bottom",strip.background = element_rect(fill="white",colour = "black"))+
        scale_color_brewer(palette = "Dark2")+
        labs(x="Probability threshold",y="Assignment success ± sd",col="Classification");h3
      
      if(filetag!=""){ggsave(paste0(dir,"Figures/pdf/",filetag,"_AssinmentSuccess~level-class_Hybrid.pdf"),h3,height = 8,width = 10)} else 
      {ggsave(paste0(dir,"Figures/pdf/AssinmentSuccess~level-class_Hybrid.pdf"),h3,height = 8,width = 10)}
      
      if(filetag!=""){ggsave(paste0(dir,"Figures/jpg/",filetag,"_AssinmentSuccess~level-class_Hybrid.jpg"),h3,height = 8,width = 10)} else 
      {ggsave(paste0(dir,"Figures/jpg/AssinmentSuccess~level-class_Hybrid.jpg"),h3,height = 8,width = 10)}

    #plot if no threshold specified
      if(addThresh){
      p4 <- ggplot(data=FinalData)+geom_line(aes(x=level,y=mprob,col=factor(nLoci)),lwd=1.25)+
        geom_line(aes(x=level,y=mprob+sdprob,col=factor(nLoci)),lty=2)+
        geom_line(aes(x=level,y=mprob-sdprob,col=factor(nLoci)),lty=2)+
        theme_bw()+facet_wrap(~class,nrow=3,scales="free_y")+theme(strip.background = element_rect(fill="white",colour = "black"))+
        scale_color_brewer(palette = "Dark2")+
        labs(x="Probability threshold",y="Assignment success ± sd",col="# Loci")+geom_vline(xintercept = Thresholds, lty=2)
      }
      
      if(!addThresh){
        p4 <- ggplot(data=FinalData)+geom_line(aes(x=level,y=mprob,col=factor(nLoci)),lwd=1.25)+
          geom_line(aes(x=level,y=mprob+sdprob,col=factor(nLoci)),lty=2)+
          geom_line(aes(x=level,y=mprob-sdprob,col=factor(nLoci)),lty=2)+
          theme_bw()+facet_wrap(~class,nrow=3,scales="free_y")+theme(strip.background = element_rect(fill="white",colour = "black"))+
          scale_color_brewer(palette = "Dark2")+
          labs(x="Probability threshold",y="Assignment success ± sd",col="# Loci")
      }
      
      if(filetag!=""){ggsave(paste0(dir,"Figures/pdf/",filetag,"_AssignmentSuccess~level-error.pdf"),p4,height = 10,width = 8)} else 
      {ggsave(paste0(dir,"Figures/pdf/AssignmentSuccess~level-error.pdf"),p4,height = 10,width = 8)}
      
      if(filetag!=""){ggsave(paste0(dir,"Figures/jpg/",filetag,"_AssignmentSuccess~level-error.jpg"),p4,height = 10,width = 8)} else 
      {ggsave(paste0(dir,"Figures/jpg/AssignmentSuccess~level-error.jpg"),p4,height = 10,width = 8)}
      
      ## combined hybrids
      
      if(addThresh){
        h4 <- ggplot(data=FinalData2)+geom_line(aes(x=level,y=mprob,col=factor(nLoci)),lwd=1.25)+
          geom_line(aes(x=level,y=mprob+sdprob,col=factor(nLoci)),lty=2)+
          geom_line(aes(x=level,y=mprob-sdprob,col=factor(nLoci)),lty=2)+
          theme_bw()+facet_grid(~class)+theme(strip.background = element_rect(fill="white",colour = "black"))+
          scale_color_brewer(palette = "Dark2")+
          labs(x="Probability threshold",y="Assignment success ± sd",col="# Loci")+geom_vline(xintercept = Thresholds, lty=2)
      }
      
      if(!addThresh){
        h4 <- ggplot(data=FinalData2)+geom_line(aes(x=level,y=mprob,col=factor(nLoci)),lwd=1.25)+
          geom_line(aes(x=level,y=mprob+sdprob,col=factor(nLoci)),lty=2)+
          geom_line(aes(x=level,y=mprob-sdprob,col=factor(nLoci)),lty=2)+
          theme_bw()+facet_grid(~class)+theme(strip.background = element_rect(fill="white",colour = "black"))+
          scale_color_brewer(palette = "Dark2")+
          labs(x="Probability threshold",y="Assignment success ± sd",col="# Loci")
      }
      
      #Save plot
      if(filetag!=""){ggsave(paste0(dir,"Figures/pdf/",filetag,"_AssignmentSuccess~level-error_Hybrid.pdf"),h4,height = 10,width = 8)} else 
      {ggsave(paste0(dir,"Figures/pdf/AssignmentSuccess~level-error_Hybrid.pdf"),h4,height = 10,width = 8)}
      
      if(filetag!=""){ggsave(paste0(dir,"Figures/jpg/",filetag,"_AssignmentSuccess~level-error_Hybrid.jpg"),h4,height = 10,width = 8)} else 
      {ggsave(paste0(dir,"Figures/jpg/AssignmentSuccess~level-error_Hybrid.jpg"),h4,height = 10,width = 8)}
      
      ## mean plot
      
      #facet labels
      FinalData$threshold <- paste0(FinalData$level*100,"%")
      
        p5 <- ggplot(filter(FinalData,level %in% Thresholds),aes(x=factor(nLoci),y=mprob,col=class,group=class))+
        geom_point(size=2.5)+geom_path(lwd=0.9)+
        geom_errorbar(aes(ymin=mprob-sdprob,ymax=mprob+sdprob),width=0.1)+
        facet_grid(~threshold)+theme_bw()+
        labs(x="Number of loci",y="Assignment success ± sd",col="Classification",group="")+scale_color_brewer(palette = "Dark2")+
        theme(legend.position="bottom",strip.background = element_rect(fill="white",colour = "black"))
        
        #Save plot
        if(filetag!=""){ggsave(paste0(dir,"Figures/pdf/",filetag,"_AssignmentSuccess~z-loci.pdf"),p5,height = 8,width = 10)} else 
        {ggsave(paste0(dir,"Figures/pdf/AssignmentSuccess~~z-loci.pdf"),p5,height = 8,width = 10)}
        
        if(filetag!=""){ggsave(paste0(dir,"Figures/jpg/",filetag,"_AssignmentSuccess~z-loci.jpg"),p5,height = 8,width = 10)} else 
        {ggsave(paste0(dir,"Figures/jpg/AssignmentSuccess~~z-loci.jpg"),p5,height = 8,width = 10)}
        
        #Combined Hybrids
        FinalData2$threshold <- paste0(FinalData2$level*100,"%")
        
        h5 <- ggplot(filter(FinalData2,level %in% Thresholds),aes(x=factor(nLoci),y=mprob,col=class,group=class))+
          geom_point(size=2.5)+geom_path(lwd=0.9)+
          geom_errorbar(aes(ymin=mprob-sdprob,ymax=mprob+sdprob),width=0.1)+
          facet_grid(~threshold)+theme_bw()+
          labs(x="Number of loci",y="Assignment success ± sd",col="Classification",group="")+scale_color_brewer(palette = "Dark2")+
          theme(legend.position="bottom",strip.background = element_rect(fill="white",colour = "black"))
      
        #Save plot
        if(filetag!=""){ggsave(paste0(dir,"Figures/pdf/",filetag,"_AssignmentSuccess~z-loci_Hybrid.pdf"),h5,height = 8,width = 10)} else 
        {ggsave(paste0(dir,"Figures/pdf/AssignmentSuccess~~z-loci_Hybrid.pdf"),h5,height = 9,width = 10)}
        
        if(filetag!=""){ggsave(paste0(dir,"Figures/jpg/",filetag,"_AssignmentSuccess~z-loci_Hybrid.jpg"),h5,height = 8,width = 10)} else 
        {ggsave(paste0(dir,"Figures/jpg/AssignmentSuccess~~z-loci_Hybrid.jpg"),h5,height = 8,width = 10)}

        
        
    ## Misclassification 'type II' error
    
      classnames <- c("Pure1","Pure2","F1","F2","BC1","BC2")
        Start <- Sys.time()
        missout <- NULL
        for (s in unique(sim_means$nLoci)){
          lsub <- filter(sim_means,nLoci == s)
          for(i in unique(sim_means$sim)){
            tempsub <- filter(lsub,sim==i)
            for(q in 50:99/100){ # probability of 50 - 99%
              tempq <- tempsub
              tempq$missclass<- classnames[apply(tempq[,classnames],1,which.max)] # what is the class of the highest NH probability
              
              tempq$missval <- 999 #place holder
              for(w in 1:nrow(tempq))
              {
                tempq[w,"missval"] <- tempq[w,classnames[which.max(tempq[w,classnames])]]
              }
              
              temp1 <- tempq[which(tempq$class!=tempq$missclass & tempq$missval>=q),] #dataset with missclassifications
              
              
              dummydf=data.frame(Var1 = c("Pure1","Pure2","F1","F2","BC1","BC2"),dummy=NA) # dataframe for dummy values
              
              for (z in classnames){
                temp2 <- filter(temp1,class == z)
                if(nrow(temp2)>0){
                  temp3 <- as.data.frame(table(temp2$missclass)/samplesize) # percentage of samples miss classed to a given class of a given type of class (i)
                  
                  temp4 <-  merge(dummydf,temp3,by="Var1",all.x = TRUE)
                  
                  tempout <- data.frame(nLoci=s,sim=i,level=q,class=z,
                                        mclass_P1=temp4[which(temp4$Var1 == "Pure1"),"Freq"],
                                        mclass_P2=temp4[which(temp4$Var1 == "Pure2"),"Freq"],
                                        mclass_F1=temp4[which(temp4$Var1 == "F1"),"Freq"],
                                        mclass_F2=temp4[which(temp4$Var1 == "F2"),"Freq"],
                                        mclass_BC1=temp4[which(temp4$Var1 == "BC1"),"Freq"],
                                        mclass_BC2=temp4[which(temp4$Var1 == "BC2"),"Freq"])
                } else
                {tempout <- data.frame(nLoci=s,sim=i,level=q,class=z,
                                       mclass_P1=NA,
                                       mclass_P2=NA,
                                       mclass_F1=NA,
                                       mclass_F2=NA,
                                       mclass_BC1=NA,
                                       mclass_BC2=NA)}
                missout <- rbind(missout,tempout)
                
              } # end of z class loop
            } # end q loop
          } #end i loop
        } # end s loop
        
        #calcluate the means among simulations
    miss_mean <- data.frame(missout%>%group_by(nLoci,level,class)%>%summarise(mprobP1 = mean(mclass_P1,na.rm=T), 
                                                                   sdprobP1 = sd(mclass_P1,na.rm=T),
                                                                   mprobP2 = mean(mclass_P2,na.rm=T), 
                                                                   sdprobP2 = sd(mclass_P2,na.rm=T),
                                                                   mprobF1 = mean(mclass_F1,na.rm=T), 
                                                                   sdprobF1 = sd(mclass_F1,na.rm=T),
                                                                   mprobF2 = mean(mclass_F2,na.rm=T), 
                                                                   sdprobF2 = sd(mclass_F2,na.rm=T),
                                                                   mprobBC1 = mean(mclass_BC1,na.rm=T), 
                                                                   sdprobBC1 = sd(mclass_BC1,na.rm=T),
                                                                   mprobBC2 = mean(mclass_BC2,na.rm=T), 
                                                                   sdprobBC2 = sd(mclass_BC2,na.rm=T))%>%ungroup())
    
    miss_mean[is.na(miss_mean)]=NA #replace NaN's with NAs
    
    #merge with the other data
    FinalData2 <- merge(miss_mean,FinalData,by=c("nLoci","level","class"))
    
    PlotData <- melt(FinalData2[c("nLoci","level","class","mprobP1","mprobP2","mprobF1","mprobF2","mprobBC1","mprobBC2")],id.vars=c("nLoci","level","class"))
    PlotDatasd <- melt(FinalData2[c("nLoci","level","class","sdprobP1","sdprobP2","sdprobF1","sdprobF2","sdprobBC1","sdprobBC2")],id.vars=c("nLoci","level","class"))
    PlotData$sd <- PlotDatasd$value
    PlotData$variable <- as.character(PlotData$variable)
    
    PlotData[which(PlotData$variable == "mprobP1"),"variable"]="Pure1"
    PlotData[which(PlotData$variable == "mprobP2"),"variable"]="Pure2"
    PlotData[which(PlotData$variable == "mprobF1"),"variable"]="F1"
    PlotData[which(PlotData$variable == "mprobF2"),"variable"]="F2"
    PlotData[which(PlotData$variable == "mprobBC1"),"variable"]="BC1"
    PlotData[which(PlotData$variable == "mprobBC2"),"variable"]="BC2"
    
    PlotData$variable=factor(PlotData$variable,levels=c("Pure1","Pure2","F1","F2","BC1","BC2"))
      
#     p6 <- ggplot(filter(PlotData,nLoci==192))+
#       geom_line(aes(x=level,y=value,col=variable),lwd=1.25)+
#       geom_line(aes(x=level,y=value+sd,col=variable),lty=2)+
#       geom_line(aes(x=level,y=value-sd,col=variable),lty=2)+
#       facet_wrap(~class,nrow=3,scales="free_y")+
#       theme_bw()+scale_color_brewer(palette = "Dark2")+
#       theme(legend.position="bottom",strip.background = element_rect(fill="white",colour = "black"))+
#       labs(x="Probability threshold",y="Proportion misassigned ± sd",col="Classification");p6
    
    #Create the plots
    
    for (i in unique(PlotData$class)){
        temp.plot <- ggplot(filter(PlotData,class==i))+
        geom_line(aes(x=level,y=value,col=variable),lwd=1.25)+
        geom_line(aes(x=level,y=value+sd,col=variable),lty=2)+
        geom_line(aes(x=level,y=value-sd,col=variable),lty=2)+
        facet_grid(~nLoci,scales="free_y")+
        theme_bw()+scale_color_brewer(palette = "Dark2")+
        theme(legend.position="bottom",strip.background = element_rect(fill="white",colour = "black"))+
        labs(x="Probability threshold",y=paste0("Proportion ",i," misassigned ± sd"),col="Classification")
        
        if(filetag!=""){ggsave(paste0(dir,"Figures/pdf/",filetag,"_",i,"_MissAssignment~z-nloci.pdf"),temp.plot,height = 6,width = 8)} else 
        {ggsave(paste0(dir,paste0("Figures/pdf/",i,"_MissAssignment~z-nloci.pdf"),temp.plot,height = 6,width = 8))}
        
        if(filetag!=""){ggsave(paste0(dir,"Figures/jpg/",filetag,"_",i,"_MissAssignment~z-nloci.jpg"),temp.plot,height = 6,width = 8)} else 
        {ggsave(paste0(dir,paste0("Figures/jpg/",i,"_MissAssignment~z-nloci.jpg"),temp.plot,height = 6,width = 8))}
        
    }
    
    ## Create summary booklet
    
    if(filetag!=""){
      pdf(file = paste0(dir,"Figures/pdf/",filetag,"_OutputBooklet.pdf"))
      print(p1);print(h1)
      print(p3);print(h3)
      print(p4);print(h4)
      print(p5);print(h5)
      for (i in unique(PlotData$class)){
        temp.plot <- ggplot(filter(PlotData,class==i))+
          geom_line(aes(x=level,y=value,col=variable),lwd=1.25)+
          geom_line(aes(x=level,y=value+sd,col=variable),lty=2)+
          geom_line(aes(x=level,y=value-sd,col=variable),lty=2)+
          facet_grid(~nLoci,scales="free_y")+
          theme_bw()+scale_color_brewer(palette = "Dark2")+
          theme(legend.position="bottom",strip.background = element_rect(fill="white",colour = "black"))+
          labs(x="Probability threshold",y=paste0("Proportion ",i," misassigned ± sd"),col="Classification")
          print(temp.plot)
          }
      dev.off()
    } else {
      pdf(file = paste0(dir,"Figures/pdf/",filetag,"_OutputBooklet.pdf"))
      print(p1);print(h1)
      print(p3);print(h3)
      print(p4);print(h4)
      print(p5);print(h5)
      for (i in unique(PlotData$class)){
        temp.plot <- ggplot(filter(PlotData,class==i))+
          geom_line(aes(x=level,y=value,col=variable),lwd=1.25)+
          geom_line(aes(x=level,y=value+sd,col=variable),lty=2)+
          geom_line(aes(x=level,y=value-sd,col=variable),lty=2)+
          facet_grid(~nLoci,scales="free_y")+
          theme_bw()+scale_color_brewer(palette = "Dark2")+
          theme(legend.position="bottom",strip.background = element_rect(fill="white",colour = "black"))+
          labs(x="Probability threshold",y=paste0("Proportion ",i," misassigned ± sd"),col="Classification")
        print(temp.plot)
      }
      dev.off()
    }

} #end function
