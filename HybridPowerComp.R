Hybridpower_comparison <-function(dir,filetag="",Thresholds=c(0.5,0.6,0.7,0.8,0.9),addThresh=FALSE){
  
## this function will estimate the power of assignment success and associated variability from New Hybrids for 6 possible classes
#  Pure1, Pure2, F1, F2, BC1, BC2. The code is based on different simulations of individuals and repeat runs through New Hybrids (NH)
# for each simulation. The output is a series of diagnostic plots and data summaries. 
  
## dir - path directory which holds the output from different runs through New Hybrids (e.g. 3 simulations with 3 replicate runs each through NH)
#        note that this directory should only hold the output folders.

## filetag - this is a name tag which will be added to the plots

## Thresholds - this is a vector of thesholds (default c(0.6,0.7,0.8,0.9)) which will to compare among different thesholds to compare among different numbers of SNPs

## addThresh - a logical vector (default: FALSE) which specifies whether the threshold values should be added to the summary plot.               


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
                if(sum(tempfile[1:200,"Pure1"],na.rm=T)<sum(tempfile[1:200,"Pure2"],na.rm=T)){
                  pure1 <- tempfile$Pure2;pure2 <- tempfile$Pure1
                  bc1 <- tempfile$BC2;bc2 <- tempfile$BC1
                  
                  tempfile$Pure1 <- pure1;tempfile$Pure2 <- pure2
                  tempfile$BC1 <- bc1;tempfile$BC2 <- bc2
                }
              
              #flag bad runs
              if(sum(tempfile$F2)>sum(sum(tempfile$Pure1),sum(tempfile$Pure2))){tempfile[,4:length(tempfile)]=NA}
              
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
      labs(y="Class probability",x="Number of SNPs")+scale_fill_manual(values=c("white","white","white"))+
      theme(strip.background = element_rect(fill="white"),legend.position="none")
   
    #save plot
  if(filetag!=""){ggsave(paste0(dir,"Figures/",filetag,"_AssignmentSuccess~simulation-nSNPs.pdf"),p1,height = 8,width = 10)}else
  {ggsave(paste0(dir,"Figures/AssignmentSuccess~simulation-nSNPs.pdf"),p1,height = 8,width = 10)}

## Look at assignment success as a function of threshold probability
    num.sim <- length(which(sim_means$sim=="S1"))/6/length(unique(sim_means$nLoci))
    
    ProbOutput <- NULL
    for (s in unique(sim_means$nLoci)){
        lsub <- filter(sim_means,nLoci == s)
        for(i in unique(sim_means$sim)){
          tempsub <- filter(lsub,sim==i)
            for(q in 50:99/100){ # probability of 50 - 99%
             
              farm.p <- length(which(tempsub$Pure1[1:200] > q))/num.sim
              wild.p <- length(which(tempsub$Pure2[201:400] > q))/num.sim
              F1.p <- length(which(tempsub$F1[401:600] > q))/num.sim
              F2.p <- length(which(tempsub$F2[601:800] > q))/num.sim
              farmBC.p <- length(which(tempsub$BC1[801:1000] > q))/num.sim
              wildBC.p <- length(which(tempsub$BC2[1001:1200] > q))/num.sim
              tempout <- data.frame(nLoci=s,sim=i,level=q,prob=c(farm.p, wild.p,F1.p,F2.p,farmBC.p,wildBC.p),
                                    class=c("Pure1","Pure2","F1","F2","BC1","BC2")) 
              ProbOutput <- rbind(ProbOutput,tempout)
              
            } #end q loop
          } #end i loop
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
      
      ## mean plot
      
      #facet labels
      FinalData$threshold <- paste0(FinalData$level*100,"%")
      
        p5 <- ggplot(filter(FinalData,level %in% Thresholds),aes(x=factor(nLoci),y=mprob,col=class,group=class))+
        geom_point(size=2.5)+geom_path(lwd=0.9)+
        geom_errorbar(aes(ymin=mprob-sdprob,ymax=mprob+sdprob),width=0.1)+
        facet_grid(~threshold)+theme_bw()+
        labs(x="Number of loci",y="Assignment success ± sd",col="Classification",group="")+scale_color_brewer(palette = "Dark2")+
        theme(legend.position="bottom",strip.background = element_rect(fill="white",colour = "black"))
      
      
#Save plots

    if(filetag!=""){ggsave(paste0(dir,"Figures/",filetag,"_AssinmentSuccess~level-class.pdf"),p3,height = 10,width = 8)} else 
    {ggsave(paste0(dir,"Figures/AssinmentSuccess~level-class.pdf"),p3,height = 10,width = 8)}
    
    if(filetag!=""){ggsave(paste0(dir,"Figures/",filetag,"_AssignmentSuccess~level-error.pdf"),p4,height = 10,width = 8)} else 
    {ggsave(paste0(dir,"Figures/AssignmentSuccess~level-error.pdf"),p4,height = 10,width = 8)}
        
    if(filetag!=""){ggsave(paste0(dir,"Figures/",filetag,"_AssignmentSuccess~z-loci.pdf"),p5,height = 10,width = 8)} else 
    {ggsave(paste0(dir,"Figures/AssignmentSuccess~~z-loci.pdf"),p5,height = 10,width = 8)}

} #end function
