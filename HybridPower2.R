Hybridpower <-function(dir,filetag="",Threshold=NULL) {
  
## this function will estimate the power of assignment success and associated variability from New Hybrids for 6 possible classes
#  Pure1, Pure2, F1, F2, BC1, BC2. The code is based on different simulations of individuals and repeat runs through New Hybrids (NH)
# for each simulation. The output is a series of diagnostic plots and data summaries. 
  
## dir - path directory which holds the output from different runs through New Hybrids (e.g. 3 simulations with 3 replicate runs each through NH)
#        note that this directory should only hold the output folders.

## filetag - this is a name tag which will be added to the plots

## threshhold - this is a theshold which will be added to the plots showing the assignment success for different levels of probability of a given 
##              class estimated by NH. Default is (NULL) so if nothing specified it will not add this to the output plots (success ~ threshold by class)


  
  #Check to make sure the packages required are there and if not install them
  packages <- c("dplyr", "tidyr", "stringr","ggplot2","reshape2")
  if (length(setdiff(packages, rownames(installed.packages()))) > 0) { 
    install.packages(setdiff(packages, rownames(installed.packages())))
  } 
  
  #load each library
  require(dplyr)
  require(ggplot2)
  require(tidyr)
  require(stringr)
  require(reshape2)


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
              #identify the simulation and repeat info
                S_ident <- gsub("_","",str_extract(pzfile,paste0("_S","[:digit:]{1}","R","[:digit:]{1}","_")))
                tempfile$sim <- substring(S_ident,1,2)
                tempfile$rep <- substring(S_ident,3,4)
                
                tempfile <- tempfile[,-grep("IndivName",colnames(tempfile))] #delete IndivName
                
              #rename the columns
                colnames(tempfile) <- c("Indv","Pure1","Pure2","F1","F2","BC1","BC2","sim","rep")
                tempfile=tempfile[,c("Indv","sim","rep","Pure1","Pure2","F1","F2","BC1","BC2")]# reorder
                
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
      sim_data <- as.data.frame(output%>%group_by(sim,Indv)%>%summarise(Pure1_sd=sd(Pure1),Pure1=mean(Pure1),
                                                                       Pure2_sd=sd(Pure2),Pure2=mean(Pure2),
                                                                       F1_sd=sd(F1),F1=mean(F1),
                                                                       F2_sd=sd(F2),F2=mean(F2),
                                                                       BC1_sd=sd(BC1),BC1=mean(BC1),
                                                                       BC2_sd=sd(BC2),BC2=mean(BC2))%>%ungroup())
    
    #pull out just the means of the replicates
      sim_means <- sim_data[,-grep("_sd",colnames(sim_data))]
    
    ## assign the classes to the data 
      nIndv <- nrow(sim_means)/6/length(unique(sim_means$sim)) #number of simulated individuals (assumes the same number for each class)
      sim_means$class=rep(rep(c("Pure1","Pure2","F1","F2","BC1","BC2"),each=nIndv),n=3)


#Compare the simulations using boxplots
    boxdata <- NULL
    for (i in unique(sim_means$class))
    {
      temp <- filter(sim_means,class==i)
      tout <- data.frame(sim=temp$sim,Indv=temp$Indv,class=i,value=temp[,i])
      boxdata <- rbind(boxdata,tout)
    }

    # Create pot
    p1=ggplot(boxdata,aes(x=sim,y=value))+geom_jitter(size=0.1)+geom_boxplot(alpha=0.8,outlier.size = 0)+theme_bw()+facet_wrap(~class,nrow=3,scales="free_y")+
    labs(y="Class probability",x="Simulation")+theme(strip.background = element_rect(fill="white"))
    
    #save plot
  if(filetag!=""){ggsave(paste0(dir,"Figures/",filetag,"_AssignmentSuccess~simulation.pdf"),p1,height = 8,width = 10)}else
  {ggsave(paste0(dir,"Figures/AssignmentSuccess~simulation.pdf"),p1,height = 8,width = 10)}

## create the New Hybrids plot
  
      NH_melt <- melt(data = sim_means[,-grep("class",colnames(sim_means))], id.vars = c("Indv","sim")) ## melt the data to allow the data to be stacked by indivudal
      colnames(NH_melt) <- c("Indv","sim", "PopProb", "CumProb") ## rename so that its prettier
      
      NH_melt$sim <- gsub("S","Simulation ",NH_melt$sim) # fix the simulation labels
      
      ## Plot colours
      col.vec <- c("red", "blue", "grey", "green", "black", "yellow", "brown")
      
      
      ## make a nice pretty little plot - and save that bad boy
      #f.name <- paste0(res.folder, res.name, ".consensus.pdf")
      #pdf(file = f.name)
      
      p2 <- ggplot(NH_melt, aes(x = Indv, y=CumProb, fill=PopProb))+
      geom_bar(stat="identity", position = "stack") + 
        scale_fill_manual(values=col.vec)+
        labs(y="Cumulative Probability",x="Individual",fill="")+
        scale_y_continuous(limits = c(0, 1), expand=c(0, 0)) + 
        scale_x_continuous(expand=c(0,0)) +
        facet_grid(sim~.)+
        theme(axis.text.x = element_text(colour = "black"), axis.title.x = element_blank(), 
              axis.text.y = element_text(colour = "black"), axis.title.y = element_text(size = 15, colour = "black"), 
              strip.background = element_rect(fill="white",colour = "black"),panel.margin=unit(1, "lines"))
      
      #save the figure
      if(filetag!=""){ggsave(paste0(dir,"Figures/",filetag,"_NewHybridsPlot~simulation.pdf"),p2,height = 8,width = 10)} else 
      {ggsave(paste0(dir,"Figures/NewHybridsPlot~simulation.pdf"),p2,height = 8,width = 10)}


## Look at assignment success as a function of threshold probability
    num.sim <- length(which(sim_means$sim=="S1"))/6
    
    ProbOutput <- NULL
    for(i in unique(sim_means$sim)){
      tempsub <- filter(sim_means,sim==i)
        for(q in 50:99/100){ # probability of 50 - 99%
         
          farm.p <- length(which(tempsub$Pure1[1:200] > q))/num.sim
          wild.p <- length(which(tempsub$Pure2[201:400] > q))/num.sim
          F1.p <- length(which(tempsub$F1[401:600] > q))/num.sim
          F2.p <- length(which(tempsub$F2[601:800] > q))/num.sim
          farmBC.p <- length(which(tempsub$BC1[801:1000] > q))/num.sim
          wildBC.p <- length(which(tempsub$BC2[1001:1200] > q))/num.sim
          tempout <- data.frame(sim=i,level=q,prob=c(farm.p, wild.p,F1.p,F2.p,farmBC.p,wildBC.p),
                                class=c("Pure1","Pure2","F1","F2","BC1","BC2")) 
          ProbOutput <- rbind(ProbOutput,tempout)
          
        }
      
    }

# get the mean and standard error for the estimates of assignment succes based on NH probabilty among simulations    
      FinalData <- data.frame(ProbOutput%>%group_by(level,class)%>%summarise(mprob = mean(prob,na.rm=T), 
                                                                            sdprob = sd(prob,na.rm=T))%>%ungroup())
      FinalData$class <- factor(FinalData$class, levels=c("Pure1","Pure2","F1","F2","BC1","BC2")) # NH class
      
    # set plotting levels
      FinalData$group <- "Pure"
      FinalData[which(FinalData$class %in% c("BC1","BC2")),"group"] <- "Back-cross"
      FinalData[which(FinalData$class %in% c("F1","F2")),"group"] <- "Generational hybrids"
      
      FinalData$group <-  factor(FinalData$group,levels=c("Pure","Generational hybrids","Back-cross"))
      FinalData$class <- factor(FinalData$class,levels=c("Pure1","Pure2","F1","F2","BC1","BC2"))
    
    #plot the class groupings
      p3 <- ggplot(FinalData,aes(x=level,y=mprob,col=class))+geom_line(lwd=1.25)+theme_bw()+facet_grid(group~.,scales="free_y")+
        theme(legend.position="bottom",strip.background = element_rect(fill="white",colour = "black"))+labs(x="Probability threshold",y="Assignment success",col="Classification");p3

    #plot if no threshold specified
      if(length(Threshold)==0){
      p4 <- ggplot(data=FinalData)+geom_line(aes(x=level,y=mprob),lwd=1.25)+
        geom_line(aes(x=level,y=mprob+sdprob),lty=2)+geom_line(aes(x=level,y=mprob-sdprob),lty=2)+
        theme_bw()+facet_wrap(~class,nrow=3,scales="free_y")+theme(strip.background = element_rect(fill="white",colour = "black"))+
        labs(x="Probability threshold",y="Assignment success ± sd");p4
      }
      
    #plot if threshold specified 
       if(length(Threshold)!=0){
        #calcualte the lower assignment success for each class
        lowerlim <- data.frame(FinalData%>%group_by(class)%>%
                                 summarise(lim=min(mprob))%>%ungroup())
        
        threshlimits <- filter(FinalData,level==Threshold)
        threshlimits$lower <- 0.5
        threshlimits$threshold <- Threshold
        threshlimits$x <- (threshlimits$threshold-0.5)/2
        threshlimits <- merge(threshlimits,lowerlim,by = "class")
        threshlimits$lim <- threshlimits$lim-threshlimits$sdprob
        
        p4 <- ggplot(data=FinalData)+
          geom_line(aes(x=level,y=mprob),lwd=1.25)+
          geom_line(aes(x=level,y=mprob+sdprob),lty=2)+
          geom_line(aes(x=level,y=mprob-sdprob),lty=2)+
          theme_bw()+facet_wrap(~class,nrow=3,scales="free_y")+
          geom_segment(data=threshlimits,aes(x=lower,xend=level,y=mprob,yend=mprob),col="black")+
          geom_segment(data=threshlimits,aes(x=level,xend=level,y=lim,yend=mprob),col="black")+
          scale_y_continuous(expand=c(0.03,0))+scale_x_continuous(expand=c(0.04,0))+
          labs(x="Probability threshold",y="Assignment success ± sd")+
          theme(strip.background = element_rect(fill="white",colour = "black"))+
          geom_text(data=threshlimits,aes(x=0.55,y = mprob-(mprob-lim)*0.1,label=round(mprob,2)));p4
      }

#Save plots

    if(filetag!=""){ggsave(paste0(dir,"Figures/",filetag,"_AssinmentSuccess~level-class.pdf"),p3,height = 10,width = 8)} else 
    {ggsave(paste0(dir,"Figures/AssinmentSuccess~level-class.pdf"),p3,height = 10,width = 8)}
    
    if(filetag!=""){ggsave(paste0(dir,"Figures/",filetag,"_AssignmentSuccess~level-error.pdf"),p4,height = 10,width = 8)} else 
    {ggsave(paste0(dir,"Figures/AssignmentSuccess~level-error.pdf"),p4,height = 10,width = 8)}

} #end function
