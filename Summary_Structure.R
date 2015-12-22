structureFun <- function(dir,n=2,PopOrder=NULL,South,nameManual=NULL){
  
  ## this function is a preliminary function which will create "clumpack" like results and summaries of the q outputs for plotting. 

  # dir is a directory where the output (_f files) of the structure can be found and results will be placed
  # n is 2 is the number of populations that should be found optimized at the moment for 2
  # order of the populatons (these are the characters in each sample ID) e.g. GOM_01 and NSC_04 would be "GOM" and "NSC"
  # south is the southern most population which will be coloured Red e.g. "GOM"
  # if the sample lables have more than one set of letters (i.e. bon_004a) but have the first three as letters then this can 
  # be used to manually grab (n) characters off of each id. This assumes that each population has the same number of letters
  # default is to not use this. 
  
  #Libraries
  require(stringr)
  require(tidyr)
  require(dplyr)
  require(reshape2)
  require(ggplot2)
  require(grid)
  require(scales)
  
  # this function will split the data apart for each row of output from Structure into the individual ID and the q values for each population
  # of a clusters
  splitFun <- function(x){
    temp <- unlist(strsplit(x," "))[which(unlist(strsplit(x," "))!="")]
    temp[c(2,6:length(temp))]
  }
  
  numFun <- function(x){as.numeric(as.character(x))} # to switch factor columns to numeric
  
  #read all the structure results files "_f"
  Files <- list.files(dir)[grep("_f",list.files(dir))]
  
  #Read each file and extract data where the number of clusters == n
  df=NULL
  for (i in Files)
    {
    
          #Read in data 
          FData <- read.table(paste0(dir,i),header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)
          
          #What is the number of populations that Structure clustered
          Pops <- FData[grep("populations assumed",FData$V1),"V1"] # find the number of populations structured
          Pops <- as.numeric(gsub("[^\\d]+", "", Pops, perl=TRUE))

          if(n==Pops)
            {
            top <- which(FData$V1=="Inferred ancestry of individuals:")+2 # index vector of the beginning of the data
            bottom <- which(FData$V1=="Estimated Allele Frequencies in each cluster")-1 # index vector of the end of the data
            FData <- as.character(FData[top:bottom,]) #Use indexes to just take the q data
            FData <- do.call(rbind,lapply(FData,splitFun)) #seperate out just the data for 
            df=rbind(df,FData)
            }# end of if (n) loop

    } #end of loop

  #Convert ouput into IDs and numbers
  df=as.data.frame(df)
  df[,1]=as.character(df[,1])
  df[,2:length(df)]=sapply(df[,2:length(df)],numFun)
  colnames(df)=c("ID",paste0("q",1:(length(df)-1)))
  df$run=rep(1:length(Files),each=nrow(FData))
  
  # based on each ID make sure the structure runs assign the individual (ID) to the same cluster order. Structure will switch
  #back and forth for some reason
  
  IDs=unique(df$ID)#vector of sample IDs
  
  df2=df # new dataframe which will be manipulatyed
  
  for (i in IDs){
    
    ind=which(df2$ID==i) # row index corresponding to the ID 'i'
    temp=df2[ind,] # subsetted data corresponding to the ID 'i'
    
    #based on the forth population which cluster is bigger
    if(temp[nrow(temp),"q1"]>temp[nrow(temp),"q2"]){offset="larger"}else{offset="smaller"} 
    
    
    # if the first cluster is the 'largest" in the last structure run then make sure all other runs fit
    if(offset == "larger"){
        for (i in 1:(nrow(temp)-1)){
          if(temp[i,"q1"]<=0.5){temp[i,c("q1","q2")]=1-temp[i,c("q1","q2")]}
        }
    }
    
    # if the first cluster is the 'smallest' in the last structure run then make sure all other runs fit
    if(offset == "smaller"){
      for (i in 1:(nrow(temp)-1)){
        if(temp[i,"q1"]>=0.5){temp[i,c("q1","q2")]=1-temp[i,c("q1","q2")]}
      }
    }
    
    #replace the structure results with the alligned resutls
   df2[ind,]=temp 
   rm(ind,temp) #clean workspace
  }
  

  #Get the average q assignment for each individual
  df3 <- as.data.frame(df2%>%group_by(ID)%>%summarise(Q1=mean(q1,na.rm=T),Q2=mean(q2,na.rm=T),SD=sd(q1,na.rm=T))%>%ungroup())
  
  #find the Populations
  if(length(nameManual)>=1){df3$Pop <- substring(df3$ID,1,3)}else{df3$Pop <- str_extract(df3$ID, "[A-Za-z]+" )}
  
  
  #standardize the q values so that the north is q1 and the south is q2
  df3a <- as.data.frame(df3%>%group_by(Pop)%>%summarise(mn=mean(Q1,na.rm=T))%>%ungroup()) #means per pop
  
  if(df3a[which(df3a$Pop==South)[1],"mn"]>=0.5){df3[,c("Q1","Q2")]=1-df3[,c("Q1","Q2")]} #switch if required
  
  #set the plotting order.
  if(length(PopOrder)>1){df3$Pop <- factor(df3$Pop,levels=PopOrder)}
  
  df4=melt(df3[,c("ID","Q1","Q2","Pop")],id.vars=c("ID","Pop")) # plot data
  
  df5=as.data.frame(df3%>%group_by(Pop)%>%summarise(q1=mean(Q1,na.rm=T),q2=mean(Q2,na.rm=T),SD=sd(Q1,na.rm=T))%>%ungroup()) #mean population data
  
  if(df5[which(df5$Pop==South),"q1"]>0.5){COLS=c("red","blue")}else{COLS=c("blue","red")} #set the red blue colouring based on the Structure clusters
  
  
  p1=ggplot(df4,aes(x=factor(ID),y=value,fill=variable,group=Pop)) + 
    geom_bar(position = "stack",stat = "identity",width=1)+
    theme_bw()+scale_y_continuous(expand = c(0,0),labels = percent)+
    facet_grid(~Pop,scales="free",space="free_x")+
    theme(axis.text.x = element_blank(),
          #axis.ticks=element_blank(),
          panel.margin = unit(0.1, "lines"),
          legend.position="none",
          axis.ticks.x=element_blank())+
    scale_fill_manual(values = COLS)+
    labs(x="Population",y="Proportion",fill="");p1

  #Return the output
  ggsave(paste0(dir,"StructurePlot.png"),p1,width=10,height=5)
  
  write.table(df3,paste0(dir,"Structure_q_values_IDs.csv"),sep=",",row.names=F,col.names=T)
  write.table(df5,paste0(dir,"Structure_q_values_Pops.csv"),sep=",",row.names=F,col.names=T)

    #End of function
}

