JobList <- function(nPops,k,nsim=500000,burnin=100000,reps=3,dir=NULL){
  
  #This function will create a jobslist which can be used by paralell structure
  
#   npops is the number of populations in your data
#   k is a number or vector of potentail clustering that STRUCTURE will look for
#     if k = 5 then it will look for 1 through 5 clusters
#     if k = c(1,3,4) it will look for 1, 3 and 4 clusters
#   nsim is the number of simulations (defaults to 500000)
#   burnin is the burnin desired for MCMC (defaults to 100000)
#   reps is the number of replicate clusterings for each level of k (defaults to 3)
#   dir is the directory you where you want the job list to be created. By default it will just return
#     to R environment
  
  
  options(scipen = 999) # this will stop R from displaying the burning and nsim values in scientic notation
  nsim=as.numeric(nsim)
  
  #Make a population vector
  if(length(k)==1)
  {
    kvec <- rep(1:k,each=reps)
  }
  
  if(length(k)>1)
  {
    kvec <- rep(k,each=reps)
  }
  
  #Population data
  temp1 <- rep(paste(as.character(1:nPops), collapse=","))
  PopVec <- as.vector(do.call("rbind", replicate(length(kvec),temp1, simplify = FALSE)))
  jobs <- paste("T",1:length(kvec),sep="")
  burn <- rep(burnin,length(kvec))
  sim <- rep(nsim,length(kvec))
  
  temp2 <- as.data.frame(cbind(jobs,PopVec,kvec,burn,sim))
  
  output <- do.call(paste,c(temp2, sep=" "))
  
  if(length(dir)>0)
  {
  write.table(output,paste(dir,nPops,"_",max(k),"_",nsim,"_",burnin,"_",reps,".txt",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE)
  }
  
  if(length(dir)<=0)
  {
  return(output)
  }
  
  
}