#load HZAR
library("hzar", lib.loc="~/R/win-library/3.2")

#setWD
setwd("~/School/MASTER_Masters/Genetic_Data/June_2014_finals/Clines/HZAR/")

#import data
outlier_cline_data <- read.csv("finalfreqs_pop_HWE_outlier_HZAR_formatted_2.csv")



###copy SI 1 from HZAR paper###
## Save all plots in a series of png files
png(width=900, height=900, res=200, filename="outlierPlot%03d.png",pointsize=8)
dev.off()

## A typical chain length.  This value is the default setting in the package.
chainLength=1e5;

## Make each model run off a separate seed
mainSeed=
  list(A=c(596,528,124,978,544,99),
       B=c(528,124,978,544,99,596),
       C=c(124,978,544,99,596,528),
       D=c(902,168,123,843,761,378))

##if you have doMC, use foreach in parallel mode to speed up computation
if(require(doMC)){
  registerdoMC()
} else {
  #use foreach in sequential mode
  registerDoSEQ();
}




#####Prepping Loci --> using 1DDist#####

####make lists for stuff!

## Blank out space in memory to hold molecular analysis
if(length(apropos("^outlier$",ignore.case=FALSE)) == 0 ||
   !is.list(outlier) ) outlier <- list()

#okay, so the list is called "outlier" - need to make a space on this list for every locus (only using the allele I want)
#make list of only the alleles I want to use
#see code_lines_for_HZAR.xlsx for all loci

#EXAMPLES#
outlier$"11_6A" <- list()
# Space to hold the observed data
outlier$"11_6A"$obs <- list()
# Space to hold the models to fit
outlier$"11_6A"$models <- list()
# Space to hold the compiled fit requests
outlier$"11_6A"$fitRs <- list()
# Space to hold the output data chains
outlier$"11_6A"$runs <- list()
# Space to hold the analysed data
outlier$"11_6A"$analysis <- list()

#EXAMPLES#
#Locus 11_6, allele A from outlier data
outlier$"11_6A"$obs <- hzar.doMolecularData1DPops(outlier_cline_data$"X1DDist",
                                                  outlier_cline_data$"X11_6.A",
                                                  outlier_cline_data$"Alleles_11_6")
#EXAMPLES#
#look at plot
hzar.plot.obsData(outlier$"11_6A"$obs)

###DO THIS FOR ALL LOCI  --> code lines all in code_lines_for_HZAR.xlsx###



#####MODELS - code is just for 1 locus, if this works will redo for all 112#####

#set the model i want to look at - check param options using ?hzar.makeCline1DFreq
#EXAMPLE - testing 4 models#
outlier$"11_6A"$models$model1 <- hzar.makeCline1DFreq(data=outlier$"11_6A"$obs, 
                                               scaling="fixed", tails="none")
outlier$"11_6A"$models$model2 <- hzar.makeCline1DFreq(data=outlier$"11_6A"$obs, 
                                                      scaling="fixed", tails="both")
outlier$"11_6A"$models$model3 <- hzar.makeCline1DFreq(data=outlier$"11_6A"$obs, 
                                                      scaling="free", tails="none")
outlier$"11_6A"$models$model4 <- hzar.makeCline1DFreq(data=outlier$"11_6A"$obs, 
                                                      scaling="free", tails="both")

#check params
print(outlier$"11_6A"$models)

#modify all models to focus on observed region 
#data collected between 0 and 4400km
outlier$"11_6A"$models <- sapply(outlier$"11_6A"$models, hzar.model.addBoxReq,
                                 -30, 4500, simplify=FALSE)

#check updates params
print(outlier$"11_6A"$models)

###compile models to prepare for fitting --> creates hzar.fitRequest from each 
#clineMetalModel object

outlier$"11_6A"$fitRs$init <- sapply(outlier$"11_6A"$models, 
                                     hzar.first.fitRequest.old.ML,
                                     obsData = outlier$"11_6A"$obs,
                                     verbose=FALSE,
                                     simplify=FALSE)

#update settings for the fitter using chainLength and mainSeed created before
outlier$"11_6A"$fitRs$init$model1$mcmcParam$chainLength <- chainLength
outlier$"11_6A"$fitRs$init$model1$mcmcParam$burnin <- chainLength %/% 10
#outlier$"11_6A"$fitRs$init$model1$mcmcParam$seed[[1]] <- mainSeed$A

outlier$"11_6A"$fitRs$init$model2$mcmcParam$chainLength <- chainLength
outlier$"11_6A"$fitRs$init$model2$mcmcParam$burnin <- chainLength %/% 10
#outlier$"11_6A"$fitRs$init$model2$mcmcParam$seed[[1]] <- mainSeed$B

outlier$"11_6A"$fitRs$init$model3$mcmcParam$chainLength <- chainLength
outlier$"11_6A"$fitRs$init$model3$mcmcParam$burnin <- chainLength %/% 10
#outlier$"11_6A"$fitRs$init$model3$mcmcParam$seed[[1]] <- mainSeed$C

outlier$"11_6A"$fitRs$init$model4$mcmcParam$chainLength <- chainLength
outlier$"11_6A"$fitRs$init$model4$mcmcParam$burnin <- chainLength %/% 10
#outlier$"11_6A"$fitRs$init$model4$mcmcParam$seed[[1]] <- mainSeed$D

#check fit request settings 
print(outlier$"11_6A"$fitRs$init)

#replicate each fit request 3 times, keeping the original seeds 
#while switching to a new seed channel
#12 total fit requests - 4 models, 3 times each
outlier$"11_6A"$fitRs$chains <- hzar.multiFitRequest(outlier$"11_6A"$fitRs$init,
                                                     each=3,
                                                    baseSeed=NULL)


##have 36 fit requests - models 4, each with 3 chain, 
#running the chain 3 times - 36 total runs? - THIS WILL TAKE A WHILE
outlier$"11_6A"$runs$doSeq <- lapply(outlier$"11_6A"$fitRs$chains,
                                     hzar.chain.doSeq, 
                                     count = 3)


#did model1 converge?
summary(do.call(mcmc.list, lapply(outlier$"11_6A"$runs$doSeq[1:3],
                                  function(x) hzar.mcmc.bindLL(x[[3]]))))
#YES

#did model2 converge?
summary(do.call(mcmc.list, lapply(outlier$"11_6A"$runs$doSeq[4:6],
                                  function(x) hzar.mcmc.bindLL(x[[3]]))))
#YES

#did model3 converge?
summary(do.call(mcmc.list, lapply(outlier$"11_6A"$runs$doSeq[7:9],
                                  function(x) hzar.mcmc.bindLL(x[[3]]))))
#YES

#did model4 converge?
summary(do.call(mcmc.list, lapply(outlier$"11_6A"$runs$doSeq[10:12],
                                  function(x) hzar.mcmc.bindLL(x[[3]]))))
#YES


#####ANALYSIS#####
#start aggregation of data for analysis

#create a data group for the null model
outlier$"11_6A"$analysis$initDGs <- list(nullModel=hzar.dataGroup.null(outlier$"11_6A"$obs))

#create a model data group for each model from the initial runs
outlier$"11_6A"$analysis$initDGs$model1 <- hzar.dataGroup.add(outlier$"11_6A"$runs$doSeq$model1)
outlier$"11_6A"$analysis$initDGs$model2 <- hzar.dataGroup.add(outlier$"11_6A"$runs$doSeq$model2)
outlier$"11_6A"$analysis$initDGs$model3 <- hzar.dataGroup.add(outlier$"11_6A"$runs$doSeq$model3)
outlier$"11_6A"$analysis$initDGs$model4 <- hzar.dataGroup.add(outlier$"11_6A"$runs$doSeq$model4)


##create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme

outlier$"11_6A"$analysis$oDG <- hzar.make.obsDataGroup(outlier$"11_6A"$analysis$initDGs)

outlier$"11_6A"$analysis$oDG <- hzar.copyModelLabels(outlier$"11_6A"$analysis$initDGs,
                                                     outlier$"11_6A"$analysis$oDG)

##convert all 36 runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object
outlier$"11_6A"$analysis$oDG <- hzar.make.obsDataGroup(lapply(outlier$"11_6A"$runs$doSeq, 
                                                               hzar.dataGroup.add),
                                                      outlier$"11_6A"$analysis$oDG);

#check to make sure there are only 5 hzar.dataGroup objects
print(summary(outlier$"11_6A"$analysis$oDG$data.groups))

#compare the 4 cline models to the null model graphically
hzar.plot.cline(outlier$"11_6A"$analysis$oDG)
hzar.plot.fzCline(outlier$"11_6A"$analysis$model.selected)

#model selection based on AICc scores
print(outlier$"11_6A"$analysis$AICcTable <- 
        hzar.AICc.hzar.obsDataGroup(outlier$"11_6A"$analysis$oDG))



#print the model with the minimum AICc score
print(outlier$"11_6A"$analysis$model.name <-
  rownames(outlier$"11_6A"$analysis$AICcTable
           )[[which.min(outlier$"11_6A"$analysis$AICcTable$AICc )]])
#[1] "model3"



#Extract the hzar.dataGroup object for the selected model
outlier$"11_6A"$analysis$model.selected <-
  outlier$"11_6A"$analysis$oDG$data.groups[[outlier$"11_6A"$analysis$model.name]]

#look at the variation in parameters for the selected model
print(hzar.getLLCutParam(outlier$"11_6A"$analysis$model.selected,
                         names(outlier$"11_6A"$analysis$model.selected$data.param)))
#center2LLLow center2LLHigh width2LLLow width2LLHigh pMin2LLLow pMin2LLHigh pMax2LLLow pMax2LLHigh
#1     769.2505      1284.908    12.13201     1369.043  0.2200383   0.3115534  0.5054706   0.7085702


####Print Params####

#print the cline width for the selected model
outlier$"11_6A"$analysis$modeldetails <- print(hzar.get.ML.cline(outlier$"11_6A"$analysis$model.selected))
print(outlier$"11_6A"$analysis$modeldetails$param.all$width)

print(outlier$"3350_66A"$analysis$modeldetails$logLike)


#plot the maximum likelihood cline for the selected model
hzar.plot.cline(outlier$"7115_30A"$analysis$model.selected)

#plot the 95% credible cline region for the selected model
hzar.plot.fzCline(outlier$"7115_30A"$analysis$model.selected)

###DONE###


####Pick loci that are actually clinal!####

#print all AICc scores
print(outlier$"11_6A"$analysis$AICcTable)

#print the min AICc score
print(min(outlier$"11_6A"$analysis$AICcTable$AICc))

#make subset of clines that have AICc values under a certain criteria, then plot







####Save data image####
save.image("~/School/MASTER_Masters/R/GeneticData/Clines_HZAR_Outlier.RData")


#####EXTRA#####
####THE FOLLOWING IS REPLACED BY HZAR.CHAIN.DOSEQ###
#####create space for initial chains
# outlier$"11_6A"$runs$init <- list()

# #run model1 for an initial chain
# outlier$"11_6A"$runs$init$model1 <- hzar.doFit(outlier$"11_6A"$fitRs$init$model1)
# #plot the trace
# plot(hzar.mcmc.bindLL(outlier$"11_6A"$runs$init$model1))
# 
# #run model2 for an initial chain
# outlier$"11_6A"$runs$init$model2 <- hzar.doFit(outlier$"11_6A"$fitRs$init$model2)
# #plot the trace
# plot(hzar.mcmc.bindLL(outlier$"11_6A"$runs$init$model2))
# 
# #run model3 for an initial chain
# outlier$"11_6A"$runs$init$model3 <- hzar.doFit(outlier$"11_6A"$fitRs$init$model3)
# #plot the trace
# plot(hzar.mcmc.bindLL(outlier$"11_6A"$runs$init$model3))
# 
# #run model4 for an initial chain
# outlier$"11_6A"$runs$init$model4 <- hzar.doFit(outlier$"11_6A"$fitRs$init$model4)
# #plot the trace
# plot(hzar.mcmc.bindLL(outlier$"11_6A"$runs$init$model4))

# ####compile a new set of fit requests using the initial chains
# outlier$"11_6A"$fitRs$chains <- lapply(outlier$"11_6A"$runs$init, 
#                                        hzar.next.fitRequest)



# ##to be thorough, randomize the initial value for each fit - 
# #just centre and width, everything else takes too long
# #centre
# 
# runif(12,-30,4500)
# 
# outlier$"11_6A"$fitRs$chains[[1]]$modelParam$init["center"]=2016.7540
# outlier$"11_6A"$fitRs$chains[[2]]$modelParam$init["center"]=3451.1281
# outlier$"11_6A"$fitRs$chains[[3]]$modelParam$init["center"]=1579.0481
# outlier$"11_6A"$fitRs$chains[[4]]$modelParam$init["center"]=4305.2671
# outlier$"11_6A"$fitRs$chains[[5]]$modelParam$init["center"]=3179.8165
# outlier$"11_6A"$fitRs$chains[[6]]$modelParam$init["center"]=906.4475
# outlier$"11_6A"$fitRs$chains[[7]]$modelParam$init["center"]=4277.4323
# outlier$"11_6A"$fitRs$chains[[8]]$modelParam$init["center"]=3999.2289
# outlier$"11_6A"$fitRs$chains[[9]]$modelParam$init["center"]=4454.4944
# outlier$"11_6A"$fitRs$chains[[10]]$modelParam$init["center"]=1149.1651
# outlier$"11_6A"$fitRs$chains[[11]]$modelParam$init["center"]=1897.6752
# outlier$"11_6A"$fitRs$chains[[12]]$modelParam$init["center"]=680.6496
#   
# #width
# runif(12, 0, 4530)
# outlier$"11_6A"$fitRs$chains[[1]]$modelParam$init["width"]=2400.3286
# outlier$"11_6A"$fitRs$chains[[2]]$modelParam$init["width"]=137.4859
# outlier$"11_6A"$fitRs$chains[[3]]$modelParam$init["width"]=2873.5229
# outlier$"11_6A"$fitRs$chains[[4]]$modelParam$init["width"]=2203.0933
# outlier$"11_6A"$fitRs$chains[[5]]$modelParam$init["width"]=3287.5375
# outlier$"11_6A"$fitRs$chains[[6]]$modelParam$init["width"]=2402.2096
# outlier$"11_6A"$fitRs$chains[[7]]$modelParam$init["width"]=3643.4791
# outlier$"11_6A"$fitRs$chains[[8]]$modelParam$init["width"]=1461.4384
# outlier$"11_6A"$fitRs$chains[[9]]$modelParam$init["width"]=507.1691
# outlier$"11_6A"$fitRs$chains[[10]]$modelParam$init["width"]=957.6766
# outlier$"11_6A"$fitRs$chains[[11]]$modelParam$init["width"]=3344.2016
# outlier$"11_6A"$fitRs$chains[[12]]$modelParam$init["width"]=2128.1165



#run a chain of 3 runs for every fit request - this will take a LONG time
#end up with 36 - 4 models, 3 fit requests each (=12) x 3 chains each (36)
#outlier$"11_6A"$runs$chains <- hzar.doChain.multi(outlier$"11_6A"$fitRs$chains,
#                                                  doPar=TRUE,
#                                                  inOrder=FALSE,
#                                                  count=3)
