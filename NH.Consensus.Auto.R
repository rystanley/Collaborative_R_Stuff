NH.consensus.auto <- function(res.folder, where.CLUMPP, res.name, plot.it="No"){
    require(plyr)
require(reshape2)
require(stringr)
require(ggplot2)
    res.files.list <- list.files(path = res.folder)
    
    num.runs <- length(res.files.list) ## CLUMPP Requires the number of runs to be specified. This is the # of times New Hybrids was run
    
    ### amalgamate the runs into one data frame
    
    out.data <- NULL
    for(i in 1:length(res.files.list)){
      name.i <- res.files.list[i]
      
      where.data.i <- paste0(res.folder, res.files.list[i])
      
      hold.data.i <- read.table(file = where.data.i, header = T)
      
      out.data <- rbind(out.data, hold.data.i)
      
    }
    
    dim(out.data)
    
    no.row.vec <- c(1:nrow(out.data))
    ind.name.vec <- noquote(paste0("(", out.data[,2], ")"))
    pop.vec <- rep(4, length= nrow(out.data))
    colon.vec <- rep(":", length = nrow(out.data))
    
    data.clumpp.format <- data.frame(no.row.vec, out.data[,1], ind.name.vec, pop.vec, colon.vec, out.data[,c(3:length(out.data))])
    
    ### output this to the clump folder
    
    CLUMPP.indv.file.PP <- paste0(where.CLUMPP, "/", res.name, ".indfile") #### individual file name plus path
    CLUMPP.ind.file <- paste0(res.name, ".indfile") ## just individual file name
    
    write.table(data.clumpp.format, file = CLUMPP.indv.file.PP, col.names = F, row.names = F, quote = FALSE, sep = "\t") ### write the combined runs out in the CLUMPP format
    
    ## CLUMPP requires a Parameter file - create here
      
    CLUMPP.data.type <- paste("DATATYPE", "0")
    
    CLUMPP.ind.file.1 <- paste("INDFILE", CLUMPP.ind.file, sep = " ")
     
    CLUMPP.outfile <- paste("OUTFILE", paste0(res.name, ".outfile"))
    
    CLUMPP.miscfile <- paste("MISCFILE", paste0(res.name, ".miscfile"))
    
    CLUMPP.K <- paste("K", length(3:length(out.data))) ## Number of "clusters" = generations/populations
    
    CLUMPP.C <- paste("C", length(unique(out.data[,1]))) ### number of individuals - calculate from the nubmer of unique values
    
    CLUMPP.R <- paste("R", i) ## number of runs - take from the number of files entered
    
    CLUMPP.M <- paste("M", 2) ### method ### use Greedy
    
    Greedy.option <- paste("GREEDY_OPTION", 2)
    
    Greedy.repeats <- paste("REPEATS", i)
    
    CLUMPP.W <- paste("W", 0) ## weight by number of individuals in each population as specified by datafile - no data file specified at this point, so leave as "No"
    
    CLUMPP.S <- paste("S", 2) ### pairwise similarity statistic to be used - leave as default for now
    
    O.Warnings <- "OVERRIDE_WARNINGS 1"
    
    Ord.run <- "ORDER_BY_RUN 0"
    
    P.permuted.data <- "PRINT_PERMUTED_DATA 1"
    
    Permuted.data.file <- paste("PERMUTED_DATAFILE", paste0(res.name, ".perm_datafile"))
    
    P.every.perm <- "PRINT_EVERY_PERM 0"
    
    E.perm.file <- paste("EVERY_PERMFILE", paste0(res.name, ".every.permfile"))
    
    P.random.input.ord <- "PRINT_RANDOM_INPUTORDER 0"
    
    params.out <- data.frame(rbind(CLUMPP.data.type, CLUMPP.ind.file.1, CLUMPP.outfile, CLUMPP.miscfile, CLUMPP.K, CLUMPP.C, CLUMPP.R, CLUMPP.M, Greedy.option, 
      Greedy.repeats, CLUMPP.W, CLUMPP.S, O.Warnings, Ord.run, P.permuted.data, Permuted.data.file, P.every.perm, E.perm.file, P.random.input.ord))
    
    
    CLUMPP.param.file <- paste0(res.name, ".paramfile")
    CLUMPP.param.file.PP <- paste0(where.CLUMPP,"/", CLUMPP.param.file) #### param file name plus path
    
    write.table(x = params.out, file = CLUMPP.param.file.PP, col.names = F, row.names = F, quote = FALSE, sep = "\t") ### write the combined runs out in the CLUMPP format
    
    ### run CLUMPP
    
    where.CLUMPP2 <- gsub(x = where.CLUMPP, pattern = " ", replacement = "\\ ", fixed = T)
    
    to.run.CLUMPP <- paste0("cd ", where.CLUMPP2, "; ./CLUMPP ", CLUMPP.param.file)
    
    system(to.run.CLUMPP)
    
    clump.res.file.list <- list.files(path = where.CLUMPP)
    
    file.to.get <- paste0(where.CLUMPP, "/", res.name, ".outfile")
   
    run.consensus <- read.table(file.to.get)
    
    run.consensus.convert.to.orig.format <- run.consensus[-c(2,3,4,5)]
    
    colnames(run.consensus.convert.to.orig.format) <- c("Indv", "Pure1", "Pure2", "F1", "F2", "BC1", "BC2")
    
    ## now export the consensus in NH format
    
    NH.format.out <- paste0(res.folder, res.name, ".consensus.csv")
    
    write.csv(x = run.consensus.convert.to.orig.format, file = NH.format.out, row.names = FALSE)

    str_detect(clump.res.file.list, pattern = res.name)
    
    clump.res.file.list.to.copy <- clump.res.file.list[which(str_detect(clump.res.file.list, pattern = res.name))]
    
    clump.res.file.from <- paste0(where.CLUMPP, "/", clump.res.file.list.to.copy)
    clump.res.file.to <- paste0(res.folder, clump.res.file.list.to.copy)
    
    file.copy(from = clump.res.file.from, to = clump.res.file.to)
    
    file.remove(clump.res.file.from)
    
    
    if(plot.it == "Yes"){
      NH_melt <- melt(data = run.consensus.convert.to.orig.format, id.vars = "Indv") ## melt the data to allow the data to be stacked by indivudal
colnames(NH_melt) <- c("Indv", "PopProb", "CumProb") ## rename so that its prettier

## lets give the plot some pretty colours
col.vec <- c("red", "blue", "grey", "green", "black", "yellow", "brown")

## to be used later if decide to break the data or if there are different numbers of individuals in each population type
# break.by <- nrow(NH_output)/6
# break.vec <- c(break.by, (break.by*2), (break.by*3), (break.by*4), (break.by*5), (break.by*6))

## make a nice pretty little plot - and save that bad boy
f.name <- paste0(res.folder, res.name, ".consensus.pdf")
pdf(file = f.name)
pretty.plot <- ggplot(NH_melt, aes(x = Indv, y=CumProb, fill=PopProb))
print(pretty.plot+geom_bar(stat="identity", position = "stack") + scale_fill_manual(values=col.vec)+ylab("Cumulative Probability")+xlab("Individual"))
dev.off()
      
      
      
    }
    
    
    
    
}