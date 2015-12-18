file.choose()
require(parallel)
## Example
#PGDSpider.GPtoX.auto(Convert.folder = "/Users/brendanwringe/Desktop/untitled folder/", where.PGDSpider = "/Users/brendanwringe/Desktop/DFO Aquaculture Interaction/Software/PGDSpider_2.0.9.0/", out.format = "NewHybrids")

### out.format must be either GENEPOP or NEWHYBRIDS - all caps

PGDSpider.GPtoX.auto <- function(Convert.folder, where.PGDSpider, out.format){
    
  where.PGDSpider <- gsub(pattern = " ", replacement = "\\ ", x = where.PGDSpider, fixed = T)
  #Convert.folder <- gsub(pattern = " ", replacement = "\\ ", x = Convert.folder, fixed = T)
  #gsub(x = where.CLUMPP, pattern = " ", replacement = "\\ ", fixed = T)
  
    conv.files.list <- list.files(path = Convert.folder)
    Convert.folder <- gsub(pattern = " ", replacement = "\\ ", x = Convert.folder, fixed = T)
    conv.files.path <- paste0("-inputfile ", Convert.folder, conv.files.list)
    
    out.file.names <- gsub(pattern = ".txt", x = conv.files.list, replacement = "NH.txt")
    out.file.path <- paste0("-outputfile ", Convert.folder, out.file.names)
    
    execute.SPIDER <- "java -Xmx1024m -Xms512m -jar PGDSpider2-cli.jar"
    
   
    if(out.format=="GENEPOP"){
      outputformat <- "-outputformat GENEPOP"
      spidcall <- paste0("-spid ", where.spid, "GENEPOP_to_GENEPOP.spid")
    }
    if(out.format=="NEWHYBRIDS"){
      outputformat <- "-outputformat NEWHYBRIDS"
      spidcall <- paste0("-spid ", where.PGDSpider, "GENEPOP_to_NEWHYBRIDS.spid")
    }
    
    
  
    
    input.format <- "-inputformat GENEPOP"
    
    goto.SPIDER <- paste0("cd ", where.PGDSpider, "; ", execute.SPIDER)
    
    run.call.SPIDER <- paste0(goto.SPIDER, " ", conv.files.path, " ", input.format, " ", out.file.path, " ",  outputformat, " ", spidcall)
    
    
    mclapply(X = run.call.SPIDER, FUN = system)
    

      
    }


