## Hybrid power analysis

library(RCurl)
SourceGitFunc <- function(url)
{
  ## URL is the raw format link from Github online   
  ## e.g. "https://raw.githubusercontent.com/rystanley/RAD_R_Functions/master/GenoPopConvert.R"
  require(RCurl)
  script <- getURL(url, ssl.verifypeer = FALSE)
  eval(parse(text = script),envir=.GlobalEnv)
}

SourceGitFunc("https://raw.githubusercontent.com/rystanley/Collaborative_R_Stuff/master/HybridPower.R")
SourceGitFunc("https://raw.githubusercontent.com/rystanley/Collaborative_R_Stuff/master/HybridPowerComp.R")


#directory where the New Hybrids runs are for a given subset of Loci
Hybridpower(dir="C:/Users/RyanStanley/OneDrive/PostDoc/DFO/Salmon/Frequency Based Sim/West/West Correct Wild Fst Top Loci/FileRun/",
            filetag="Top192WestNL",Threshold=0.75)

#directory where the New Hybrids runs are compared among different subsets of Loci (i.e. 48, 96, 144, 192 and 240)
Hybridpower_comparison(dir="C:/Users/RyanStanley/OneDrive/PostDoc/DFO/Salmon/Frequency Based Sim/West/West Correct Wild Fst Top Loci/FileRun/",
                       filetag="WestNL")