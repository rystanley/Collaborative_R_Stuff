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

Hybridpower(dir="C:/Users/RyanStanley/OneDrive/PostDoc/DFO/Salmon/Frequency Based Sim/West/West Correct Wild Fst Top Loci/FileRun/",
            filetag="Top192",Threshold=0.7)