SourceGitFunc <- function(url)
{
## URL is the raw format link from Github online 
## e.g. "https://raw.githubusercontent.com/rystanley/RAD_R_Functions/master/GenoPopConvert.R"
require(RCurl)
script <- getURL(url, ssl.verifypeer = FALSE)
eval(parse(text = script),envir=.GlobalEnv)
}