SourceGitFunc <- function(url)
{
## URL is the raw format link from Github online 
## e.g. "https://raw.githubusercontent.com/rystanley/Rad_R_Functions/master/GenePopConvert.R"
require(RCurl)
script <- getURL(url)
eval(parse(text = script))
}