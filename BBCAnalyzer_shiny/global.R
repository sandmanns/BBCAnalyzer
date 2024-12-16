###1) Install all required packages by executing lines 7 to 20 (one after another)
###
###   Important: if RStudio asks you to restart: click "Yes"
###              if RStudio asks you to update all (a), some (s) or none (n)
###              of the installed packages: type "n"

install.packages("shiny")
install.packages("shinyjs")
install.packages("png")
install.packages("XML")
install.packages("RCurl")

source("https://bioconductor.org/biocLite.R")
biocLite("SummarizedExperiment")
biocLite("VariantAnnotation")
biocLite("Rsamtools")
biocLite("GenomicRanges")
biocLite("IRanges")
biocLite("Biostrings")
biocLite("BSgenome")


###2) Please define the directory, where the downloaded folder 
###   "BBCAnalyzer_local" can be found and execute lines 31, 33 and 34:
###
###   Important: for Windows-User: use double-backslash
###              (e.g. C:\\home\\BBCAnalyzer_local\\)
###              for Linux-User: use single-forwardslash
###              (e.g. /home/BBCAnalyzer_local/)

bbcanalyzer_dir<-file.path("/home/BBCAnalyzer_local/")

library(shiny)
runApp(appDir=bbcanalyzer_dir)


###If you want to re-open the web application in the same R session, 
###execute only line 32
###
###If you want to re-open the web application in a new R session,
###execute lines 29, 31 and 32