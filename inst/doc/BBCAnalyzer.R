### R code from vignette source 'BBCAnalyzer.Rnw'

###################################################
### code chunk number 1: BBCAnalyzer.Rnw:21-22
###################################################
options(width=100)


###################################################
### code chunk number 2: 3installing (eval = FALSE)
###################################################
## BiocManager::install("BBCAnalyzer")


###################################################
### code chunk number 3: 4loading
###################################################
library(BBCAnalyzer)


###################################################
### code chunk number 4: 5
###################################################
sample_file <- system.file("extdata", "SampleNames.txt", package = "BBCAnalyzer")
samples <- read.table(sample_file)
samples


###################################################
### code chunk number 5: 6
###################################################
target_file <- system.file("extdata", "targetRegions.txt", package = "BBCAnalyzer")
targetRegions <- read.table(target_file)
targetRegions


###################################################
### code chunk number 6: 7
###################################################
vcf_file <- system.file("extdata", "Example_IonTorrent.vcf", package = "BBCAnalyzer")
vcf_IT <- read.table(vcf_file)
vcf_IT


###################################################
### code chunk number 7: 8
###################################################
#Example1
library("BSgenome.Hsapiens.UCSC.hg19")
ref_genome<-BSgenome.Hsapiens.UCSC.hg19

output1<-analyzeBases(sample_names=sample_file,
            bam_input=system.file("extdata",package="BBCAnalyzer"),
            target_regions=target_file,
            vcf_input="",
            output=system.file("extdata",package="BBCAnalyzer"),
            output_pictures="",
            known_file=system.file("extdata","dbsnp_138.part.vcf.gz",package="BBCAnalyzer"),
            genome=ref_genome,
            MQ_threshold=60,
            BQ_threshold=50,
            frequency_threshold=0.01,
            qual_lower_bound=58,
            qual_upper_bound=63,
            marks=c(0.01),
            relative=TRUE,
            per_sample=TRUE)


###################################################
### code chunk number 8: 9
###################################################
head(output1[[1]][[1]])


###################################################
### code chunk number 9: 9
###################################################
head(output1[[2]][[1]])


###################################################
### code chunk number 10: 9
###################################################
head(output1[[1]][[2]])


###################################################
### code chunk number 11: 9
###################################################
head(output1[[2]][[2]])


###################################################
### code chunk number 12: 9
###################################################
head(output1[[3]][[1]])


###################################################
### code chunk number 13: 9
###################################################
head(output1[[3]][[2]])


###################################################
### code chunk number 14: 9
###################################################
head(output1[[4]][[1]])


###################################################
### code chunk number 15: 9
###################################################
head(output1[[4]][[2]])


###################################################
### code chunk number 16: Example2 (eval = FALSE)
###################################################
## analyzeBasesPlotOnly(sample_names=sample_file,
##             vcf_input=system.file("extdata",package="BBCAnalyzer"),
##             output="",
##             known_file=system.file("extdata","dbsnp_138.part.vcf.gz",package="BBCAnalyzer"),            
##             output_list=output_list,
##             qual_lower_bound=58,
##             qual_upper_bound=63,
##             marks=c(0.2,0.4,0.6,0.8,1),
##             relative=FALSE,
##             per_sample=FALSE)


###################################################
### code chunk number 17: Example3 (eval = FALSE)
###################################################
## output<-analyzeBases(
##             sample_names=system.file("extdata","SampleNames_small.txt",package="BBCAnalyzer"),
##             bam_input=system.file("extdata",package="BBCAnalyzer"),
##             target_regions=system.file("extdata","targetRegions_small.txt",package="BBCAnalyzer"),
##             vcf_input="",
##             output=system.file("extdata",package="BBCAnalyzer"),
##             output_pictures=system.file("extdata",package="BBCAnalyzer"),
##             known_file="", 
##             genome=ref_genome,
##             MQ_threshold=60,
##             BQ_threshold=50,
##             frequency_threshold=0.01,
##             qual_lower_bound=58,
##             qual_upper_bound=63,
##             marks=c(0.01),
##             relative=TRUE,
##             per_sample=TRUE)


