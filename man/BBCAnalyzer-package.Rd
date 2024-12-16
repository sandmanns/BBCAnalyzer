\name{BBCAnalyzer-package}
\alias{BBCAnalyzer}
\alias{Bases By Cigar Analyzer}
\alias{BasesByCigar Analyzer}
\alias{BasesByCigarAnalyzer}
\alias{Bases Coverage}
\alias{Cigar Analysis}
\docType{package}
\title{
\packageTitle{BBCAnalyzer}
}
\description{
\packageDescription{BBCAnalyzer}
}
\details{

\packageDESCRIPTION{BBCAnalyzer} 
In the use case of medical diagnostics, a tool performing detailed analyses of those locations where mutations may be expected -- but not always called -- appears to be most useful. Low allele frequency and bad base quality do often explain a lacking call. Yet, this information is not included in a VCF-file and difficult to obtain from other existing tool. Furthermore, with regards to the comparison of different sequencing techniques, it seems helpful to have a tool for visualizing the background at a selection of locations where e.g. one technique calls a variant but another technique does not.

BBCAnalyzer (Bases By CIGAR Analyzer) is a tool for visualizing the number of counted bases, deletions and insertions at any given position in any genome in comparison to the reference bases. Relative frequencies, base qualities, known mutations or polymorphisms and called variants may be included into the plots as well.

\packageIndices{BBCAnalyzer}
The package contains two functions - \code{analyzeBases} and \code{analyzeBasesPlotOnly}. The major use of BBCAnalyzer is documented in the description of the function \code{analyzeBases}. The function \code{analyzeBasesPlotOnly} serves as an extension.
}
\author{
\packageAuthor{BBCAnalyzer}

Maintainer: \packageMaintainer{BBCAnalyzer}
}
\references{
More information on the bam format can be found at: \url{http://samtools.github.io/hts-specs/SAMv1.pdf}
}
\keyword{ package }
\seealso{
\code{\link{analyzeBases}}, \code{\link{analyzeBasesPlotOnly}}
}
\examples{
library("BSgenome.Hsapiens.UCSC.hg19")
ref_genome<-BSgenome.Hsapiens.UCSC.hg19
\donttest{
output<-analyzeBases(sample_names=system.file("extdata","SampleNames_small.txt",package="BBCAnalyzer"),
                     bam_input=system.file("extdata",package="BBCAnalyzer"),
                     target_regions=system.file("extdata","targetRegions_small.txt",package="BBCAnalyzer"),
                     vcf_input="",
                     output=system.file("extdata",package="BBCAnalyzer"),
                     output_pictures=system.file("extdata",package="BBCAnalyzer"),
                     known_file="",
                     genome=ref_genome,
                     MQ_threshold=60,
                     BQ_threshold=50,
                     frequency_threshold=0.01,
                     qual_lower_bound=58,
                     qual_upper_bound=63,
                     marks=c(0.01),
                     relative=TRUE,
                     per_sample=TRUE)
             
} 
}
