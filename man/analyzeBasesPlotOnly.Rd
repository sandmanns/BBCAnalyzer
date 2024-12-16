\name{analyzeBasesPlotOnly}
\alias{analyzeBasesPlotOnly}
\alias{analyze Bases Plot Only}
\title{
Plots the number of reads at all previously defined positions
}
\description{
To allow for a quick change in the way the analysis of the bases, deletions and insertions at defined positions in sequence alignment data is visualized, the function \code{analyzeBasesPlotOnly} may be used. It solely performs the last step of the whole analysis pipeline -- the plotting of the results.
}
\usage{
analyzeBasesPlotOnly(sample_names, vcf_input, output, known_file, 
                     output_list, qual_lower_bound, qual_upper_bound, 
                     marks, relative, per_sample)
}
\arguments{
  \item{sample_names}{
The file containing the names of the samples to be analyzed.
}
  \item{vcf_input}{
The folder containing the vcf files or a VcfFileList of the samples to be analyzed or a multiple-sample vcf file. The argument may be left blank.
}
  \item{output}{
The folder to write the output (plots) into. The argument may be left blank. In this case, the plots are returned to the workspace.
}
  \item{known_file}{
The name of a tabix file containing e.g. known polymorphisms or mutations (dbSNP). The argument may be left blank.
}
  \item{output_list}{
The name of the list that is returned by the function \code{analyzeBases}.
}
  \item{qual_lower_bound}{
The lower bound for the mean quality that shall be color-coded in the plots. All bases with a mean quality below \code{qual_lower_bound} are colored with the lightest color defined for the corresponding base. If the bases shall not be color-coded according to their mean quality, \code{qual_lower_bound} has to be identical compared to \code{qual_upper_bound}.
}
  \item{qual_upper_bound}{
The upper bound for the mean quality that shall be color-coded in the plots. All bases with a mean quality above \code{qual_upper_bound} are colored with the darkest color defined for the corresponding base. If the bases shall not be color-coded according to their mean quality, \code{qual_upper_bound} has to be identical compared to \code{qual_lower_bound}.
}
  \item{marks}{
A vector of floats [0,1] defining the levels at which marks shall be drawn in the plot.
}
  \item{relative}{
A boolean defining whether the relative (\code{true}) or the absolute (\code{false}) number of reads shall be plotted.
}
  \item{per_sample}{
A boolean defining whether one plot per sample (\code{true}) or one plot per target (\code{false}) shall be created.
}
}
\details{
About the input:

Names of the samples to be analyzed have to be provided by a file (\code{sample_names}). There has to be one sample name per line without the ".bam"-suffix.

If vcf files shall be considered, the corresponding files of the samples to be analyzed have to be provided in a folder (\code{vcf_input}). There has to be one file per sample. The names of the files have to match the sample names provided by \code{sample_names}. If a vcf file is not available for one or more samples, empty vcf files may be used instead.

The value that gets returned by the function \code{analyzeBases} has to be provided (\code{output_list}). Otherwise, \code{analyzeBasesPlotOnly} will not be able to use the previously generated output.


About plotting the results:

The absolute number of the detected bases, deletions and insertions for each sample at each targeted position is plotted if \code{relative==FALSE}. Otherwise the relative frequencies of the detected bases, deletions and insertions for each sample at each targeted position get plotted. The bars are colored according to the base (adenine: green; cytosine: blue; guanine: yellow; thymine: red; deletion: black; insertion: lilac edge) and their mean quality (high mean quality: dark color; low mean quality: light color). The reference bases (using the defined package \code{ref_genome}) are plotted on the negative y-axis below each position. If a file containing known variants or mutations is provided (\code{known_file}), more than one reference base is plotted at the corresponding position. For each position to be analyzed, lines are drawn at the heights of the ratios defined in \code{marks}. Every targeted position is labelled according to the chromosome and the position. The function copes with different inserted bases at one position (stacked bars) and insertions $>$1bp, even if these are not covered by all samples. If the analysis shall consider vcf files as well (\code{vcf_input} not blank), the expected number of detected bases, deletions and insertions -- according to the vcf file -- is added to the plot using dashed lines.

Plot per sample: One barplot per sample is created. The output is saved as [Sample].png.

Plot per target: One barplot per targeted position is created. The output is saved as chr[number];[position].png.
}

\value{
No value is returned.
}

\references{
db SNP -- Short Genetic variations: \url{http://www.ncbi.nlm.nih.gov/SNP/}
}
\author{
Sarah Sandmann <sarah.sandmann@uni-muenster.de>
}

\seealso{
\code{\link{BBCAnalyzer}}, \code{\link{analyzeBases}}
}
\examples{
library("BSgenome.Hsapiens.UCSC.hg19")
ref_genome<-BSgenome.Hsapiens.UCSC.hg19

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
                     
analyzeBasesPlotOnly(sample_names=system.file("extdata","SampleNames_small.txt",package="BBCAnalyzer"),
                     vcf_input="",
                     output=system.file("extdata",package="BBCAnalyzer"),
                     known_file="",
                     output_list=output,
                     qual_lower_bound=58,
                     qual_upper_bound=63,
                     marks=c(0.25,0.5,0.75,1),
                     relative=FALSE,
                     per_sample=TRUE)                     
}
