\name{analyzeBases}
\alias{analyzeBases}
\alias{analyze Bases}
\title{
Analyze the bases at all previously defined positions
}
\description{
\bold{BBCAnalyzer} performs an analysis of the bases, deletions and insertions at defined positions in sequence alignment data and visualizes the results. The function \code{analyzeBases} performs the whole process of analyzing the data and plotting the results.
}
\usage{
analyzeBases(sample_names, bam_input, target_regions, vcf_input, output,
             output_pictures, known_file, genome, MQ_threshold, BQ_threshold, 
             frequency_threshold, qual_lower_bound, qual_upper_bound, 
             marks, relative, per_sample)
}

\arguments{
  \item{sample_names}{
The txt-file containing the names of the samples to be analyzed.
}
  \item{bam_input}{
The folder containing the alignment data in bam- and bai format or BamFileList.
}
  \item{target_regions}{
The txt-file containing the target regions to be analyzed.
}
  \item{vcf_input}{
The folder containing the vcf files or a VcfFileList of the samples to be analyzed or a multiple-sample vcf file. The argument may be left blank.
}
  \item{output}{
The folder to write the output into.
}
  \item{output_pictures}{
The folder to write the plots into. The argument may be left blank. In this case, the plots are returned to the workspace.
}
  \item{known_file}{
The name of a tabix file containing e.g. known polymorphisms or mutations (dbSNP). The argument may be left blank.
}
  \item{genome}{
A BSgenome object defining reference genome that shall be used for comparison (e.g. BSgenome.Hsapiens.UCSC.hg19).
}
  \item{MQ_threshold}{
A PHRED-scaled value to be used as a mapping quality threshold. All reads with a mapping quality below this threshold are excluded from analysis. Every base in an excluded read gets marked in an output file.
}
  \item{BQ_threshold}{
A PHRED-scaled value to be used as a base quality threshold. All bases with a base quality below this threshold are excluded from analysis. Every excluded base gets marked in an output file. The number of excluded bases per position gets counted.
}
  \item{frequency_threshold}{
A frequency to be used as a threshold for variants to be reported ([0,1]).
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

The bam- and the corresponding bai files of the samples to be analyzed have to be provided in a folder (\code{bam_input}). The names of the files have to match the sample names provided by \code{sample_names}.

The target regions have to be provided by a file (\code{target_regions}). The file may either contain regions (chromosome, tab, first base of a region, tab, last base of a region) or positions (chromosome, tab, position). A mixture of both is not supported. Yet, a region may cover only one base, i.e. the first and last base of a region may be identical.

If vcf files shall be considered, the corresponding files of the samples to be analyzed have to be provided in a folder (\code{vcf_input}). There has to be one file per sample. The names of the files have to match the sample names provided by \code{sample_names}. If a vcf file is not available for one or more samples, empty vcf files may be used instead.

 
About the analysis of the data:

Determine target: If the \code{target_regions} file contains regions to be analyzed, the different positions covered by the regions are determined. If the file already contains single positions, the program directly proceeds with the next step. It is not necessary that the regions or positions are ordered. If a known insert is supposed to be analyzed, the position of the base succeeding the insert has to be given.

Analyze reads: The reads at every targeted position get analyzed. By the help of the CIGAR string the bases, deletions and insertions are determined. The output is saved as [Sample].bases.txt. For every base - also the inserted ones - the base quality is determined. The output is saved as [Sample].quality.txt. Reads with a mapping quality below \code{MQ_threshold} get excluded from the analysis. Instead of a base, "MQ" is noted in the base file. Instead of a quality value, "-2" is noted in the quality file. The function copes with uncovered positions ("NotCovered" in the base- and the quality file) and insertions >1bp (repeated analysis of the position).

Analyze frequency: The number of detected bases, deletions and insertions at every position is summed up. Additionally, the mean quality of the detected bases - including the insertions and for the inserted bases only - is calculated. Bases with a quality below \code{BQ_threshold} are excluded and counted separately. The output is saved as [Sample].frequency.txt. If the analysis shall consider vcf files as well (\code{vcf_input} not blank), the alternate alleles and the genotypes - as far as they are available for the positions analyzed - are written out as well. Furthermore, for every called variant it is noted whether it is an insert. The function copes with up to two different alternate alleles per position. 

Report variants: The ratios of the detected bases, deletions and insertions (additionally) at every position are determined. According to the determined ratios, up to six different calls get reported. If \code{frequency_threshold} is set, minor variants with ratios below this threshold do not get reported. The output is saved as [Sample].calling.txt. The function copes with insertions >1bp even if the position at which an insert is detected is not covered by all samples being analyzed. If the analysis shall consider vcf files as well (\code{vcf_input} not blank), the call - taking the reference allele, the alternate allele(s) and the genotype into account - is written out as well. For every called variant it is noted whether it is an insert and whether the genotype is heterozygous.


About plotting the results:

The absolute number of the detected bases, deletions and insertions for each sample at each targeted position is plotted if \code{relative==FALSE}. Otherwise the relative frequencies of the detected bases, deletions and insertions for each sample at each targeted position get plotted. The bars are colored according to the base (adenine: green; cytosine: blue; guanine: yellow; thymine: red; deletion: black; insertion: lilac edge) and their mean quality (high mean quality: dark color; low mean quality: light color). The reference bases (using the defined package \code{ref_genome}) are plotted on the negative y-axis below each position. If a file containing known variants or mutations is provided (\code{known_file}), more than one reference base is plotted at the corresponding position. For each position to be analyzed, lines are drawn at the heights of the ratios defined in \code{marks}. Every targeted position is labelled according to the chromosome and the position. The function copes with different inserted bases at one position (stacked bars) and insertions $>$1bp, even if these are not covered by all samples. If the analysis shall consider vcf files as well (\code{vcf_input} not blank), the expected number of detected bases, deletions and insertions -- according to the vcf file -- is added to the plot using dashed lines.

Plot per sample: One barplot per sample is created. The output is saved as [Sample].png.

Plot per target: One barplot per targeted position is created. The output is saved as chr[number];[position].png.
}
\value{
A list is returned:

\item{bases }{A list containing the bases for all samples and positions being analyzed (identical compared to [Sample].bases.txt).}
\item{quality}{A list containing the corresponding qualities for all samples and positions being analyzed (identical compared to [Sample].quality.txt).}
\item{frequency }{A list containing the summed up number of bases for all samples and positions being analyzed (identical compared to [Sample].frequency.txt).}
\item{calling }{A list containing the relative frequencies and the potential calls for all samples and positions being analyzed (identical compared to [Sample].calling.txt).}

}
\references{
More information on the bam format can be found at: \url{http://samtools.github.io/hts-specs/SAMv1.pdf}

db SNP -- Short Genetic variations: \url{http://www.ncbi.nlm.nih.gov/SNP/}
}
\author{
Sarah Sandmann <sarah.sandmann@uni-muenster.de>
}

\seealso{
\code{\link{BBCAnalyzer}}, \code{\link{analyzeBasesPlotOnly}}
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
