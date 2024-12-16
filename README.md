# BBCAnalyzer: a visual approach to facilitate variant calling

<p align="center">
    <img height="600" src="https://uni-muenster.sciebo.de/s/nc98cPy3z7IsFJK/download">
</p>


Sandmann, S., de Graaf, A.O. & Dugas, M. BBCAnalyzer: a visual approach to facilitate variant calling. BMC Bioinformatics 18, 133 (2017). https://doi.org/10.1186/s12859-017-1549-4


The analysis of alignment data is an essential step in the process of analyzing sequencing data. Mutations like single nucleotide base changes, deletions and insertions may be identified in this step.

Usually, the available programs -- like GATK \cite{gatk} or SAMtools \cite{samtools} -- take a bam file, containing the alignment data in a compressed form, as input and return a vcf file containing the called variants and a selection of parameters characterizing them. The overall depth as well as the allele frequencies are returned as two of these parameters.

However, the programs do usually not report the unfiltered number of reads covering a certain mutated base, but they apply a set of partially complex filtration steps. Thresholds -- internally or user-defined -- are usually used to reject additional variants with a minor ratio. Furthermore, to our knowledge it is not possible to perform a comparable analysis of positions where no variant is called.

The integrative genomics viewer (IGV) \cite{igv} provides a possibility for investigating the different number of bases, deletions and insertions at any position in the genome. Yet, it is time consuming to load different samples into the program and to look at the different positions of interest. Moreover, there is no way of automatically summing up and visualizing the base counts in IGV. 

With regards to the comparison of different sequencing techniques, it appears useful to have a tool that is able to visualize the background at a selection of locations where e.g. one technique calls a variant but another technique does not. Furthermore, it seems helpful to have a possibility for quickly analyzing positions of expected mutations, which have not been called.

This package provides a possibility for visualizing the number of counted bases, deletions and insertions at any given position in any genome in comparison to the reference bases. Additionally,

* markers for the relative base frequencies, 
* the mean qualities of the detected bases,
* known mutations or polymorphisms (e.g. based on dbSNP \cite{dbsnp}) and
* called variants in the data

may equally be included into the plots.


## Requirements
To run the BBCAnalyzer, you need R (Version 4.1.0 or higher).

## Installation
BBCAnalyzer is available at Bioconductor. The latest version can be installed via:

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("BBCAnalyzer")
```

The latest version available at github can be installed via:

```
if (!requireNamespace("devtools", quietly=TRUE))
  install.packages("devtools")
devtools::install_github("sandmanns/BBCAnalyzer")
```

## Important Note

The BBCAnalyzer available at Bioconductor just contains 'common' R functions. If you prefer working with an R Shiny GUI, download the BBCAnalyzer_shiny folder.


## Detailed documentation
For detailed documentation, please check out the manual and the vignette available within this repsitory or on the bioconductor website (https://bioconductor.org/packages/release/bioc/html/BBCAnalyzer.html).

In case of errors or feature requests, do not hesitate to open an issue or contact Sarah Sandmann (sarah.sandmann@uni-muenster.de).
