library(shiny)
library(shinyjs)
# Define UI for application that draws a histogram
shinyUI(fluidPage(
    # Application title
    titlePanel("BBCAnalyzer"),
    shinyjs::useShinyjs(),
    fluidRow(
        column(3,wellPanel(
            textAreaInput('sample_names', 'Samples to analyze', cols=1,rows=1000,height="100","Example_454\nExample_IonTorrent"),
            textInput('bam_input', 'Define folder containing bam- and bai-files', "/home/BBCAnalyzer_local/input/"),
            textAreaInput('target_regions', 'Target regions to analyze', cols=1,rows=1000,height="100","20\t31024054\n21\t36206893\n4\t106157698"),
            textInput('vcf_input', 'Define folder containing vcf files (optional)', ""),
            textInput('output', 'Define output folder', "/home/BBCAnalyzer_local/output/"),
            textInput('known_file', 'Define tabix file containing known variants (optional)', ""),
            selectInput('genome', 'Select reference genome for analysis',choices = c("BSgenome.Alyrata.JGI.v1","BSgenome.Amellifera.BeeBase.assembly4","BSgenome.Amellifera.UCSC.apiMel2",
                                                                                     "BSgenome.Amellifera.UCSC.apiMel2.masked","BSgenome.Athaliana.TAIR.04232008","BSgenome.Athaliana.TAIR.TAIR9",
                                                                                     "BSgenome.Btaurus.UCSC.bosTau3","BSgenome.Btaurus.UCSC.bosTau3.masked","BSgenome.Btaurus.UCSC.bosTau4",
                                                                                     "BSgenome.Btaurus.UCSC.bosTau4.masked","BSgenome.Btaurus.UCSC.bosTau6","BSgenome.Btaurus.UCSC.bosTau6.masked",
                                                                                     "BSgenome.Btaurus.UCSC.bosTau8","BSgenome.Celegans.UCSC.ce10","BSgenome.Celegans.UCSC.ce11","BSgenome.Celegans.UCSC.ce2",
                                                                                     "BSgenome.Celegans.UCSC.ce6","BSgenome.Cfamiliaris.UCSC.canFam2","BSgenome.Cfamiliaris.UCSC.canFam2.masked",
                                                                                     "BSgenome.Cfamiliaris.UCSC.canFam3","BSgenome.Cfamiliaris.UCSC.canFam3.masked","BSgenome.Dmelanogaster.UCSC.dm2",
                                                                                     "BSgenome.Dmelanogaster.UCSC.dm2.masked","BSgenome.Dmelanogaster.UCSC.dm3","BSgenome.Dmelanogaster.UCSC.dm3.masked",
                                                                                     "BSgenome.Dmelanogaster.UCSC.dm6","BSgenome.Drerio.UCSC.danRer10","BSgenome.Drerio.UCSC.danRer5",
                                                                                     "BSgenome.Drerio.UCSC.danRer5.masked","BSgenome.Drerio.UCSC.danRer6","BSgenome.Drerio.UCSC.danRer6.masked",
                                                                                     "BSgenome.Drerio.UCSC.danRer7","BSgenome.Drerio.UCSC.danRer7.masked","BSgenome.Ecoli.NCBI.20080805",
                                                                                     "BSgenome.Gaculeatus.UCSC.gasAcu1","BSgenome.Gaculeatus.UCSC.gasAcu1.masked","BSgenome.Ggallus.UCSC.galGal3",
                                                                                     "BSgenome.Ggallus.UCSC.galGal3.masked","BSgenome.Ggallus.UCSC.galGal4","BSgenome.Ggallus.UCSC.galGal4.masked",
                                                                                     "BSgenome.Hsapiens.1000genomes.hs37d5","BSgenome.Hsapiens.NCBI.GRCh38","BSgenome.Hsapiens.UCSC.hg17",
                                                                                     "BSgenome.Hsapiens.UCSC.hg17.masked","BSgenome.Hsapiens.UCSC.hg18","BSgenome.Hsapiens.UCSC.hg18.masked",
                                                                                     "BSgenome.Hsapiens.UCSC.hg19","BSgenome.Hsapiens.UCSC.hg19.masked","BSgenome.Hsapiens.UCSC.hg38",
                                                                                     "BSgenome.Hsapiens.UCSC.hg38.masked","BSgenome.Mfascicularis.NCBI.5.0","BSgenome.Mfuro.UCSC.musFur1",
                                                                                     "BSgenome.Mmulatta.UCSC.rheMac2","BSgenome.Mmulatta.UCSC.rheMac2.masked","BSgenome.Mmulatta.UCSC.rheMac3",
                                                                                     "BSgenome.Mmulatta.UCSC.rheMac3.masked","BSgenome.Mmusculus.UCSC.mm10","BSgenome.Mmusculus.UCSC.mm10.masked",
                                                                                     "BSgenome.Mmusculus.UCSC.mm8","BSgenome.Mmusculus.UCSC.mm8.masked","BSgenome.Mmusculus.UCSC.mm9",
                                                                                     "BSgenome.Mmusculus.UCSC.mm9.masked","BSgenome.Osativa.MSU.MSU7","BSgenome.Ptroglodytes.UCSC.panTro2",
                                                                                     "BSgenome.Ptroglodytes.UCSC.panTro2.masked","BSgenome.Ptroglodytes.UCSC.panTro3",
                                                                                     "BSgenome.Ptroglodytes.UCSC.panTro3.masked","BSgenome.Rnorvegicus.UCSC.rn4","BSgenome.Rnorvegicus.UCSC.rn4.masked",
                                                                                     "BSgenome.Rnorvegicus.UCSC.rn5","BSgenome.Rnorvegicus.UCSC.rn5.masked","BSgenome.Rnorvegicus.UCSC.rn6",
                                                                                     "BSgenome.Scerevisiae.UCSC.sacCer1","BSgenome.Scerevisiae.UCSC.sacCer2","BSgenome.Scerevisiae.UCSC.sacCer3",
                                                                                     "BSgenome.Sscrofa.UCSC.susScr3","BSgenome.Sscrofa.UCSC.susScr3.masked","BSgenome.Tgondii.ToxoDB.7.0",
                                                                                     "BSgenome.Tguttata.UCSC.taeGut1","BSgenome.Tguttata.UCSC.taeGut1.masked","BSgenome.Tguttata.UCSC.taeGut2",
                                                                                     "BSgenome.Vvinifera.URGI.IGGP12Xv0","BSgenome.Vvinifera.URGI.IGGP12Xv2"),selected="BSgenome.Hsapiens.UCSC.hg19"),
            sliderInput('MQ_threshold', 'Mapping quality threshold', 60,min = 0, max = 99),
            sliderInput('BQ_threshold', 'Base quality threshold', 50,min = 0, max = 99),
            sliderInput('frequency_threshold', 'Frequency threshold for variant reporting', 0.01,min = 0, max = 1),
            sliderInput('qual_bounds', 'Lower- and upper mean quality bound for color-coding',dragRange = T,min=0,max=99,value = c(58,63)),
            checkboxGroupInput('marks','Select levels at which marks shall be drawn',
                               choices=c(0.01,0.02,seq(0.05,1,0.05)),inline=T,width=150),
            radioButtons('relative_choice',"Plot number of reads",choices = c("Relative","Absolute"),selected="Relative"),
            radioButtons('per_sample_choice',"Create one plot per",choices = c("Sample","Position"),selected = "Sample"),
            actionButton("do", "Start complete analysis with BBCAnalyzer")
            
        )),
        column(3,wellPanel(
            textAreaInput('sample_names2', 'Samples to analyze', cols=1,rows=1000,height="100","Example_454\nExample_IonTorrent"),
            radioButtons('vcf_input2', 'Consider vcf file information (only possible if evaluated in complete analysis)',
                         choices=c("Yes","No"),selected="No"),
            textInput('output2', 'Define output folder', "/home/BBCAnalyzer_local/output/"),
            textInput('known_file2', 'Define tabix file containing known variants (optional)', ""),
            sliderInput('qual_bounds2', 'Lower- and upper mean quality bound for color-coding',dragRange = T,min=0,max=99,value = c(58,63)),
            checkboxGroupInput('marks2','Select levels at which marks shall be drawn',
                               choices=c(0.01,0.02,seq(0.05,1,0.05)),inline=T,width=150),
            radioButtons('relative_choice2',"Plot number of reads",choices = c("Relative","Absolute"),selected="Relative"),
            radioButtons('per_sample_choice2',"Create one plot per",choices = c("Sample","Position"),selected = "Sample"),
            actionButton("do2", "Create plots only")
        )),
        
        column(6,wellPanel(
            div(id = "text"),
            plotOutput("plot")
        ))
    )
)
)
