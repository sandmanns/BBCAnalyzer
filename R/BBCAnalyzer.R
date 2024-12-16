checkBamInput<-function(sample_names,bam_input){
    if(is.character(bam_input)==TRUE){
        for(i in 1:length(sample_names[,1])){
            if(file.exists(paste(bam_input,"/",sample_names[i,1],".bam",sep=""))==FALSE){
                message("File: ",
                        paste(bam_input,"/",sample_names[i,1],".bam",sep=""),
                        " is missing")
            }
            if(file.exists(paste(bam_input,"/",sample_names[i,1],".bai",sep=""))==FALSE){
                message("File: ",
                        paste(bam_input,"/",sample_names[i,1],".bai",sep=""),
                        " is missing")
            }
        }
        return()
    }
    if(is.character(bam_input)==FALSE){
        for(i in 1:length(sample_names[,1])){
            if(substr(basename(path(bam_input)[i]),1,nchar(basename(path(bam_input)[i]))-4)
               !=sample_names[i,1]){
                message("Bam file name does not match sample name (",
                        substr(basename(path(bam_input)[i]),1,nchar(basename(path(bam_input)[i]))-4),
                        "!=",sample_names[i,1],")")
            }
            if(file.exists(path(bam_input[i]))==FALSE){
                message("File: ",path(bam_input[i])," is missing")
            }
            if(file.exists(paste(substr(path(bam_input)[i],1,nchar(path(bam_input)[i])-1),"i",sep=""))==FALSE){
                message("File: ",
                        paste(substr(path(bam_input)[i],1,nchar(path(bam_input)[i])-1),"i",sep=""),
                        " is missing")
            }
        }
    }
    return()
}

checkVcfInput<-function(sample_names,vcf_input){
    if(is.character(vcf_input)==TRUE){
        if(substr(vcf_input,nchar(vcf_input)-3,nchar(vcf_input))!=".vcf"){
            for(i in 1:length(sample_names[,1])){
                if(file.exists(paste(vcf_input,"/",sample_names[i,1],".vcf",sep=""))==FALSE){
                    message("File: ",paste(vcf_input,"/",sample_names[i,1],".vcf",sep="")," is missing")
                }
            }
            return()
        }
        if(substr(vcf_input,nchar(vcf_input)-3,nchar(vcf_input))==".vcf"){
            vcf<-readVcf(vcf_input,"hg19")
            hdr<-header(vcf)
            for(i in 1:length(sample_names[,1])){
                if(sample_names[i,1]!=samples(hdr)[i]){
                    message("VCF file for sample ",sample_names[i,1],
                            " may be missing, or defined order of sample names 
                            does not match order of samples in the vcf file.")
                }
            }
            return()
        }
    }
    if(is.character(vcf_input)==FALSE){
        for(i in 1:length(sample_names[,1])){
            if(substr(basename(path(vcf_input)[i]),1,nchar(basename(path(vcf_input)[i]))-4)
               !=sample_names[i,1]){
                message("Vcf file name does not match sample name (",
                        substr(basename(path(vcf_input)[i]),1,nchar(basename(path(vcf_input)[i]))-4),
                        "!=",sample_names[i,1],")")
            }
            if(file.exists(path(vcf_input[i]))==FALSE){
                message("File: ",path(vcf_input[i])," is missing")
            }
        }
    }
    return()
}

#Pipeline
analyzeBases<-function(sample_names,
                       bam_input,
                       target_regions,
                       vcf_input,
                       output,
                       output_pictures,
                       known_file,
                       genome,
                       MQ_threshold,
                       BQ_threshold,
                       frequency_threshold,
                       qual_lower_bound,
                       qual_upper_bound,
                       marks,
                       relative,
                       per_sample){
    samples<-read.table(sample_names)
    targetRegions<-read.table(target_regions)
    dbsnp_file<-known_file
    target<-determineTarget(targetRegions)
    checkBamInput(samples,bam_input)
    bq<-analyzeReads(samples,target,bam_input,output,MQ_threshold)
    if(nchar(vcf_input)==0){
        frequency<-frequencyAnalysis(samples,output,bq,BQ_threshold,genome)
        calling<-reportVariants(samples,output,frequency,frequency_threshold,
                                target)
        if(per_sample==TRUE&&relative==FALSE){
            plotCountsPerSample(samples,output_pictures,calling,dbsnp_file,
                                frequency,qual_lower_bound,qual_upper_bound,
                                marks)
        }
        if(per_sample==TRUE&&relative==TRUE){
            plotCountsPerSampleRelative(samples,output_pictures,calling,
                                        dbsnp_file,frequency,qual_lower_bound,
                                        qual_upper_bound,marks)
        }
        if(per_sample==FALSE&&relative==FALSE){
            plotCountsPerTarget(samples,output_pictures,calling,dbsnp_file,
                                frequency,qual_lower_bound,qual_upper_bound,
                                marks)
        }
        if(per_sample==FALSE&&relative==TRUE){
            plotCountsPerTargetRelative(samples,output_pictures,calling,
                                        dbsnp_file,frequency,qual_lower_bound,
                                        qual_upper_bound,marks)
        }
    }
    if(nchar(vcf_input)>0){
        checkVcfInput(samples,vcf_input)
        frequency<-frequencyAnalysisVcf(samples,vcf_input,output,bq,
                                        BQ_threshold,genome)
        calling<-reportVariantsVcf(samples,output,frequency,frequency_threshold,
                                   target)
        if(per_sample==TRUE&&relative==FALSE){
            plotCountsPerSampleVcf(samples,output_pictures,calling,dbsnp_file,
                                   frequency,qual_lower_bound,qual_upper_bound,
                                   marks)
        }
        if(per_sample==TRUE&&relative==TRUE){
            plotCountsPerSampleRelativeVcf(samples,output_pictures,calling,
                                           dbsnp_file,frequency,
                                           qual_lower_bound,qual_upper_bound,
                                           marks)
        }
        if(per_sample==FALSE&&relative==FALSE){
            plotCountsPerTargetVcf(samples,output_pictures,calling,dbsnp_file,
                                   frequency,qual_lower_bound,qual_upper_bound,
                                   marks)
        }
        if(per_sample==FALSE&&relative==TRUE){
            plotCountsPerTargetRelativeVcf(samples,output_pictures,calling,
                                           dbsnp_file,frequency,
                                           qual_lower_bound,qual_upper_bound,
                                           marks)
        }
    }
    output<-list(bases=bq[[1]],quality=bq[[2]],frequency=frequency,
                 calling=calling)
    return(output)
}

analyzeBasesPlotOnly<-function(sample_names,
                               vcf_input,
                               output,
                               known_file,
                               output_list,
                               qual_lower_bound,
                               qual_upper_bound,
                               marks,
                               relative,
                               per_sample){
    samples<-read.table(sample_names)
    dbsnp_file<-known_file
    calling<-output_list[[4]]
    frequency<-output_list[[3]]
    if(nchar(vcf_input)==0){
        if(per_sample==TRUE&&relative==FALSE){
            plotCountsPerSample(samples,output,calling,dbsnp_file,frequency,
                                qual_lower_bound,qual_upper_bound,marks)
        }
        if(per_sample==TRUE&&relative==TRUE){
            plotCountsPerSampleRelative(samples,output,calling,dbsnp_file,
                                        frequency,qual_lower_bound,
                                        qual_upper_bound,marks)
        }
        if(per_sample==FALSE&&relative==FALSE){
            plotCountsPerTarget(samples,output,calling,dbsnp_file,frequency,
                                qual_lower_bound,qual_upper_bound,marks)
        }
        if(per_sample==FALSE&&relative==TRUE){
            plotCountsPerTargetRelative(samples,output,calling,dbsnp_file,
                                        frequency,qual_lower_bound,
                                        qual_upper_bound,marks)
        }
    }
    if(nchar(vcf_input)>0){
        checkVcfInput(samples,vcf_input)
        if(per_sample==TRUE&&relative==FALSE){
            plotCountsPerSampleVcf(samples,output,calling,dbsnp_file,frequency,
                                   qual_lower_bound,qual_upper_bound,marks)
        }
        if(per_sample==TRUE&&relative==TRUE){
            plotCountsPerSampleRelativeVcf(samples,output,calling,dbsnp_file,
                                           frequency,qual_lower_bound,
                                           qual_upper_bound,marks)
        }
        if(per_sample==FALSE&&relative==FALSE){
            plotCountsPerTargetVcf(samples,output,calling,dbsnp_file,frequency,
                                   qual_lower_bound,qual_upper_bound,marks)
        }
        if(per_sample==FALSE&&relative==TRUE){
            plotCountsPerTargetRelativeVcf(samples,output,calling,dbsnp_file,
                                           frequency,qual_lower_bound,
                                           qual_upper_bound,marks)
        }
    }
    return()
}
