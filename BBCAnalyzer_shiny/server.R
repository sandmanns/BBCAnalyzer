library(shiny)
shinyServer(function(input, output, session) {
    library('SummarizedExperiment')
    library('VariantAnnotation')
    library('Rsamtools')
    library('grDevices')
    library('GenomicRanges')
    library('IRanges')
    library('Biostrings')
    library('png')
    
    checkBamInput<-function(sample_names,bam_input){
        if(is.character(bam_input)==TRUE){
            for(i in 1:length(sample_names[,1])){
                if(file.exists(paste(bam_input,"/",sample_names[i,1],".bam",sep=""))==FALSE){
                    shinyjs::html("text", paste("File: ",bam_input,"/",sample_names[i,1],".bam",
                                                " is missing<br>",sep=""), add = TRUE)
                    return(FALSE)
                }
                if(file.exists(paste(bam_input,"/",sample_names[i,1],".bai",sep=""))==FALSE){
                    shinyjs::html("text", paste("File: ",bam_input,"/",sample_names[i,1],".bai",
                                                " is missing<br>",sep=""), add = TRUE)
                    return(FALSE)
                }
            }
            return(TRUE)
        }
    }
    
    checkVcfInput<-function(sample_names,vcf_input){
        for(i in 1:length(sample_names[,1])){
            if(file.exists(paste(vcf_input,"/",sample_names[i,1],".vcf",sep=""))==FALSE){
                shinyjs::html("text", paste("File: ",vcf_input,"/",sample_names[i,1],".vcf",
                                            " is missing<br>",sep=""), add = TRUE)
                return(FALSE)
            }
        }
        return(TRUE)
    }
    
    #Determine single target bases from target regions:
    determineTarget<-function(targetRegions){
        shinyjs::html("text", paste("Determine Target<br>",sep=""), add = TRUE)
        if(length(targetRegions[1,])==3){
            targetLine<-1
            target2Line<-1
            
            target<-data.frame(chr=NA,pos=NA)
            
            while(targetLine<=length(targetRegions[,1])){
                if(targetLine==1){
                    target[targetLine,1]<-as.character(targetRegions[targetLine,1])
                    target[targetLine,2]<-targetRegions[targetLine,2]
                }
                if(targetLine>1){
                    target[target2Line,1]<-as.character(targetRegions[targetLine,1])
                    target[target2Line,2]<-targetRegions[targetLine,2]
                }
                target2Line<-target2Line+1
                
                start<-as.numeric(targetRegions[targetLine,2])+1
                end<-as.numeric(targetRegions[targetLine,3])
                
                if(start-1!=end){
                    while(end>start){
                        target[target2Line,1]<-as.character(targetRegions[targetLine,1])
                        target[target2Line,2]<-start
                        
                        start<-start+1    
                        target2Line<-target2Line+1
                    }  
                    target[target2Line,1]<-as.character(targetRegions[targetLine,1])
                    target[target2Line,2]<-end 
                    target2Line<-target2Line+1
                }
                targetLine<-targetLine+1
            }
        }
        if(length(targetRegions[1,])==2){
            target<-targetRegions
        }
        return(target)
    }
    
    #Analyze reads at every base in target:
    #use mapping quality<60 as threshold for reads with a bad alignment
    analyzeReads<-function(samples,target,directory,directory2,mapq_threshold){
        shinyjs::html("text", paste("Analyze Reads<br>",sep=""), add = TRUE)
        bases<-list()
        quality<-list()
        
        for(j in 1:length(samples[,1])){
            bases[[j]]<-data.frame()
            quality[[j]]<-data.frame()
        }
        counter<-0
        for(j in 1:length(samples[,1])){
            shinyjs::html("text", paste("&nbsp&nbsp&nbspSample ",j," out of ",length(samples[,1]),"<br>",sep=""), add = TRUE)
            bases_insert<-data.frame()
            quality_insert<-data.frame()
            counter<-counter+1
            long_inserts<-0
            for(i in 1:length(target[,1])){
                shinyjs::html("text", paste("&nbsp&nbsp&nbsp&nbsp&nbsp&nbspPosition ",i," out of ",length(target[,1]),"<br>",sep=""), add = TRUE)
                chr<-as.character(target[i,1])
                start<-target[i,2]
                
                what<-c("rname","pos","strand","seq","cigar","qual","mapq")
                which<-GRanges(seqnames=chr,ranges=IRanges(c(start),c(start)))
                param<-ScanBamParam(what=what,which=which)
                
                if(is.character(directory)==TRUE){
                    bamFile<-file.path(paste(directory,"/",samples[j,1],".bam",sep=""))
                    bam<-scanBam(bamFile,
                                 index=file.path(paste(directory,"/",sep=""),samples[j,1]),
                                 param=param)
                }
                if(is.character(directory)==FALSE){
                    bam<-scanBam(path(directory)[j],
                                 index=substr(path(directory)[j],1,nchar(path(directory)[j])-4),
                                 param=param)
                }            
                bam_info_rname<-bam[[1]]$rname
                bam_info_pos<-bam[[1]]$pos
                bam_info_strand<-bam[[1]]$strand
                bam_info_seq<-bam[[1]]$seq
                bam_info_cigar<-bam[[1]]$cigar
                bam_info_qual<-PhredQuality(bam[[1]]$qual)
                bam_info_mapq<-bam[[1]]$mapq
                count_long_inserts<-1
                
                if(length(bam_info_rname)>0){  
                    progress_analysis<-seq(0,length(bam_info_rname),length.out=101)
                    progress_analysis<-round(progress_analysis[2:101])
                    shinyjs::html("text", paste("&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp|",sep=""), add = TRUE)
                    for(k in 1:length(bam_info_rname)){ 
                        if(k==progress_analysis[1]){
                            shinyjs::html("text", paste("*",sep=""), add = TRUE)
                            if((length(progress_analysis)%%10)==0&&length(progress_analysis)!=100){
                                shinyjs::html("text", paste("|",sep=""), add = TRUE)
                            }
                            if(length(progress_analysis)==1){
                                shinyjs::html("text", paste("|<br>",sep=""), add = TRUE)
                            }
                            progress_analysis<-progress_analysis[2:length(progress_analysis)]
                        }
                        dif<-0
                        dif2<-0
                        match<-FALSE
                        found_snp<-FALSE
                        found_indel<-FALSE
                        del_bam<-FALSE
                        mapq<-TRUE
                        insert_length<-1
                        
                        if(bam_info_mapq[k]<mapq_threshold){
                            mapq<-FALSE
                        }
                        if(vcountPattern("I",bam_info_cigar[k])==0&&vcountPattern("D",bam_info_cigar[k])==0&&vcountPattern("N",bam_info_cigar[k])==0&&vcountPattern("S",bam_info_cigar[k])==0&&vcountPattern("H",bam_info_cigar[k])==0&&vcountPattern("P",bam_info_cigar[k])==0&&vcountPattern("X",bam_info_cigar[k])==0){
                            match<-TRUE        
                        }
                        
                        if(vcountPattern("I",bam_info_cigar[k])>0||vcountPattern("D",bam_info_cigar[k])>0||vcountPattern("N",bam_info_cigar[k])>0||vcountPattern("S",bam_info_cigar[k])>0||vcountPattern("H",bam_info_cigar[k])>0||vcountPattern("P",bam_info_cigar[k])>0||vcountPattern("X",bam_info_cigar[k])>0){
                            pointer<-bam_info_pos[k]        
                            
                            while(pointer<=start&&found_snp==FALSE&&found_indel==FALSE&&nchar(bam_info_cigar[k])>=2){
                                shifted<-FALSE
                                #1-9 bases
                                #"M"
                                if(pointer<=start&&found_snp==FALSE&&found_indel==FALSE&&nchar(bam_info_cigar[k])>=2&&shifted==FALSE&&vcountPattern(substr(bam_info_cigar[k],2,2),"1,2,3,4,5,6,7,8,9,0")==0&&vcountPattern(substr(bam_info_cigar[k],2,2),"M")>0){
                                    add<-as.integer(substr(bam_info_cigar[k],1,1))
                                    pointer<-pointer+add
                                    bam_info_cigar[k]<-substr(bam_info_cigar[k],3,1000)
                                    shifted<-TRUE
                                    if(pointer>start){
                                        found_snp=TRUE
                                    }
                                }
                                #"I"
                                if(pointer<=start&&found_snp==FALSE&&found_indel==FALSE&&nchar(bam_info_cigar[k])>=2&&shifted==FALSE&&vcountPattern(substr(bam_info_cigar[k],2,2),"1,2,3,4,5,6,7,8,9,0")==0&&vcountPattern(substr(bam_info_cigar[k],2,2),"I")>0){
                                    add<-as.integer(substr(bam_info_cigar[k],1,1))
                                    dif<-dif+add
                                    bam_info_cigar[k]<-substr(bam_info_cigar[k],3,1000)
                                    shifted<-TRUE
                                    if((pointer)==start){
                                        insert_length<-pointer+add-start
                                        dif<-dif-add
                                        found_indel<-TRUE
                                    }
                                }
                                #"D"
                                if(pointer<=start&&found_snp==FALSE&&found_indel==FALSE&&nchar(bam_info_cigar[k])>=2&&shifted==FALSE&&vcountPattern(substr(bam_info_cigar[k],2,2),"1,2,3,4,5,6,7,8,9,0")==0&&vcountPattern(substr(bam_info_cigar[k],2,2),"D")>0){
                                    add<-as.integer(substr(bam_info_cigar[k],1,1))
                                    dif<-dif-add
                                    dif2<-dif2-add
                                    pointer<-pointer+add
                                    bam_info_cigar[k]<-substr(bam_info_cigar[k],3,1000)
                                    shifted<-TRUE
                                    if(pointer>start){
                                        found_indel<-TRUE
                                        del_bam<-TRUE
                                    }
                                }
                                #"S"
                                if(pointer<=start&&found_snp==FALSE&&found_indel==FALSE&&nchar(bam_info_cigar[k])>=2&&shifted==FALSE&&vcountPattern(substr(bam_info_cigar[k],2,2),"1,2,3,4,5,6,7,8,9,0")==0&&vcountPattern(substr(bam_info_cigar[k],2,2),"S")==1){
                                    add<-as.integer(substr(bam_info_cigar[k],1,1))
                                    dif<-dif+add
                                    dif2<-dif2+add
                                    pointer<-pointer+add
                                    bam_info_cigar[k]<-substr(bam_info_cigar[k],3,1000)
                                    shifted<-TRUE
                                }
                                
                                #10-99 bases
                                #"M"
                                if(pointer<=start&&found_snp==FALSE&&found_indel==FALSE&&nchar(bam_info_cigar[k])>=3&&shifted==FALSE&&vcountPattern(substr(bam_info_cigar[k],3,3),"1,2,3,4,5,6,7,8,9,0")==0&&vcountPattern(substr(bam_info_cigar[k],3,3),"M")>0){
                                    add<-as.integer(substr(bam_info_cigar[k],1,2))
                                    pointer<-pointer+add
                                    bam_info_cigar[k]<-substr(bam_info_cigar[k],4,1000)
                                    shifted<-TRUE
                                    if(pointer>start){
                                        found_snp=TRUE
                                    }
                                }
                                #"I"
                                if(pointer<=start&&found_snp==FALSE&&found_indel==FALSE&&nchar(bam_info_cigar[k])>=3&&shifted==FALSE&&vcountPattern(substr(bam_info_cigar[k],3,3),"1,2,3,4,5,6,7,8,9,0")==0&&vcountPattern(substr(bam_info_cigar[k],3,3),"I")>0){
                                    add<-as.integer(substr(bam_info_cigar[k],1,2))
                                    dif<-dif+add
                                    bam_info_cigar[k]<-substr(bam_info_cigar[k],4,1000)
                                    shifted<-TRUE
                                    if((pointer)==start){
                                        insert_length<-pointer+add-start
                                        dif<-dif-add
                                        found_indel<-TRUE
                                    }
                                }
                                #"D"
                                if(pointer<=start&&found_snp==FALSE&&found_indel==FALSE&&nchar(bam_info_cigar[k])>=3&&shifted==FALSE&&vcountPattern(substr(bam_info_cigar[k],3,3),"1,2,3,4,5,6,7,8,9,0")==0&&vcountPattern(substr(bam_info_cigar[k],3,3),"D")>0){
                                    add<-as.integer(substr(bam_info_cigar[k],1,2))
                                    dif<-dif-add
                                    dif2<-dif2-add
                                    pointer<-pointer+add
                                    bam_info_cigar[k]<-substr(bam_info_cigar[k],4,1000)
                                    shifted<-TRUE
                                    if(pointer>start){
                                        found_indel<-TRUE
                                        del_bam<-TRUE
                                    }
                                }
                                #"S"
                                if(pointer<=start&&found_snp==FALSE&&found_indel==FALSE&&nchar(bam_info_cigar[k])>=3&&shifted==FALSE&&vcountPattern(substr(bam_info_cigar[k],3,3),"1,2,3,4,5,6,7,8,9,0")==0&&vcountPattern(substr(bam_info_cigar[k],3,3),"S")==1){
                                    add<-as.integer(substr(bam_info_cigar[k],1,2))
                                    dif<-dif+add
                                    dif2<-dif2+add
                                    pointer<-pointer+add
                                    bam_info_cigar[k]<-substr(bam_info_cigar[k],4,1000)
                                    shifted<-TRUE
                                }
                                
                                #100-999 bases
                                #"M"
                                if(pointer<=start&&found_snp==FALSE&&found_indel==FALSE&&nchar(bam_info_cigar[k])>=4&&shifted==FALSE&&vcountPattern(substr(bam_info_cigar[k],4,4),"1,2,3,4,5,6,7,8,9,0")==0&&vcountPattern(substr(bam_info_cigar[k],4,4),"M")>0){
                                    add<-as.integer(substr(bam_info_cigar[k],1,3))
                                    pointer<-pointer+add
                                    bam_info_cigar[k]<-substr(bam_info_cigar[k],5,1000)
                                    shifted<-TRUE
                                    if(pointer>start){
                                        found_snp=TRUE
                                    }
                                }
                                #"I"
                                if(pointer<=start&&found_snp==FALSE&&found_indel==FALSE&&nchar(bam_info_cigar[k])>=4&&shifted==FALSE&&vcountPattern(substr(bam_info_cigar[k],4,4),"1,2,3,4,5,6,7,8,9,0")==0&&vcountPattern(substr(bam_info_cigar[k],4,4),"I")>0){
                                    add<-as.integer(substr(bam_info_cigar[k],1,3))
                                    dif<-dif+add
                                    bam_info_cigar[k]<-substr(bam_info_cigar[k],5,1000)
                                    shifted<-TRUE
                                    if((pointer)==start){
                                        insert_length<-pointer+add-start
                                        dif<-dif-add
                                        found_indel<-TRUE
                                    }
                                }
                                #"D"
                                if(pointer<=start&&found_snp==FALSE&&found_indel==FALSE&&nchar(bam_info_cigar[k])>=4&&shifted==FALSE&&vcountPattern(substr(bam_info_cigar[k],4,4),"1,2,3,4,5,6,7,8,9,0")==0&&vcountPattern(substr(bam_info_cigar[k],4,4),"D")>0){
                                    add<-as.integer(substr(bam_info_cigar[k],1,3))
                                    dif<-dif-add
                                    dif2<-dif2-add
                                    pointer<-pointer+add
                                    bam_info_cigar[k]<-substr(bam_info_cigar[k],5,1000)
                                    shifted<-TRUE
                                    if(pointer>start){
                                        found_indel<-TRUE
                                        del_bam<-TRUE
                                    }
                                }
                                #"S"
                                if(pointer<=start&&found_snp==FALSE&&found_indel==FALSE&&nchar(bam_info_cigar[k])>=4&&shifted==FALSE&&vcountPattern(substr(bam_info_cigar[k],4,4),"1,2,3,4,5,6,7,8,9,0")==0&&vcountPattern(substr(bam_info_cigar[k],4,4),"S")==1){
                                    add<-as.integer(substr(bam_info_cigar[k],1,3))
                                    dif<-dif+add
                                    dif2<-dif2+add
                                    pointer<-pointer+add
                                    bam_info_cigar[k]<-substr(bam_info_cigar[k],5,1000)
                                    shifted<-TRUE
                                }
                                
                                #1000-9999 bases
                                #"M"
                                if(pointer<=start&&found_snp==FALSE&&found_indel==FALSE&&nchar(bam_info_cigar[k])>=5&&shifted==FALSE&&vcountPattern(substr(bam_info_cigar[k],5,5),"1,2,3,4,5,6,7,8,9,0")==0&&vcountPattern(substr(bam_info_cigar[k],5,5),"M")>0){
                                    add<-as.integer(substr(bam_info_cigar[k],1,4))
                                    pointer<-pointer+add
                                    bam_info_cigar[k]<-substr(bam_info_cigar[k],6,1000)
                                    shifted<-TRUE
                                    if(pointer>start){
                                        found_snp=TRUE
                                    }
                                }
                                #"I"
                                if(pointer<=start&&found_snp==FALSE&&found_indel==FALSE&&nchar(bam_info_cigar[k])>=5&&shifted==FALSE&&vcountPattern(substr(bam_info_cigar[k],5,5),"1,2,3,4,5,6,7,8,9,0")==0&&vcountPattern(substr(bam_info_cigar[k],5,5),"I")>0){
                                    add<-as.integer(substr(bam_info_cigar[k],1,4))
                                    dif<-dif+add
                                    bam_info_cigar[k]<-substr(bam_info_cigar[k],6,1000)
                                    shifted<-TRUE
                                    if((pointer)==start){
                                        insert_length<-pointer+add-start
                                        dif<-dif-add
                                        found_indel<-TRUE
                                    }
                                }
                                #"D"
                                if(pointer<=start&&found_snp==FALSE&&found_indel==FALSE&&nchar(bam_info_cigar[k])>=5&&shifted==FALSE&&vcountPattern(substr(bam_info_cigar[k],5,5),"1,2,3,4,5,6,7,8,9,0")==0&&vcountPattern(substr(bam_info_cigar[k],5,5),"D")>0){
                                    add<-as.integer(substr(bam_info_cigar[k],1,4))
                                    dif<-dif-add
                                    dif2<-dif2-add
                                    pointer<-pointer+add
                                    bam_info_cigar[k]<-substr(bam_info_cigar[k],6,1000)
                                    shifted<-TRUE
                                    if(pointer>start){
                                        found_indel<-TRUE
                                        del_bam<-TRUE
                                    }
                                }
                                #"S"
                                if(pointer<=start&&found_snp==FALSE&&found_indel==FALSE&&nchar(bam_info_cigar[k])>=5&&shifted==FALSE&&vcountPattern(substr(bam_info_cigar[k],5,5),"1,2,3,4,5,6,7,8,9,0")==0&&vcountPattern(substr(bam_info_cigar[k],5,5),"S")==1){
                                    add<-as.integer(substr(bam_info_cigar[k],1,4))
                                    dif<-dif+add
                                    dif2<-dif2+add
                                    pointer<-pointer+add
                                    bam_info_cigar[k]<-substr(bam_info_cigar[k],6,1000)
                                    shifted<-TRUE
                                }
                            }
                        }
                        
                        qual<-encoding(bam_info_qual)[unlist(strsplit(as.character(bam_info_qual[[k]]),split=""))]+33
                        
                        pos<-start+dif-bam_info_pos[k]+1
                        pos_qual<-start+dif-bam_info_pos[k]+1
                        
                        if(mapq==FALSE){
                            bases[[counter]][k,i+long_inserts]<-"MQ"
                            quality[[counter]][k,i+long_inserts]<--2
                        }
                        if(found_snp==TRUE||match==TRUE){
                            bases[[counter]][k,i+long_inserts]<-substr(bam_info_seq[k],pos,pos)
                            quality[[counter]][k,i+long_inserts]<-qual[pos_qual]
                        }
                        if(found_indel==TRUE){
                            if(del_bam==TRUE){
                                bases[[counter]][k,i+long_inserts]<-"D"
                                quality[[counter]][k,i+long_inserts]<--1
                            }
                            if(del_bam==FALSE){
                                bases[[counter]][k,i+long_inserts]<-paste("I, ",substr(bam_info_seq[k],pos,pos),sep="")
                                quality[[counter]][k,i+long_inserts]<-qual[pos_qual]
                                if(insert_length>1){
                                    for(m in 2:insert_length){
                                        bases_insert[count_long_inserts,m-1]<-paste("I, ",substr(bam_info_seq[k],pos+m-1,pos+m-1),sep="")
                                        quality_insert[count_long_inserts,m-1]<-qual[pos_qual+m-1]
                                    }
                                    count_long_inserts<-count_long_inserts+1
                                }
                            }
                        }
                    }
                }
                
                if(length(bam_info_rname)==0){
                    bases[[counter]][1,i+long_inserts]<-"NotCovered"
                    quality[[counter]][1,i+long_inserts]<-"NotCovered"
                }
                
                names(bases[[counter]])[i+long_inserts]<-paste(target[i,1],";",target[i,2],sep="")
                names(quality[[counter]])[i+long_inserts]<-paste(target[i,1],";",target[i,2],sep="")
                
                if(count_long_inserts>1){
                    for(n in 1:length(bases_insert[1,])){
                        long_inserts<-long_inserts+1
                        counter_inserts_bq<-1
                        for(m in 1:length(bases[[counter]][,i+long_inserts-n])){
                            if(!is.na(bases[[counter]][m,i+long_inserts-n])&&substr(bases[[counter]][m,i+long_inserts-n],1,1)=="I"){
                                bases[[counter]][m,i+long_inserts]<-bases_insert[counter_inserts_bq,n]
                                quality[[counter]][m,i+long_inserts]<-quality_insert[counter_inserts_bq,n]
                            }
                            if(!is.na(bases[[counter]][m,i+long_inserts-n])&&substr(bases[[counter]][m,i+long_inserts-n],1,1)!="I"){
                                bases[[counter]][m,i+long_inserts]<-bases[[counter]][m,i+long_inserts-n]
                                quality[[counter]][m,i+long_inserts]<-quality[[counter]][m,i+long_inserts-n]
                            }
                        }
                        names(bases[[counter]])[i+long_inserts]<-paste(target[i,1],";",target[i,2],sep="")
                        names(quality[[counter]])[i+long_inserts]<-paste(target[i,1],";",target[i,2],sep="")
                    }
                }
            }
            write.table(bases[[counter]],paste(directory2,"/",samples[counter,1],".bases.txt",sep=""),row.names=FALSE,sep="\t",quote=FALSE)
            write.table(quality[[counter]],paste(directory2,"/",samples[counter,1],".quality.txt",sep=""),row.names=FALSE,sep="\t",quote=FALSE)
        }
        bases_quality<-list(bases,quality) 
        return(bases_quality)
    }
    
    processVcfFile<-function(vcf_input,j){
        vcf_temp<-read.table(vcf_input,header=FALSE)
        vcf_temp2<-vcf_temp[,c(1:9,9+j)]
        vcf<-vcf_temp2[vcf_temp2[,10]!="./.",]
        return(vcf)
    }
    
    #Frequency analysis:
    #use (no) quality threshold
    frequencyAnalysis<-function(samples,
                                directory2,
                                bases_quality,
                                quality_threshold,
                                genome){
        shinyjs::html("text", paste("Analyze Frequency<br>",sep=""), add = TRUE)
        frequency_weight<-list()
        bases<-bases_quality[[1]]
        quality<-bases_quality[[2]]
        
        for(j in 1:length(samples[,1])){
            frequency_weight[[j]]<-data.frame(chr=NA,pos=NA,ref=NA,a=NA,a_qual=NA,
                                              c=NA,c_qual=NA,g=NA,g_qual=NA,t=NA,
                                              t_qual=NA,Del=NA,Ins=NA,Ins_a=NA,
                                              Ins_a_qual=NA,Ins_c=NA,Ins_c_qual=NA,
                                              Ins_g=NA,Ins_g_qual=NA,Ins_t=NA,
                                              Ins_t_qual=NA,excluded=NA)
        }
        for(j in 1:length(samples[,1])){
            shinyjs::html("text", paste("&nbsp&nbsp&nbspSample ",j,"<br>",sep=""), add = TRUE)    
            for(i in 1:length(bases_quality[[1]][[j]][1,])){
                if(vcountPattern(".",names(bases_quality[[1]][[j]])[i])==1){
                    start<-as.numeric(strsplit(strsplit(names(bases_quality[[1]][[j]])[i],split="\\.")[[1]][1],split=";")[[1]][2])
                }
                if(vcountPattern(".",names(bases_quality[[1]][[j]])[i])==0){
                    start<-as.numeric(strsplit(names(bases_quality[[1]][[j]])[i],split=";")[[1]][2])
                }
                chr<-paste("chr",as.character(strsplit(names(bases_quality[[1]][[j]])[i],split=";")[[1]][1]),sep="")
                ref<-genome[[chr]][start]
                
                frequency_weight[[j]][i,1]<-chr
                frequency_weight[[j]][i,2]<-start
                frequency_weight[[j]][i,3]<-as.character(ref)
                
                a_qual<-c()
                c_qual<-c()
                g_qual<-c()
                t_qual<-c()
                ins_a_qual<-c()
                ins_c_qual<-c()
                ins_g_qual<-c()
                ins_t_qual<-c()
                
                if(bases[[j]][1,i]!="NotCovered"){
                    for(k in 1:sum(!is.na(bases[[j]][,i]))){
                        if(!is.na(bases[[j]][k,i])&&quality[[j]][k,i]>quality_threshold){
                            if(bases[[j]][k,i]=="A"){
                                frequency_weight[[j]][i,4]<-sum(frequency_weight[[j]][i,4],1,na.rm=TRUE)
                                a_qual[frequency_weight[[j]][i,4]]<-quality[[j]][k,i]
                            }
                            if(bases[[j]][k,i]=="C"){
                                frequency_weight[[j]][i,6]<-sum(frequency_weight[[j]][i,6],1,na.rm=TRUE)
                                c_qual[frequency_weight[[j]][i,6]]<-quality[[j]][k,i]
                            }
                            if(bases[[j]][k,i]=="G"){
                                frequency_weight[[j]][i,8]<-sum(frequency_weight[[j]][i,8],1,na.rm=TRUE)
                                g_qual[frequency_weight[[j]][i,8]]<-quality[[j]][k,i]
                            }
                            if(bases[[j]][k,i]=="T"){
                                frequency_weight[[j]][i,10]<-sum(frequency_weight[[j]][i,10],1,na.rm=TRUE)
                                t_qual[frequency_weight[[j]][i,10]]<-quality[[j]][k,i]
                            }
                            if(substring(bases[[j]][k,i],1,1)=="I"){
                                frequency_weight[[j]][i,13]<-sum(frequency_weight[[j]][i,13],1,na.rm=TRUE)
                                if(substring(bases[[j]][k,i],4,4)=="A"){
                                    frequency_weight[[j]][i,4]<-sum(frequency_weight[[j]][i,4],1,na.rm=TRUE)
                                    a_qual[frequency_weight[[j]][i,4]]<-quality[[j]][k,i]
                                    frequency_weight[[j]][i,14]<-sum(frequency_weight[[j]][i,14],1,na.rm=TRUE)
                                    ins_a_qual[frequency_weight[[j]][i,14]]<-quality[[j]][k,i]
                                }
                                if(substring(bases[[j]][k,i],4,4)=="C"){
                                    frequency_weight[[j]][i,6]<-sum(frequency_weight[[j]][i,6],1,na.rm=TRUE)
                                    c_qual[frequency_weight[[j]][i,6]]<-quality[[j]][k,i]
                                    frequency_weight[[j]][i,16]<-sum(frequency_weight[[j]][i,16],1,na.rm=TRUE)
                                    ins_c_qual[frequency_weight[[j]][i,16]]<-quality[[j]][k,i]
                                }
                                if(substring(bases[[j]][k,i],4,4)=="G"){
                                    frequency_weight[[j]][i,8]<-sum(frequency_weight[[j]][i,8],1,na.rm=TRUE)
                                    g_qual[frequency_weight[[j]][i,8]]<-quality[[j]][k,i]
                                    frequency_weight[[j]][i,18]<-sum(frequency_weight[[j]][i,18],1,na.rm=TRUE)
                                    ins_g_qual[frequency_weight[[j]][i,18]]<-quality[[j]][k,i]
                                }
                                if(substring(bases[[j]][k,i],4,4)=="T"){
                                    frequency_weight[[j]][i,10]<-sum(frequency_weight[[j]][i,10],1,na.rm=TRUE)
                                    t_qual[frequency_weight[[j]][i,10]]<-quality[[j]][k,i]
                                    frequency_weight[[j]][i,20]<-sum(frequency_weight[[j]][i,20],1,na.rm=TRUE)
                                    ins_t_qual[frequency_weight[[j]][i,20]]<-quality[[j]][k,i]
                                }
                            }
                        }
                        if(!is.na(bases[[j]][k,i])&&quality[[j]][k,i]<=quality_threshold&&quality[[j]][k,i]!=-1){
                            frequency_weight[[j]][i,22]<-sum(frequency_weight[[j]][i,22],1,na.rm=TRUE)
                        }
                        if(!is.na(bases[[j]][k,i])&&quality[[j]][k,i]==-1&&bases[[j]][k,i]=="D"){
                            frequency_weight[[j]][i,12]<-sum(frequency_weight[[j]][i,12],1,na.rm=TRUE)
                        }
                    }  
                    
                    if(!is.na(frequency_weight[[j]][i,4])){
                        frequency_weight[[j]][i,5]<-mean(a_qual)
                    }
                    if(!is.na(frequency_weight[[j]][i,6])){
                        frequency_weight[[j]][i,7]<-mean(c_qual)
                    }
                    if(!is.na(frequency_weight[[j]][i,8])){
                        frequency_weight[[j]][i,9]<-mean(g_qual)
                    }
                    if(!is.na(frequency_weight[[j]][i,10])){
                        frequency_weight[[j]][i,11]<-mean(t_qual)
                    }
                    if(!is.na(frequency_weight[[j]][i,14])){
                        frequency_weight[[j]][i,15]<-mean(ins_a_qual)
                    }
                    if(!is.na(frequency_weight[[j]][i,16])){
                        frequency_weight[[j]][i,17]<-mean(ins_c_qual)
                    }
                    if(!is.na(frequency_weight[[j]][i,18])){
                        frequency_weight[[j]][i,19]<-mean(ins_g_qual)
                    }
                    if(!is.na(frequency_weight[[j]][i,20])){
                        frequency_weight[[j]][i,21]<-mean(ins_t_qual)
                    }
                }
            }
            write.table(frequency_weight[[j]],
                        paste(directory2,"/",samples[j,1],".frequency.txt",sep=""),
                        row.names=FALSE,sep="\t",quote=FALSE)
        }
        return(frequency_weight)
    }
    
    
    #Frequency analysis:
    #use (no) quality threshold
    frequencyAnalysisVcf<-function(samples,
                                   vcf_input,
                                   directory2,
                                   bases_quality,
                                   quality_threshold,
                                   genome){
        shinyjs::html("text", paste("Analyze Frequency<br>",sep=""), add = TRUE)  
        frequency_weight<-list()
        bases<-bases_quality[[1]]
        quality<-bases_quality[[2]]
        
        for(j in 1:length(samples[,1])){
            frequency_weight[[j]]<-data.frame(chr=NA,pos=NA,ref=NA,a=NA,a_qual=NA,
                                              c=NA,c_qual=NA,g=NA,g_qual=NA,t=NA,
                                              t_qual=NA,Del=NA,Ins=NA,Ins_a=NA,
                                              Ins_a_qual=NA,Ins_c=NA,Ins_c_qual=NA,
                                              Ins_g=NA,Ins_g_qual=NA,Ins_t=NA,
                                              Ins_t_qual=NA,excluded=NA,
                                              vcf_call_1=NA,vcf_call_2=NA,GT=NA,
                                              Insert=NA)
        }
        for(j in 1:length(samples[,1])){
            shinyjs::html("text", paste("&nbsp&nbsp&nbspSample ",j,"<br>",sep=""), add = TRUE)  
            for(i in 1:length(bases_quality[[1]][[j]][1,])){
                if(vcountPattern(".",names(bases_quality[[1]][[j]])[i])==1){
                    start<-as.numeric(strsplit(strsplit(names(bases_quality[[1]][[j]])[i],split="\\.")[[1]][1],split=";")[[1]][2])
                }
                if(vcountPattern(".",names(bases_quality[[1]][[j]])[i])==0){
                    start<-as.numeric(strsplit(names(bases_quality[[1]][[j]])[i],split=";")[[1]][2])
                }
                chr<-paste("chr",as.character(strsplit(names(bases_quality[[1]][[j]])[i],split=";")[[1]][1]),sep="")
                ref<-genome[[chr]][start]
                
                frequency_weight[[j]][i,1]<-chr
                frequency_weight[[j]][i,2]<-start
                frequency_weight[[j]][i,3]<-as.character(ref)
                
                a_qual<-c()
                c_qual<-c()
                g_qual<-c()
                t_qual<-c()
                ins_a_qual<-c()
                ins_c_qual<-c()
                ins_g_qual<-c()
                ins_t_qual<-c()
                
                if(bases[[j]][1,i]!="NotCovered"){
                    for(k in 1:sum(!is.na(bases[[j]][,i]))){
                        if(!is.na(bases[[j]][k,i])&&quality[[j]][k,i]>quality_threshold){
                            if(bases[[j]][k,i]=="A"){
                                frequency_weight[[j]][i,4]<-sum(frequency_weight[[j]][i,4],1,na.rm=TRUE)
                                a_qual[frequency_weight[[j]][i,4]]<-quality[[j]][k,i]
                            }
                            if(bases[[j]][k,i]=="C"){
                                frequency_weight[[j]][i,6]<-sum(frequency_weight[[j]][i,6],1,na.rm=TRUE)
                                c_qual[frequency_weight[[j]][i,6]]<-quality[[j]][k,i]
                            }
                            if(bases[[j]][k,i]=="G"){
                                frequency_weight[[j]][i,8]<-sum(frequency_weight[[j]][i,8],1,na.rm=TRUE)
                                g_qual[frequency_weight[[j]][i,8]]<-quality[[j]][k,i]
                            }
                            if(bases[[j]][k,i]=="T"){
                                frequency_weight[[j]][i,10]<-sum(frequency_weight[[j]][i,10],1,na.rm=TRUE)
                                t_qual[frequency_weight[[j]][i,10]]<-quality[[j]][k,i]
                            }
                            if(substring(bases[[j]][k,i],1,1)=="I"){
                                frequency_weight[[j]][i,13]<-sum(frequency_weight[[j]][i,13],1,na.rm=TRUE)
                                if(substring(bases[[j]][k,i],4,4)=="A"){
                                    frequency_weight[[j]][i,4]<-sum(frequency_weight[[j]][i,4],1,na.rm=TRUE)
                                    a_qual[frequency_weight[[j]][i,4]]<-quality[[j]][k,i]
                                    frequency_weight[[j]][i,14]<-sum(frequency_weight[[j]][i,14],1,na.rm=TRUE)
                                    ins_a_qual[frequency_weight[[j]][i,14]]<-quality[[j]][k,i]
                                }
                                if(substring(bases[[j]][k,i],4,4)=="C"){
                                    frequency_weight[[j]][i,6]<-sum(frequency_weight[[j]][i,6],1,na.rm=TRUE)
                                    c_qual[frequency_weight[[j]][i,6]]<-quality[[j]][k,i]
                                    frequency_weight[[j]][i,16]<-sum(frequency_weight[[j]][i,16],1,na.rm=TRUE)
                                    ins_c_qual[frequency_weight[[j]][i,16]]<-quality[[j]][k,i]
                                }
                                if(substring(bases[[j]][k,i],4,4)=="G"){
                                    frequency_weight[[j]][i,8]<-sum(frequency_weight[[j]][i,8],1,na.rm=TRUE)
                                    g_qual[frequency_weight[[j]][i,8]]<-quality[[j]][k,i]
                                    frequency_weight[[j]][i,18]<-sum(frequency_weight[[j]][i,18],1,na.rm=TRUE)
                                    ins_g_qual[frequency_weight[[j]][i,18]]<-quality[[j]][k,i]
                                }
                                if(substring(bases[[j]][k,i],4,4)=="T"){
                                    frequency_weight[[j]][i,10]<-sum(frequency_weight[[j]][i,10],1,na.rm=TRUE)
                                    t_qual[frequency_weight[[j]][i,10]]<-quality[[j]][k,i]
                                    frequency_weight[[j]][i,20]<-sum(frequency_weight[[j]][i,20],1,na.rm=TRUE)
                                    ins_t_qual[frequency_weight[[j]][i,20]]<-quality[[j]][k,i]
                                }
                            }
                        }
                        if(!is.na(bases[[j]][k,i])&&quality[[j]][k,i]<=quality_threshold&&quality[[j]][k,i]!=-1){
                            frequency_weight[[j]][i,22]<-sum(frequency_weight[[j]][i,22],1,na.rm=TRUE)
                        }
                        if(!is.na(bases[[j]][k,i])&&quality[[j]][k,i]==-1&&bases[[j]][k,i]=="D"){
                            frequency_weight[[j]][i,12]<-sum(frequency_weight[[j]][i,12],1,na.rm=TRUE)
                        }
                    }  
                    
                    if(!is.na(frequency_weight[[j]][i,4])){
                        frequency_weight[[j]][i,5]<-mean(a_qual)
                    }
                    if(!is.na(frequency_weight[[j]][i,6])){
                        frequency_weight[[j]][i,7]<-mean(c_qual)
                    }
                    if(!is.na(frequency_weight[[j]][i,8])){
                        frequency_weight[[j]][i,9]<-mean(g_qual)
                    }
                    if(!is.na(frequency_weight[[j]][i,10])){
                        frequency_weight[[j]][i,11]<-mean(t_qual)
                    }
                    if(!is.na(frequency_weight[[j]][i,14])){
                        frequency_weight[[j]][i,15]<-mean(ins_a_qual)
                    }
                    if(!is.na(frequency_weight[[j]][i,16])){
                        frequency_weight[[j]][i,17]<-mean(ins_c_qual)
                    }
                    if(!is.na(frequency_weight[[j]][i,18])){
                        frequency_weight[[j]][i,19]<-mean(ins_g_qual)
                    }
                    if(!is.na(frequency_weight[[j]][i,20])){
                        frequency_weight[[j]][i,21]<-mean(ins_t_qual)
                    }
                }
                
                
                vcf<-read.table(paste(vcf_input,"/",samples[j,1],".vcf",sep=""),header=FALSE)
                if(length(vcf[,1])>0){
                    for(k in 1:length(vcf[,1])){
                        del<-FALSE
                        ins<-FALSE
                        del_length<-c()
                        ins_length<-c()
                        #deletion
                        if(nchar(as.character(vcf[k,4]))>1){
                            del<-TRUE
                            if(vcountPattern(",",as.character(vcf[k,5]))==1){
                                del_length[1]<-nchar(as.character(vcf[k,4]))-nchar(strsplit(as.character(vcf[k,5]),split=",")[[1]][1])
                                del_length[2]<-nchar(as.character(vcf[k,4]))-nchar(strsplit(as.character(vcf[k,5]),split=",")[[1]][2])
                            }
                            if(vcountPattern(",",as.character(vcf[k,5]))==0){
                                del_length[1]<-nchar(as.character(vcf[k,4]))-nchar(as.character(vcf[k,5]))
                            }
                        }
                        #insert
                        if((nchar(as.character(vcf[k,5]))>1
                            &&vcountPattern(",",as.character(vcf[k,5]))==0)
                           ||(vcountPattern(",",as.character(vcf[k,5]))==1
                              &&nchar(strsplit(as.character(vcf[k,5]),split=",")[[1]][1])>1
                              &&nchar(strsplit(as.character(vcf[k,5]),split=",")[[1]][2])>1)){
                            ins<-TRUE
                            if(vcountPattern(",",as.character(vcf[k,5]))==1){
                                ins_length[1]<-nchar(as.character(vcf[k,5]))-nchar(strsplit(as.character(vcf[k,4]),split=",")[[1]][1])
                                ins_length[2]<-nchar(as.character(vcf[k,5]))-nchar(strsplit(as.character(vcf[k,4]),split=",")[[1]][2])
                            }
                            if(vcountPattern(",",as.character(vcf[k,5]))==0){
                                ins_length[1]<-nchar(as.character(vcf[k,5]))-nchar(as.character(vcf[k,4]))
                            }
                        }
                        
                        #snv
                        if(k<=length(vcf[,1])&&del==FALSE&&ins==FALSE
                           &&as.character(vcf[k,1])==substr(chr,4,nchar(chr))
                           &&vcf[k,2]==start){
                            if(vcountPattern(as.character(vcf[k,5]),pattern=",")==0){
                                frequency_weight[[j]][i,23]<-as.character(vcf[k,5])
                            }
                            if(vcountPattern(as.character(vcf[k,5]),pattern=",")==1){
                                frequency_weight[[j]][i,23]<-strsplit(as.character(vcf[k,5]),split=",")[[1]][1]
                                frequency_weight[[j]][i,23]<-strsplit(as.character(vcf[k,5]),split=",")[[1]][2]
                            }
                            frequency_weight[[j]][i,25]<-strsplit(as.character(vcf[k,10]),split=":")[[1]][1]
                            k<-length(vcf[,1])+1
                        }   
                        #deletion
                        if(k<length(vcf[,1])&&del==TRUE&&is.na(del_length[2])
                           &&as.character(vcf[k,1])==substr(chr,4,nchar(chr))
                           &&vcf[k,2]<=start-1&&vcf[k,2]>=start-del_length[1]){
                            frequency_weight[[j]][i,23]<-""
                            frequency_weight[[j]][i,25]<-strsplit(as.character(vcf[k,10]),split=":")[[1]][1]
                            k<-length(vcf[,1])+1
                        }
                        if(k<length(vcf[,1])&&del==TRUE
                           &&!is.na(del_length[2])
                           &&as.character(vcf[k,1])==substr(chr,4,nchar(chr))
                           &&vcf[k,2]<=start-1
                           &&(vcf[k,2]>=start-del_length[1]||vcf[k,2]>=start-del_length[2])){
                            frequency_weight[[j]][i,23]<-substr(strsplit(as.character(vcf[k,5]),split=",")[[1]][1],
                                                                nchar(as.character(vcf[k,4]))-(start-vcf[k,2]-1),
                                                                nchar(as.character(vcf[k,4]))-(start-vcf[k,2]-1))
                            frequency_weight[[j]][i,24]<-substr(strsplit(as.character(vcf[k,5]),split=",")[[1]][2],
                                                                nchar(as.character(vcf[k,4]))-(start-vcf[k,2]-1),
                                                                nchar(as.character(vcf[k,4]))-(start-vcf[k,2]-1))
                            frequency_weight[[j]][i,25]<-strsplit(as.character(vcf[k,10]),split=":")[[1]][1]
                            k<-length(vcf[,1])+1
                        }
                        #insert
                        if(k<length(vcf[,1])&&ins==TRUE&&is.na(ins_length[2])
                           &&as.character(vcf[k,1])==substr(chr,4,nchar(chr))
                           &&vcf[k,2]==start-1){
                            frequency_weight[[j]][i,23]<-substr(as.character(vcf[k,5]),2,nchar(as.character(vcf[k,5])))
                            frequency_weight[[j]][i,25]<-strsplit(as.character(vcf[k,10]),split=":")[[1]][1]
                            frequency_weight[[j]][i,26]<-TRUE
                            k<-length(vcf[,1])+1
                        }
                        if(k<length(vcf[,1])&&ins==TRUE&&!is.na(ins_length[2])
                           &&as.character(vcf[k,1])==substr(chr,4,nchar(chr))
                           &&vcf[k,2]==start-1){
                            frequency_weight[[j]][i,23]<-substr(strsplit(as.character(vcf[k,5]),split=",")[[1]][1],
                                                                2,nchar(as.character(vcf[k,5])))
                            frequency_weight[[j]][i,24]<-substr(strsplit(as.character(vcf[k,5]),split=",")[[1]][2],
                                                                2,nchar(as.character(vcf[k,5])))
                            frequency_weight[[j]][i,25]<-strsplit(as.character(vcf[k,10]),split=":")[[1]][1]
                            frequency_weight[[j]][i,26]<-TRUE
                            k<-length(vcf[,1])+1
                        }
                    } 
                } 
            }
            write.table(frequency_weight[[j]],
                        paste(directory2,"/",samples[j,1],".frequency.txt",sep=""),
                        row.names=FALSE,sep="\t",quote=FALSE)
        }
        return(frequency_weight)
    }
    
    #Variant calling:
    #define threshold for frequency of variant
    reportVariants<-function(samples,
                             directory2,
                             frequency_weight,
                             frequency_threshold,
                             target){
        shinyjs::html("text", paste("Report Variants<br>",sep=""), add = TRUE)
        calling<-list()
        for(j in 1:length(samples[,1])){
            calling[[j]]<-data.frame(chr=NA,pos=NA,ref=NA,a_freq=NA,c_freq=NA,
                                     g_freq=NA,t_freq=NA,del_freq=NA,ins_freq=NA,
                                     call1=NA,call2=NA,call3=NA,call4=NA,call5=NA,
                                     call6=NA)
        }
        
        max.variants<-0  
        line<-rep(1,length(samples[,1]))
        start.line<-1
        for(i in 1:length(target[,1])){
            maxvar<-0
            for(j in 1:length(samples[,1])){
                maxvarsample<-0
                while(line[j]<=length(frequency_weight[[j]][,1])
                      &&paste("chr",as.character(target[start.line,1]),sep="")==frequency_weight[[j]][line[j],1]
                      &&target[start.line,2]==frequency_weight[[j]][line[j],2]){
                    maxvarsample<-maxvarsample+1
                    line[j]<-line[j]+1
                }
                if(maxvar<maxvarsample){
                    maxvar<-maxvarsample
                }
            }
            max.variants<-max.variants+maxvar
            start.line<-start.line+1
        } 
        
        for(j in 1:length(samples[,1])){
            shinyjs::html("text", paste("&nbsp&nbsp&nbspSample ",j,"<br>",sep=""), add = TRUE)  
            long_insert<-rep(0,length(samples[,1])) 
            
            for(i in 1:max.variants){
                not_covered_long_insert<-0
                if(i>1){
                    for(k in 1:length(samples[,1])){
                        if(!is.na(frequency_weight[[k]][i-long_insert[k],1])
                           &&frequency_weight[[k]][i-long_insert[k]-1,1]==frequency_weight[[k]][i-long_insert[k],1]
                           &&frequency_weight[[k]][i-long_insert[k]-1,2]==frequency_weight[[k]][i-long_insert[k],2]){
                            not_covered_long_insert<-not_covered_long_insert+1
                        }
                    } 
                    for(k in 1:length(samples[,1])){
                        if(not_covered_long_insert!=0
                           &&not_covered_long_insert!=length(samples[,1])
                           &&!is.na(frequency_weight[[k]][i-long_insert[k],1])
                           &&!(frequency_weight[[k]][i-long_insert[k]-1,1]==frequency_weight[[k]][i-long_insert[k],1]
                               &&frequency_weight[[k]][i-long_insert[k]-1,2]==frequency_weight[[k]][i-long_insert[k],2])){
                            long_insert[k]<-long_insert[k]+1
                        }
                        if(not_covered_long_insert!=0&&not_covered_long_insert!=length(samples[,1])
                           &&is.na(frequency_weight[[k]][i-long_insert[k],1])){
                            long_insert[k]<-long_insert[k]+1
                        } 
                    }
                }
                
                calling[[j]][i,1]<-frequency_weight[[j]][i-long_insert[j],1]
                calling[[j]][i,2]<-frequency_weight[[j]][i-long_insert[j],2]
                calling[[j]][i,3]<-frequency_weight[[j]][i-long_insert[j],3]
                
                if(!is.na(frequency_weight[[j]][i-long_insert[j],4])
                   ||!is.na(frequency_weight[[j]][i-long_insert[j],6])
                   ||!is.na(frequency_weight[[j]][i-long_insert[j],8])
                   ||!is.na(frequency_weight[[j]][i-long_insert[j],10]
                            ||!is.na(frequency_weight[[j]][i-long_insert[j],12]))){
                    calling[[j]][i,4]<-frequency_weight[[j]][i-long_insert[j],4]/sum(frequency_weight[[j]][i-long_insert[j],c(4,6,8,10,12)],na.rm=TRUE)
                    calling[[j]][i,5]<-frequency_weight[[j]][i-long_insert[j],6]/sum(frequency_weight[[j]][i-long_insert[j],c(4,6,8,10,12)],na.rm=TRUE)
                    calling[[j]][i,6]<-frequency_weight[[j]][i-long_insert[j],8]/sum(frequency_weight[[j]][i-long_insert[j],c(4,6,8,10,12)],na.rm=TRUE)
                    calling[[j]][i,7]<-frequency_weight[[j]][i-long_insert[j],10]/sum(frequency_weight[[j]][i-long_insert[j],c(4,6,8,10,12)],na.rm=TRUE)
                    calling[[j]][i,8]<-frequency_weight[[j]][i-long_insert[j],12]/sum(frequency_weight[[j]][i-long_insert[j],c(4,6,8,10,12)],na.rm=TRUE)
                    calling[[j]][i,9]<-frequency_weight[[j]][i-long_insert[j],13]/sum(frequency_weight[[j]][i-long_insert[j],c(4,6,8,10,12)],na.rm=TRUE) 
                }
                
                freqs<-c(calling[[j]][i,4],calling[[j]][i,5],calling[[j]][i,6],
                         calling[[j]][i,7],calling[[j]][i,8],calling[[j]][i,9])
                acgt<-c("A","C","G","T","Del","Ins")
                k<-10
                
                while(sum(freqs,na.rm=TRUE)>0){
                    if(!is.na(calling[[j]][i,3+which.max(freqs)])
                       &&freqs[which.max(freqs)]>frequency_threshold){
                        calling[[j]][i,k]<-acgt[which.max(freqs)]
                        freqs[which.max(freqs)]<-0
                        k<-k+1
                    }
                    if(!is.na(calling[[j]][i,3+which.max(freqs)])
                       &&freqs[which.max(freqs)]<=frequency_threshold){
                        freqs[which.max(freqs)]<-0
                    }
                }
            }
            write.table(calling[[j]],
                        paste(directory2,"/",samples[j,1],".calling.txt",sep=""),
                        row.names=FALSE,sep="\t",quote=FALSE)
        }
        return(calling)
    }
    
    #Variant calling:
    #define threshold for frequency of variant 
    reportVariantsVcf<-function(samples,
                                directory2,
                                frequency_weight,
                                frequency_threshold,
                                target){
        shinyjs::html("text", paste("Report Variants","<br>",sep=""), add = TRUE)  
        calling<-list()
        for(j in 1:length(samples[,1])){
            calling[[j]]<-data.frame(chr=NA,pos=NA,ref=NA,a_freq=NA,c_freq=NA,
                                     g_freq=NA,t_freq=NA,del_freq=NA,ins_freq=NA,
                                     call1=NA,call2=NA,call3=NA,call4=NA,call5=NA,
                                     call6=NA,vcf_call_1=NA,vcf_call_2=NA,Insert=NA,
                                     heterozygous=NA)
        }
        
        max.variants<-0  
        line<-rep(1,length(samples[,1]))
        start.line<-1
        for(i in 1:length(target[,1])){
            maxvar<-0
            for(j in 1:length(samples[,1])){
                maxvarsample<-0
                while(line[j]<=length(frequency_weight[[j]][,1])
                      &&paste("chr",as.character(target[start.line,1]),sep="")==frequency_weight[[j]][line[j],1]
                      &&target[start.line,2]==frequency_weight[[j]][line[j],2]){
                    maxvarsample<-maxvarsample+1
                    line[j]<-line[j]+1
                }
                if(maxvar<maxvarsample){
                    maxvar<-maxvarsample
                }
            }
            max.variants<-max.variants+maxvar
            start.line<-start.line+1
        } 
        
        for(j in 1:length(samples[,1])){
            shinyjs::html("text", paste("&nbsp&nbsp&nbspSample ",j,"<br>",sep=""), add = TRUE)  
            long_insert<-rep(0,length(samples[,1])) 
            
            for(i in 1:max.variants){
                not_covered_long_insert<-0
                if(i>1){
                    for(k in 1:length(samples[,1])){
                        if(!is.na(frequency_weight[[k]][i-long_insert[k],1])
                           &&frequency_weight[[k]][i-long_insert[k]-1,1]==frequency_weight[[k]][i-long_insert[k],1]
                           &&frequency_weight[[k]][i-long_insert[k]-1,2]==frequency_weight[[k]][i-long_insert[k],2]){
                            not_covered_long_insert<-not_covered_long_insert+1
                        }
                    } 
                    for(k in 1:length(samples[,1])){
                        if(not_covered_long_insert!=0
                           &&not_covered_long_insert!=length(samples[,1])
                           &&!is.na(frequency_weight[[k]][i-long_insert[k],1])
                           &&!(frequency_weight[[k]][i-long_insert[k]-1,1]==frequency_weight[[k]][i-long_insert[k],1]
                               &&frequency_weight[[k]][i-long_insert[k]-1,2]==frequency_weight[[k]][i-long_insert[k],2])){
                            long_insert[k]<-long_insert[k]+1
                        }
                        if(not_covered_long_insert!=0
                           &&not_covered_long_insert!=length(samples[,1])
                           &&is.na(frequency_weight[[k]][i-long_insert[k],1])){
                            long_insert[k]<-long_insert[k]+1
                        } 
                    }
                }
                
                calling[[j]][i,1]<-frequency_weight[[j]][i-long_insert[j],1]
                calling[[j]][i,2]<-frequency_weight[[j]][i-long_insert[j],2]
                calling[[j]][i,3]<-frequency_weight[[j]][i-long_insert[j],3]
                
                if(!is.na(frequency_weight[[j]][i-long_insert[j],4])
                   ||!is.na(frequency_weight[[j]][i-long_insert[j],6])
                   ||!is.na(frequency_weight[[j]][i-long_insert[j],8])
                   ||!is.na(frequency_weight[[j]][i-long_insert[j],10]
                            ||!is.na(frequency_weight[[j]][i-long_insert[j],12]))){
                    calling[[j]][i,4]<-frequency_weight[[j]][i-long_insert[j],4]/sum(frequency_weight[[j]][i-long_insert[j],c(4,6,8,10,12)],na.rm=TRUE)
                    calling[[j]][i,5]<-frequency_weight[[j]][i-long_insert[j],6]/sum(frequency_weight[[j]][i-long_insert[j],c(4,6,8,10,12)],na.rm=TRUE)
                    calling[[j]][i,6]<-frequency_weight[[j]][i-long_insert[j],8]/sum(frequency_weight[[j]][i-long_insert[j],c(4,6,8,10,12)],na.rm=TRUE)
                    calling[[j]][i,7]<-frequency_weight[[j]][i-long_insert[j],10]/sum(frequency_weight[[j]][i-long_insert[j],c(4,6,8,10,12)],na.rm=TRUE)
                    calling[[j]][i,8]<-frequency_weight[[j]][i-long_insert[j],12]/sum(frequency_weight[[j]][i-long_insert[j],c(4,6,8,10,12)],na.rm=TRUE)
                    calling[[j]][i,9]<-frequency_weight[[j]][i-long_insert[j],13]/sum(frequency_weight[[j]][i-long_insert[j],c(4,6,8,10,12)],na.rm=TRUE) 
                }
                
                calling[[j]][i,18]<-frequency_weight[[j]][i-long_insert[j],26]
                
                freqs<-c(calling[[j]][i,4],calling[[j]][i,5],calling[[j]][i,6],
                         calling[[j]][i,7],calling[[j]][i,8],calling[[j]][i,9])
                acgt<-c("A","C","G","T","Del","Ins")
                k<-10
                
                while(sum(freqs,na.rm=TRUE)>0){
                    if(!is.na(calling[[j]][i,3+which.max(freqs)])
                       &&freqs[which.max(freqs)]>frequency_threshold){
                        calling[[j]][i,k]<-acgt[which.max(freqs)]
                        freqs[which.max(freqs)]<-0
                        k<-k+1
                    }
                    if(!is.na(calling[[j]][i,3+which.max(freqs)])
                       &&freqs[which.max(freqs)]<=frequency_threshold){
                        freqs[which.max(freqs)]<-0
                    }
                }
                
                #Call 
                if(!is.na(frequency_weight[[j]][i-long_insert[j],25])){
                    if(frequency_weight[[j]][i-long_insert[j],25]=="0/1"){
                        calling[[j]][i,19]<-TRUE
                        if(!is.na(frequency_weight[[j]][i-long_insert[j],23])
                           &&nchar(frequency_weight[[j]][i-long_insert[j],23])<2){
                            calling[[j]][i,16]<-calling[[j]][i,3]
                            calling[[j]][i,17]<-frequency_weight[[j]][i-long_insert[j],23]
                        }
                        if(!is.na(frequency_weight[[j]][i-long_insert[j],23])
                           &&nchar(frequency_weight[[j]][i-long_insert[j],23])>=2){
                            calling[[j]][i,16]<-calling[[j]][i,3]
                            for(n in 1:nchar(frequency_weight[[j]][i-long_insert[j],23])){
                                if(i-nchar(frequency_weight[[j]][i-long_insert[j],23])+n>0
                                   &&calling[[j]][i-nchar(frequency_weight[[j]][i-long_insert[j],23])+n,1]==calling[[j]][i,1]
                                   &&calling[[j]][i-nchar(frequency_weight[[j]][i-long_insert[j],23])+n,2]==calling[[j]][i,2]){
                                    calling[[j]][i-nchar(frequency_weight[[j]][i-long_insert[j],23])+n,17]<-substr(frequency_weight[[j]][i-long_insert[j],23],n,n) 
                                }
                            }
                        }
                    }
                    if(frequency_weight[[j]][i-long_insert[j],25]=="1/1"){
                        if(!is.na(frequency_weight[[j]][i-long_insert[j],24])&&nchar(frequency_weight[[j]][i-long_insert[j],23])<2){
                            calling[[j]][i,16]<-frequency_weight[[j]][i-long_insert[j],23]
                            calling[[j]][i,17]<-frequency_weight[[j]][i-long_insert[j],23]
                        }
                        if(!is.na(frequency_weight[[j]][i-long_insert[j],23])
                           &&nchar(frequency_weight[[j]][i-long_insert[j],23])>=2){
                            for(n in 1:nchar(frequency_weight[[j]][i-long_insert[j],23])){
                                if(i-nchar(frequency_weight[[j]][i-long_insert[j],23])+n>0
                                   &&calling[[j]][i-nchar(frequency_weight[[j]][i-long_insert[j],23]),1]==calling[[j]][i,1]
                                   &&calling[[j]][i-nchar(frequency_weight[[j]][i-long_insert[j],23]),2]==calling[[j]][i,2]){
                                    calling[[j]][i-nchar(frequency_weight[[j]][i-long_insert[j],23])+n,16]<-substr(frequency_weight[[j]][i-long_insert[j],23],n,n)
                                    calling[[j]][i-nchar(frequency_weight[[j]][i-long_insert[j],23])+n,17]<-substr(frequency_weight[[j]][i-long_insert[j],23],n,n) 
                                }
                            }
                        }
                    }
                    if(frequency_weight[[j]][i-long_insert[j],25]=="1/2"){
                        calling[[j]][i,19]<-TRUE
                        if(!is.na(frequency_weight[[j]][i-long_insert[j],23])
                           &&!is.na(frequency_weight[[j]][i-long_insert[j],24])
                           &&nchar(frequency_weight[[j]][i-long_insert[j],23])<2
                           &&nchar(frequency_weight[[j]][i-long_insert[j],24])<2){
                            calling[[j]][i,16]<-frequency_weight[[j]][i-long_insert[j],23]
                            calling[[j]][i,17]<-frequency_weight[[j]][i-long_insert[j],24]
                        }
                        if(!is.na(frequency_weight[[j]][i-long_insert[j],23])
                           &&!is.na(frequency_weight[[j]][i-long_insert[j],24])
                           &&nchar(frequency_weight[[j]][i-long_insert[j],23])>=2
                           &&nchar(frequency_weight[[j]][i-long_insert[j],24])>=2){
                            max.insert<-max(nchar(frequency_weight[[j]][i-long_insert[j],23]),nchar(frequency_weight[[j]][i-long_insert[j],24]))
                            for(n in 1:max.insert){
                                if(i-max.insert+n>0
                                   &&calling[[j]][i-max.insert+n,1]==calling[[j]][i,1]
                                   &&calling[[j]][i-max.insert+n,2]==calling[[j]][i,2]){
                                    calling[[j]][i-nchar(frequency_weight[[j]][i-long_insert[j],23])+n,16]<-substr(frequency_weight[[j]][i-long_insert[j],23],n,n)
                                    calling[[j]][i-nchar(frequency_weight[[j]][i-long_insert[j],24])+n,17]<-substr(frequency_weight[[j]][i-long_insert[j],24],n,n) 
                                }
                            }
                        }
                    }
                    if(frequency_weight[[j]][i-long_insert[j],25]=="2/2"){
                        if(!is.na(frequency_weight[[j]][i-long_insert[j],24])
                           &&nchar(frequency_weight[[j]][i-long_insert[j],24])<2){
                            calling[[j]][i,16]<-frequency_weight[[j]][i-long_insert[j],24]
                            calling[[j]][i,17]<-frequency_weight[[j]][i-long_insert[j],24]
                        }
                        if(!is.na(frequency_weight[[j]][i-long_insert[j],24])
                           &&nchar(frequency_weight[[j]][i-long_insert[j],24])>=2){
                            for(n in 1:nchar(frequency_weight[[j]][i-long_insert[j],24])){
                                if(i-nchar(frequency_weight[[j]][i-long_insert[j],24])+n>0
                                   &&calling[[j]][i-nchar(frequency_weight[[j]][i-long_insert[j],24]),1]==calling[[j]][i,1]
                                   &&calling[[j]][i-nchar(frequency_weight[[j]][i-long_insert[j],24]),2]==calling[[j]][i,2]){
                                    calling[[j]][i-nchar(frequency_weight[[j]][i-long_insert[j],24])+n,16]<-substr(frequency_weight[[j]][i-long_insert[j],24],n,n)
                                    calling[[j]][i-nchar(frequency_weight[[j]][i-long_insert[j],24])+n,17]<-substr(frequency_weight[[j]][i-long_insert[j],24],n,n) 
                                }
                            }
                        }
                    }  
                }  
            }
            write.table(calling[[j]],
                        paste(directory2,"/",samples[j,1],".calling.txt",sep=""),
                        row.names=FALSE,sep="\t",quote=FALSE)
        }
        return(calling)
    }
    
    #Plot counts
    plotCountsPerSample_abs<-function(calling,frequency_weight,long_insert,i,j,k,
                                      plot_abs){
        if(k==1){
            plot_abs[i,k]<-frequency_weight[[j]][i-long_insert,3+k]
            if(is.na(plot_abs[i,k])){
                plot_abs[i,k]<-0
            }
        }
        if(k==3){
            plot_abs[i,2]<-frequency_weight[[j]][i-long_insert,3+k]
            if(is.na(plot_abs[i,2])){
                plot_abs[i,2]<-0
            }
        }
        if(k==5){
            plot_abs[i,3]<-frequency_weight[[j]][i-long_insert,3+k]
            if(is.na(plot_abs[i,3])){
                plot_abs[i,3]<-0
            }
        }
        if(k==7){
            plot_abs[i,4]<-frequency_weight[[j]][i-long_insert,3+k]
            if(is.na(plot_abs[i,4])){
                plot_abs[i,4]<-0
            }
        }
        if(k==9){
            plot_abs[i,5]<-frequency_weight[[j]][i-long_insert,3+k]
            if(is.na(plot_abs[i,5])){
                plot_abs[i,5]<-0
            }
        }
        if(k==10){
            plot_abs[i,6]<-frequency_weight[[j]][i-long_insert,3+k]
            if(is.na(plot_abs[i,6])){
                plot_abs[i,6]<-0
            }
        }
        return(plot_abs)
    }
    
    plotCountsPerSample_col<-function(calling,frequency_weight,long_insert,
                                      qual_lower_bound,qual_upper_bound,
                                      A_col,C_col,G_col,T_col,i,j,k,plot_col){
        if(k==1&&!is.na(frequency_weight[[j]][i-long_insert,4+k])){
            if(round(frequency_weight[[j]][i-long_insert,4+k],digits=1)<qual_lower_bound){
                plot_col[i,1]<-"chartreuse"
            }
            if(round(frequency_weight[[j]][i-long_insert,4+k],digits=1)>qual_upper_bound){
                plot_col[i,1]<-"forestgreen"
            }
            if(qual_lower_bound<=round(frequency_weight[[j]][i-long_insert,4+k],digits=1)
               &&round(frequency_weight[[j]][i-long_insert,4+k],digits=1)<=qual_upper_bound){
                plot_col[i,1]<-A_col[round(frequency_weight[[j]][i-long_insert,4+k],digits=1)-qual_lower_bound+1]
            }
        }
        if(k==3&&!is.na(frequency_weight[[j]][i-long_insert,4+k])){
            if(round(frequency_weight[[j]][i-long_insert,4+k],digits=1)<qual_lower_bound){
                plot_col[i,2]<-"cyan"
            }
            if(round(frequency_weight[[j]][i-long_insert,4+k],digits=1)>qual_upper_bound){
                plot_col[i,2]<-"blue4"
            }
            if(qual_lower_bound<=round(frequency_weight[[j]][i-long_insert,4+k],digits=1)
               &&round(frequency_weight[[j]][i-long_insert,4+k],digits=1)<=qual_upper_bound){
                plot_col[i,2]<-C_col[round(frequency_weight[[j]][i-long_insert,4+k],digits=1)-qual_lower_bound+1]
            }
        }
        if(k==5&&!is.na(frequency_weight[[j]][i-long_insert,4+k])){
            if(round(frequency_weight[[j]][i-long_insert,4+k],digits=1)<qual_lower_bound){
                plot_col[i,3]<-"yellow"
            }
            if(round(frequency_weight[[j]][i-long_insert,4+k],digits=1)>qual_upper_bound){
                plot_col[i,3]<-"darkgoldenrod3"
            }
            if(qual_lower_bound<=round(frequency_weight[[j]][i-long_insert,4+k],digits=1)
               &&round(frequency_weight[[j]][i-long_insert,4+k],digits=1)<=qual_upper_bound){
                plot_col[i,3]<-G_col[round(frequency_weight[[j]][i-long_insert,4+k],digits=1)-qual_lower_bound+1]
            }
        }
        if(k==7&&!is.na(frequency_weight[[j]][i-long_insert,4+k])){
            if(round(frequency_weight[[j]][i-long_insert,4+k],digits=1)<qual_lower_bound){
                plot_col[i,4]<-"firebrick1"
            }
            if(round(frequency_weight[[j]][i-long_insert,4+k],digits=1)>qual_upper_bound){
                plot_col[i,4]<-"darkred"
            }
            if(qual_lower_bound<=round(frequency_weight[[j]][i-long_insert,4+k],digits=1)
               &&round(frequency_weight[[j]][i-long_insert,4+k],digits=1)<=qual_upper_bound){
                plot_col[i,4]<-T_col[round(frequency_weight[[j]][i-long_insert,4+k],digits=1)-qual_lower_bound+1]
            }
        }
        if(k==9&&!is.na(frequency_weight[[j]][i-long_insert,3+k])){
            plot_col[i,5]<-"black"
        }
        if(k==10&&!is.na(frequency_weight[[j]][i-long_insert,3+k])){
            plot_col[i,6]<-"purple"
        }
        return(plot_col)
    }
    
    plotCountsPerSample_insert<-function(calling,frequency_weight,long_insert,i,j,k,
                                         plot_insert){
        if(k==1){
            plot_insert[i,k]<-frequency_weight[[j]][i-long_insert,13+k]
            if(is.na(plot_insert[i,k])){
                plot_insert[i,k]<-0
            }
        }
        if(k==3){
            plot_insert[i,2]<-frequency_weight[[j]][i-long_insert,13+k]
            if(is.na(plot_insert[i,2])){
                plot_insert[i,2]<-0
            }
        }
        if(k==5){
            plot_insert[i,3]<-frequency_weight[[j]][i-long_insert,13+k]
            if(is.na(plot_insert[i,3])){
                plot_insert[i,3]<-0
            }
        }
        if(k==7){
            plot_insert[i,4]<-frequency_weight[[j]][i-long_insert,13+k]
            if(is.na(plot_insert[i,4])){
                plot_insert[i,4]<-0
            }
        } 
        return(plot_insert)
    }
    
    plotCountsPerSample_col_insert<-function(calling,frequency_weight,long_insert,
                                             qual_lower_bound,qual_upper_bound,
                                             A_col,C_col,G_col,T_col,i,j,k,plot_col_insert){
        if(k==1&&!is.na(frequency_weight[[j]][(i-long_insert),(14+k)])){
            if(round(frequency_weight[[j]][(i-long_insert),(14+k)])<qual_lower_bound){
                plot_col_insert[i,1]<-"chartreuse"
            }
            if(round(frequency_weight[[j]][(i-long_insert),(14+k)])>qual_upper_bound){
                plot_col_insert[i,1]<-"forestgreen"
            }
            if(qual_lower_bound<=round(frequency_weight[[j]][(i-long_insert),(14+k)])
               &&round(frequency_weight[[j]][(i-long_insert),(14+k)])<=qual_upper_bound){
                plot_col_insert[i,1]<-A_col[round(frequency_weight[[j]][(i-long_insert),(14+k)])-qual_lower_bound+1]
            }
        }
        if(k==3&&!is.na(frequency_weight[[j]][(i-long_insert),(14+k)])){
            if(round(frequency_weight[[j]][(i-long_insert),(14+k)])<qual_lower_bound){
                plot_col_insert[i,1:2]<-"cyan"
            }
            if(round(frequency_weight[[j]][(i-long_insert),(14+k)])>qual_upper_bound){
                plot_col_insert[i,1:2]<-"blue4"
            }
            if(qual_lower_bound<=round(frequency_weight[[j]][(i-long_insert),(14+k)])
               &&round(frequency_weight[[j]][(i-long_insert),(14+k)])<=qual_upper_bound){
                plot_col_insert[i,1:2]<-C_col[round(frequency_weight[[j]][(i-long_insert),(14+k)])-qual_lower_bound+1]
            }
        }
        if(k==5&&!is.na(frequency_weight[[j]][(i-long_insert),(14+k)])){
            if(round(frequency_weight[[j]][(i-long_insert),(14+k)])<qual_lower_bound){
                plot_col_insert[i,1:3]<-"yellow"
            }
            if(round(frequency_weight[[j]][(i-long_insert),(14+k)])>qual_upper_bound){
                plot_col_insert[i,1:3]<-"darkgoldenrod3"
            }
            if(qual_lower_bound<=round(frequency_weight[[j]][(i-long_insert),(14+k)])
               &&round(frequency_weight[[j]][(i-long_insert),(14+k)])<=qual_upper_bound){
                plot_col_insert[i,1:3]<-G_col[round(frequency_weight[[j]][(i-long_insert),(14+k)])-qual_lower_bound+1]
            }
        }
        if(k==7&&!is.na(frequency_weight[[j]][(i-long_insert),(14+k)])){
            if(round(frequency_weight[[j]][(i-long_insert),(14+k)])<qual_lower_bound){
                plot_col_insert[i,1:4]<-"firebrick1"
            }
            if(round(frequency_weight[[j]][(i-long_insert),(14+k)])>qual_upper_bound){
                plot_col_insert[i,1:4]<-"darkred"
            }
            if(qual_lower_bound<=round(frequency_weight[[j]][(i-long_insert),(14+k)])
               &&round(frequency_weight[[j]][(i-long_insert),(14+k)])<=qual_upper_bound){
                plot_col_insert[i,1:4]<-T_col[round(frequency_weight[[j]][(i-long_insert),(14+k)])-qual_lower_bound+1]
            }
        }  
        return(plot_col_insert)
    }
    
    plotCountsPerSample_ref<-function(calling,dbsnp_file,vcf.comp,j,i,plot_ref){
        if(dbsnp_file!=""&&length(rowRanges(vcf.comp))==1
           &&length(rowRanges(vcf.comp)$REF)>0
           &&length(rowRanges(vcf.comp)$ALT)>0
           &&nchar(as.character(rowRanges(vcf.comp)$REF))==1
           &&nchar(as.character(rowRanges(vcf.comp)$ALT[[1]]))==1){
            ref<-as.character(rowRanges(vcf.comp)$REF)
            alt<-as.character(rowRanges(vcf.comp)$ALT[[1]][1])
            alt2<-""
            if(length(vcf.comp$ALT)==2){
                alt2<-(as.character(rowRanges(vcf.comp)$ALT[[2]][1]))
            }
            if(nchar(alt2)==0){
                if(ref=="A"){
                    plot_ref[i,1]<--1
                }
                if(ref=="C"){
                    plot_ref[i,2]<--1
                }
                if(ref=="G"){
                    plot_ref[i,3]<--1
                }
                if(ref=="T"){
                    plot_ref[i,4]<--1
                }
                if(alt=="A"){
                    plot_ref[i,1]<--1
                }
                if(alt=="C"){
                    plot_ref[i,2]<--1
                }
                if(alt=="G"){
                    plot_ref[i,3]<--1
                }
                if(alt=="T"){
                    plot_ref[i,4]<--1
                }
            }
            if(nchar(alt2)>0){
                if(alt2=="A"){
                    plot_ref[i,1]<--1
                }
                if(alt2=="C"){
                    plot_ref[i,2]<--1
                }
                if(alt2=="G"){
                    plot_ref[i,3]<--1
                }
                if(alt2=="T"){
                    plot_ref[i,4]<--1
                }
            }
        }
        if(dbsnp_file==""||length(rowRanges(vcf.comp))==0
           ||nchar(as.character(rowRanges(vcf.comp)$REF))!=1
           ||nchar(as.character(rowRanges(vcf.comp)$ALT[[1]]))!=1){
            if(calling[[j]][i,3]=="A"){
                plot_ref[i,1]<--1
            }
            if(calling[[j]][i,3]=="C"){
                plot_ref[i,2]<--1
            }
            if(calling[[j]][i,3]=="G"){
                plot_ref[i,3]<--1
            }
            if(calling[[j]][i,3]=="T"){
                plot_ref[i,4]<--1
            }
        } 
        return(plot_ref)
    }
    
    plotCountsPerSample_plot<-function(directory2,samples,plot_abs,plot_col,abs_max,
                                       calling,plot_insert,plot_col_insert,A_col,
                                       C_col,G_col,T_col,plot_ref,marks,
                                       qual_lower_bound,qual_upper_bound,j){
        if(directory2!=""){
            png(filename=paste(directory2,"/",samples[j,1],".png",sep=""),width=1600,height=800)
        }
        par(mar=c(16,6,2,1)) 
        barplot(t(plot_abs[[j]]),beside=TRUE,col=t(plot_col[[j]]),
                xlim=c(-1.2*length(plot_abs[[j]][,1]),7*length(plot_abs[[j]][,1])),
                ylim=c((-1*0.1*abs_max),abs_max),ylab="Number of reads",
                mgp=c(4.5,1,-2),main=samples[j,1],
                names.arg=paste(calling[[1]][,1],calling[[1]][,2],sep=";"),
                las=2,cex.axis=2,cex.names=2,cex.lab=2,cex.main=2)
        par(new=TRUE)
        barplot((t(plot_ref[[j]]))*0.1*abs_max,beside=TRUE,
                col=c("forestgreen","blue4","darkgoldenrod3","darkred","black","purple"),
                xlim=c(-1.2*length(plot_abs[[j]][,1]),7*length(plot_abs[[j]][,1])),
                ylim=c((-1*0.1*abs_max),abs_max),yaxt='n')
        for(i in 1:length(plot_abs[[j]][,1])){
            for(k in 4:1){
                rect(xleft=(7*i-0.8),xright=(7*i-0.15),ybottom=0,ytop=plot_insert[[j]][i,k],border = "white",
                     col=plot_col_insert[[j]][i,k])
            }
            if(plot_insert[[j]][i,4]!=0)
                rect(xleft=(7*i-1),xright=(7*i),ybottom=0,ytop=plot_insert[[j]][i,4],
                     col=NA,border="purple",lwd=5)
        }
        counter<-0
        for(i in 1:length(plot_abs[[j]][,1])){
            if(length(marks)>0){
                for(k in 1:length(marks)){
                    lines(x=c((1+counter*7),(7+counter*7)),
                          y=c((sum(plot_abs[[j]][i,1:5],na.rm=TRUE)*as.numeric(marks[k])),
                              (sum(plot_abs[[j]][i,1:5],na.rm=TRUE)*as.numeric(marks[k]))))
                }
                counter<-counter+1
            }
        }
        counter<-0
        for(i in 1:length(plot_abs[[j]][,1])){
            lines(x=c((1+counter*7),(7+counter*7)),
                  y=c(0,0))
            counter<-counter+1
        }
        text(x=-0.5*length(plot_abs[[j]][,1]),y=c(abs_max*0.8*(-0.1)),"Reference",cex=2)
        steps<-seq(1,length(plot_abs[[j]][,1])*0.65,
                   by=(length(plot_abs[[j]][,1])*0.65-1)/length(A_col))-1.2*length(plot_abs[[j]][,1])
        if(length(plot_abs[[j]][,1])==1){
            steps<-seq(1*0.65,1,by=(1-1*0.65)/length(A_col))-1.75*1
        }
        if(length(plot_abs[[j]][,1])==2){
            steps<-seq(0.7,2*0.65,by=(2*0.65-0.7)/length(A_col))-1.2*2
        }
        rect(xleft=steps[1:length(steps)-1],xright=steps[2:length(steps)],
             ybottom=abs_max*0.88,ytop=abs_max*0.93,lwd=0,col=A_col,border=NA)
        rect(xleft=steps[1:length(steps)-1],xright=steps[2:length(steps)],
             ybottom=abs_max*0.81,ytop=abs_max*0.86,lwd=0,col=C_col,border=NA)
        rect(xleft=steps[1:length(steps)-1],xright=steps[2:length(steps)],
             ybottom=abs_max*0.74,ytop=abs_max*0.79,lwd=0,col=G_col,border=NA)
        rect(xleft=steps[1:length(steps)-1],xright=steps[2:length(steps)],
             ybottom=abs_max*0.67,ytop=abs_max*0.72,lwd=0,col=T_col,border=NA)
        text(x=length(plot_abs[[j]][,1])*0.75-1.2*length(plot_abs[[j]][,1]),
             y=c(abs_max*0.90),"A",cex=2)
        text(x=length(plot_abs[[j]][,1])*0.75-1.2*length(plot_abs[[j]][,1]),
             y=c(abs_max*0.83),"C",cex=2)
        text(x=length(plot_abs[[j]][,1])*0.75-1.2*length(plot_abs[[j]][,1]),
             y=c(abs_max*0.76),"G",cex=2)
        text(x=length(plot_abs[[j]][,1])*0.75-1.2*length(plot_abs[[j]][,1]),
             y=c(abs_max*0.69),"T",cex=2)
        text(x=length(plot_abs[[j]][,1])*0.78-1.2*length(plot_abs[[j]][,1]),
             y=c(abs_max*0.62),"Del",cex=2)
        text(x=length(plot_abs[[j]][,1])*0.78-1.2*length(plot_abs[[j]][,1]),
             y=c(abs_max*0.55),"Ins",cex=2)
        text(x=length(plot_abs[[j]][,1])*1.15-1.2*length(plot_abs[[j]][,1]),
             y=c(abs_max*0.97),"Mean quality",cex=2)
        if(length(plot_abs[[j]][,1])>2){
            rect(xleft=1-1.2*length(plot_abs[[j]][,1]),
                 xright=length(plot_abs[[j]][,1])*0.65-1.2*length(plot_abs[[j]][,1]),
                 ybottom=abs_max*0.60,ytop=abs_max*0.65,lwd=0,col="black")
            rect(xleft=1-1.2*length(plot_abs[[j]][,1]),
                 xright=length(plot_abs[[j]][,1])*0.65-1.2*length(plot_abs[[j]][,1]),
                 ybottom=abs_max*0.53,ytop=abs_max*0.58,lwd=4,
                 border="purple",col="white")
            text(x=1-1.2*length(plot_abs[[j]][,1]),y=c(abs_max*0.97),
                 qual_lower_bound,cex=2)
            text(x=length(plot_abs[[j]][,1])*0.65-1.2*length(plot_abs[[j]][,1]),
                 y=c(abs_max*0.97),qual_upper_bound,cex=2)
        }
        if(length(plot_abs[[j]][,1])==1){
            rect(xleft=-1.1,xright=-0.75,ybottom=abs_max*0.60,
                 ytop=abs_max*0.65,lwd=0,col="black")
            rect(xleft=-1.1,xright=-0.75,ybottom=abs_max*0.53,
                 ytop=abs_max*0.58,lwd=4,border="purple",col="white")
            text(x=-1.1,y=c(abs_max*0.97),qual_lower_bound,cex=2)
            text(x=-0.75,y=c(abs_max*0.97),qual_upper_bound,cex=2)
        }
        if(length(plot_abs[[j]][,1])==2){
            rect(xleft=-1.7,xright=-1.1,ybottom=abs_max*0.60,ytop=abs_max*0.65,
                 lwd=0,col="black")
            rect(xleft=-1.7,xright=-1.1,ybottom=abs_max*0.53,ytop=abs_max*0.58,
                 lwd=4,border="purple",col="white")
            text(x=-1.7,y=c(abs_max*0.97),qual_lower_bound,cex=2)
            text(x=-1.1,y=c(abs_max*0.97),qual_upper_bound,cex=2)
        }
        if(directory2!=""){
            dev.off()
        }
        return()
    }
    
    
    plotCountsPerSampleRelative_col<-function(calling,frequency_weight,long_insert,
                                              qual_lower_bound,qual_upper_bound,
                                              A_col,C_col,G_col,T_col,i,j,k,plot_col){
        if(k==4&&!is.na(frequency_weight[[j]][i-long_insert,1+k])){
            if(round(frequency_weight[[j]][i-long_insert,1+k],digits=1)<qual_lower_bound){
                plot_col[i,1]<-"chartreuse"
            }
            if(round(frequency_weight[[j]][i-long_insert,1+k],digits=1)>qual_upper_bound){
                plot_col[i,1]<-"forestgreen"
            }
            if(qual_lower_bound<=round(frequency_weight[[j]][i-long_insert,1+k],digits=1)&&round(frequency_weight[[j]][i-long_insert,1+k],digits=1)<=qual_upper_bound){
                plot_col[i,1]<-A_col[round(frequency_weight[[j]][i-long_insert,1+k],digits=1)-qual_lower_bound+1]
            }
        }
        if(k==5&&!is.na(frequency_weight[[j]][i-long_insert,2+k])){
            if(round(frequency_weight[[j]][i-long_insert,2+k],digits=1)<qual_lower_bound){
                plot_col[i,2]<-"cyan"
            }
            if(round(frequency_weight[[j]][i-long_insert,2+k],digits=1)>qual_upper_bound){
                plot_col[i,2]<-"blue4"
            }
            if(qual_lower_bound<=round(frequency_weight[[j]][i-long_insert,2+k],digits=1)&&round(frequency_weight[[j]][i-long_insert,2+k],digits=1)<=qual_upper_bound){
                plot_col[i,2]<-C_col[round(frequency_weight[[j]][i-long_insert,2+k],digits=1)-qual_lower_bound+1]
            }
        }
        if(k==6&&!is.na(frequency_weight[[j]][i-long_insert,3+k])){
            if(round(frequency_weight[[j]][i-long_insert,3+k],digits=1)<qual_lower_bound){
                plot_col[i,3]<-"yellow"
            }
            if(round(frequency_weight[[j]][i-long_insert,3+k],digits=1)>qual_upper_bound){
                plot_col[i,3]<-"darkgoldenrod3"
            }
            if(qual_lower_bound<=round(frequency_weight[[j]][i-long_insert,3+k],digits=1)&&round(frequency_weight[[j]][i-long_insert,3+k],digits=1)<=qual_upper_bound){
                plot_col[i,3]<-G_col[round(frequency_weight[[j]][i-long_insert,3+k],digits=1)-qual_lower_bound+1]
            }
        }
        if(k==7&&!is.na(frequency_weight[[j]][i-long_insert,4+k])){
            if(round(frequency_weight[[j]][i-long_insert,4+k],digits=1)<qual_lower_bound){
                plot_col[i,4]<-"firebrick1"
            }
            if(round(frequency_weight[[j]][i-long_insert,4+k],digits=1)>qual_upper_bound){
                plot_col[i,4]<-"darkred"
            }
            if(qual_lower_bound<=round(frequency_weight[[j]][i-long_insert,4+k],digits=1)&&round(frequency_weight[[j]][i-long_insert,4+k],digits=1)<=qual_upper_bound){
                plot_col[i,4]<-T_col[round(frequency_weight[[j]][i-long_insert,4+k],digits=1)-qual_lower_bound+1]
            }
        }
        if(k==8&&!is.na(frequency_weight[[j]][i-long_insert,4+k])){
            plot_col[i,5]<-"black"
        }
        if(k==9&&!is.na(frequency_weight[[j]][i-long_insert,4+k])){
            plot_col[i,6]<-"purple"
        }
        return(plot_col)
    }
    
    plotCountsPerSampleRelative_insert<-function(calling,frequency_weight,
                                                 long_insert,i,j,k,plot_insert){
        if(k==1&&!(is.na(frequency_weight[[j]][i-long_insert,4])
                   &&is.na(frequency_weight[[j]][i-long_insert,6])
                   &&is.na(frequency_weight[[j]][i-long_insert,8])
                   &&is.na(frequency_weight[[j]][i-long_insert,10])
                   &&is.na(frequency_weight[[j]][i-long_insert,12]))){
            plot_insert[i,k]<-frequency_weight[[j]][i-long_insert,13+k]/sum(frequency_weight[[j]][i-long_insert,c(4,6,8,10,12)],na.rm=TRUE)
            if(is.na(plot_insert[i,k])){
                plot_insert[i,k]<-0
            }
        }
        if(k==3&&!(is.na(frequency_weight[[j]][i-long_insert,4])
                   &&is.na(frequency_weight[[j]][i-long_insert,6])
                   &&is.na(frequency_weight[[j]][i-long_insert,8])
                   &&is.na(frequency_weight[[j]][i-long_insert,10])
                   &&is.na(frequency_weight[[j]][i-long_insert,12]))){
            plot_insert[i,2]<-frequency_weight[[j]][i-long_insert,13+k]/sum(frequency_weight[[j]][i-long_insert,c(4,6,8,10,12)],na.rm=TRUE)
            if(is.na(plot_insert[i,2])){
                plot_insert[i,2]<-0
            }
        }
        if(k==5&&!(is.na(frequency_weight[[j]][i-long_insert,4])
                   &&is.na(frequency_weight[[j]][i-long_insert,6])
                   &&is.na(frequency_weight[[j]][i-long_insert,8])
                   &&is.na(frequency_weight[[j]][i-long_insert,10])
                   &&is.na(frequency_weight[[j]][i-long_insert,12]))){
            plot_insert[i,3]<-frequency_weight[[j]][i-long_insert,13+k]/sum(frequency_weight[[j]][i-long_insert,c(4,6,8,10,12)],na.rm=TRUE)
            if(is.na(plot_insert[i,3])){
                plot_insert[i,3]<-0
            }
        }
        if(k==7&&!(is.na(frequency_weight[[j]][i-long_insert,4])
                   &&is.na(frequency_weight[[j]][i-long_insert,6])
                   &&is.na(frequency_weight[[j]][i-long_insert,8])
                   &&is.na(frequency_weight[[j]][i-long_insert,10])
                   &&is.na(frequency_weight[[j]][i-long_insert,12]))){
            plot_insert[i,4]<-frequency_weight[[j]][i-long_insert,13+k]/sum(frequency_weight[[j]][i-long_insert,c(4,6,8,10,12)],na.rm=TRUE)
            if(is.na(plot_insert[i,4])){
                plot_insert[i,4]<-0
            }
        }
        return(plot_insert)
    }
    
    plotCountsPerSampleRelative_col_insert<-function(calling,frequency_weight,
                                                     long_insert,qual_lower_bound,
                                                     qual_upper_bound,
                                                     A_col,C_col,G_col,T_col,i,j,k,
                                                     plot_col_insert){
        if(k==1&&!is.na(frequency_weight[[j]][(i-long_insert),(14+k)])){
            if(round(frequency_weight[[j]][(i-long_insert),(14+k)])<qual_lower_bound){
                plot_col_insert[i,1]<-"chartreuse"
            }
            if(round(frequency_weight[[j]][(i-long_insert),(14+k)])>qual_upper_bound){
                plot_col_insert[i,1]<-"forestgreen"
            }
            if(qual_lower_bound<=round(frequency_weight[[j]][(i-long_insert),(14+k)])
               &&round(frequency_weight[[j]][(i-long_insert),(14+k)])<=qual_upper_bound){
                plot_col_insert[i,1]<-A_col[round(frequency_weight[[j]][(i-long_insert),(14+k)])-qual_lower_bound+1]
            }
        }
        if(k==3&&!is.na(frequency_weight[[j]][(i-long_insert),(14+k)])){
            if(round(frequency_weight[[j]][(i-long_insert),(14+k)])<qual_lower_bound){
                plot_col_insert[i,1:2]<-"cyan"
            }
            if(round(frequency_weight[[j]][(i-long_insert),(14+k)])>qual_upper_bound){
                plot_col_insert[i,1:2]<-"blue4"
            }
            if(qual_lower_bound<=round(frequency_weight[[j]][(i-long_insert),(14+k)])
               &&round(frequency_weight[[j]][(i-long_insert),(14+k)])<=qual_upper_bound){
                plot_col_insert[i,1:2]<-C_col[round(frequency_weight[[j]][(i-long_insert),(14+k)])-qual_lower_bound+1]
            }
        }
        if(k==5&&!is.na(frequency_weight[[j]][(i-long_insert),(14+k)])){
            if(round(frequency_weight[[j]][(i-long_insert),(14+k)])<qual_lower_bound){
                plot_col_insert[i,1:3]<-"yellow"
            }
            if(round(frequency_weight[[j]][(i-long_insert),(14+k)])>qual_upper_bound){
                plot_col_insert[i,1:3]<-"darkgoldenrod3"
            }
            if(qual_lower_bound<=round(frequency_weight[[j]][(i-long_insert),(14+k)])
               &&round(frequency_weight[[j]][(i-long_insert),(14+k)])<=qual_upper_bound){
                plot_col_insert[i,1:3]<-G_col[round(frequency_weight[[j]][(i-long_insert),(14+k)])-qual_lower_bound+1]
            }
        }
        if(k==7&&!is.na(frequency_weight[[j]][(i-long_insert),(14+k)])){
            if(round(frequency_weight[[j]][(i-long_insert),(14+k)])<qual_lower_bound){
                plot_col_insert[i,1:4]<-"firebrick1"
            }
            if(round(frequency_weight[[j]][(i-long_insert),(14+k)])>qual_upper_bound){
                plot_col_insert[i,1:4]<-"darkred"
            }
            if(qual_lower_bound<=round(frequency_weight[[j]][(i-long_insert),(14+k)])
               &&round(frequency_weight[[j]][(i-long_insert),(14+k)])<=qual_upper_bound){
                plot_col_insert[i,1:4]<-T_col[round(frequency_weight[[j]][(i-long_insert),(14+k)])-qual_lower_bound+1]
            }
        }
        return(plot_col_insert)
    }
    
    plotCountsPerSampleRelative_plot<-function(directory2,samples,plot_rel,plot_col,
                                               rel_max,calling,plot_insert,
                                               plot_col_insert,A_col,C_col,G_col,
                                               T_col,plot_ref,marks,
                                               qual_lower_bound,qual_upper_bound,j){
        if(directory2!=""){
            png(filename=paste(directory2,"/",samples[j,1],".png",sep=""),width=1600,height=800)
        }
        par(mar=c(16,6,2,1)) 
        barplot(t(plot_rel[[j]]),beside=TRUE,col=t(plot_col[[j]]),
                xlim=c(-1.2*length(plot_rel[[j]][,1]),7*length(plot_rel[[j]][,1])),
                ylim=c((-1*0.1*rel_max),rel_max),ylab="Relative frequency",
                mgp=c(4.5,1,-2),main=samples[j,1],names.arg=
                    paste(calling[[1]][,1],calling[[1]][,2],sep=";"),las=2,
                cex.axis=2,cex.names=2,cex.lab=2,cex.main=2)
        par(new=TRUE)
        barplot((t(plot_ref[[j]]))*0.1*rel_max,beside=TRUE,
                col=c("forestgreen","blue4","darkgoldenrod3","darkred","black","purple"),
                xlim=c(-1.2*length(plot_rel[[j]][,1]),7*length(plot_rel[[j]][,1])),
                ylim=c((-1*0.1*rel_max),rel_max),yaxt='n')
        for(i in 1:length(plot_rel[[j]][,1])){
            for(k in 4:1){
                rect(xleft=(7*i-0.8),xright=(7*i-0.15),ybottom=0,
                     ytop=plot_insert[[j]][i,k],border="white",col=plot_col_insert[[j]][i,k])
            }
            if(plot_insert[[j]][i,4]!=0)
                rect(xleft=(7*i-1),xright=(7*i),ybottom=0,
                     ytop=plot_insert[[j]][i,4],col=NA,border="purple",lwd=5)
        }
        counter<-0
        for(i in 1:length(plot_rel[[j]][,1])){
            if(length(marks)>0){
                for(k in 1:length(marks)){
                    lines(x=c((1+counter*7),(7+counter*7)),y=c(marks[k],marks[k]))
                }
                counter<-counter+1 
            }
        }
        counter<-0
        for(i in 1:length(plot_rel[[j]][,1])){
            lines(x=c((1+counter*7),(7+counter*7)),
                  y=c(0,0))
            counter<-counter+1
        }
        text(x=-0.5*length(plot_rel[[j]][,1]),y=c(rel_max*0.8*(-0.1)),"Reference",cex=2)
        steps<-seq(1,length(plot_rel[[j]][,1])*0.65,
                   by=(length(plot_rel[[j]][,1])*0.65-1)/length(A_col))-1.2*length(plot_rel[[j]][,1])
        if(length(plot_rel[[j]][,1])==1){
            steps<-seq(1*0.65,1,by=(1-1*0.65)/length(A_col))-1.75*1
        }
        if(length(plot_rel[[j]][,1])==2){
            steps<-seq(0.7,2*0.65,by=(2*0.65-0.7)/length(A_col))-1.2*2
        }
        rect(xleft=steps[1:length(steps)-1],xright=steps[2:length(steps)],
             ybottom=rel_max*0.88,ytop=rel_max*0.93,lwd=0,col=A_col,border=NA)
        rect(xleft=steps[1:length(steps)-1],xright=steps[2:length(steps)],
             ybottom=rel_max*0.81,ytop=rel_max*0.86,lwd=0,col=C_col,border=NA)
        rect(xleft=steps[1:length(steps)-1],xright=steps[2:length(steps)],
             ybottom=rel_max*0.74,ytop=rel_max*0.79,lwd=0,col=G_col,border=NA)
        rect(xleft=steps[1:length(steps)-1],xright=steps[2:length(steps)],
             ybottom=rel_max*0.67,ytop=rel_max*0.72,lwd=0,col=T_col,border=NA)
        text(x=length(plot_rel[[j]][,1])*0.75-1.2*length(plot_rel[[j]][,1]),
             y=c(rel_max*0.90),"A",cex=2)
        text(x=length(plot_rel[[j]][,1])*0.75-1.2*length(plot_rel[[j]][,1]),
             y=c(rel_max*0.83),"C",cex=2)
        text(x=length(plot_rel[[j]][,1])*0.75-1.2*length(plot_rel[[j]][,1]),
             y=c(rel_max*0.76),"G",cex=2)
        text(x=length(plot_rel[[j]][,1])*0.75-1.2*length(plot_rel[[j]][,1]),
             y=c(rel_max*0.69),"T",cex=2)
        text(x=length(plot_rel[[j]][,1])*0.78-1.2*length(plot_rel[[j]][,1]),
             y=c(rel_max*0.62),"Del",cex=2)
        text(x=length(plot_rel[[j]][,1])*0.78-1.2*length(plot_rel[[j]][,1]),
             y=c(rel_max*0.55),"Ins",cex=2)
        text(x=length(plot_rel[[j]][,1])*1.15-1.2*length(plot_rel[[j]][,1]),
             y=c(rel_max*0.97),"Mean quality",cex=2)
        if(length(plot_rel[[j]][,1])>2){
            rect(xleft=1-1.2*length(plot_rel[[j]][,1]),
                 xright=length(plot_rel[[j]][,1])*0.65-1.2*length(plot_rel[[j]][,1]),
                 ybottom=rel_max*0.60,ytop=rel_max*0.65,lwd=0,col="black")
            rect(xleft=1-1.2*length(plot_rel[[j]][,1]),
                 xright=length(plot_rel[[j]][,1])*0.65-1.2*length(plot_rel[[j]][,1]),
                 ybottom=rel_max*0.53,ytop=rel_max*0.58,lwd=4,border="purple",
                 col="white")
            text(x=1-1.2*length(plot_rel[[j]][,1]),y=c(rel_max*0.97),
                 qual_lower_bound,cex=2)
            text(x=length(plot_rel[[j]][,1])*0.65-1.2*length(plot_rel[[j]][,1]),
                 y=c(rel_max*0.97),qual_upper_bound,cex=2)
        }
        if(length(plot_rel[[j]][,1])==1){
            rect(xleft=-1.1,xright=-0.75,ybottom=rel_max*0.60,ytop=rel_max*0.65,
                 lwd=0,col="black")
            rect(xleft=-1.1,xright=-0.75,ybottom=rel_max*0.53,ytop=rel_max*0.58,
                 lwd=4,border="purple",col="white")
            text(x=-1.1,y=c(rel_max*0.97),qual_lower_bound,cex=2)
            text(x=-0.75,y=c(rel_max*0.97),qual_upper_bound,cex=2)
        }
        if(length(plot_rel[[j]][,1])==2){
            rect(xleft=-1.7,xright=-1.1,ybottom=rel_max*0.60,ytop=rel_max*0.65,
                 lwd=0,col="black")
            rect(xleft=-1.7,xright=-1.1,ybottom=rel_max*0.53,ytop=rel_max*0.58,
                 lwd=4,border="purple",col="white")
            text(x=-1.7,y=c(rel_max*0.97),qual_lower_bound,cex=2)
            text(x=-1.1,y=c(rel_max*0.97),qual_upper_bound,cex=2)
        }
        if(directory2!=""){
            dev.off()
        }
    }
    
    
    plotCountsPerTarget_abs<-function(samples,frequency_weight,long_insert,i,j,k,
                                      plot_abs){
        if(k==1){
            plot_abs[i,k]<-frequency_weight[[i]][j-long_insert[i],3+k]
            if(is.na(plot_abs[i,k])){
                plot_abs[i,k]<-0
            }
        }
        if(k==3){
            plot_abs[i,2]<-frequency_weight[[i]][j-long_insert[i],3+k]
            if(is.na(plot_abs[i,2])){
                plot_abs[i,2]<-0
            }
        }
        if(k==5){
            plot_abs[i,3]<-frequency_weight[[i]][j-long_insert[i],3+k]
            if(is.na(plot_abs[i,3])){
                plot_abs[i,3]<-0
            }
        }
        if(k==7){
            plot_abs[i,4]<-frequency_weight[[i]][j-long_insert[i],3+k]
            if(is.na(plot_abs[i,4])){
                plot_abs[i,4]<-0
            }
        }
        if(k==9){
            plot_abs[i,5]<-frequency_weight[[i]][j-long_insert[i],3+k]
            if(is.na(plot_abs[i,5])){
                plot_abs[i,5]<-0
            }
        }
        if(k==10){
            plot_abs[i,6]<-frequency_weight[[i]][j-long_insert[i],3+k]
            if(is.na(plot_abs[i,6])){
                plot_abs[i,6]<-0
            }
        }
        return(plot_abs)
    }
    
    plotCountsPerTarget_col<-function(samples,frequency_weight,long_insert,
                                      qual_lower_bound,qual_upper_bound,
                                      A_col,C_col,G_col,T_col,i,j,k,plot_col){
        if(k==1&&!is.na(frequency_weight[[i]][j-long_insert[i],4+k])){
            if(round(frequency_weight[[i]][j-long_insert[i],4+k])<qual_lower_bound){
                plot_col[i,1]<-"chartreuse"
            }
            if(round(frequency_weight[[i]][j-long_insert[i],4+k])>qual_upper_bound){
                plot_col[i,1]<-"forestgreen"
            }
            if(qual_lower_bound<=round(frequency_weight[[i]][j-long_insert[i],4+k])
               &&round(frequency_weight[[i]][j-long_insert[i],4+k])<=qual_upper_bound){
                plot_col[i,1]<-A_col[round(frequency_weight[[i]][j-long_insert[i],4+k])-qual_lower_bound+1]
            }
        }
        if(k==3&&!is.na(frequency_weight[[i]][j-long_insert[i],4+k])){
            if(round(frequency_weight[[i]][j-long_insert[i],4+k])<qual_lower_bound){
                plot_col[i,2]<-"cyan"
            }
            if(round(frequency_weight[[i]][j-long_insert[i],4+k])>qual_upper_bound){
                plot_col[i,2]<-"blue4"
            }
            if(qual_lower_bound<=round(frequency_weight[[i]][j-long_insert[i],4+k])
               &&round(frequency_weight[[i]][j-long_insert[i],4+k])<=qual_upper_bound){
                plot_col[i,2]<-C_col[round(frequency_weight[[i]][j-long_insert[i],4+k])-qual_lower_bound+1]
            }
        }
        if(k==5&&!is.na(frequency_weight[[i]][j-long_insert[i],4+k])){
            if(round(frequency_weight[[i]][j-long_insert[i],4+k])<qual_lower_bound){
                plot_col[i,3]<-"yellow"
            }
            if(round(frequency_weight[[i]][j-long_insert[i],4+k])>qual_upper_bound){
                plot_col[i,3]<-"darkgoldenrod3"
            }
            if(qual_lower_bound<=round(frequency_weight[[i]][j-long_insert[i],4+k])
               &&round(frequency_weight[[i]][j-long_insert[i],4+k])<=qual_upper_bound){
                plot_col[i,3]<-G_col[round(frequency_weight[[i]][j-long_insert[i],4+k])-qual_lower_bound+1]
            }
        }
        if(k==7&&!is.na(frequency_weight[[i]][j-long_insert[i],4+k])){
            if(round(frequency_weight[[i]][j-long_insert[i],4+k])<qual_lower_bound){
                plot_col[i,4]<-"firebrick1"
            }
            if(round(frequency_weight[[i]][j-long_insert[i],4+k])>qual_upper_bound){
                plot_col[i,4]<-"darkred"
            }
            if(qual_lower_bound<=round(frequency_weight[[i]][j-long_insert[i],4+k])
               &&round(frequency_weight[[i]][j-long_insert[i],4+k])<=qual_upper_bound){
                plot_col[i,4]<-T_col[round(frequency_weight[[i]][j-long_insert[i],4+k])-qual_lower_bound+1]
            }
        }
        if(k==9&&!is.na(frequency_weight[[i]][j-long_insert[i],3+k])){
            plot_col[i,5]<-"black"
        }
        if(k==10&&!is.na(frequency_weight[[i]][j-long_insert[i],3+k])){
            plot_col[i,6]<-"purple"
        } 
        return(plot_col)
    }
    
    plotCountsPerTarget_insert<-function(samples,frequency_weight,long_insert,i,j,k,
                                         plot_insert){
        if(k==1){
            plot_insert[i,k]<-frequency_weight[[i]][j-long_insert[i],13+k]
            if(is.na(plot_insert[i,k])){
                plot_insert[i,k]<-0
            }
        }
        if(k==3){
            plot_insert[i,2]<-frequency_weight[[i]][j-long_insert[i],13+k]
            if(is.na(plot_insert[i,2])){
                plot_insert[i,2]<-0
            }
        }
        if(k==5){
            plot_insert[i,3]<-frequency_weight[[i]][j-long_insert[i],13+k]
            if(is.na(plot_insert[i,3])){
                plot_insert[i,3]<-0
            }
        }
        if(k==7){
            plot_insert[i,4]<-frequency_weight[[i]][j-long_insert[i],13+k]
            if(is.na(plot_insert[i,4])){
                plot_insert[i,4]<-0
            }
        }
        return(plot_insert)
    }
    
    plotCountsPerTarget_col_insert<-function(samples,frequency_weight,long_insert,
                                             qual_lower_bound,qual_upper_bound,
                                             A_col,C_col,G_col,T_col,i,j,k,
                                             plot_col_insert){
        if(k==1&&!is.na(frequency_weight[[i]][j-long_insert[i],14+k])){
            if(round(frequency_weight[[i]][j-long_insert[i],14+k])<qual_lower_bound){
                plot_col_insert[i,1]<-"chartreuse"
            }
            if(round(frequency_weight[[i]][j-long_insert[i],14+k])>qual_upper_bound){
                plot_col_insert[i,1]<-"forestgreen"
            }
            if(qual_lower_bound<=round(frequency_weight[[i]][j-long_insert[i],14+k])
               &&round(frequency_weight[[i]][j-long_insert[i],14+k])<=qual_upper_bound){
                plot_col_insert[i,1]<-A_col[round(frequency_weight[[i]][j-long_insert[i],14+k])-qual_lower_bound+1]
            }
        }
        if(k==3&&!is.na(frequency_weight[[i]][j-long_insert[i],14+k])){
            if(round(frequency_weight[[i]][j-long_insert[i],14+k])<qual_lower_bound){
                plot_col_insert[i,1:2]<-"cyan"
            }
            if(round(frequency_weight[[i]][j-long_insert[i],14+k])>qual_upper_bound){
                plot_col_insert[i,1:2]<-"blue4"
            }
            if(qual_lower_bound<=round(frequency_weight[[i]][j-long_insert[i],14+k])
               &&round(frequency_weight[[i]][j-long_insert[i],14+k])<=qual_upper_bound){
                plot_col_insert[i,1:2]<-C_col[round(frequency_weight[[i]][j-long_insert[i],14+k])-qual_lower_bound+1]
            }
        }
        if(k==5&&!is.na(frequency_weight[[i]][j-long_insert[i],14+k])){
            if(round(frequency_weight[[i]][j-long_insert[i],14+k])<qual_lower_bound){
                plot_col_insert[i,1:3]<-"yellow"
            }
            if(round(frequency_weight[[i]][j-long_insert[i],14+k])>qual_upper_bound){
                plot_col_insert[i,1:3]<-"darkgoldenrod3"
            }
            if(qual_lower_bound<=round(frequency_weight[[i]][j-long_insert[i],14+k])
               &&round(frequency_weight[[i]][j-long_insert[i],14+k])<=qual_upper_bound){
                plot_col_insert[i,1:3]<-G_col[round(frequency_weight[[i]][j-long_insert[i],14+k])-qual_lower_bound+1]
            }
        }
        if(k==7&&!is.na(frequency_weight[[i]][j-long_insert[i],14+k])){
            if(round(frequency_weight[[i]][j-long_insert[i],14+k])<qual_lower_bound){
                plot_col_insert[i,1:4]<-"firebrick1"
            }
            if(round(frequency_weight[[i]][j-long_insert[i],14+k])>qual_upper_bound){
                plot_col_insert[i,1:4]<-"darkred"
            }
            if(qual_lower_bound<=round(frequency_weight[[i]][j-long_insert[i],14+k])
               &&round(frequency_weight[[i]][j-long_insert[i],14+k])<=qual_upper_bound){
                plot_col_insert[i,1:4]<-T_col[round(frequency_weight[[i]][j-long_insert[i],14+k])-qual_lower_bound+1]
            }
        } 
        return(plot_col_insert)
    }
    
    plotCountsPerTarget_ref<-function(samples,calling,dbsnp_file,vcf.comp,j,i,
                                      plot_ref){
        if(dbsnp_file!=""&&length(rowRanges(vcf.comp))==1
           &&length(rowRanges(vcf.comp)$REF)>0
           &&length(rowRanges(vcf.comp)$ALT)>0
           &&nchar(as.character(rowRanges(vcf.comp)$REF))==1
           &&nchar(as.character(rowRanges(vcf.comp)$ALT[[1]]))==1){
            ref<-as.character(rowRanges(vcf.comp)$REF)
            alt<-as.character(rowRanges(vcf.comp)$ALT[[1]][1])
            alt2<-""
            if(length(vcf.comp$ALT)==2){
                alt2<-(as.character(rowRanges(vcf.comp)$ALT[[2]][1]))
            }
            if(nchar(alt2)==0){
                if(ref=="A"){
                    plot_ref[i,1]<--1
                }
                if(ref=="C"){
                    plot_ref[i,2]<--1
                }
                if(ref=="G"){
                    plot_ref[i,3]<--1
                }
                if(ref=="T"){
                    plot_ref[i,4]<--1
                }
                if(alt=="A"){
                    plot_ref[i,1]<--1
                }
                if(alt=="C"){
                    plot_ref[i,2]<--1
                }
                if(alt=="G"){
                    plot_ref[i,3]<--1
                }
                if(alt=="T"){
                    plot_ref[i,4]<--1
                }
            }
            if(nchar(alt2)>0){
                if(alt2=="A"){
                    plot_ref[i,1]<--1
                }
                if(alt2=="C"){
                    plot_ref[i,2]<--1
                }
                if(alt2=="G"){
                    plot_ref[i,3]<--1
                }
                if(alt2=="T"){
                    plot_ref[i,4]<--1
                }
            }
        }
        if(dbsnp_file==""||length(rowRanges(vcf.comp))==0
           ||nchar(as.character(rowRanges(vcf.comp)$REF))!=1
           ||nchar(as.character(rowRanges(vcf.comp)$ALT[[1]]))!=1){
            if(calling[[i]][j,3]=="A"){
                plot_ref[i,1]<--1
            }
            if(calling[[i]][j,3]=="C"){
                plot_ref[i,2]<--1
            }
            if(calling[[i]][j,3]=="G"){
                plot_ref[i,3]<--1
            }
            if(calling[[i]][j,3]=="T"){
                plot_ref[i,4]<--1
            }
        }  
        return(plot_ref)
    }
    
    plotCountsPerTarget_plot<-function(directory2,samples,plot_abs,plot_col,abs_max,
                                       calling,plot_insert,plot_col_insert,A_col,
                                       C_col,G_col,T_col,plot_ref,marks,
                                       qual_lower_bound,qual_upper_bound,counter,
                                       j){
        if(directory2!=""){
            if(j>1&&calling[[1]][j-1,1]==calling[[1]][j,1]
               &&calling[[1]][j-1,2]==calling[[1]][j,2]){
                counter<-counter+1
                png(filename=paste(directory2,"/",calling[[1]][j,1],";",
                                   calling[[1]][j,2],"_",counter,".png",sep=""),
                    width=1600,height=600)
            }
            if(j==1||calling[[1]][j-1,1]!=calling[[1]][j,1]
               ||calling[[1]][j-1,2]!=calling[[1]][j,2]){
                png(filename=paste(directory2,"/",calling[[1]][j,1],";",
                                   calling[[1]][j,2],".png",sep=""),
                    width=1600,height=600)
                counter<-0
            }
        }
        par(mar=c(16,6,2,1)) 
        if(counter==0){
            barplot(t(plot_abs[[j]]),beside=TRUE,col=t(plot_col[[j]]),
                    xlim=c(-1.2*length(plot_abs[[j]][,1]),7*length(plot_abs[[j]][,1])),
                    ylim=c((-1*0.1*abs_max),abs_max),ylab="Number of reads",
                    main=paste(calling[[1]][j,1],";",calling[[1]][j,2],sep=""),
                    names.arg=t(samples[1]),las=2,cex.axis=2,cex.names=2,cex.lab=2,
                    cex.main=2)
        }
        if(counter>0){
            barplot(t(plot_abs[[j]]),beside=TRUE,col=t(plot_col[[j]]),
                    xlim=c(-1.2*length(plot_abs[[j]][,1]),7*length(plot_abs[[j]][,1])),
                    ylim=c((-1*0.1*abs_max),abs_max),ylab="Number of reads",
                    main=paste(calling[[1]][j,1],";",calling[[1]][j,2],"_",counter,sep=""),
                    names.arg=t(samples[1]),las=2,cex.axis=2,cex.names=2,
                    cex.lab=2,cex.main=2)
        }
        par(new=TRUE)
        barplot((t(plot_ref[[j]]))*0.1*abs_max,beside=TRUE,
                col=c("forestgreen","blue4","darkgoldenrod3","darkred","black","purple"),
                xlim=c(-1.2*length(plot_abs[[j]][,1]),7*length(plot_abs[[j]][,1])),
                ylim=c((-1*0.1*abs_max),abs_max),yaxt='n')
        for(i in 1:length(plot_abs[[j]][,1])){
            for(k in 4:1){
                rect(xleft=(7*i-0.8),xright=(7*i-0.15),ybottom=0,
                     ytop=plot_insert[[j]][i,k],border="white",col=plot_col_insert[[j]][i,k])
            }
            if(plot_insert[[j]][i,4]!=0)
                rect(xleft=(7*i-1),xright=(7*i),ybottom=0,
                     ytop=plot_insert[[j]][i,4],col=NA,border="purple",lwd=5)
        }
        counter2<-0
        for(i in 1:length(plot_abs[[j]][,1])){
            if(length(marks)>0){
                for(k in 1:length(marks)){
                    lines(x=c((1+counter2*7),(7+counter2*7)),
                          y=c((sum(plot_abs[[j]][i,1:5],na.rm=TRUE)*as.numeric(marks[k])),
                              (sum(plot_abs[[j]][i,1:5],na.rm=TRUE)*as.numeric(marks[k]))))
                }
                counter2<-counter2+1
            }
        }
        counter2<-0
        for(i in 1:length(plot_abs[[j]][,1])){
            lines(x=c((1+counter2*7),(7+counter2*7)),
                  y=c(0,0))
            counter2<-counter2+1
        }
        text(x=-0.5*length(plot_abs[[j]][,1]),y=c(abs_max*0.8*(-0.1)),"Reference",cex=2)
        steps<-seq(1,length(plot_abs[[j]][,1])*0.65,
                   by=(length(plot_abs[[j]][,1])*0.65-1)/length(A_col))-1.2*length(plot_abs[[j]][,1])
        if(length(plot_abs[[j]][,1])==1){
            steps<-seq(1*0.65,1,by=(1-1*0.65)/length(A_col))-1.75*1
        }
        if(length(plot_abs[[j]][,1])==2){
            steps<-seq(0.7,2*0.65,by=(2*0.65-0.7)/length(A_col))-1.2*2
        }   
        rect(xleft=steps[1:length(steps)-1],xright=steps[2:length(steps)],
             ybottom=abs_max*0.88,ytop=abs_max*0.93,lwd=0,col=A_col,border=NA)
        rect(xleft=steps[1:length(steps)-1],xright=steps[2:length(steps)],
             ybottom=abs_max*0.81,ytop=abs_max*0.86,lwd=0,col=C_col,border=NA)
        rect(xleft=steps[1:length(steps)-1],xright=steps[2:length(steps)],
             ybottom=abs_max*0.74,ytop=abs_max*0.79,lwd=0,col=G_col,border=NA)
        rect(xleft=steps[1:length(steps)-1],xright=steps[2:length(steps)],
             ybottom=abs_max*0.67,ytop=abs_max*0.72,lwd=0,col=T_col,border=NA)
        text(x=length(plot_abs[[j]][,1])*0.75-1.2*length(plot_abs[[j]][,1]),
             y=c(abs_max*0.90),"A",cex=2)
        text(x=length(plot_abs[[j]][,1])*0.75-1.2*length(plot_abs[[j]][,1]),
             y=c(abs_max*0.83),"C",cex=2)
        text(x=length(plot_abs[[j]][,1])*0.75-1.2*length(plot_abs[[j]][,1]),
             y=c(abs_max*0.76),"G",cex=2)
        text(x=length(plot_abs[[j]][,1])*0.75-1.2*length(plot_abs[[j]][,1]),
             y=c(abs_max*0.69),"T",cex=2)
        text(x=length(plot_abs[[j]][,1])*0.78-1.2*length(plot_abs[[j]][,1]),
             y=c(abs_max*0.62),"Del",cex=2)
        text(x=length(plot_abs[[j]][,1])*0.78-1.2*length(plot_abs[[j]][,1]),
             y=c(abs_max*0.55),"Ins",cex=2)
        text(x=length(plot_abs[[j]][,1])*1.15-1.2*length(plot_abs[[j]][,1]),
             y=c(abs_max*0.97),"Mean quality",cex=2)
        if(length(plot_abs[[j]][,1])>2){
            rect(xleft=1-1.2*length(plot_abs[[j]][,1]),
                 xright=length(plot_abs[[j]][,1])*0.65-1.2*length(plot_abs[[j]][,1]),
                 ybottom=abs_max*0.60,ytop=abs_max*0.65,lwd=0,col="black")
            rect(xleft=1-1.2*length(plot_abs[[j]][,1]),
                 xright=length(plot_abs[[j]][,1])*0.65-1.2*length(plot_abs[[j]][,1]),
                 ybottom=abs_max*0.53,ytop=abs_max*0.58,lwd=4,
                 border="purple",col="white")
            text(x=1-1.2*length(plot_abs[[j]][,1]),y=c(abs_max*0.97),
                 qual_lower_bound,cex=2)
            text(x=length(plot_abs[[j]][,1])*0.65-1.2*length(plot_abs[[j]][,1]),
                 y=c(abs_max*0.97),qual_upper_bound,cex=2)
        }
        if(length(plot_abs[[j]][,1])==1){
            rect(xleft=-1.1,xright=-0.75,ybottom=abs_max*0.60,ytop=abs_max*0.65,
                 lwd=0,col="black")
            rect(xleft=-1.1,xright=-0.75,ybottom=abs_max*0.53,ytop=abs_max*0.58,
                 lwd=4,border="purple",col="white")
            text(x=-1.1,y=c(abs_max*0.97),qual_lower_bound,cex=2)
            text(x=-0.75,y=c(abs_max*0.97),qual_upper_bound,cex=2)
        }
        if(length(plot_abs[[j]][,1])==2){
            rect(xleft=-1.7,xright=-1.1,ybottom=abs_max*0.60,ytop=abs_max*0.65,
                 lwd=0,col="black")
            rect(xleft=-1.7,xright=-1.1,ybottom=abs_max*0.53,ytop=abs_max*0.58,
                 lwd=4,border="purple",col="white")
            text(x=-1.7,y=c(abs_max*0.97),qual_lower_bound,cex=2)
            text(x=-1.1,y=c(abs_max*0.97),qual_upper_bound,cex=2)
        }
        if(directory2!=""){
            dev.off()
        }
        return(counter)
    }
    
    
    plotCountsPerTargetRelative_col<-function(samples,frequency_weight,long_insert,
                                              qual_lower_bound,qual_upper_bound,
                                              A_col,C_col,G_col,T_col,i,j,k,plot_col){
        if(k==4&&!is.na(frequency_weight[[i]][j-long_insert[i],1+k])){
            if(round(frequency_weight[[i]][j-long_insert[i],1+k])<qual_lower_bound){
                plot_col[i,1]<-"chartreuse"
            }
            if(round(frequency_weight[[i]][j-long_insert[i],1+k])>qual_upper_bound){
                plot_col[i,1]<-"forestgreen"
            }
            if(qual_lower_bound<=round(frequency_weight[[i]][j-long_insert[i],1+k])&&round(frequency_weight[[i]][j-long_insert[i],1+k])<=qual_upper_bound){
                plot_col[i,1]<-A_col[round(frequency_weight[[i]][j-long_insert[i],1+k])-qual_lower_bound+1]
            }
        }
        if(k==5&&!is.na(frequency_weight[[i]][j-long_insert[i],2+k])){
            if(round(frequency_weight[[i]][j-long_insert[i],2+k])<qual_lower_bound){
                plot_col[i,2]<-"cyan"
            }
            if(round(frequency_weight[[i]][j-long_insert[i],2+k])>qual_upper_bound){
                plot_col[i,2]<-"blue4"
            }
            if(qual_lower_bound<=round(frequency_weight[[i]][j-long_insert[i],2+k])&&round(frequency_weight[[i]][j-long_insert[i],2+k])<=qual_upper_bound){
                plot_col[i,2]<-C_col[round(frequency_weight[[i]][j-long_insert[i],2+k])-qual_lower_bound+1]
            }
        }
        if(k==6&&!is.na(frequency_weight[[i]][j-long_insert[i],3+k])){
            if(round(frequency_weight[[i]][j-long_insert[i],3+k])<qual_lower_bound){
                plot_col[i,3]<-"yellow"
            }
            if(round(frequency_weight[[i]][j-long_insert[i],3+k])>qual_upper_bound){
                plot_col[i,3]<-"darkgoldenrod3"
            }
            if(qual_lower_bound<=round(frequency_weight[[i]][j-long_insert[i],3+k])&&round(frequency_weight[[i]][j-long_insert[i],3+k])<=qual_upper_bound){
                plot_col[i,3]<-G_col[round(frequency_weight[[i]][j-long_insert[i],3+k])-qual_lower_bound+1]
            }
        }
        if(k==7&&!is.na(frequency_weight[[i]][j-long_insert[i],4+k])){
            if(round(frequency_weight[[i]][j-long_insert[i],4+k])<qual_lower_bound){
                plot_col[i,4]<-"firebrick1"
            }
            if(round(frequency_weight[[i]][j-long_insert[i],4+k])>qual_upper_bound){
                plot_col[i,4]<-"darkred"
            }
            if(qual_lower_bound<=round(frequency_weight[[i]][j-long_insert[i],4+k])&&round(frequency_weight[[i]][j-long_insert[i],4+k])<=qual_upper_bound){
                plot_col[i,4]<-T_col[round(frequency_weight[[i]][j-long_insert[i],4+k])-qual_lower_bound+1]
            }
        }
        if(k==8&&!is.na(frequency_weight[[i]][j-long_insert[i],4+k])){
            plot_col[i,5]<-"black"
        }
        if(k==9&&!is.na(frequency_weight[[i]][j-long_insert[i],4+k])){
            plot_col[i,6]<-"purple"
        }
        return(plot_col)
    }
    
    plotCountsPerTargetRelative_insert<-function(samples,frequency_weight,
                                                 long_insert,i,j,k,plot_insert){
        if(k==1&&!(is.na(frequency_weight[[i]][j-long_insert[i],4])
                   &&is.na(frequency_weight[[i]][j-long_insert[i],6])
                   &&is.na(frequency_weight[[i]][j-long_insert[i],8])
                   &&is.na(frequency_weight[[i]][j-long_insert[i],10])
                   &&is.na(frequency_weight[[i]][j-long_insert[i],12]))){
            plot_insert[i,k]<-frequency_weight[[i]][j-long_insert[i],13+k]/sum(frequency_weight[[i]][j-long_insert[i],c(4,6,8,10,12)],na.rm=TRUE)
            if(is.na(plot_insert[i,k])){
                plot_insert[i,k]<-0
            }
        }
        if(k==3&&!(is.na(frequency_weight[[i]][j-long_insert[i],4])
                   &&is.na(frequency_weight[[i]][j-long_insert[i],6])
                   &&is.na(frequency_weight[[i]][j-long_insert[i],8])
                   &&is.na(frequency_weight[[i]][j-long_insert[i],10])
                   &&is.na(frequency_weight[[i]][j-long_insert[i],12]))){
            plot_insert[i,2]<-frequency_weight[[i]][j-long_insert[i],13+k]/sum(frequency_weight[[i]][j-long_insert[i],c(4,6,8,10,12)],na.rm=TRUE)
            if(is.na(plot_insert[i,2])){
                plot_insert[i,2]<-0
            }
        }
        if(k==5&&!(is.na(frequency_weight[[i]][j-long_insert[i],4])
                   &&is.na(frequency_weight[[i]][j-long_insert[i],6])
                   &&is.na(frequency_weight[[i]][j-long_insert[i],8])
                   &&is.na(frequency_weight[[i]][j-long_insert[i],10])
                   &&is.na(frequency_weight[[i]][j-long_insert[i],12]))){
            plot_insert[i,3]<-frequency_weight[[i]][j-long_insert[i],13+k]/sum(frequency_weight[[i]][j-long_insert[i],c(4,6,8,10,12)],na.rm=TRUE)
            if(is.na(plot_insert[i,3])){
                plot_insert[i,3]<-0
            }
        }
        if(k==7&&!(is.na(frequency_weight[[i]][j-long_insert[i],4])
                   &&is.na(frequency_weight[[i]][j-long_insert[i],6])
                   &&is.na(frequency_weight[[i]][j-long_insert[i],8])
                   &&is.na(frequency_weight[[i]][j-long_insert[i],10])
                   &&is.na(frequency_weight[[i]][j-long_insert[i],12]))){
            plot_insert[i,4]<-frequency_weight[[i]][j-long_insert[i],13+k]/sum(frequency_weight[[i]][j-long_insert[i],c(4,6,8,10,12)],na.rm=TRUE)
            if(is.na(plot_insert[i,4])){
                plot_insert[i,4]<-0
            }
        }
        return(plot_insert)
    }
    
    plotCountsPerTargetRelative_col_insert<-function(samples,frequency_weight,
                                                     long_insert,qual_lower_bound,
                                                     qual_upper_bound,
                                                     A_col,C_col,G_col,T_col,i,j,k,
                                                     plot_col_insert){
        if(k==1&&!is.na(frequency_weight[[i]][j-long_insert[i],14+k])){
            if(round(frequency_weight[[i]][j-long_insert[i],14+k])<qual_lower_bound){
                plot_col_insert[i,1]<-"chartreuse"
            }
            if(round(frequency_weight[[i]][j-long_insert[i],14+k])>qual_upper_bound){
                plot_col_insert[i,1]<-"forestgreen"
            }
            if(qual_lower_bound<=round(frequency_weight[[i]][j-long_insert[i],14+k])
               &&round(frequency_weight[[i]][j-long_insert[i],14+k])<=qual_upper_bound){
                plot_col_insert[i,1]<-A_col[round(frequency_weight[[i]][j-long_insert[i],14+k])-qual_lower_bound+1]
            }
        }
        if(k==3&&!is.na(frequency_weight[[i]][j-long_insert[i],14+k])){
            if(round(frequency_weight[[i]][j-long_insert[i],14+k])<qual_lower_bound){
                plot_col_insert[i,1:2]<-"cyan"
            }
            if(round(frequency_weight[[i]][j-long_insert[i],14+k])>qual_upper_bound){
                plot_col_insert[i,1:2]<-"blue4"
            }
            if(qual_lower_bound<=round(frequency_weight[[i]][j-long_insert[i],14+k])
               &&round(frequency_weight[[i]][j-long_insert[i],14+k])<=qual_upper_bound){
                plot_col_insert[i,1:2]<-C_col[round(frequency_weight[[i]][j-long_insert[i],14+k])-qual_lower_bound+1]
            }
        }
        if(k==5&&!is.na(frequency_weight[[i]][j-long_insert[i],14+k])){
            if(round(frequency_weight[[i]][j-long_insert[i],14+k])<qual_lower_bound){
                plot_col_insert[i,1:3]<-"yellow"
            }
            if(round(frequency_weight[[i]][j-long_insert[i],14+k])>qual_upper_bound){
                plot_col_insert[i,1:3]<-"darkgoldenrod3"
            }
            if(qual_lower_bound<=round(frequency_weight[[i]][j-long_insert[i],14+k])
               &&round(frequency_weight[[i]][j-long_insert[i],14+k])<=qual_upper_bound){
                plot_col_insert[i,1:3]<-G_col[round(frequency_weight[[i]][j-long_insert[i],14+k])-qual_lower_bound+1]
            }
        }
        if(k==7&&!is.na(frequency_weight[[i]][j-long_insert[i],14+k])){
            if(round(frequency_weight[[i]][j-long_insert[i],14+k])<qual_lower_bound){
                plot_col_insert[i,1:4]<-"firebrick1"
            }
            if(round(frequency_weight[[i]][j-long_insert[i],14+k])>qual_upper_bound){
                plot_col_insert[i,1:4]<-"darkred"
            }
            if(qual_lower_bound<=round(frequency_weight[[i]][j-long_insert[i],14+k])
               &&round(frequency_weight[[i]][j-long_insert[i],14+k])<=qual_upper_bound){
                plot_col_insert[i,1:4]<-T_col[round(frequency_weight[[i]][j-long_insert[i],14+k])-qual_lower_bound+1]
            }
        }  
        return(plot_col_insert)
    }
    
    plotCountsPerTargetRelative_plot<-function(directory2,samples,plot_rel,plot_col,
                                               rel_max,calling,plot_insert,
                                               plot_col_insert,A_col,C_col,G_col,
                                               T_col,plot_ref,marks,
                                               qual_lower_bound,qual_upper_bound,
                                               counter,j){
        if(directory2!=""){
            if(j>1&&calling[[1]][j-1,1]==calling[[1]][j,1]
               &&calling[[1]][j-1,2]==calling[[1]][j,2]){
                counter<-counter+1
                png(filename=paste(directory2,"/",calling[[1]][j,1],";",calling[[1]][j,2],"_",counter,".png",sep=""),
                    width=1600,height=600)
            }
            if(j==1||calling[[1]][j-1,1]!=calling[[1]][j,1]||calling[[1]][j-1,2]!=calling[[1]][j,2]){
                png(filename=paste(directory2,"/",calling[[1]][j,1],";",calling[[1]][j,2],".png",sep=""),
                    width=1600,height=600)
                counter<-0
            }
        }
        par(mar=c(16,6,2,1)) 
        if(counter==0){
            barplot(t(plot_rel[[j]]),beside=TRUE,col=t(plot_col[[j]]),
                    xlim=c(-1.2*length(plot_rel[[j]][,1]),7*length(plot_rel[[j]][,1])),
                    ylim=c((-1*0.1*rel_max),rel_max),ylab="Relative frequency",
                    main=paste(calling[[1]][j,1],";",calling[[1]][j,2],sep=""),
                    names.arg=t(samples[1]),las=2,cex.axis=2,cex.names=2,cex.lab=2,
                    cex.main=2)
        }
        if(counter>0){
            barplot(t(plot_rel[[j]]),beside=TRUE,col=t(plot_col[[j]]),
                    xlim=c(-1.2*length(plot_rel[[j]][,1]),7*length(plot_rel[[j]][,1])),
                    ylim=c((-1*0.1*rel_max),rel_max),ylab="Relative frequency",
                    main=paste(calling[[1]][j,1],";",calling[[1]][j,2],"_",counter,sep=""),
                    names.arg=t(samples[1]),las=2,cex.axis=2,cex.names=2,cex.lab=2,
                    cex.main=2)
        }
        par(new=TRUE)
        barplot((t(plot_ref[[j]]))*0.1*rel_max,beside=TRUE,
                col=c("forestgreen","blue4","darkgoldenrod3","darkred","black","purple"),
                xlim=c(-1.2*length(plot_rel[[j]][,1]),7*length(plot_rel[[j]][,1])),
                ylim=c((-1*0.1*rel_max),rel_max),yaxt='n')
        for(i in 1:length(plot_rel[[j]][,1])){
            for(k in 4:1){
                rect(xleft=(7*i-0.8),xright=(7*i-0.15),ybottom=0,
                     ytop=plot_insert[[j]][i,k],border="white",col=plot_col_insert[[j]][i,k])
            }
            if(plot_insert[[j]][i,4]!=0)
                rect(xleft=(7*i-1),xright=(7*i),ybottom=0,
                     ytop=plot_insert[[j]][i,4],col=NA,border="purple",lwd=5)
        }
        counter2<-0
        for(i in 1:length(plot_rel[[j]][,1])){
            if(length(marks)>0){
                for(k in 1:length(marks)){
                    lines(x=c((1+counter2*7),(7+counter2*7)),y=c(marks[k],marks[k]))
                }
                counter2<-counter2+1
            }
        }
        counter2<-0
        for(i in 1:length(plot_rel[[j]][,1])){
            lines(x=c((1+counter2*7),(7+counter2*7)),
                  y=c(0,0))
            counter2<-counter2+1
        }
        text(x=-0.5*length(plot_rel[[j]][,1]),y=c(rel_max*0.8*(-0.1)),"Reference",cex=2)
        steps<-seq(1,length(plot_rel[[j]][,1])*0.65,
                   by=(length(plot_rel[[j]][,1])*0.65-1)/length(A_col))-1.2*length(plot_rel[[j]][,1])
        if(length(plot_rel[[j]][,1])==1){
            steps<-seq(1*0.65,1,by=(1-1*0.65)/length(A_col))-1.75*1
        }
        if(length(plot_rel[[j]][,1])==2){
            steps<-seq(0.7,2*0.65,by=(2*0.65-0.7)/length(A_col))-1.2*2
        }    
        rect(xleft=steps[1:length(steps)-1],xright=steps[2:length(steps)],
             ybottom=rel_max*0.88,ytop=rel_max*0.93,lwd=0,col=A_col,border=NA)
        rect(xleft=steps[1:length(steps)-1],xright=steps[2:length(steps)],
             ybottom=rel_max*0.81,ytop=rel_max*0.86,lwd=0,col=C_col,border=NA)
        rect(xleft=steps[1:length(steps)-1],xright=steps[2:length(steps)],
             ybottom=rel_max*0.74,ytop=rel_max*0.79,lwd=0,col=G_col,border=NA)
        rect(xleft=steps[1:length(steps)-1],xright=steps[2:length(steps)],
             ybottom=rel_max*0.67,ytop=rel_max*0.72,lwd=0,col=T_col,border=NA)
        text(x=length(plot_rel[[j]][,1])*0.75-1.2*length(plot_rel[[j]][,1]),
             y=c(rel_max*0.90),"A",cex=2)
        text(x=length(plot_rel[[j]][,1])*0.75-1.2*length(plot_rel[[j]][,1]),
             y=c(rel_max*0.83),"C",cex=2)
        text(x=length(plot_rel[[j]][,1])*0.75-1.2*length(plot_rel[[j]][,1]),
             y=c(rel_max*0.76),"G",cex=2)
        text(x=length(plot_rel[[j]][,1])*0.75-1.2*length(plot_rel[[j]][,1]),
             y=c(rel_max*0.69),"T",cex=2)
        text(x=length(plot_rel[[j]][,1])*0.78-1.2*length(plot_rel[[j]][,1]),
             y=c(rel_max*0.62),"Del",cex=2)
        text(x=length(plot_rel[[j]][,1])*0.78-1.2*length(plot_rel[[j]][,1]),
             y=c(rel_max*0.55),"Ins",cex=2)
        text(x=length(plot_rel[[j]][,1])*1.15-1.2*length(plot_rel[[j]][,1]),
             y=c(rel_max*0.97),"Mean quality",cex=2)
        if(length(plot_rel[[j]][,1])>2){
            rect(xleft=1-1.2*length(plot_rel[[j]][,1]),
                 xright=length(plot_rel[[j]][,1])*0.65-1.2*length(plot_rel[[j]][,1]),
                 ybottom=rel_max*0.60,ytop=rel_max*0.65,lwd=0,col="black")
            rect(xleft=1-1.2*length(plot_rel[[j]][,1]),
                 xright=length(plot_rel[[j]][,1])*0.65-1.2*length(plot_rel[[j]][,1]),
                 ybottom=rel_max*0.53,ytop=rel_max*0.58,lwd=4,border="purple",
                 col="white")
            text(x=1-1.2*length(plot_rel[[j]][,1]),y=c(rel_max*0.97),
                 qual_lower_bound,cex=2)
            text(x=length(plot_rel[[j]][,1])*0.65-1.2*length(plot_rel[[j]][,1]),
                 y=c(rel_max*0.97),qual_upper_bound,cex=2)
        }
        if(length(plot_rel[[j]][,1])==1){
            rect(xleft=-1.1,xright=-0.75,ybottom=rel_max*0.60,ytop=rel_max*0.65,
                 lwd=0,col="black")
            rect(xleft=-1.1,xright=-0.75,ybottom=rel_max*0.53,ytop=rel_max*0.58,
                 lwd=4,border="purple",col="white")
            text(x=-1.1,y=c(rel_max*0.97),qual_lower_bound,cex=2)
            text(x=-0.75,y=c(rel_max*0.97),qual_upper_bound,cex=2)
        }
        if(length(plot_rel[[j]][,1])==2){
            rect(xleft=-1.7,xright=-1.1,ybottom=rel_max*0.60,ytop=rel_max*0.65,
                 lwd=0,col="black")
            rect(xleft=-1.7,xright=-1.1,ybottom=rel_max*0.53,ytop=rel_max*0.58,
                 lwd=4,border="purple",col="white")
            text(x=-1.7,y=c(rel_max*0.97),qual_lower_bound,cex=2)
            text(x=-1.1,y=c(rel_max*0.97),qual_upper_bound,cex=2)
        }
        if(directory2!=""){
            dev.off()
        }
        return(counter)
    }
    
    
    plotCountsPerSampleVcf_plot<-function(directory2,samples,plot_abs,plot_col,
                                          abs_max,calling,plot_insert,
                                          plot_col_insert,A_col,C_col,G_col,T_col,
                                          plot_ref,marks,qual_lower_bound,
                                          qual_upper_bound,j){
        if(directory2!=""){
            png(filename=paste(directory2,"/",samples[j,1],".png",sep=""),width=1600,
                height=800)   
        }
        par(mar=c(16,6,2,1)) 
        barplot(t(plot_abs[[j]]),beside=TRUE,col=t(plot_col[[j]]),
                xlim=c(-1.2*length(plot_abs[[j]][,1]),7*length(plot_abs[[j]][,1])),
                ylim=c((-1*0.1*abs_max),abs_max),ylab="Number of reads",
                mgp=c(4.5,1,-2),main=samples[j,1],
                names.arg=paste(calling[[1]][,1],calling[[1]][,2],sep=";"),
                las=2,cex.axis=2,cex.names=2,cex.lab=2,cex.main=2)
        par(new=TRUE)
        barplot((t(plot_ref[[j]]))*0.1*abs_max,beside=TRUE,
                col=c("forestgreen","blue4","darkgoldenrod3","darkred","black","purple"),
                xlim=c(-1.2*length(plot_abs[[j]][,1]),7*length(plot_abs[[j]][,1])),
                ylim=c((-1*0.1*abs_max),abs_max),yaxt='n')
        for(i in 1:length(plot_abs[[j]][,1])){
            for(k in 4:1){
                rect(xleft=(7*i-0.8),xright=(7*i-0.15),ybottom=0,
                     ytop=plot_insert[[j]][i,k],border="white",col=plot_col_insert[[j]][i,k])
            }
            if(plot_insert[[j]][i,4]!=0)
                rect(xleft=(7*i-1),xright=(7*i),ybottom=0,
                     ytop=plot_insert[[j]][i,4],col=NA,border="purple",lwd=5)
        }
        counter<-0
        for(i in 1:length(plot_abs[[j]][,1])){
            if(length(marks)){
                for(k in 1:length(marks)){
                    lines(x=c((1+counter*7),(7+counter*7)),
                          y=c((sum(plot_abs[[j]][i,1:5],na.rm=TRUE)*as.numeric(marks[k])),
                              (sum(plot_abs[[j]][i,1:5],na.rm=TRUE)*as.numeric(marks[k]))))
                }
                counter<-counter+1 
            }
        }  
        counter<-0
        for(i in 1:length(plot_abs[[j]][,1])){
            lines(x=c((1+counter*7),(7+counter*7)),
                  y=c(0,0))
            counter<-counter+1
        }
        for(i in 1:length(calling[[j]][,1])){
            if(is.na(calling[[j]][i,16])&&sum(plot_abs[[j]][i,1:5],na.rm=TRUE)!=0){
                if(calling[[j]][i,3]=="A"){
                    rect(xleft=1+(i-1)*7,xright=2+(i-1)*7,ybottom=0,
                         ytop=sum(plot_abs[[j]][i,1:5],na.rm=TRUE),border="grey50",
                         lty=3,lwd=4)
                }
                if(calling[[j]][i,3]=="C"){
                    rect(xleft=2+(i-1)*7,xright=3+(i-1)*7,ybottom=0,
                         ytop=sum(plot_abs[[j]][i,1:5],na.rm=TRUE),border="grey50",
                         lty=3,lwd=4)
                }
                if(calling[[j]][i,3]=="G"){
                    rect(xleft=3+(i-1)*7,xright=4+(i-1)*7,ybottom=0,
                         ytop=sum(plot_abs[[j]][i,1:5],na.rm=TRUE),border="grey50",
                         lty=3,lwd=4)
                }
                if(calling[[j]][i,3]=="T"){
                    rect(xleft=4+(i-1)*7,xright=5+(i-1)*7,ybottom=0,
                         ytop=sum(plot_abs[[j]][i,1:5],na.rm=TRUE),border="grey50",
                         lty=3,lwd=4)
                }
            }
            if(!is.na(calling[[j]][i,16])){
                if(calling[[j]][i,16]==calling[[j]][i,17]&&is.na(calling[[j]][i,18])){
                    if(calling[[j]][i,16]=="A"){
                        rect(xleft=1+(i-1)*7,xright=2+(i-1)*7,ybottom=0,
                             ytop=sum(plot_abs[[j]][i,1:5],na.rm=TRUE),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[j]][i,16]=="C"){
                        rect(xleft=2+(i-1)*7,xright=3+(i-1)*7,ybottom=0,
                             ytop=sum(plot_abs[[j]][i,1:5],na.rm=TRUE),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[j]][i,16]=="G"){
                        rect(xleft=3+(i-1)*7,xright=4+(i-1)*7,ybottom=0,
                             ytop=sum(plot_abs[[j]][i,1:5],na.rm=TRUE),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[j]][i,16]=="T"){
                        rect(xleft=4+(i-1)*7,xright=5+(i-1)*7,ybottom=0,
                             ytop=sum(plot_abs[[j]][i,1:5],na.rm=TRUE),
                             border="grey50",lty=3,lwd=4)
                    } 
                    if(calling[[j]][i,16]==""){
                        rect(xleft=5+(i-1)*7,xright=6+(i-1)*7,ybottom=0,
                             ytop=sum(plot_abs[[j]][i,1:5],na.rm=TRUE),
                             border="grey50",lty=3,lwd=4)
                    } 
                }
                if(calling[[j]][i,16]!=calling[[j]][i,17]
                   &&is.na(calling[[j]][i,18])){
                    if(calling[[j]][i,16]=="A"){
                        rect(xleft=1+(i-1)*7,xright=2+(i-1)*7,ybottom=0,
                             ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[j]][i,16]=="C"){
                        rect(xleft=2+(i-1)*7,xright=3+(i-1)*7,ybottom=0,
                             ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[j]][i,16]=="G"){
                        rect(xleft=3+(i-1)*7,xright=4+(i-1)*7,ybottom=0,
                             ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[j]][i,16]=="T"){
                        rect(xleft=4+(i-1)*7,xright=5+(i-1)*7,ybottom=0,
                             ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[j]][i,16]==""){
                        rect(xleft=5+(i-1)*7,xright=6+(i-1)*7,ybottom=0,
                             ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[j]][i,17]=="A"){
                        rect(xleft=1+(i-1)*7,xright=2+(i-1)*7,ybottom=0,
                             ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[j]][i,17]=="C"){
                        rect(xleft=2+(i-1)*7,xright=3+(i-1)*7,ybottom=0,
                             ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[j]][i,17]=="G"){
                        rect(xleft=3+(i-1)*7,xright=4+(i-1)*7,ybottom=0,
                             ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[j]][i,17]=="T"){
                        rect(xleft=4+(i-1)*7,xright=5+(i-1)*7,ybottom=0,
                             ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[j]][i,17]==""){
                        rect(xleft=5+(i-1)*7,xright=6+(i-1)*7,ybottom=0,
                             ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                             border="grey50",lty=3,lwd=4)
                    }
                }
                if(is.na(calling[[j]][i,19])&&!is.na(calling[[j]][i,18])){
                    rect(xleft=6+(i-1)*7,xright=7+(i-1)*7,ybottom=0,
                         ytop=sum(plot_abs[[j]][i,1:5],na.rm=TRUE),
                         border="grey50",lty=3,lwd=4)
                    if(calling[[j]][i,16]=="A"){
                        rect(xleft=1+(i-1)*7,xright=2+(i-1)*7,ybottom=0,
                             ytop=sum(plot_abs[[j]][i,1:5],na.rm=TRUE),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[j]][i,16]=="C"){
                        rect(xleft=2+(i-1)*7,xright=3+(i-1)*7,ybottom=0,
                             ytop=sum(plot_abs[[j]][i,1:5],na.rm=TRUE),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[j]][i,16]=="G"){
                        rect(xleft=3+(i-1)*7,xright=4+(i-1)*7,ybottom=0,
                             ytop=sum(plot_abs[[j]][i,1:5],na.rm=TRUE),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[j]][i,16]=="T"){
                        rect(xleft=4+(i-1)*7,xright=5+(i-1)*7,ybottom=0,
                             ytop=sum(plot_abs[[j]][i,1:5],na.rm=TRUE),
                             border="grey50",lty=3,lwd=4)
                    } 
                }
                if(!is.na(calling[[j]][i,19])&&!is.na(calling[[j]][i,18])){
                    rect(xleft=6+(i-1)*7,xright=7+(i-1)*7,ybottom=0,
                         ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                         border="grey50",lty=3,lwd=4)
                    if(calling[[j]][i,16]==calling[[j]][i,17]){
                        if(calling[[j]][i,16]=="A"){
                            rect(xleft=1+(i-1)*7,xright=2+(i-1)*7,ybottom=0,
                                 ytop=sum(plot_abs[[j]][i,1:5],na.rm=TRUE),
                                 border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[j]][i,16]=="C"){
                            rect(xleft=2+(i-1)*7,xright=3+(i-1)*7,ybottom=0,
                                 ytop=sum(plot_abs[[j]][i,1:5],na.rm=TRUE),
                                 border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[j]][i,16]=="G"){
                            rect(xleft=3+(i-1)*7,xright=4+(i-1)*7,ybottom=0,
                                 ytop=sum(plot_abs[[j]][i,1:5],na.rm=TRUE),
                                 border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[j]][i,16]=="T"){
                            rect(xleft=4+(i-1)*7,xright=5+(i-1)*7,ybottom=0,
                                 ytop=sum(plot_abs[[j]][i,1:5],na.rm=TRUE),
                                 border="grey50",lty=3,lwd=4)
                        } 
                    }
                    if(calling[[j]][i,16]!=calling[[j]][i,17]){
                        if(calling[[j]][i,16]=="A"){
                            rect(xleft=1+(i-1)*7,xright=2+(i-1)*7,ybottom=0,
                                 ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                                 border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[j]][i,16]=="C"){
                            rect(xleft=2+(i-1)*7,xright=3+(i-1)*7,ybottom=0,
                                 ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                                 border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[j]][i,16]=="G"){
                            rect(xleft=3+(i-1)*7,xright=4+(i-1)*7,ybottom=0,
                                 ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                                 border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[j]][i,16]=="T"){
                            rect(xleft=4+(i-1)*7,xright=5+(i-1)*7,ybottom=0,
                                 ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                                 border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[j]][i,16]==""){
                            rect(xleft=5+(i-1)*7,xright=6+(i-1)*7,ybottom=0,
                                 ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                                 border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[j]][i,17]=="A"){
                            rect(xleft=1+(i-1)*7,xright=2+(i-1)*7,ybottom=0,
                                 ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                                 border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[j]][i,17]=="C"){
                            rect(xleft=2+(i-1)*7,xright=3+(i-1)*7,ybottom=0,
                                 ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                                 border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[j]][i,17]=="G"){
                            rect(xleft=3+(i-1)*7,xright=4+(i-1)*7,ybottom=0,
                                 ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                                 border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[j]][i,17]=="T"){
                            rect(xleft=4+(i-1)*7,xright=5+(i-1)*7,ybottom=0,
                                 ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                                 border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[j]][i,17]==""){
                            rect(xleft=5+(i-1)*7,xright=6+(i-1)*7,ybottom=0,
                                 ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                                 border="grey50",lty=3,lwd=4)
                        }
                    }
                }
            }
        }
        
        text(x=-0.5*length(plot_abs[[j]][,1]),y=c(abs_max*0.8*(-0.1)),"Reference",cex=2)
        steps<-seq(1,length(plot_abs[[j]][,1])*0.65,
                   by=(length(plot_abs[[j]][,1])*0.65-1)/length(A_col))-1.2*length(plot_abs[[j]][,1])
        if(length(plot_abs[[j]][,1])==1){
            steps<-seq(1*0.65,1,by=(1-1*0.65)/length(A_col))-1.75*1
        }
        if(length(plot_abs[[j]][,1])==2){
            steps<-seq(0.7,2*0.65,by=(2*0.65-0.7)/length(A_col))-1.2*2
        }
        rect(xleft=steps[1:length(steps)-1],xright=steps[2:length(steps)],
             ybottom=abs_max*0.88,ytop=abs_max*0.93,lwd=0,col=A_col,border=NA)
        rect(xleft=steps[1:length(steps)-1],xright=steps[2:length(steps)],
             ybottom=abs_max*0.81,ytop=abs_max*0.86,lwd=0,col=C_col,border=NA)
        rect(xleft=steps[1:length(steps)-1],xright=steps[2:length(steps)],
             ybottom=abs_max*0.74,ytop=abs_max*0.79,lwd=0,col=G_col,border=NA)
        rect(xleft=steps[1:length(steps)-1],xright=steps[2:length(steps)],
             ybottom=abs_max*0.67,ytop=abs_max*0.72,lwd=0,col=T_col,border=NA)
        text(x=length(plot_abs[[j]][,1])*0.75-1.2*length(plot_abs[[j]][,1]),
             y=c(abs_max*0.90),"A",cex=2)
        text(x=length(plot_abs[[j]][,1])*0.75-1.2*length(plot_abs[[j]][,1]),
             y=c(abs_max*0.83),"C",cex=2)
        text(x=length(plot_abs[[j]][,1])*0.75-1.2*length(plot_abs[[j]][,1]),
             y=c(abs_max*0.76),"G",cex=2)
        text(x=length(plot_abs[[j]][,1])*0.75-1.2*length(plot_abs[[j]][,1]),
             y=c(abs_max*0.69),"T",cex=2)
        text(x=length(plot_abs[[j]][,1])*0.78-1.2*length(plot_abs[[j]][,1]),
             y=c(abs_max*0.62),"Del",cex=2)
        text(x=length(plot_abs[[j]][,1])*0.78-1.2*length(plot_abs[[j]][,1]),
             y=c(abs_max*0.55),"Ins",cex=2)
        text(x=length(plot_abs[[j]][,1])*0.98-1.2*length(plot_abs[[j]][,1]),
             y=c(abs_max*0.48),"Expected",cex=2)
        text(x=length(plot_abs[[j]][,1])*1.15-1.2*length(plot_abs[[j]][,1]),
             y=c(abs_max*0.97),"Mean quality",cex=2)
        if(length(plot_abs[[j]][,1])>2){
            rect(xleft=1-1.2*length(plot_abs[[j]][,1]),
                 xright=length(plot_abs[[j]][,1])*0.65-1.2*length(plot_abs[[j]][,1]),
                 ybottom=abs_max*0.60,ytop=abs_max*0.65,lwd=0,col="black")
            rect(xleft=1-1.2*length(plot_abs[[j]][,1]),
                 xright=length(plot_abs[[j]][,1])*0.65-1.2*length(plot_abs[[j]][,1]),
                 ybottom=abs_max*0.53,ytop=abs_max*0.58,lwd=4,border="purple",
                 col="white")
            rect(xleft=1-1.2*length(plot_abs[[j]][,1]),
                 xright=length(plot_abs[[j]][,1])*0.65-1.2*length(plot_abs[[j]][,1]),
                 ybottom=abs_max*0.46,ytop=abs_max*0.51,lwd=2,border="grey50",
                 col="white",lty=3)
            text(x=1-1.2*length(plot_abs[[j]][,1]),y=c(abs_max*0.97),
                 qual_lower_bound,cex=2)
            text(x=length(plot_abs[[j]][,1])*0.65-1.2*length(plot_abs[[j]][,1]),
                 y=c(abs_max*0.97),qual_upper_bound,cex=2)
        }
        if(length(plot_abs[[j]][,1])==1){
            rect(xleft=-1.1,xright=-0.75,ybottom=abs_max*0.60,ytop=abs_max*0.65,
                 lwd=0,col="black")
            rect(xleft=-1.1,xright=-0.75,ybottom=abs_max*0.53,ytop=abs_max*0.58,
                 lwd=4,border="purple",col="white")
            rect(xleft=-1.1,xright=-0.75,ybottom=abs_max*0.46,ytop=abs_max*0.51,
                 lwd=2,border="grey50",col="white",lty=3)
            text(x=-1.1,y=c(abs_max*0.97),qual_lower_bound,cex=2)
            text(x=-0.75,y=c(abs_max*0.97),qual_upper_bound,cex=2)
        }
        if(length(plot_abs[[j]][,1])==2){
            rect(xleft=-1.7,xright=-1.1,ybottom=abs_max*0.60,ytop=abs_max*0.65,
                 lwd=0,col="black")
            rect(xleft=-1.7,xright=-1.1,ybottom=abs_max*0.53,ytop=abs_max*0.58,
                 lwd=4,border="purple",col="white")
            rect(xleft=-1.7,xright=-1.1,ybottom=abs_max*0.46,ytop=abs_max*0.51,
                 lwd=2,border="grey50",col="white",lty=3)
            text(x=-1.7,y=c(abs_max*0.97),qual_lower_bound,cex=2)
            text(x=-1.1,y=c(abs_max*0.97),qual_upper_bound,cex=2)
        }
        if(directory2!=""){
            dev.off()   
        }
        return()
    }
    
    plotCountsPerSampleRelativeVcf_plot<-function(directory2,samples,plot_rel,
                                                  plot_col,rel_max,calling,
                                                  plot_insert,plot_col_insert,
                                                  A_col,C_col,G_col,T_col,plot_ref,
                                                  marks,qual_lower_bound,
                                                  qual_upper_bound,j){
        if(directory2!=""){
            png(filename=paste(directory2,"/",samples[j,1],".png",sep=""),width=1600,
                height=800)
        }
        par(mar=c(16,6,2,1)) 
        barplot(t(plot_rel[[j]]),beside=TRUE,col=t(plot_col[[j]]),
                xlim=c(-1.2*length(plot_rel[[j]][,1]),7*length(plot_rel[[j]][,1])),
                ylim=c((-1*0.1*rel_max),rel_max),ylab="Relative frequency",
                mgp=c(4.5,1,-2),main=samples[j,1],
                names.arg=paste(calling[[1]][,1],calling[[1]][,2],sep=";"),
                las=2,cex.axis=2,cex.names=2,cex.lab=2,cex.main=2)
        par(new=TRUE)
        barplot((t(plot_ref[[j]]))*0.1*rel_max,beside=TRUE,
                col=c("forestgreen","blue4","darkgoldenrod3","darkred","black","purple"),
                xlim=c(-1.2*length(plot_rel[[j]][,1]),7*length(plot_rel[[j]][,1])),
                ylim=c((-1*0.1*rel_max),rel_max),yaxt='n')
        for(i in 1:length(plot_rel[[j]][,1])){
            for(k in 4:1){
                rect(xleft=(7*i-0.8),xright=(7*i-0.15),ybottom=0,
                     ytop=plot_insert[[j]][i,k],border="white",col=plot_col_insert[[j]][i,k])
            }
            if(plot_insert[[j]][i,4]!=0)
                rect(xleft=(7*i-1),xright=(7*i),ybottom=0,
                     ytop=plot_insert[[j]][i,4],col=NA,border="purple",lwd=5)
        }
        counter<-0
        for(i in 1:length(plot_rel[[j]][,1])){
            if(length(marks)>0){
                for(k in 1:length(marks)){
                    lines(x=c((1+counter*7),(7+counter*7)),y=c(marks[k],marks[k]))
                }
                counter<-counter+1 
            }
        }
        counter<-0
        for(i in 1:length(plot_rel[[j]][,1])){
            lines(x=c((1+counter*7),(7+counter*7)),
                  y=c(0,0))
            counter<-counter+1
        }
        for(i in 1:length(calling[[j]][,1])){
            if(is.na(calling[[j]][i,16])&&sum(plot_rel[[j]][i,1:5],na.rm=TRUE)!=0){
                if(calling[[j]][i,3]=="A"){
                    rect(xleft=1+(i-1)*7,xright=2+(i-1)*7,ybottom=0,ytop=1,
                         border="grey50",lty=3,lwd=4)
                }
                if(calling[[j]][i,3]=="C"){
                    rect(xleft=2+(i-1)*7,xright=3+(i-1)*7,ybottom=0,ytop=1,
                         border="grey50",lty=3,lwd=4)
                }
                if(calling[[j]][i,3]=="G"){
                    rect(xleft=3+(i-1)*7,xright=4+(i-1)*7,ybottom=0,ytop=1,
                         border="grey50",lty=3,lwd=4)
                }
                if(calling[[j]][i,3]=="T"){
                    rect(xleft=4+(i-1)*7,xright=5+(i-1)*7,ybottom=0,ytop=1,
                         border="grey50",lty=3,lwd=4)
                }
            }
            if(!is.na(calling[[j]][i,16])){
                if(calling[[j]][i,16]==calling[[j]][i,17]&&is.na(calling[[j]][i,18])){
                    if(calling[[j]][i,16]=="A"){
                        rect(xleft=1+(i-1)*7,xright=2+(i-1)*7,ybottom=0,ytop=1,
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[j]][i,16]=="C"){
                        rect(xleft=2+(i-1)*7,xright=3+(i-1)*7,ybottom=0,ytop=1,
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[j]][i,16]=="G"){
                        rect(xleft=3+(i-1)*7,xright=4+(i-1)*7,ybottom=0,ytop=1,
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[j]][i,16]=="T"){
                        rect(xleft=4+(i-1)*7,xright=5+(i-1)*7,ybottom=0,ytop=1,
                             border="grey50",lty=3,lwd=4)
                    } 
                    if(calling[[j]][i,16]==""){
                        rect(xleft=5+(i-1)*7,xright=6+(i-1)*7,ybottom=0,ytop=1,
                             border="grey50",lty=3,lwd=4)
                    } 
                }
                if(calling[[j]][i,16]!=calling[[j]][i,17]&&is.na(calling[[j]][i,18])){
                    if(calling[[j]][i,16]=="A"){
                        rect(xleft=1+(i-1)*7,xright=2+(i-1)*7,ybottom=0,ytop=(1/2),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[j]][i,16]=="C"){
                        rect(xleft=2+(i-1)*7,xright=3+(i-1)*7,ybottom=0,ytop=(1/2),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[j]][i,16]=="G"){
                        rect(xleft=3+(i-1)*7,xright=4+(i-1)*7,ybottom=0,ytop=(1/2),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[j]][i,16]=="T"){
                        rect(xleft=4+(i-1)*7,xright=5+(i-1)*7,ybottom=0,ytop=(1/2),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[j]][i,16]==""){
                        rect(xleft=5+(i-1)*7,xright=6+(i-1)*7,ybottom=0,ytop=(1/2),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[j]][i,17]=="A"){
                        rect(xleft=1+(i-1)*7,xright=2+(i-1)*7,ybottom=0,ytop=(1/2),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[j]][i,17]=="C"){
                        rect(xleft=2+(i-1)*7,xright=3+(i-1)*7,ybottom=0,ytop=(1/2),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[j]][i,17]=="G"){
                        rect(xleft=3+(i-1)*7,xright=4+(i-1)*7,ybottom=0,ytop=(1/2),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[j]][i,17]=="T"){
                        rect(xleft=4+(i-1)*7,xright=5+(i-1)*7,ybottom=0,ytop=(1/2),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[j]][i,17]==""){
                        rect(xleft=5+(i-1)*7,xright=6+(i-1)*7,ybottom=0,ytop=(1/2),
                             border="grey50",lty=3,lwd=4)
                    }
                }
                if(is.na(calling[[j]][i,19])&&!is.na(calling[[j]][i,18])){
                    rect(xleft=6+(i-1)*7,xright=7+(i-1)*7,ybottom=0,ytop=1,
                         border="grey50",lty=3,lwd=4)
                    if(calling[[j]][i,16]=="A"){
                        rect(xleft=1+(i-1)*7,xright=2+(i-1)*7,ybottom=0,ytop=1,
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[j]][i,16]=="C"){
                        rect(xleft=2+(i-1)*7,xright=3+(i-1)*7,ybottom=0,ytop=1,
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[j]][i,16]=="G"){
                        rect(xleft=3+(i-1)*7,xright=4+(i-1)*7,ybottom=0,ytop=1,
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[j]][i,16]=="T"){
                        rect(xleft=4+(i-1)*7,xright=5+(i-1)*7,ybottom=0,ytop=1,
                             border="grey50",lty=3,lwd=4)
                    } 
                }
                if(!is.na(calling[[j]][i,19])&&!is.na(calling[[j]][i,18])){
                    rect(xleft=6+(i-1)*7,xright=7+(i-1)*7,ybottom=0,ytop=(1/2),
                         border="grey50",lty=3,lwd=4)
                    if(calling[[j]][i,16]==calling[[j]][i,17]){
                        if(calling[[j]][i,16]=="A"){
                            rect(xleft=1+(i-1)*7,xright=2+(i-1)*7,ybottom=0,ytop=1,
                                 border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[j]][i,16]=="C"){
                            rect(xleft=2+(i-1)*7,xright=3+(i-1)*7,ybottom=0,ytop=1,
                                 border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[j]][i,16]=="G"){
                            rect(xleft=3+(i-1)*7,xright=4+(i-1)*7,ybottom=0,ytop=1,
                                 border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[j]][i,16]=="T"){
                            rect(xleft=4+(i-1)*7,xright=5+(i-1)*7,ybottom=0,ytop=1,
                                 border="grey50",lty=3,lwd=4)
                        } 
                    }
                    if(calling[[j]][i,16]!=calling[[j]][i,17]){
                        if(calling[[j]][i,16]=="A"){
                            rect(xleft=1+(i-1)*7,xright=2+(i-1)*7,ybottom=0,
                                 ytop=(1/2),border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[j]][i,16]=="C"){
                            rect(xleft=2+(i-1)*7,xright=3+(i-1)*7,ybottom=0,
                                 ytop=(1/2),border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[j]][i,16]=="G"){
                            rect(xleft=3+(i-1)*7,xright=4+(i-1)*7,ybottom=0,
                                 ytop=(1/2),border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[j]][i,16]=="T"){
                            rect(xleft=4+(i-1)*7,xright=5+(i-1)*7,ybottom=0,
                                 ytop=(1/2),border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[j]][i,16]==""){
                            rect(xleft=5+(i-1)*7,xright=6+(i-1)*7,ybottom=0,
                                 ytop=(1/2),border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[j]][i,17]=="A"){
                            rect(xleft=1+(i-1)*7,xright=2+(i-1)*7,ybottom=0,
                                 ytop=(1/2),border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[j]][i,17]=="C"){
                            rect(xleft=2+(i-1)*7,xright=3+(i-1)*7,ybottom=0,
                                 ytop=(1/2),border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[j]][i,17]=="G"){
                            rect(xleft=3+(i-1)*7,xright=4+(i-1)*7,ybottom=0,
                                 ytop=(1/2),border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[j]][i,17]=="T"){
                            rect(xleft=4+(i-1)*7,xright=5+(i-1)*7,ybottom=0,
                                 ytop=(1/2),border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[j]][i,17]==""){
                            rect(xleft=5+(i-1)*7,xright=6+(i-1)*7,ybottom=0,
                                 ytop=(1/2),border="grey50",lty=3,lwd=4)
                        }
                    }
                }
            }
        }
        
        text(x=-0.5*length(plot_rel[[j]][,1]),y=c(rel_max*0.8*(-0.1)),"Reference",cex=2)
        steps<-seq(1,length(plot_rel[[j]][,1])*0.65,
                   by=(length(plot_rel[[j]][,1])*0.65-1)/length(A_col))-1.2*length(plot_rel[[j]][,1])
        if(length(plot_rel[[j]][,1])==1){
            steps<-seq(1*0.65,1,by=(1-1*0.65)/length(A_col))-1.75*1
        }
        if(length(plot_rel[[j]][,1])==2){
            steps<-seq(0.7,2*0.65,by=(2*0.65-0.7)/length(A_col))-1.2*2
        }
        rect(xleft=steps[1:length(steps)-1],xright=steps[2:length(steps)],
             ybottom=rel_max*0.88,ytop=rel_max*0.93,lwd=0,col=A_col,border=NA)
        rect(xleft=steps[1:length(steps)-1],xright=steps[2:length(steps)],
             ybottom=rel_max*0.81,ytop=rel_max*0.86,lwd=0,col=C_col,border=NA)
        rect(xleft=steps[1:length(steps)-1],xright=steps[2:length(steps)],
             ybottom=rel_max*0.74,ytop=rel_max*0.79,lwd=0,col=G_col,border=NA)
        rect(xleft=steps[1:length(steps)-1],xright=steps[2:length(steps)],
             ybottom=rel_max*0.67,ytop=rel_max*0.72,lwd=0,col=T_col,border=NA)
        text(x=length(plot_rel[[j]][,1])*0.75-1.2*length(plot_rel[[j]][,1]),
             y=c(rel_max*0.90),"A",cex=2)
        text(x=length(plot_rel[[j]][,1])*0.75-1.2*length(plot_rel[[j]][,1]),
             y=c(rel_max*0.83),"C",cex=2)
        text(x=length(plot_rel[[j]][,1])*0.75-1.2*length(plot_rel[[j]][,1]),
             y=c(rel_max*0.76),"G",cex=2)
        text(x=length(plot_rel[[j]][,1])*0.75-1.2*length(plot_rel[[j]][,1]),
             y=c(rel_max*0.69),"T",cex=2)
        text(x=length(plot_rel[[j]][,1])*0.78-1.2*length(plot_rel[[j]][,1]),
             y=c(rel_max*0.62),"Del",cex=2)
        text(x=length(plot_rel[[j]][,1])*0.78-1.2*length(plot_rel[[j]][,1]),
             y=c(rel_max*0.55),"Ins",cex=2)
        text(x=length(plot_rel[[j]][,1])*0.98-1.2*length(plot_rel[[j]][,1]),
             y=c(rel_max*0.48),"Expected",cex=2)
        text(x=length(plot_rel[[j]][,1])*1.15-1.2*length(plot_rel[[j]][,1]),
             y=c(rel_max*0.97),"Mean quality",cex=2)
        if(length(plot_rel[[j]][,1])>2){
            rect(xleft=1-1.2*length(plot_rel[[j]][,1]),
                 xright=length(plot_rel[[j]][,1])*0.65-1.2*length(plot_rel[[j]][,1]),
                 ybottom=rel_max*0.60,ytop=rel_max*0.65,lwd=0,col="black")
            rect(xleft=1-1.2*length(plot_rel[[j]][,1]),
                 xright=length(plot_rel[[j]][,1])*0.65-1.2*length(plot_rel[[j]][,1]),
                 ybottom=rel_max*0.53,ytop=rel_max*0.58,lwd=4,border="purple",
                 col="white")
            rect(xleft=1-1.2*length(plot_rel[[j]][,1]),
                 xright=length(plot_rel[[j]][,1])*0.65-1.2*length(plot_rel[[j]][,1]),
                 ybottom=rel_max*0.46,ytop=rel_max*0.51,lwd=2,border="grey50",
                 col="white",lty=3)
            text(x=1-1.2*length(plot_rel[[j]][,1]),y=c(rel_max*0.97),
                 qual_lower_bound,cex=2)
            text(x=length(plot_rel[[j]][,1])*0.65-1.2*length(plot_rel[[j]][,1]),
                 y=c(rel_max*0.97),qual_upper_bound,cex=2)
        }
        if(length(plot_rel[[j]][,1])==1){
            rect(xleft=-1.1,xright=-0.75,ybottom=rel_max*0.60,ytop=rel_max*0.65,
                 lwd=0,col="black")
            rect(xleft=-1.1,xright=-0.75,ybottom=rel_max*0.53,ytop=rel_max*0.58,
                 lwd=4,border="purple",col="white")
            rect(xleft=-1.1,xright=-0.75,ybottom=rel_max*0.46,ytop=rel_max*0.51,
                 lwd=2,border="grey50",col="white",lty=3)
            text(x=-1.1,y=c(rel_max*0.97),qual_lower_bound,cex=2)
            text(x=-0.75,y=c(rel_max*0.97),qual_upper_bound,cex=2)
        }
        if(length(plot_rel[[j]][,1])==2){
            rect(xleft=-1.7,xright=-1.1,ybottom=rel_max*0.60,ytop=rel_max*0.65,
                 lwd=0,col="black")
            rect(xleft=-1.7,xright=-1.1,ybottom=rel_max*0.53,ytop=rel_max*0.58,
                 lwd=4,border="purple",col="white")
            rect(xleft=-1.7,xright=-1.1,ybottom=rel_max*0.46,ytop=rel_max*0.51,
                 lwd=2,border="grey50",col="white",lty=3)
            text(x=-1.7,y=c(rel_max*0.97),qual_lower_bound,cex=2)
            text(x=-1.1,y=c(rel_max*0.97),qual_upper_bound,cex=2)
        }
        if(directory2!=""){
            dev.off()
        }
        return()
    }
    
    plotCountsPerTargetVcf_plot<-function(directory2,samples,plot_abs,plot_col,
                                          abs_max,calling,plot_insert,
                                          plot_col_insert,A_col,C_col,G_col,T_col,
                                          plot_ref,marks,qual_lower_bound,
                                          qual_upper_bound,counter,j){
        if(directory2!=""){
            if(j>1&&calling[[1]][j-1,1]==calling[[1]][j,1]
               &&calling[[1]][j-1,2]==calling[[1]][j,2]){
                counter<-counter+1
                png(filename=paste(directory2,"/",calling[[1]][j,1],";",calling[[1]][j,2],"_",counter,".png",sep=""),
                    width=1600,height=600)
            }
            if(j==1||calling[[1]][j-1,1]!=calling[[1]][j,1]||calling[[1]][j-1,2]!=calling[[1]][j,2]){
                png(filename=paste(directory2,"/",calling[[1]][j,1],";",calling[[1]][j,2],".png",sep=""),
                    width=1600,height=600)
                counter<-0
            }
        }
        par(mar=c(16,6,2,1)) 
        if(counter==0){
            barplot(t(plot_abs[[j]]),beside=TRUE,col=t(plot_col[[j]]),
                    xlim=c(-1.2*length(plot_abs[[j]][,1]),7*length(plot_abs[[j]][,1])),
                    ylim=c((-1*0.1*abs_max),abs_max),ylab="Number of reads",
                    main=paste(calling[[1]][j,1],";",calling[[1]][j,2],sep=""),
                    names.arg=t(samples[1]),las=2,cex.axis=2,cex.names=2,cex.lab=2,
                    cex.main=2)
        }
        if(counter>0){
            barplot(t(plot_abs[[j]]),beside=TRUE,col=t(plot_col[[j]]),
                    xlim=c(-1.2*length(plot_abs[[j]][,1]),7*length(plot_abs[[j]][,1])),
                    ylim=c((-1*0.1*abs_max),abs_max),ylab="Number of reads",
                    main=paste(calling[[1]][j,1],";",calling[[1]][j,2],"_",counter,sep=""),
                    names.arg=t(samples[1]),las=2,cex.axis=2,cex.names=2,cex.lab=2,
                    cex.main=2)
        }
        par(new=TRUE)
        barplot((t(plot_ref[[j]]))*0.1*abs_max,beside=TRUE,
                col=c("forestgreen","blue4","darkgoldenrod3","darkred","black","purple"),
                xlim=c(-1.2*length(plot_abs[[j]][,1]),7*length(plot_abs[[j]][,1])),
                ylim=c((-1*0.1*abs_max),abs_max),yaxt='n')
        for(i in 1:length(plot_abs[[j]][,1])){
            for(k in 4:1){
                rect(xleft=(7*i-0.8),xright=(7*i-0.15),ybottom=0,
                     ytop=plot_insert[[j]][i,k],border=NA,col=plot_col_insert[[j]][i,k])
            }
            if(plot_insert[[j]][i,4]!=0)
                rect(xleft=(7*i-1),xright=(7*i),ybottom=0,
                     ytop=plot_insert[[j]][i,4],col=NA,border="purple",lwd=5)
        }
        counter2<-0
        for(i in 1:length(plot_abs[[j]][,1])){
            if(length(marks)>0){
                for(k in 1:length(marks)){
                    lines(x=c((1+counter2*7),(7+counter2*7)),
                          y=c((sum(plot_abs[[j]][i,1:5],na.rm=TRUE)*as.numeric(marks[k])),
                              (sum(plot_abs[[j]][i,1:5],na.rm=TRUE)*as.numeric(marks[k]))))
                }
                counter2<-counter2+1
            }
        }
        counter2<-0
        for(i in 1:length(plot_abs[[j]][,1])){
            lines(x=c((1+counter2*7),(7+counter2*7)),
                  y=c(0,0))
            counter2<-counter2+1
        }
        for(i in 1:length(calling)){
            if(is.na(calling[[i]][j,16])&&sum(plot_abs[[j]][i,1:5],na.rm=TRUE)!=0){
                if(calling[[i]][j,3]=="A"){
                    rect(xleft=1+(i-1)*7,xright=2+(i-1)*7,ybottom=0,
                         ytop=sum(plot_abs[[j]][i,1:5],na.rm=TRUE),border="grey50",
                         lty=3,lwd=4)
                }
                if(calling[[i]][j,3]=="C"){
                    rect(xleft=2+(i-1)*7,xright=3+(i-1)*7,ybottom=0,
                         ytop=sum(plot_abs[[j]][i,1:5],na.rm=TRUE),border="grey50",
                         lty=3,lwd=4)
                }
                if(calling[[i]][j,3]=="G"){
                    rect(xleft=3+(i-1)*7,xright=4+(i-1)*7,ybottom=0,
                         ytop=sum(plot_abs[[j]][i,1:5],na.rm=TRUE),border="grey50",
                         lty=3,lwd=4)
                }
                if(calling[[i]][j,3]=="T"){
                    rect(xleft=4+(i-1)*7,xright=5+(i-1)*7,ybottom=0,
                         ytop=sum(plot_abs[[j]][i,1:5],na.rm=TRUE),border="grey50",
                         lty=3,lwd=4)
                }
            }
            if(!is.na(calling[[i]][j,16])){
                if(calling[[i]][j,16]==calling[[i]][j,17]&&is.na(calling[[i]][j,18])){
                    if(calling[[i]][j,16]=="A"){
                        rect(xleft=1+(i-1)*7,xright=2+(i-1)*7,ybottom=0,
                             ytop=sum(plot_abs[[j]][i,1:5],na.rm=TRUE),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[i]][j,16]=="C"){
                        rect(xleft=2+(i-1)*7,xright=3+(i-1)*7,ybottom=0,
                             ytop=sum(plot_abs[[j]][i,1:5],na.rm=TRUE),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[i]][j,16]=="G"){
                        rect(xleft=3+(i-1)*7,xright=4+(i-1)*7,ybottom=0,
                             ytop=sum(plot_abs[[j]][i,1:5],na.rm=TRUE),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[i]][j,16]=="T"){
                        rect(xleft=4+(i-1)*7,xright=5+(i-1)*7,ybottom=0,
                             ytop=sum(plot_abs[[j]][i,1:5],na.rm=TRUE),
                             border="grey50",lty=3,lwd=4)
                    } 
                    if(calling[[i]][j,16]==""){
                        rect(xleft=5+(i-1)*7,xright=6+(i-1)*7,ybottom=0,
                             ytop=sum(plot_abs[[j]][i,1:5],na.rm=TRUE),
                             border="grey50",lty=3,lwd=4)
                    } 
                }
                if(calling[[i]][j,16]!=calling[[i]][j,17]&&is.na(calling[[i]][j,18])){
                    if(calling[[i]][j,16]=="A"){
                        rect(xleft=1+(i-1)*7,xright=2+(i-1)*7,ybottom=0,
                             ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[i]][j,16]=="C"){
                        rect(xleft=2+(i-1)*7,xright=3+(i-1)*7,ybottom=0,
                             ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[i]][j,16]=="G"){
                        rect(xleft=3+(i-1)*7,xright=4+(i-1)*7,ybottom=0,
                             ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[i]][j,16]=="T"){
                        rect(xleft=4+(i-1)*7,xright=5+(i-1)*7,ybottom=0,
                             ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[i]][j,16]==""){
                        rect(xleft=5+(i-1)*7,xright=6+(i-1)*7,ybottom=0,
                             ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[i]][j,17]=="A"){
                        rect(xleft=1+(i-1)*7,xright=2+(i-1)*7,ybottom=0,
                             ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[i]][j,17]=="C"){
                        rect(xleft=2+(i-1)*7,xright=3+(i-1)*7,ybottom=0,
                             ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[i]][j,17]=="G"){
                        rect(xleft=3+(i-1)*7,xright=4+(i-1)*7,ybottom=0,
                             ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[i]][j,17]=="T"){
                        rect(xleft=4+(i-1)*7,xright=5+(i-1)*7,ybottom=0,
                             ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[i]][j,17]==""){
                        rect(xleft=5+(i-1)*7,xright=6+(i-1)*7,ybottom=0,
                             ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                             border="grey50",lty=3,lwd=4)
                    }
                }
                if(is.na(calling[[i]][j,19])&&!is.na(calling[[i]][j,18])){
                    rect(xleft=6+(i-1)*7,xright=7+(i-1)*7,ybottom=0,
                         ytop=sum(plot_abs[[j]][i,1:5],na.rm=TRUE),
                         border="grey50",lty=3,lwd=4)
                    if(calling[[i]][j,16]=="A"){
                        rect(xleft=1+(i-1)*7,xright=2+(i-1)*7,ybottom=0,
                             ytop=sum(plot_abs[[j]][i,1:5],na.rm=TRUE),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[i]][j,16]=="C"){
                        rect(xleft=2+(i-1)*7,xright=3+(i-1)*7,ybottom=0,
                             ytop=sum(plot_abs[[j]][i,1:5],na.rm=TRUE),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[i]][j,16]=="G"){
                        rect(xleft=3+(i-1)*7,xright=4+(i-1)*7,ybottom=0,
                             ytop=sum(plot_abs[[j]][i,1:5],na.rm=TRUE),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[i]][j,16]=="T"){
                        rect(xleft=4+(i-1)*7,xright=5+(i-1)*7,ybottom=0,
                             ytop=sum(plot_abs[[j]][i,1:5],na.rm=TRUE),
                             border="grey50",lty=3,lwd=4)
                    } 
                }
                if(!is.na(calling[[i]][j,19])&&!is.na(calling[[i]][j,18])){
                    rect(xleft=6+(i-1)*7,xright=7+(i-1)*7,ybottom=0,
                         ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                         border="grey50",lty=3,lwd=4)
                    if(calling[[i]][j,16]==calling[[i]][j,17]){
                        if(calling[[i]][j,16]=="A"){
                            rect(xleft=1+(i-1)*7,xright=2+(i-1)*7,ybottom=0,
                                 ytop=sum(plot_abs[[j]][i,1:5],na.rm=TRUE),
                                 border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[i]][j,16]=="C"){
                            rect(xleft=2+(i-1)*7,xright=3+(i-1)*7,ybottom=0,
                                 ytop=sum(plot_abs[[j]][i,1:5],na.rm=TRUE),
                                 border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[i]][j,16]=="G"){
                            rect(xleft=3+(i-1)*7,xright=4+(i-1)*7,ybottom=0,
                                 ytop=sum(plot_abs[[j]][i,1:5],na.rm=TRUE),
                                 border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[i]][j,16]=="T"){
                            rect(xleft=4+(i-1)*7,xright=5+(i-1)*7,ybottom=0,
                                 ytop=sum(plot_abs[[j]][i,1:5],na.rm=TRUE),
                                 border="grey50",lty=3,lwd=4)
                        } 
                    }
                    if(calling[[i]][j,16]!=calling[[i]][j,17]){
                        if(calling[[i]][j,16]=="A"){
                            rect(xleft=1+(i-1)*7,xright=2+(i-1)*7,ybottom=0,
                                 ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                                 border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[i]][j,16]=="C"){
                            rect(xleft=2+(i-1)*7,xright=3+(i-1)*7,ybottom=0,
                                 ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                                 border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[i]][j,16]=="G"){
                            rect(xleft=3+(i-1)*7,xright=4+(i-1)*7,ybottom=0,
                                 ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                                 border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[i]][j,16]=="T"){
                            rect(xleft=4+(i-1)*7,xright=5+(i-1)*7,ybottom=0,
                                 ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                                 border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[i]][j,16]==""){
                            rect(xleft=5+(i-1)*7,xright=6+(i-1)*7,ybottom=0,
                                 ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                                 border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[i]][j,17]=="A"){
                            rect(xleft=1+(i-1)*7,xright=2+(i-1)*7,ybottom=0,
                                 ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                                 border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[i]][j,17]=="C"){
                            rect(xleft=2+(i-1)*7,xright=3+(i-1)*7,ybottom=0,
                                 ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                                 border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[i]][j,17]=="G"){
                            rect(xleft=3+(i-1)*7,xright=4+(i-1)*7,ybottom=0,
                                 ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                                 border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[i]][j,17]=="T"){
                            rect(xleft=4+(i-1)*7,xright=5+(i-1)*7,ybottom=0,
                                 ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                                 border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[i]][j,17]==""){
                            rect(xleft=5+(i-1)*7,xright=6+(i-1)*7,ybottom=0,
                                 ytop=(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)/2),
                                 border="grey50",lty=3,lwd=4)
                        }
                    }
                }
            }
        }    
        
        text(x=-0.5*length(plot_abs[[j]][,1]),
             y=c(abs_max*0.8*(-0.1)),"Reference",cex=2)
        steps<-seq(1,length(plot_abs[[j]][,1])*0.65,
                   by=(length(plot_abs[[j]][,1])*0.65-1)/length(A_col))-1.2*length(plot_abs[[j]][,1])
        if(length(plot_abs[[j]][,1])==1){
            steps<-seq(1*0.65,1,by=(1-1*0.65)/length(A_col))-1.75*1
        }
        if(length(plot_abs[[j]][,1])==2){
            steps<-seq(0.7,2*0.65,by=(2*0.65-0.7)/length(A_col))-1.2*2
        }
        rect(xleft=steps[1:length(steps)-1],xright=steps[2:length(steps)],
             ybottom=abs_max*0.88,ytop=abs_max*0.93,lwd=0,col=A_col,border=NA)
        rect(xleft=steps[1:length(steps)-1],xright=steps[2:length(steps)],
             ybottom=abs_max*0.81,ytop=abs_max*0.86,lwd=0,col=C_col,border=NA)
        rect(xleft=steps[1:length(steps)-1],xright=steps[2:length(steps)],
             ybottom=abs_max*0.74,ytop=abs_max*0.79,lwd=0,col=G_col,border=NA)
        rect(xleft=steps[1:length(steps)-1],xright=steps[2:length(steps)],
             ybottom=abs_max*0.67,ytop=abs_max*0.72,lwd=0,col=T_col,border=NA)
        text(x=length(plot_abs[[j]][,1])*0.75-1.2*length(plot_abs[[j]][,1]),
             y=c(abs_max*0.90),"A",cex=2)
        text(x=length(plot_abs[[j]][,1])*0.75-1.2*length(plot_abs[[j]][,1]),
             y=c(abs_max*0.83),"C",cex=2)
        text(x=length(plot_abs[[j]][,1])*0.75-1.2*length(plot_abs[[j]][,1]),
             y=c(abs_max*0.76),"G",cex=2)
        text(x=length(plot_abs[[j]][,1])*0.75-1.2*length(plot_abs[[j]][,1]),
             y=c(abs_max*0.69),"T",cex=2)
        text(x=length(plot_abs[[j]][,1])*0.78-1.2*length(plot_abs[[j]][,1]),
             y=c(abs_max*0.62),"Del",cex=2)
        text(x=length(plot_abs[[j]][,1])*0.78-1.2*length(plot_abs[[j]][,1]),
             y=c(abs_max*0.55),"Ins",cex=2)
        text(x=length(plot_abs[[j]][,1])*0.98-1.2*length(plot_abs[[j]][,1]),
             y=c(abs_max*0.48),"Expected",cex=2)
        text(x=length(plot_abs[[j]][,1])*1.15-1.2*length(plot_abs[[j]][,1]),
             y=c(abs_max*0.97),"Mean quality",cex=2)
        if(length(plot_abs[[j]][,1])>2){
            rect(xleft=1-1.2*length(plot_abs[[j]][,1]),
                 xright=length(plot_abs[[j]][,1])*0.65-1.2*length(plot_abs[[j]][,1]),
                 ybottom=abs_max*0.60,ytop=abs_max*0.65,lwd=0,col="black")
            rect(xleft=1-1.2*length(plot_abs[[j]][,1]),
                 xright=length(plot_abs[[j]][,1])*0.65-1.2*length(plot_abs[[j]][,1]),
                 ybottom=abs_max*0.53,ytop=abs_max*0.58,lwd=4,border="purple",
                 col="white")
            rect(xleft=1-1.2*length(plot_abs[[j]][,1]),
                 xright=length(plot_abs[[j]][,1])*0.65-1.2*length(plot_abs[[j]][,1]),
                 ybottom=abs_max*0.46,ytop=abs_max*0.51,lwd=2,border="grey50",
                 col="white",lty=3)
            text(x=1-1.2*length(plot_abs[[j]][,1]),y=c(abs_max*0.97),
                 qual_lower_bound,cex=2)
            text(x=length(plot_abs[[j]][,1])*0.65-1.2*length(plot_abs[[j]][,1]),
                 y=c(abs_max*0.97),qual_upper_bound,cex=2)
        }
        if(length(plot_abs[[j]][,1])==1){
            rect(xleft=-1.1,xright=-0.75,ybottom=abs_max*0.60,ytop=abs_max*0.65,
                 lwd=0,col="black")
            rect(xleft=-1.1,xright=-0.75,ybottom=abs_max*0.53,ytop=abs_max*0.58,
                 lwd=4,border="purple",col="white")
            rect(xleft=-1.1,xright=-0.75,ybottom=abs_max*0.46,ytop=abs_max*0.51,
                 lwd=2,border="grey50",col="white",lty=3)
            text(x=-1.1,y=c(abs_max*0.97),qual_lower_bound,cex=2)
            text(x=-0.75,y=c(abs_max*0.97),qual_upper_bound,cex=2)
        }
        if(length(plot_abs[[j]][,1])==2){
            rect(xleft=-1.7,xright=-1.1,ybottom=abs_max*0.60,ytop=abs_max*0.65,
                 lwd=0,col="black")
            rect(xleft=-1.7,xright=-1.1,ybottom=abs_max*0.53,ytop=abs_max*0.58,
                 lwd=4,border="purple",col="white")
            rect(xleft=-1.7,xright=-1.1,ybottom=abs_max*0.46,ytop=abs_max*0.51,
                 lwd=2,border="grey50",col="white",lty=3)
            text(x=-1.7,y=c(abs_max*0.97),qual_lower_bound,cex=2)
            text(x=-1.1,y=c(abs_max*0.97),qual_upper_bound,cex=2)
        }
        if(directory2!=""){
            dev.off()
        }
        return(counter)
    }
    
    plotCountsPerTargetRelativeVcf_plot<-function(directory2,samples,plot_rel,
                                                  plot_col,rel_max,calling,
                                                  plot_insert,plot_col_insert,
                                                  A_col,C_col,G_col,T_col,plot_ref,
                                                  marks,qual_lower_bound,
                                                  qual_upper_bound,counter,j){
        if(directory2!=""){
            if(j>1&&calling[[1]][j-1,1]==calling[[1]][j,1]
               &&calling[[1]][j-1,2]==calling[[1]][j,2]){
                counter<-counter+1
                png(filename=paste(directory2,"/",calling[[1]][j,1],";",calling[[1]][j,2],"_",counter,".png",sep=""),
                    width=1600,height=600)
            }
            if(j==1||calling[[1]][j-1,1]!=calling[[1]][j,1]
               ||calling[[1]][j-1,2]!=calling[[1]][j,2]){
                png(filename=paste(directory2,"/",calling[[1]][j,1],";",calling[[1]][j,2],".png",sep=""),
                    width=1600,height=600)
                counter<-0
            }
        }
        par(mar=c(16,6,2,1)) 
        if(counter==0){
            barplot(t(plot_rel[[j]]),beside=TRUE,col=t(plot_col[[j]]),
                    xlim=c(-1.2*length(plot_rel[[j]][,1]),7*length(plot_rel[[j]][,1])),
                    ylim=c((-1*0.1*rel_max),rel_max),ylab="Relative frequency",
                    main=paste(calling[[1]][j,1],";",calling[[1]][j,2],sep=""),
                    names.arg=t(samples[1]),las=2,cex.axis=2,cex.names=2,cex.lab=2,
                    cex.main=2)
        }
        if(counter>0){
            barplot(t(plot_rel[[j]]),beside=TRUE,col=t(plot_col[[j]]),
                    xlim=c(-1.2*length(plot_rel[[j]][,1]),7*length(plot_rel[[j]][,1])),
                    ylim=c((-1*0.1*rel_max),rel_max),ylab="Relative frequency",
                    main=paste(calling[[1]][j,1],";",calling[[1]][j,2],"_",counter,sep=""),
                    names.arg=t(samples[1]),las=2,cex.axis=2,cex.names=2,cex.lab=2,
                    cex.main=2)
        }
        par(new=TRUE)
        barplot((t(plot_ref[[j]]))*0.1*rel_max,beside=TRUE,
                col=c("forestgreen","blue4","darkgoldenrod3","darkred","black","purple"),
                xlim=c(-1.2*length(plot_rel[[j]][,1]),7*length(plot_rel[[j]][,1])),
                ylim=c((-1*0.1*rel_max),rel_max),yaxt='n')
        for(i in 1:length(plot_rel[[j]][,1])){
            for(k in 4:1){
                rect(xleft=(7*i-0.8),xright=(7*i-0.15),ybottom=0,
                     ytop=plot_insert[[j]][i,k],border="white",col=plot_col_insert[[j]][i,k])
            }
            if(plot_insert[[j]][i,4]!=0)
                rect(xleft=(7*i-1),xright=(7*i),ybottom=0,
                     ytop=plot_insert[[j]][i,4],col=NA,border="purple",lwd=5)
        }
        counter2<-0
        for(i in 1:length(plot_rel[[j]][,1])){
            if(length(marks)>0){
                for(k in 1:length(marks)){
                    lines(x=c((1+counter2*7),(7+counter2*7)),y=c(marks[k],marks[k]))
                }
                counter2<-counter2+1
            }
        }
        counter2<-0
        for(i in 1:length(plot_rel[[j]][,1])){
            lines(x=c((1+counter2*7),(7+counter2*7)),
                  y=c(0,0))
            counter2<-counter2+1
        }
        for(i in 1:length(calling)){
            if(is.na(calling[[i]][j,16])&&sum(plot_rel[[j]][i,1:5],na.rm=TRUE)!=0){
                if(calling[[i]][j,3]=="A"){
                    rect(xleft=1+(i-1)*7,xright=2+(i-1)*7,ybottom=0,ytop=1,
                         border="grey50",lty=3,lwd=4)
                }
                if(calling[[i]][j,3]=="C"){
                    rect(xleft=2+(i-1)*7,xright=3+(i-1)*7,ybottom=0,ytop=1,
                         border="grey50",lty=3,lwd=4)
                }
                if(calling[[i]][j,3]=="G"){
                    rect(xleft=3+(i-1)*7,xright=4+(i-1)*7,ybottom=0,ytop=1,
                         border="grey50",lty=3,lwd=4)
                }
                if(calling[[i]][j,3]=="T"){
                    rect(xleft=4+(i-1)*7,xright=5+(i-1)*7,ybottom=0,ytop=1,
                         border="grey50",lty=3,lwd=4)
                }
            }
            if(!is.na(calling[[i]][j,16])){
                if(calling[[i]][j,16]==calling[[i]][j,17]
                   &&is.na(calling[[i]][j,18])){
                    if(calling[[i]][j,16]=="A"){
                        rect(xleft=1+(i-1)*7,xright=2+(i-1)*7,ybottom=0,ytop=1,
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[i]][j,16]=="C"){
                        rect(xleft=2+(i-1)*7,xright=3+(i-1)*7,ybottom=0,ytop=1,
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[i]][j,16]=="G"){
                        rect(xleft=3+(i-1)*7,xright=4+(i-1)*7,ybottom=0,ytop=1,
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[i]][j,16]=="T"){
                        rect(xleft=4+(i-1)*7,xright=5+(i-1)*7,ybottom=0,ytop=1,
                             border="grey50",lty=3,lwd=4)
                    } 
                    if(calling[[i]][j,16]==""){
                        rect(xleft=5+(i-1)*7,xright=6+(i-1)*7,ybottom=0,ytop=1,
                             border="grey50",lty=3,lwd=4)
                    } 
                }
                if(calling[[i]][j,16]!=calling[[i]][j,17]&&is.na(calling[[i]][j,18])){
                    if(calling[[i]][j,16]=="A"){
                        rect(xleft=1+(i-1)*7,xright=2+(i-1)*7,ybottom=0,ytop=(1/2),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[i]][j,16]=="C"){
                        rect(xleft=2+(i-1)*7,xright=3+(i-1)*7,ybottom=0,ytop=(1/2),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[i]][j,16]=="G"){
                        rect(xleft=3+(i-1)*7,xright=4+(i-1)*7,ybottom=0,ytop=(1/2),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[i]][j,16]=="T"){
                        rect(xleft=4+(i-1)*7,xright=5+(i-1)*7,ybottom=0,ytop=(1/2),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[i]][j,16]==""){
                        rect(xleft=5+(i-1)*7,xright=6+(i-1)*7,ybottom=0,ytop=(1/2),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[i]][j,17]=="A"){
                        rect(xleft=1+(i-1)*7,xright=2+(i-1)*7,ybottom=0,ytop=(1/2),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[i]][j,17]=="C"){
                        rect(xleft=2+(i-1)*7,xright=3+(i-1)*7,ybottom=0,ytop=(1/2),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[i]][j,17]=="G"){
                        rect(xleft=3+(i-1)*7,xright=4+(i-1)*7,ybottom=0,ytop=(1/2),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[i]][j,17]=="T"){
                        rect(xleft=4+(i-1)*7,xright=5+(i-1)*7,ybottom=0,ytop=(1/2),
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[i]][j,17]==""){
                        rect(xleft=5+(i-1)*7,xright=6+(i-1)*7,ybottom=0,ytop=(1/2),
                             border="grey50",lty=3,lwd=4)
                    }
                }
                if(is.na(calling[[i]][j,19])&&!is.na(calling[[i]][j,18])){
                    rect(xleft=6+(i-1)*7,xright=7+(i-1)*7,ybottom=0,ytop=1,
                         border="grey50",lty=3,lwd=4)
                    if(calling[[i]][j,16]=="A"){
                        rect(xleft=1+(i-1)*7,xright=2+(i-1)*7,ybottom=0,ytop=1,
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[i]][j,16]=="C"){
                        rect(xleft=2+(i-1)*7,xright=3+(i-1)*7,ybottom=0,ytop=1,
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[i]][j,16]=="G"){
                        rect(xleft=3+(i-1)*7,xright=4+(i-1)*7,ybottom=0,ytop=1,
                             border="grey50",lty=3,lwd=4)
                    }
                    if(calling[[i]][j,16]=="T"){
                        rect(xleft=4+(i-1)*7,xright=5+(i-1)*7,ybottom=0,ytop=1,
                             border="grey50",lty=3,lwd=4)
                    } 
                }
                if(!is.na(calling[[i]][j,19])&&!is.na(calling[[i]][j,18])){
                    rect(xleft=6+(i-1)*7,xright=7+(i-1)*7,ybottom=0,ytop=(1/2),
                         border="grey50",lty=3,lwd=4)
                    if(calling[[i]][j,16]==calling[[i]][j,17]){
                        if(calling[[i]][j,16]=="A"){
                            rect(xleft=1+(i-1)*7,xright=2+(i-1)*7,ybottom=0,ytop=1,
                                 border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[i]][j,16]=="C"){
                            rect(xleft=2+(i-1)*7,xright=3+(i-1)*7,ybottom=0,ytop=1,
                                 border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[i]][j,16]=="G"){
                            rect(xleft=3+(i-1)*7,xright=4+(i-1)*7,ybottom=0,ytop=1,
                                 border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[i]][j,16]=="T"){
                            rect(xleft=4+(i-1)*7,xright=5+(i-1)*7,ybottom=0,ytop=1,
                                 border="grey50",lty=3,lwd=4)
                        } 
                    }
                    if(calling[[i]][j,16]!=calling[[i]][j,17]){
                        if(calling[[i]][j,16]=="A"){
                            rect(xleft=1+(i-1)*7,xright=2+(i-1)*7,ybottom=0,
                                 ytop=(1/2),border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[i]][j,16]=="C"){
                            rect(xleft=2+(i-1)*7,xright=3+(i-1)*7,ybottom=0,
                                 ytop=(1/2),border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[i]][j,16]=="G"){
                            rect(xleft=3+(i-1)*7,xright=4+(i-1)*7,ybottom=0,
                                 ytop=(1/2),border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[i]][j,16]=="T"){
                            rect(xleft=4+(i-1)*7,xright=5+(i-1)*7,ybottom=0,
                                 ytop=(1/2),border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[i]][j,16]==""){
                            rect(xleft=5+(i-1)*7,xright=6+(i-1)*7,ybottom=0,
                                 ytop=(1/2),border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[i]][j,17]=="A"){
                            rect(xleft=1+(i-1)*7,xright=2+(i-1)*7,ybottom=0,
                                 ytop=(1/2),border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[i]][j,17]=="C"){
                            rect(xleft=2+(i-1)*7,xright=3+(i-1)*7,ybottom=0,
                                 ytop=(1/2),border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[i]][j,17]=="G"){
                            rect(xleft=3+(i-1)*7,xright=4+(i-1)*7,ybottom=0,
                                 ytop=(1/2),border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[i]][j,17]=="T"){
                            rect(xleft=4+(i-1)*7,xright=5+(i-1)*7,ybottom=0,
                                 ytop=(1/2),border="grey50",lty=3,lwd=4)
                        }
                        if(calling[[i]][j,17]==""){
                            rect(xleft=5+(i-1)*7,xright=6+(i-1)*7,ybottom=0,
                                 ytop=(1/2),border="grey50",lty=3,lwd=4)
                        }
                    }
                }
            }
        }    
        
        text(x=-0.5*length(plot_rel[[j]][,1]),
             y=c(rel_max*0.8*(-0.1)),"Reference",cex=2)
        steps<-seq(1,length(plot_rel[[j]][,1])*0.65,
                   by=(length(plot_rel[[j]][,1])*0.65-1)/length(A_col))-1.2*length(plot_rel[[j]][,1])
        if(length(plot_rel[[j]][,1])==1){
            steps<-seq(1*0.65,1,by=(1-1*0.65)/length(A_col))-1.75*1
        }
        if(length(plot_rel[[j]][,1])==2){
            steps<-seq(0.7,2*0.65,by=(2*0.65-0.7)/length(A_col))-1.2*2
        }
        rect(xleft=steps[1:length(steps)-1],xright=steps[2:length(steps)],
             ybottom=rel_max*0.88,ytop=rel_max*0.93,lwd=0,col=A_col,border=NA)
        rect(xleft=steps[1:length(steps)-1],xright=steps[2:length(steps)],
             ybottom=rel_max*0.81,ytop=rel_max*0.86,lwd=0,col=C_col,border=NA)
        rect(xleft=steps[1:length(steps)-1],xright=steps[2:length(steps)],
             ybottom=rel_max*0.74,ytop=rel_max*0.79,lwd=0,col=G_col,border=NA)
        rect(xleft=steps[1:length(steps)-1],xright=steps[2:length(steps)],
             ybottom=rel_max*0.67,ytop=rel_max*0.72,lwd=0,col=T_col,border=NA)
        text(x=length(plot_rel[[j]][,1])*0.75-1.2*length(plot_rel[[j]][,1]),
             y=c(rel_max*0.90),"A",cex=2)
        text(x=length(plot_rel[[j]][,1])*0.75-1.2*length(plot_rel[[j]][,1]),
             y=c(rel_max*0.83),"C",cex=2)
        text(x=length(plot_rel[[j]][,1])*0.75-1.2*length(plot_rel[[j]][,1]),
             y=c(rel_max*0.76),"G",cex=2)
        text(x=length(plot_rel[[j]][,1])*0.75-1.2*length(plot_rel[[j]][,1]),
             y=c(rel_max*0.69),"T",cex=2)
        text(x=length(plot_rel[[j]][,1])*0.78-1.2*length(plot_rel[[j]][,1]),
             y=c(rel_max*0.62),"Del",cex=2)
        text(x=length(plot_rel[[j]][,1])*0.78-1.2*length(plot_rel[[j]][,1]),
             y=c(rel_max*0.55),"Ins",cex=2)
        text(x=length(plot_rel[[j]][,1])*0.98-1.2*length(plot_rel[[j]][,1]),
             y=c(rel_max*0.48),"Expected",cex=2)
        text(x=length(plot_rel[[j]][,1])*1.15-1.2*length(plot_rel[[j]][,1]),
             y=c(rel_max*0.97),"Mean quality",cex=2)
        if(length(plot_rel[[j]][,1])>2){
            rect(xleft=1-1.2*length(plot_rel[[j]][,1]),
                 xright=length(plot_rel[[j]][,1])*0.65-1.2*length(plot_rel[[j]][,1]),
                 ybottom=rel_max*0.60,ytop=rel_max*0.65,lwd=0,col="black")
            rect(xleft=1-1.2*length(plot_rel[[j]][,1]),
                 xright=length(plot_rel[[j]][,1])*0.65-1.2*length(plot_rel[[j]][,1]),
                 ybottom=rel_max*0.53,ytop=rel_max*0.58,lwd=4,border="purple",col="white")
            rect(xleft=1-1.2*length(plot_rel[[j]][,1]),
                 xright=length(plot_rel[[j]][,1])*0.65-1.2*length(plot_rel[[j]][,1]),
                 ybottom=rel_max*0.46,ytop=rel_max*0.51,lwd=2,border="grey50",col="white",lty=3)
            text(x=1-1.2*length(plot_rel[[j]][,1]),
                 y=c(rel_max*0.97),qual_lower_bound,cex=2)
            text(x=length(plot_rel[[j]][,1])*0.65-1.2*length(plot_rel[[j]][,1]),
                 y=c(rel_max*0.97),qual_upper_bound,cex=2)
        }
        if(length(plot_rel[[j]][,1])==1){
            rect(xleft=-1.1,xright=-0.75,ybottom=rel_max*0.60,ytop=rel_max*0.65,
                 lwd=0,col="black")
            rect(xleft=-1.1,xright=-0.75,ybottom=rel_max*0.53,ytop=rel_max*0.58,
                 lwd=4,border="purple",col="white")
            rect(xleft=-1.1,xright=-0.75,ybottom=rel_max*0.46,ytop=rel_max*0.51,
                 lwd=2,border="grey50",col="white",lty=3)
            text(x=-1.1,y=c(rel_max*0.97),qual_lower_bound,cex=2)
            text(x=-0.75,y=c(rel_max*0.97),qual_upper_bound,cex=2)
        }
        if(length(plot_rel[[j]][,1])==2){
            rect(xleft=-1.7,xright=-1.1,ybottom=rel_max*0.60,ytop=rel_max*0.65,
                 lwd=0,col="black")
            rect(xleft=-1.7,xright=-1.1,ybottom=rel_max*0.53,ytop=rel_max*0.58,
                 lwd=4,border="purple",col="white")
            rect(xleft=-1.7,xright=-1.1,ybottom=rel_max*0.46,ytop=rel_max*0.51,
                 lwd=2,border="grey50",col="white",lty=3)
            text(x=-1.7,y=c(rel_max*0.97),qual_lower_bound,cex=2)
            text(x=-1.1,y=c(rel_max*0.97),qual_upper_bound,cex=2)
        }
        if(directory2!=""){
            dev.off()
        }
        return(counter)
    }
    
    
    
    #plot absolute (number of reads), plot frequency threshold, 
    #quality by color (intensity) (define upper and lower bound for color-coding)
    plotCountsPerSample<-function(samples,
                                  directory2,
                                  calling,
                                  dbsnp_file,
                                  frequency_weight,
                                  qual_lower_bound,
                                  qual_upper_bound,
                                  marks){
        plot_abs<-list()
        plot_ref<-list()
        plot_col<-list()
        plot_insert<-list()
        plot_col_insert<-list()
        A_Pal<-colorRampPalette(c('chartreuse','forestgreen'))
        A_col<-A_Pal(10*(qual_upper_bound-qual_lower_bound)+1)
        C_Pal<-colorRampPalette(c('cyan','blue4'))
        C_col<-C_Pal(10*(qual_upper_bound-qual_lower_bound)+1)
        G_Pal<-colorRampPalette(c('yellow','darkgoldenrod3'))
        G_col<-G_Pal(10*(qual_upper_bound-qual_lower_bound)+1)
        T_Pal<-colorRampPalette(c('firebrick1','darkred'))
        T_col<-T_Pal(10*(qual_upper_bound-qual_lower_bound)+1)
        if((qual_upper_bound-qual_lower_bound)==0){
            A_col<-"forestgreen"
            C_col<-"blue4"
            G_col<-"darkgoldenrod3"
            T_col<-"darkred"
        }   
        for(j in 1:length(samples[,1])){
            plot_abs[[j]]<-matrix(rep(0,6*length(calling[[j]][,1])),ncol=6)
            plot_ref[[j]]<-matrix(rep(0,6*length(calling[[j]][,1])),ncol=6)
            plot_col[[j]]<-matrix(rep(0,6*length(calling[[j]][,1])),ncol=6)
            plot_insert[[j]]<-matrix(rep(0,4*length(calling[[j]][,1])),ncol=4)
            plot_col_insert[[j]]<-matrix(rep(NA,4*length(calling[[j]][,1])),ncol=4)
        }
        abs_max<-0
        for(j in 1:length(samples[,1])){
            long_insert<-0
            for(i in 1:length(calling[[j]][,1])){
                if(!is.na(frequency_weight[[j]][i-long_insert,1])
                   &&paste(calling[[j]][i,1],calling[[j]][i,2])
                   !=paste(frequency_weight[[j]][i-long_insert,1],frequency_weight[[j]][i-long_insert,2])){
                    long_insert<-long_insert+1
                }
                if(is.na(frequency_weight[[j]][i-long_insert,1])){
                    long_insert<-long_insert+1
                }
                for(k in c(1,3,5,7,9,10)){
                    plot_abs[[j]]<-plotCountsPerSample_abs(calling,
                                                           frequency_weight,
                                                           long_insert,
                                                           i,j,k,plot_abs[[j]])
                    plot_col[[j]]<-plotCountsPerSample_col(calling,
                                                           frequency_weight,
                                                           long_insert,
                                                           qual_lower_bound,
                                                           qual_upper_bound,
                                                           A_col,C_col,G_col,T_col,
                                                           i,j,k,plot_col[[j]])
                    
                }
                if(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)>abs_max){
                    abs_max<-sum(plot_abs[[j]][i,1:5],na.rm=TRUE)
                }
                
                for(k in c(1,3,5,7)){
                    plot_insert[[j]]<-plotCountsPerSample_insert(calling,
                                                                 frequency_weight,
                                                                 long_insert,
                                                                 i,j,k,plot_insert[[j]])
                    plot_col_insert[[j]]<-plotCountsPerSample_col_insert(calling,
                                                                         frequency_weight,
                                                                         long_insert,
                                                                         qual_lower_bound,
                                                                         qual_upper_bound,
                                                                         A_col,C_col,
                                                                         G_col,T_col,
                                                                         i,j,k,plot_col_insert[[j]])
                }
                
                plot_insert[[j]][i,2]<-sum(plot_insert[[j]][i,1],plot_insert[[j]][i,2])
                plot_insert[[j]][i,3]<-sum(plot_insert[[j]][i,2],plot_insert[[j]][i,3])
                plot_insert[[j]][i,4]<-sum(plot_insert[[j]][i,3],plot_insert[[j]][i,4])
                
                chr<-substr(calling[[j]][i,1],4,nchar(calling[[j]][i,1]))
                pos<-calling[[j]][i,2]
                if(dbsnp_file!=""){
                    vcfFile<-system.file(dbsnp_file, package="VariantAnnotation")
                    param = ScanVcfParam(fixed=character(),info=character(), geno=character(), 
                                         samples=character(),which=GRanges(seqnames=chr,ranges=IRanges(pos,pos)))
                    vcf.comp<-readVcf(dbsnp_file,"hg19",param=param)
                }
                plot_ref[[j]]<-plotCountsPerSample_ref(calling,
                                                       dbsnp_file,
                                                       vcf.comp,
                                                       j,i,plot_ref[[j]])  
            }
        }
        if(abs_max==0){
            abs_max<-1
        }
        for(j in 1:length(plot_abs)){
            plotCountsPerSample_plot(directory2,samples,plot_abs,plot_col,abs_max,
                                     calling,plot_insert,plot_col_insert,A_col,
                                     C_col,G_col,T_col,plot_ref,marks,
                                     qual_lower_bound,qual_upper_bound,j)
        }
        return()
    }
    
    
    #plot relative (number of reads), plot frequency threshold, 
    #quality by color (intensity) (define upper and lower bound for color-coding)
    plotCountsPerSampleRelative<-function(samples,
                                          directory2,
                                          calling,
                                          dbsnp_file,
                                          frequency_weight,
                                          qual_lower_bound,
                                          qual_upper_bound,
                                          marks){
        plot_rel<-list()
        plot_ref<-list()
        plot_col<-list()
        plot_insert<-list()
        plot_col_insert<-list()
        A_Pal<-colorRampPalette(c('chartreuse','forestgreen'))
        A_col<-A_Pal(10*(qual_upper_bound-qual_lower_bound)+1)
        C_Pal<-colorRampPalette(c('cyan','blue4'))
        C_col<-C_Pal(10*(qual_upper_bound-qual_lower_bound)+1)
        G_Pal<-colorRampPalette(c('yellow','darkgoldenrod3'))
        G_col<-G_Pal(10*(qual_upper_bound-qual_lower_bound)+1)
        T_Pal<-colorRampPalette(c('firebrick1','darkred'))
        T_col<-T_Pal(10*(qual_upper_bound-qual_lower_bound)+1)
        if((qual_upper_bound-qual_lower_bound)==0){
            A_col<-"forestgreen"
            C_col<-"blue4"
            G_col<-"darkgoldenrod3"
            T_col<-"darkred"
        }
        
        for(j in 1:length(samples[,1])){
            plot_rel[[j]]<-matrix(rep(0,6*length(calling[[j]][,1])),ncol=6)
            plot_ref[[j]]<-matrix(rep(0,6*length(calling[[j]][,1])),ncol=6)
            plot_col[[j]]<-matrix(rep(0,6*length(calling[[j]][,1])),ncol=6)
            plot_insert[[j]]<-matrix(rep(0,4*length(calling[[j]][,1])),ncol=4)
            plot_col_insert[[j]]<-matrix(rep(NA,4*length(calling[[j]][,1])),ncol=4)
        }
        rel_max<-1
        for(j in 1:length(samples[,1])){
            long_insert<-0
            for(i in 1:length(calling[[j]][,1])){
                if(!is.na(frequency_weight[[j]][i-long_insert,1])
                   &&paste(calling[[j]][i,1],calling[[j]][i,2])
                   !=paste(frequency_weight[[j]][i-long_insert,1],frequency_weight[[j]][i-long_insert,2])){
                    long_insert<-long_insert+1
                }
                if(is.na(frequency_weight[[j]][i-long_insert,1])){
                    long_insert<-long_insert+1
                }
                for(k in 4:9){
                    plot_rel[[j]][i,k-3]<-calling[[j]][i,k]
                    if(is.na(plot_rel[[j]][i,k-3])){
                        plot_rel[[j]][i,k-3]<-0
                    }
                    plot_col[[j]]<-plotCountsPerSampleRelative_col(calling,
                                                                   frequency_weight,
                                                                   long_insert,
                                                                   qual_lower_bound,
                                                                   qual_upper_bound,
                                                                   A_col,C_col,
                                                                   G_col,T_col,
                                                                   i,j,k,plot_col[[j]])   
                }
                
                for(k in c(1,3,5,7)){
                    plot_insert[[j]]<-plotCountsPerSampleRelative_insert(calling,
                                                                         frequency_weight,
                                                                         long_insert,
                                                                         i,j,k,
                                                                         plot_insert[[j]])
                    plot_col_insert[[j]]<-plotCountsPerSampleRelative_col_insert(calling,
                                                                                 frequency_weight,
                                                                                 long_insert,
                                                                                 qual_lower_bound,
                                                                                 qual_upper_bound,
                                                                                 A_col,C_col,
                                                                                 G_col,T_col,
                                                                                 i,j,k,
                                                                                 plot_col_insert[[j]])
                }
                plot_insert[[j]][i,2]<-sum(plot_insert[[j]][i,1],plot_insert[[j]][i,2])
                plot_insert[[j]][i,3]<-sum(plot_insert[[j]][i,2],plot_insert[[j]][i,3])
                plot_insert[[j]][i,4]<-sum(plot_insert[[j]][i,3],plot_insert[[j]][i,4])
                
                chr<-substr(calling[[j]][i,1],4,nchar(calling[[j]][i,1]))
                pos<-calling[[j]][i,2]
                if(dbsnp_file!=""){
                    vcfFile<-system.file(dbsnp_file, package="VariantAnnotation")
                    param = ScanVcfParam(fixed=character(),info=character(), geno=character(), 
                                         samples=character(),which=GRanges(seqnames=chr,ranges=IRanges(pos,pos)))
                    vcf.comp<-readVcf(dbsnp_file,"hg19",param=param)
                }
                plot_ref[[j]]<-plotCountsPerSample_ref(calling,
                                                       dbsnp_file,
                                                       vcf.comp,
                                                       j,i,plot_ref[[j]])  
            }
        }
        
        for(j in 1:length(plot_rel)){
            plotCountsPerSampleRelative_plot(directory2,samples,plot_rel,plot_col,
                                             rel_max,calling,plot_insert,
                                             plot_col_insert,A_col,C_col,G_col,
                                             T_col,plot_ref,marks,qual_lower_bound,
                                             qual_upper_bound,j)
        }
        #return()
    }
    
    #Plot counts per Target
    plotCountsPerTarget<-function(samples,
                                  directory2,
                                  calling,
                                  dbsnp_file,
                                  frequency_weight,
                                  qual_lower_bound,
                                  qual_upper_bound,
                                  marks){
        plot_abs<-list()
        plot_ref<-list()
        plot_col<-list()
        plot_insert<-list()
        plot_col_insert<-list()
        A_Pal<-colorRampPalette(c('chartreuse','forestgreen'))
        A_col<-A_Pal(10*(qual_upper_bound-qual_lower_bound)+1)
        C_Pal<-colorRampPalette(c('cyan','blue4'))
        C_col<-C_Pal(10*(qual_upper_bound-qual_lower_bound)+1)
        G_Pal<-colorRampPalette(c('yellow','darkgoldenrod3'))
        G_col<-G_Pal(10*(qual_upper_bound-qual_lower_bound)+1)
        T_Pal<-colorRampPalette(c('firebrick1','darkred'))
        T_col<-T_Pal(10*(qual_upper_bound-qual_lower_bound)+1)  
        if((qual_upper_bound-qual_lower_bound)==0){
            A_col<-"forestgreen"
            C_col<-"blue4"
            G_col<-"darkgoldenrod3"
            T_col<-"darkred"
        }
        
        for(j in 1:length(calling[[1]][,1])){
            plot_abs[[j]]<-matrix(rep(0,6*length(samples[,1])),ncol=6)
            plot_ref[[j]]<-matrix(rep(0,6*length(samples[,1])),ncol=6)
            plot_col[[j]]<-matrix(rep(0,6*length(samples[,1])),ncol=6)
            plot_insert[[j]]<-matrix(rep(0,4*length(samples[,1])),ncol=4)
            plot_col_insert[[j]]<-matrix(rep(NA,4*length(samples[,1])),ncol=4)
        }
        abs_max<-0
        long_insert<-rep(0,length(samples[,1]))
        
        for(j in 1:length(calling[[1]][,1])){
            for(i in 1:length(samples[,1])){
                if(!is.na(frequency_weight[[i]][j-long_insert[i],1])
                   &&paste(calling[[i]][j,1],calling[[i]][j,2])
                   !=paste(frequency_weight[[i]][j-long_insert[i],1],frequency_weight[[i]][j-long_insert[i],2])){
                    long_insert[i]<-long_insert[i]+1
                }
                if(is.na(frequency_weight[[i]][j-long_insert[i],1])){
                    long_insert[i]<-long_insert[i]+1
                }
                for(k in c(1,3,5,7,9,10)){
                    plot_abs[[j]]<-plotCountsPerTarget_abs(samples,
                                                           frequency_weight,
                                                           long_insert,
                                                           i,j,k,plot_abs[[j]])
                    plot_col[[j]]<-plotCountsPerTarget_col(samples,
                                                           frequency_weight,
                                                           long_insert,
                                                           qual_lower_bound,
                                                           qual_upper_bound,
                                                           A_col,C_col,G_col,T_col,
                                                           i,j,k,plot_col[[j]])
                }
                if(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)>abs_max){
                    abs_max<-sum(plot_abs[[j]][i,1:5],na.rm=TRUE)
                }
                for(k in c(1,3,5,7)){
                    plot_insert[[j]]<-plotCountsPerTarget_insert(samples,
                                                                 frequency_weight,
                                                                 long_insert,
                                                                 i,j,k,plot_insert[[j]])
                    plot_col_insert[[j]]<-plotCountsPerTarget_col_insert(samples,
                                                                         frequency_weight,
                                                                         long_insert,
                                                                         qual_lower_bound,
                                                                         qual_upper_bound,
                                                                         A_col,C_col,
                                                                         G_col,T_col,
                                                                         i,j,k,
                                                                         plot_col_insert[[j]])
                }
                plot_insert[[j]][i,2]<-sum(plot_insert[[j]][i,1],plot_insert[[j]][i,2])
                plot_insert[[j]][i,3]<-sum(plot_insert[[j]][i,2],plot_insert[[j]][i,3])
                plot_insert[[j]][i,4]<-sum(plot_insert[[j]][i,3],plot_insert[[j]][i,4])
                
                chr<-substr(calling[[i]][j,1],4,nchar(calling[[i]][j,1]))
                pos<-calling[[i]][j,2]
                if(dbsnp_file!=""){
                    vcfFile<-system.file(dbsnp_file, package="VariantAnnotation")
                    param = ScanVcfParam(fixed=character(),info=character(), geno=character(), 
                                         samples=character(),which=GRanges(seqnames=chr,ranges=IRanges(pos,pos)))
                    vcf.comp<-readVcf(dbsnp_file,"hg19",param=param)
                }
                plot_ref[[j]]<-plotCountsPerTarget_ref(samples,
                                                       calling,
                                                       dbsnp_file,
                                                       vcf.comp,
                                                       j,i,plot_ref[[j]])
            }
        }
        if(abs_max==0){
            abs_max<-1
        }
        counter2<-0
        for(j in 1:length(calling[[1]][,1])){
            counter<-counter2
            counter2<-plotCountsPerTarget_plot(directory2,samples,plot_abs,plot_col,
                                               abs_max,calling,plot_insert,
                                               plot_col_insert,A_col,C_col,G_col,
                                               T_col,plot_ref,marks,
                                               qual_lower_bound,qual_upper_bound,
                                               counter,j)
        }
        return()
    }
    
    #Plot counts per Target
    plotCountsPerTargetRelative<-function(samples,
                                          directory2,
                                          calling,
                                          dbsnp_file,
                                          frequency_weight,
                                          qual_lower_bound,
                                          qual_upper_bound,
                                          marks){
        plot_rel<-list()
        plot_ref<-list()
        plot_col<-list()
        plot_insert<-list()
        plot_col_insert<-list()
        A_Pal<-colorRampPalette(c('chartreuse','forestgreen'))
        A_col<-A_Pal(10*(qual_upper_bound-qual_lower_bound)+1)
        C_Pal<-colorRampPalette(c('cyan','blue4'))
        C_col<-C_Pal(10*(qual_upper_bound-qual_lower_bound)+1)
        G_Pal<-colorRampPalette(c('yellow','darkgoldenrod3'))
        G_col<-G_Pal(10*(qual_upper_bound-qual_lower_bound)+1)
        T_Pal<-colorRampPalette(c('firebrick1','darkred'))
        T_col<-T_Pal(10*(qual_upper_bound-qual_lower_bound)+1)
        if((qual_upper_bound-qual_lower_bound)==0){
            A_col<-"forestgreen"
            C_col<-"blue4"
            G_col<-"darkgoldenrod3"
            T_col<-"darkred"
        }
        
        for(j in 1:length(calling[[1]][,1])){
            plot_rel[[j]]<-matrix(rep(0,6*length(samples[,1])),ncol=6)
            plot_ref[[j]]<-matrix(rep(0,6*length(samples[,1])),ncol=6)
            plot_col[[j]]<-matrix(rep(0,6*length(samples[,1])),ncol=6)
            plot_insert[[j]]<-matrix(rep(0,4*length(samples[,1])),ncol=4)
            plot_col_insert[[j]]<-matrix(rep(NA,4*length(samples[,1])),ncol=4)
        }
        rel_max<-1
        long_insert<-rep(0,length(samples[,1]))
        
        for(j in 1:length(calling[[1]][,1])){
            for(i in 1:length(samples[,1])){
                if(!is.na(frequency_weight[[i]][j-long_insert[i],1])
                   &&paste(calling[[i]][j,1],calling[[i]][j,2])
                   !=paste(frequency_weight[[i]][j-long_insert[i],1],
                           frequency_weight[[i]][j-long_insert[i],2])){
                    long_insert[i]<-long_insert[i]+1
                }
                if(is.na(frequency_weight[[i]][j-long_insert[i],1])){
                    long_insert[i]<-long_insert[i]+1
                }
                for(k in 4:9){
                    plot_rel[[j]][i,k-3]<-calling[[i]][j,k]
                    if(is.na(plot_rel[[j]][i,k-3])){
                        plot_rel[[j]][i,k-3]<-0
                    }
                    plot_col[[j]]<-plotCountsPerTargetRelative_col(samples,
                                                                   frequency_weight,
                                                                   long_insert,
                                                                   qual_lower_bound,
                                                                   qual_upper_bound,
                                                                   A_col,C_col,
                                                                   G_col,T_col,
                                                                   i,j,k,plot_col[[j]])    
                }            
                for(k in c(1,3,5,7)){
                    plot_insert[[j]]<-plotCountsPerTargetRelative_insert(samples,
                                                                         frequency_weight,
                                                                         long_insert,
                                                                         i,j,k,
                                                                         plot_insert[[j]])
                    plot_col_insert[[j]]<-plotCountsPerTargetRelative_col_insert(samples,
                                                                                 frequency_weight,
                                                                                 long_insert,
                                                                                 qual_lower_bound,
                                                                                 qual_upper_bound,
                                                                                 A_col,C_col,
                                                                                 G_col,T_col,
                                                                                 i,j,k,
                                                                                 plot_col_insert[[j]]) 
                }
                plot_insert[[j]][i,2]<-sum(plot_insert[[j]][i,1],plot_insert[[j]][i,2])
                plot_insert[[j]][i,3]<-sum(plot_insert[[j]][i,2],plot_insert[[j]][i,3])
                plot_insert[[j]][i,4]<-sum(plot_insert[[j]][i,3],plot_insert[[j]][i,4])
                
                chr<-substr(calling[[i]][j,1],4,nchar(calling[[i]][j,1]))
                pos<-calling[[i]][j,2]
                if(dbsnp_file!=""){
                    vcfFile<-system.file(dbsnp_file, package="VariantAnnotation")
                    param = ScanVcfParam(fixed=character(),info=character(), geno=character(), 
                                         samples=character(),which=GRanges(seqnames=chr,ranges=IRanges(pos,pos)))
                    vcf.comp<-readVcf(dbsnp_file,"hg19",param=param)
                }
                plot_ref[[j]]<-plotCountsPerTarget_ref(samples,
                                                       calling,
                                                       dbsnp_file,
                                                       vcf.comp,
                                                       j,i,plot_ref[[j]])  
            }
        }
        counter2<-0
        for(j in 1:length(calling[[1]][,1])){
            counter<-counter2
            counter2<-plotCountsPerTargetRelative_plot(directory2,samples,plot_rel,plot_col,
                                                       rel_max,calling,plot_insert,
                                                       plot_col_insert,A_col,C_col,G_col,
                                                       T_col,plot_ref,marks,qual_lower_bound,
                                                       qual_upper_bound,counter,j)  
        }
        return()
    }
    
    #plot absolute (number of reads), plot frequency threshold, quality by color 
    #(intensity) (define upper and lower bound for color-coding)
    plotCountsPerSampleVcf<-function(samples,
                                     directory2,
                                     calling,
                                     dbsnp_file,
                                     frequency_weight,
                                     qual_lower_bound,
                                     qual_upper_bound,
                                     marks){
        plot_abs<-list()
        plot_ref<-list()
        plot_col<-list()
        plot_insert<-list()
        plot_col_insert<-list()
        A_Pal<-colorRampPalette(c('chartreuse','forestgreen'))
        A_col<-A_Pal(10*(qual_upper_bound-qual_lower_bound)+1)
        C_Pal<-colorRampPalette(c('cyan','blue4'))
        C_col<-C_Pal(10*(qual_upper_bound-qual_lower_bound)+1)
        G_Pal<-colorRampPalette(c('yellow','darkgoldenrod3'))
        G_col<-G_Pal(10*(qual_upper_bound-qual_lower_bound)+1)
        T_Pal<-colorRampPalette(c('firebrick1','darkred'))
        T_col<-T_Pal(10*(qual_upper_bound-qual_lower_bound)+1)
        if((qual_upper_bound-qual_lower_bound)==0){
            A_col<-"forestgreen"
            C_col<-"blue4"
            G_col<-"darkgoldenrod3"
            T_col<-"darkred"
        }
        
        for(j in 1:length(samples[,1])){
            plot_abs[[j]]<-matrix(rep(0,6*length(calling[[j]][,1])),ncol=6)
            plot_ref[[j]]<-matrix(rep(0,6*length(calling[[j]][,1])),ncol=6)
            plot_col[[j]]<-matrix(rep(0,6*length(calling[[j]][,1])),ncol=6)
            plot_insert[[j]]<-matrix(rep(0,4*length(calling[[j]][,1])),ncol=4)
            plot_col_insert[[j]]<-matrix(rep(NA,4*length(calling[[j]][,1])),ncol=4)
        }
        abs_max<-0
        for(j in 1:length(samples[,1])){
            long_insert<-0
            for(i in 1:length(calling[[j]][,1])){
                if(!is.na(frequency_weight[[j]][i-long_insert,1])
                   &&paste(calling[[j]][i,1],calling[[j]][i,2])
                   !=paste(frequency_weight[[j]][i-long_insert,1],frequency_weight[[j]][i-long_insert,2])){
                    long_insert<-long_insert+1
                }
                if(is.na(frequency_weight[[j]][i-long_insert,1])){
                    long_insert<-long_insert+1
                }
                for(k in c(1,3,5,7,9,10)){
                    plot_abs[[j]]<-plotCountsPerSample_abs(calling,
                                                           frequency_weight,
                                                           long_insert,
                                                           i,j,k,plot_abs[[j]])
                    plot_col[[j]]<-plotCountsPerSample_col(calling,
                                                           frequency_weight,
                                                           long_insert,
                                                           qual_lower_bound,
                                                           qual_upper_bound,
                                                           A_col,C_col,G_col,T_col,
                                                           i,j,k,plot_col[[j]])   
                }
                if(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)>abs_max){
                    abs_max<-sum(plot_abs[[j]][i,1:5],na.rm=TRUE)
                }
                
                for(k in c(1,3,5,7)){
                    plot_insert[[j]]<-plotCountsPerSample_insert(calling,
                                                                 frequency_weight,
                                                                 long_insert,
                                                                 i,j,k,plot_insert[[j]])
                    plot_col_insert[[j]]<-plotCountsPerSample_col_insert(calling,
                                                                         frequency_weight,
                                                                         long_insert,
                                                                         qual_lower_bound,
                                                                         qual_upper_bound,
                                                                         A_col,C_col,
                                                                         G_col,T_col,
                                                                         i,j,k,
                                                                         plot_col_insert[[j]]) 
                }
                plot_insert[[j]][i,2]<-sum(plot_insert[[j]][i,1],plot_insert[[j]][i,2])
                plot_insert[[j]][i,3]<-sum(plot_insert[[j]][i,2],plot_insert[[j]][i,3])
                plot_insert[[j]][i,4]<-sum(plot_insert[[j]][i,3],plot_insert[[j]][i,4])
                
                chr<-substr(calling[[j]][i,1],4,nchar(calling[[j]][i,1]))
                pos<-calling[[j]][i,2]
                if(dbsnp_file!=""){
                    vcfFile<-system.file(dbsnp_file, package="VariantAnnotation")
                    param = ScanVcfParam(fixed=character(),info=character(), geno=character(), 
                                         samples=character(),which=GRanges(seqnames=chr,ranges=IRanges(pos,pos)))
                    vcf.comp<-readVcf(dbsnp_file,"hg19",param=param)
                }
                plot_ref<-plotCountsPerSample_ref(calling,
                                                  dbsnp_file,
                                                  vcf.comp,
                                                  j,i,plot_ref[[j]])
                
            }
        }
        if(abs_max==0){
            abs_max<-1
        }
        for(j in 1:length(plot_abs)){
            plotCountsPerSampleVcf_plot(directory2,samples,plot_abs,plot_col,
                                        abs_max,calling,plot_insert,plot_col_insert,
                                        A_col,C_col,G_col,T_col,plot_ref,marks,
                                        qual_lower_bound,qual_upper_bound,j)
        }
        return()
    }
    
    #plot relative (number of reads), plot frequency threshold, quality by color 
    #(intensity) (define upper and lower bound for color-coding)
    plotCountsPerSampleRelativeVcf<-function(samples,
                                             directory2,
                                             calling,
                                             dbsnp_file,
                                             frequency_weight,
                                             qual_lower_bound,
                                             qual_upper_bound,
                                             marks){
        plot_rel<-list()
        plot_ref<-list()
        plot_col<-list()
        plot_insert<-list()
        plot_col_insert<-list()
        A_Pal<-colorRampPalette(c('chartreuse','forestgreen'))
        A_col<-A_Pal(10*(qual_upper_bound-qual_lower_bound)+1)
        C_Pal<-colorRampPalette(c('cyan','blue4'))
        C_col<-C_Pal(10*(qual_upper_bound-qual_lower_bound)+1)
        G_Pal<-colorRampPalette(c('yellow','darkgoldenrod3'))
        G_col<-G_Pal(10*(qual_upper_bound-qual_lower_bound)+1)
        T_Pal<-colorRampPalette(c('firebrick1','darkred'))
        T_col<-T_Pal(10*(qual_upper_bound-qual_lower_bound)+1)
        if((qual_upper_bound-qual_lower_bound)==0){
            A_col<-"forestgreen"
            C_col<-"blue4"
            G_col<-"darkgoldenrod3"
            T_col<-"darkred"
        }
        
        for(j in 1:length(samples[,1])){
            plot_rel[[j]]<-matrix(rep(0,6*length(calling[[j]][,1])),ncol=6)
            plot_ref[[j]]<-matrix(rep(0,6*length(calling[[j]][,1])),ncol=6)
            plot_col[[j]]<-matrix(rep(0,6*length(calling[[j]][,1])),ncol=6)
            plot_insert[[j]]<-matrix(rep(0,4*length(calling[[j]][,1])),ncol=4)
            plot_col_insert[[j]]<-matrix(rep(NA,4*length(calling[[j]][,1])),ncol=4)
        }
        rel_max<-1
        for(j in 1:length(samples[,1])){
            long_insert<-0
            for(i in 1:length(calling[[j]][,1])){
                if(!is.na(frequency_weight[[j]][i-long_insert,1])
                   &&paste(calling[[j]][i,1],calling[[j]][i,2])
                   !=paste(frequency_weight[[j]][i-long_insert,1],frequency_weight[[j]][i-long_insert,2])){
                    long_insert<-long_insert+1
                }
                if(is.na(frequency_weight[[j]][i-long_insert,1])){
                    long_insert<-long_insert+1
                }
                for(k in 4:9){
                    plot_rel[[j]][i,k-3]<-calling[[j]][i,k]
                    if(is.na(plot_rel[[j]][i,k-3])){
                        plot_rel[[j]][i,k-3]<-0
                    }
                    
                    plot_col[[j]]<-plotCountsPerSampleRelative_col(calling,
                                                                   frequency_weight,
                                                                   long_insert,
                                                                   qual_lower_bound,
                                                                   qual_upper_bound,
                                                                   A_col,C_col,
                                                                   G_col,T_col,
                                                                   i,j,k,plot_col[[j]])
                }
                for(k in c(1,3,5,7)){
                    plot_insert[[j]]<-plotCountsPerSampleRelative_insert(calling,
                                                                         frequency_weight,
                                                                         long_insert,
                                                                         i,j,k,
                                                                         plot_insert[[j]])
                    plot_col_insert[[j]]<-plotCountsPerSampleRelative_col_insert(calling,
                                                                                 frequency_weight,
                                                                                 long_insert,
                                                                                 qual_lower_bound,
                                                                                 qual_upper_bound,
                                                                                 A_col,C_col,
                                                                                 G_col,T_col,
                                                                                 i,j,k,
                                                                                 plot_col_insert[[j]])
                }
                plot_insert[[j]][i,2]<-sum(plot_insert[[j]][i,1],plot_insert[[j]][i,2])
                plot_insert[[j]][i,3]<-sum(plot_insert[[j]][i,2],plot_insert[[j]][i,3])
                plot_insert[[j]][i,4]<-sum(plot_insert[[j]][i,3],plot_insert[[j]][i,4])
                
                chr<-substr(calling[[j]][i,1],4,nchar(calling[[j]][i,1]))
                pos<-calling[[j]][i,2]
                if(dbsnp_file!=""){
                    vcfFile<-system.file(dbsnp_file, package="VariantAnnotation")
                    param = ScanVcfParam(fixed=character(),info=character(), geno=character(), 
                                         samples=character(),which=GRanges(seqnames=chr,ranges=IRanges(pos,pos)))
                    vcf.comp<-readVcf(dbsnp_file,"hg19",param=param)
                }
                plot_ref[[j]]<-plotCountsPerSample_ref(calling,dbsnp_file,vcf.comp,j,i,
                                                       plot_ref[[j]])
            }
        }
        for(j in 1:length(plot_rel)){
            plotCountsPerSampleRelativeVcf_plot(directory2,samples,plot_rel,
                                                plot_col,rel_max,calling,
                                                plot_insert,plot_col_insert,A_col,
                                                C_col,G_col,T_col,plot_ref,marks,
                                                qual_lower_bound,qual_upper_bound,j)
        }
        return()
    }
    
    #Plot counts per Target
    plotCountsPerTargetVcf<-function(samples,
                                     directory2,
                                     calling,
                                     dbsnp_file,
                                     frequency_weight,
                                     qual_lower_bound,
                                     qual_upper_bound,
                                     marks){
        plot_abs<-list()
        plot_ref<-list()
        plot_col<-list()
        plot_insert<-list()
        plot_col_insert<-list()
        A_Pal<-colorRampPalette(c('chartreuse','forestgreen'))
        A_col<-A_Pal(10*(qual_upper_bound-qual_lower_bound)+1)
        C_Pal<-colorRampPalette(c('cyan','blue4'))
        C_col<-C_Pal(10*(qual_upper_bound-qual_lower_bound)+1)
        G_Pal<-colorRampPalette(c('yellow','darkgoldenrod3'))
        G_col<-G_Pal(10*(qual_upper_bound-qual_lower_bound)+1)
        T_Pal<-colorRampPalette(c('firebrick1','darkred'))
        T_col<-T_Pal(10*(qual_upper_bound-qual_lower_bound)+1)
        if((qual_upper_bound-qual_lower_bound)==0){
            A_col<-"forestgreen"
            C_col<-"blue4"
            G_col<-"darkgoldenrod3"
            T_col<-"darkred"
        }
        
        for(j in 1:length(calling[[1]][,1])){
            plot_abs[[j]]<-matrix(rep(0,6*length(samples[,1])),ncol=6)
            plot_ref[[j]]<-matrix(rep(0,6*length(samples[,1])),ncol=6)
            plot_col[[j]]<-matrix(rep(0,6*length(samples[,1])),ncol=6)
            plot_insert[[j]]<-matrix(rep(0,4*length(samples[,1])),ncol=4)
            plot_col_insert[[j]]<-matrix(rep(NA,4*length(samples[,1])),ncol=4)
        }
        abs_max<-0
        long_insert<-rep(0,length(samples[,1]))
        
        for(j in 1:length(calling[[1]][,1])){
            for(i in 1:length(samples[,1])){
                if(!is.na(frequency_weight[[i]][j-long_insert[i],1])
                   &&paste(calling[[i]][j,1],calling[[i]][j,2])
                   !=paste(frequency_weight[[i]][j-long_insert[i],1],frequency_weight[[i]][j-long_insert[i],2])){
                    long_insert[i]<-long_insert[i]+1
                }
                if(is.na(frequency_weight[[i]][j-long_insert[i],1])){
                    long_insert[i]<-long_insert[i]+1
                }
                for(k in c(1,3,5,7,9,10)){
                    plot_abs[[j]]<-plotCountsPerTarget_abs(samples,
                                                           frequency_weight,
                                                           long_insert,
                                                           i,j,k,plot_abs[[j]])
                    plot_col[[j]]<-plotCountsPerTarget_col(samples,
                                                           frequency_weight,
                                                           long_insert,
                                                           qual_lower_bound,
                                                           qual_upper_bound,
                                                           A_col,C_col,G_col,T_col,
                                                           i,j,k,plot_col[[j]])  
                }
                if(sum(plot_abs[[j]][i,1:5],na.rm=TRUE)>abs_max){
                    abs_max<-sum(plot_abs[[j]][i,1:5],na.rm=TRUE)
                }
                for(k in c(1,3,5,7)){
                    plot_insert[[j]]<-plotCountsPerTarget_insert(samples,
                                                                 frequency_weight,
                                                                 long_insert,
                                                                 i,j,k,plot_insert[[j]])
                    plot_col_insert[[j]]<-plotCountsPerTarget_col_insert(samples,
                                                                         frequency_weight,
                                                                         long_insert,
                                                                         qual_lower_bound,
                                                                         qual_upper_bound,
                                                                         A_col,C_col,
                                                                         G_col,T_col,
                                                                         i,j,k,
                                                                         plot_col_insert[[j]])
                }
                plot_insert[[j]][i,2]<-sum(plot_insert[[j]][i,1],plot_insert[[j]][i,2])
                plot_insert[[j]][i,3]<-sum(plot_insert[[j]][i,2],plot_insert[[j]][i,3])
                plot_insert[[j]][i,4]<-sum(plot_insert[[j]][i,3],plot_insert[[j]][i,4])
                
                chr<-substr(calling[[i]][j,1],4,nchar(calling[[i]][j,1]))
                pos<-calling[[i]][j,2]
                if(dbsnp_file!=""){
                    vcfFile<-system.file(dbsnp_file, package="VariantAnnotation")
                    param = ScanVcfParam(fixed=character(),info=character(), geno=character(), 
                                         samples=character(),which=GRanges(seqnames=chr,ranges=IRanges(pos,pos)))
                    vcf.comp<-readVcf(dbsnp_file,"hg19",param=param)
                }
                plot_ref[[j]]<-plotCountsPerTarget_ref(samples,calling,
                                                       dbsnp_file,vcf.comp,j,i,
                                                       plot_ref[[j]])
            }
        }
        if(abs_max==0){
            abs_max<-1
        }
        counter2<-0
        for(j in 1:length(calling[[1]][,1])){
            counter<-counter2
            counter2<-plotCountsPerTargetVcf_plot(directory2,samples,plot_abs,plot_col,
                                                  abs_max,calling,plot_insert,plot_col_insert,
                                                  A_col,C_col,G_col,T_col,plot_ref,marks,
                                                  qual_lower_bound,qual_upper_bound,counter,j)
        }
        return()
    }
    
    #Plot counts per Target
    plotCountsPerTargetRelativeVcf<-function(samples,
                                             directory2,
                                             calling,
                                             dbsnp_file,
                                             frequency_weight,
                                             qual_lower_bound,
                                             qual_upper_bound,
                                             marks){
        plot_rel<-list()
        plot_ref<-list()
        plot_col<-list()
        plot_insert<-list()
        plot_col_insert<-list()
        A_Pal<-colorRampPalette(c('chartreuse','forestgreen'))
        A_col<-A_Pal(10*(qual_upper_bound-qual_lower_bound)+1)
        C_Pal<-colorRampPalette(c('cyan','blue4'))
        C_col<-C_Pal(10*(qual_upper_bound-qual_lower_bound)+1)
        G_Pal<-colorRampPalette(c('yellow','darkgoldenrod3'))
        G_col<-G_Pal(10*(qual_upper_bound-qual_lower_bound)+1)
        T_Pal<-colorRampPalette(c('firebrick1','darkred'))
        T_col<-T_Pal(10*(qual_upper_bound-qual_lower_bound)+1)
        if((qual_upper_bound-qual_lower_bound)==0){
            A_col<-"forestgreen"
            C_col<-"blue4"
            G_col<-"darkgoldenrod3"
            T_col<-"darkred"
        }
        
        for(j in 1:length(calling[[1]][,1])){
            plot_rel[[j]]<-matrix(rep(0,6*length(samples[,1])),ncol=6)
            plot_ref[[j]]<-matrix(rep(0,6*length(samples[,1])),ncol=6)
            plot_col[[j]]<-matrix(rep(0,6*length(samples[,1])),ncol=6)
            plot_insert[[j]]<-matrix(rep(0,4*length(samples[,1])),ncol=4)
            plot_col_insert[[j]]<-matrix(rep(NA,4*length(samples[,1])),ncol=4)
        }
        rel_max<-1
        long_insert<-rep(0,length(samples[,1]))
        
        for(j in 1:length(calling[[1]][,1])){
            for(i in 1:length(samples[,1])){
                if(!is.na(frequency_weight[[i]][j-long_insert[i],1])
                   &&paste(calling[[i]][j,1],calling[[i]][j,2])
                   !=paste(frequency_weight[[i]][j-long_insert[i],1],frequency_weight[[i]][j-long_insert[i],2])){
                    long_insert[i]<-long_insert[i]+1
                }
                if(is.na(frequency_weight[[i]][j-long_insert[i],1])){
                    long_insert[i]<-long_insert[i]+1
                }
                for(k in 4:9){
                    plot_rel[[j]][i,k-3]<-calling[[i]][j,k]
                    if(is.na(plot_rel[[j]][i,k-3])){
                        plot_rel[[j]][i,k-3]<-0
                    }
                    plot_col[[j]]<-plotCountsPerTargetRelative_col(samples,
                                                                   frequency_weight,
                                                                   long_insert,
                                                                   qual_lower_bound,
                                                                   qual_upper_bound,
                                                                   A_col,C_col,
                                                                   G_col,T_col,
                                                                   i,j,k,plot_col[[j]]) 
                }
                for(k in c(1,3,5,7)){
                    plot_insert[[j]]<-plotCountsPerTargetRelative_insert(samples,
                                                                         frequency_weight,
                                                                         long_insert,
                                                                         i,j,k,
                                                                         plot_col[[j]])
                    plot_col_insert[[j]]<-plotCountsPerTargetRelative_col_insert(samples,
                                                                                 frequency_weight,
                                                                                 long_insert,
                                                                                 qual_lower_bound,
                                                                                 qual_upper_bound,
                                                                                 A_col,C_col,
                                                                                 G_col,T_col,
                                                                                 i,j,k,
                                                                                 plot_col_insert[[j]]) 
                }
                plot_insert[[j]][i,2]<-sum(plot_insert[[j]][i,1],plot_insert[[j]][i,2])
                plot_insert[[j]][i,3]<-sum(plot_insert[[j]][i,2],plot_insert[[j]][i,3])
                plot_insert[[j]][i,4]<-sum(plot_insert[[j]][i,3],plot_insert[[j]][i,4])
                
                chr<-substr(calling[[i]][j,1],4,nchar(calling[[i]][j,1]))
                pos<-calling[[i]][j,2]
                if(dbsnp_file!=""){
                    vcfFile<-system.file(dbsnp_file, package="VariantAnnotation")
                    param = ScanVcfParam(fixed=character(),info=character(), geno=character(), 
                                         samples=character(),which=GRanges(seqnames=chr,ranges=IRanges(pos,pos)))
                    vcf.comp<-readVcf(dbsnp_file,"hg19",param=param)
                }
                plot_ref[[j]]<-plotCountsPerTarget_ref(samples,calling,
                                                       dbsnp_file,vcf.comp,j,i,
                                                       plot_ref[[j]]) 
            }
        }
        counter2<-0
        for(j in 1:length(calling[[1]][,1])){
            counter<-counter2
            counter2<-plotCountsPerTargetRelativeVcf_plot(directory2,samples,plot_rel,
                                                          plot_col,rel_max,calling,
                                                          plot_insert,plot_col_insert,A_col,
                                                          C_col,G_col,T_col,plot_ref,marks,
                                                          qual_lower_bound,qual_upper_bound,
                                                          counter,j)
        }
        return()
    }
    
    observeEvent(input$do, {
        shinyjs::html("text", "Complete analysis with BBCAnalyzer starts...<br><br>")
        if(file.exists(input$bam_input)==F){
            shinyjs::html("text", paste("Folder containing bam- and bai-files does not exist","<br>",sep=""), add = TRUE) 
            return()
        }
        if(file.exists(input$output)==F){
            shinyjs::html("text", paste("Output folder does not exist","<br>",sep=""), add = TRUE) 
            return()
        }
        if(input$vcf_input!=""&&file.exists(input$vcf_input)==F){
            shinyjs::html("text", paste("Folder containing vcf files does not exist","<br>",sep=""), add = TRUE) 
            return()
        }
        if(input$known_file!=""&&file.exists(input$known_file)==F){
            shinyjs::html("text", paste("Tabix file containing known variants does not exist","<br>",sep=""), add = TRUE) 
            return()
        }
        #Pipeline
        library("BSgenome")
        genome_choices<-available.genomes()
        genome_selected<-genome_choices[input$genome==genome_choices]
        genome_installed<-installed.genomes()
        genome_selected2<-genome_installed[genome_selected==genome_installed]
        if(length(genome_selected2)==0){
            shinyjs::html("text", paste("Selected genome not yet installed","<br>",sep=""), add = TRUE)
            shinyjs::html("text", paste("Currently installing genome: ",genome_selected," (might take some time...)","<br>",sep=""), add = TRUE)
            source("http://bioconductor.org/biocLite.R")
            biocLite(genome_selected)   
        }
        library(genome_selected,character.only = T)
        ref_genome<-getBSgenome(genome_selected)
        
        samples_pre<-strsplit(input$sample_names,split="\\n")
        samples<-as.data.frame(samples_pre[[1]])
        
        targetRegions_pre<-strsplit(input$target_regions,split="\\n")
        if(vcountPattern(pattern="\t",targetRegions_pre[[1]][1])==1){
            targetRegions<-data.frame(chr=NA,start=NA)
        }
        if(vcountPattern(pattern="\t",targetRegions_pre[[1]][1])==2){
            targetRegions<-data.frame(chr=NA,start=NA,end=NA)
        }
        for(i in 1:length(targetRegions_pre[[1]])){
            if(i==1){
                targetRegions[1,]<-strsplit(targetRegions_pre[[1]][1],split="\\t")[[1]]
            }
            if(i>1){
                targetRegions<-rbind(targetRegions,strsplit(targetRegions_pre[[1]][i],split="\\t")[[1]])
            }
        }
        for(i in 1:length(targetRegions[1,])){
            targetRegions[,1]<-as.numeric(targetRegions[,1])
        }
        target<-determineTarget(targetRegions)
        target[,1]<-as.numeric(target[,1])
        target[,2]<-as.numeric(target[,2])
        bam_available<-checkBamInput(samples,input$bam_input)
        if(bam_available==FALSE){
            return()
        }
        bq<-analyzeReads(samples,target,input$bam_input,input$output,input$MQ_threshold)
        
        if(input$relative_choice=="Relative"){
            relative<-TRUE
        }
        if(input$relative_choice=="Absolute"){
            relative<-FALSE
        }
        if(input$per_sample_choice=="Sample"){
            per_sample<-TRUE
        }
        if(input$per_sample_choice=="Position"){
            per_sample<-FALSE
        }
        qual_lower_bound<-input$qual_bounds[1]
        qual_upper_bound<-input$qual_bounds[2]
        
        if(nchar(input$vcf_input)==0){
            frequency<-frequencyAnalysis(samples,input$output,bq,input$BQ_threshold,ref_genome)
            calling<-reportVariants(samples,input$output,frequency,input$frequency_threshold,
                                    target)
            if(per_sample==TRUE&&relative==FALSE){
                plotCountsPerSample(samples,input$output,calling,input$known_file,
                                    frequency,qual_lower_bound,qual_upper_bound,
                                    input$marks)
                output$plot <- renderImage({outfile <- tempfile(fileext='.png')
                
                png(outfile, width=800, height=400*length(samples[,1]))
                par(mar=rep(0,4))
                layout(matrix(1:length(samples[,1]),ncol=1,byrow=T))
                for(i in 1:length(samples[,1])){
                    plot1<-readPNG(paste(input$output,"/",samples[i,1],".png",sep=""))
                    plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
                    rasterImage(plot1,0,0,1,1)
                }
                dev.off()
                
                list(src = outfile,alt = "Live plotting failed")
                }, deleteFile = TRUE)
            }
            if(per_sample==TRUE&&relative==TRUE){
                plotCountsPerSampleRelative(samples,input$output,calling,
                                            input$known_file,frequency,qual_lower_bound,
                                            qual_upper_bound,input$marks)
                
                output$plot <- renderImage({outfile <- tempfile(fileext='.png')
                
                png(outfile, width=800, height=400*length(samples[,1]))
                par(mar=rep(0,4))
                layout(matrix(1:length(samples[,1]),ncol=1,byrow=T))
                for(i in 1:length(samples[,1])){
                    plot1<-readPNG(paste(input$output,"/",samples[i,1],".png",sep=""))
                    plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
                    rasterImage(plot1,0,0,1,1)
                }
                dev.off()
                
                list(src = outfile,alt = "Live plotting failed")
                }, deleteFile = TRUE)
            }
            if(per_sample==FALSE&&relative==FALSE){
                plotCountsPerTarget(samples,input$output,calling,input$known_file,
                                    frequency,qual_lower_bound,qual_upper_bound,
                                    input$marks)
                output$plot <- renderImage({outfile <- tempfile(fileext='.png')
                
                png(outfile, width=800, height=300*length(target[,1]))
                par(mar=rep(0,4))
                layout(matrix(1:length(target[,1]),ncol=1,byrow=T))
                for(i in 1:length(target[,1])){
                    plot1<-readPNG(paste(input$output,"/","chr",target[i,1],";",target[i,2],".png",sep=""))
                    plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
                    rasterImage(plot1,0,0,1,1)
                }
                dev.off()
                
                list(src = outfile,alt = "Live plotting failed")
                }, deleteFile = TRUE)
            }
            if(per_sample==FALSE&&relative==TRUE){
                plotCountsPerTargetRelative(samples,input$output,calling,
                                            input$known_file,frequency,qual_lower_bound,
                                            qual_upper_bound,input$marks)
                output$plot <- renderImage({outfile <- tempfile(fileext='.png')
                
                png(outfile, width=800, height=300*length(target[,1]))
                par(mar=rep(0,4))
                layout(matrix(1:length(target[,1]),ncol=1,byrow=T))
                for(i in 1:length(target[,1])){
                    plot1<-readPNG(paste(input$output,"/","chr",target[i,1],";",target[i,2],".png",sep=""))
                    plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
                    rasterImage(plot1,0,0,1,1)
                }
                dev.off()
                
                list(src = outfile,alt = "Live plotting failed")
                }, deleteFile = TRUE)
            }
        }
        if(nchar(input$vcf_input)>0){
            vcf_available<-checkVcfInput(samples,input$vcf_input)
            if(vcf_available==FALSE){
                return()
            }
            frequency<-frequencyAnalysisVcf(samples,input$vcf_input,input$output,bq,
                                            input$BQ_threshold,ref_genome)
            calling<-reportVariantsVcf(samples,input$output,frequency,input$frequency_threshold,
                                       target)
            if(per_sample==TRUE&&relative==FALSE){
                plotCountsPerSampleVcf(samples,input$output,calling,input$known_file,
                                       frequency,qual_lower_bound,qual_upper_bound,
                                       input$marks)
                output$plot <- renderImage({outfile <- tempfile(fileext='.png')
                
                png(outfile, width=800, height=400*length(samples[,1]))
                par(mar=rep(0,4))
                layout(matrix(1:length(samples[,1]),ncol=1,byrow=T))
                for(i in 1:length(samples[,1])){
                    plot1<-readPNG(paste(input$output,"/",samples[i,1],".png",sep=""))
                    plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
                    rasterImage(plot1,0,0,1,1)
                }
                dev.off()
                
                list(src = outfile,alt = "Live plotting failed")
                }, deleteFile = TRUE)
            }
            if(per_sample==TRUE&&relative==TRUE){
                plotCountsPerSampleRelativeVcf(samples,input$output,calling,
                                               input$known_file,frequency,
                                               qual_lower_bound,qual_upper_bound,
                                               input$marks)
                output$plot <- renderImage({outfile <- tempfile(fileext='.png')
                
                png(outfile, width=800, height=400*length(samples[,1]))
                par(mar=rep(0,4))
                layout(matrix(1:length(samples[,1]),ncol=1,byrow=T))
                for(i in 1:length(samples[,1])){
                    plot1<-readPNG(paste(input$output,"/",samples[i,1],".png",sep=""))
                    plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
                    rasterImage(plot1,0,0,1,1)
                }
                dev.off()
                
                list(src = outfile,alt = "Live plotting failed")
                }, deleteFile = TRUE)
            }
            if(per_sample==FALSE&&relative==FALSE){
                plotCountsPerTargetVcf(samples,input$output,calling,input$known_file,
                                       frequency,qual_lower_bound,qual_upper_bound,
                                       input$marks)
                output$plot <- renderImage({outfile <- tempfile(fileext='.png')
                
                png(outfile, width=800, height=300*length(target[,1]))
                par(mar=rep(0,4))
                layout(matrix(1:length(target[,1]),ncol=1,byrow=T))
                for(i in 1:length(target[,1])){
                    plot1<-readPNG(paste(input$output,"/","chr",target[i,1],";",target[i,2],".png",sep=""))
                    plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
                    rasterImage(plot1,0,0,1,1)
                }
                dev.off()
                
                list(src = outfile,alt = "Live plotting failed")
                }, deleteFile = TRUE)
            }
            if(per_sample==FALSE&&relative==TRUE){
                plotCountsPerTargetRelativeVcf(samples,input$output,calling,
                                               input$known_file,frequency,
                                               qual_lower_bound,qual_upper_bound,
                                               input$marks)
                output$plot <- renderImage({outfile <- tempfile(fileext='.png')
                
                png(outfile, width=800, height=300*length(target[,1]))
                par(mar=rep(0,4))
                layout(matrix(1:length(target[,1]),ncol=1,byrow=T))
                for(i in 1:length(target[,1])){
                    plot1<-readPNG(paste(input$output,"/","chr",target[i,1],";",target[i,2],".png",sep=""))
                    plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
                    rasterImage(plot1,0,0,1,1)
                }
                dev.off()
                
                list(src = outfile,alt = "Live plotting failed")
                }, deleteFile = TRUE)
            }
        }
        shinyjs::html("text", paste("<br><br>Complete analysis with BBCAnalyzer successful!",sep=""), add = TRUE)
    })
    
    observeEvent(input$do2, {
        shinyjs::html("text", paste("Creating plots with BBCAnalyzer...<br><br>",sep=""), add = FALSE)
        if(file.exists(input$output2)==F){
            shinyjs::html("text", paste("Output folder does not exist","<br>",sep=""), add = TRUE) 
            return()
        }
        samples_pre<-strsplit(input$sample_names2,split="\\n")
        samples<-as.data.frame(samples_pre[[1]])
        
        calling<-list()
        for(i in 1:length(samples[,1])){
            if(file.exists(paste(input$output2,samples[i,1],".calling.txt",sep=""))){
                test<-read.table(paste(input$output2,samples[i,1],".calling.txt",sep=""),stringsAsFactors = F,header=T)
                if(length(test[1,])==15){
                    calling[[i]]<-read.table(paste(input$output2,samples[i,1],".calling.txt",sep=""),stringsAsFactors = F,header=T,
                                             colClasses = c("character","numeric","character",rep("numeric",6),rep("character",6)))  
                }
                if(length(test[1,])==19){
                    calling[[i]]<-read.table(paste(input$output2,samples[i,1],".calling.txt",sep=""),stringsAsFactors = F,header=T,
                                             colClasses = c("character","numeric","character",rep("numeric",6),rep("character",10)))  
                }
            }
            if(file.exists(paste(input$output2,samples[i,1],".calling.txt",sep=""))==FALSE){
                shinyjs::html("text", paste("Analysis of sample: ",samples[i,1]," wit BBCAnalyzer not yet complete <br>",
                                            "Please perform complete analysis of sample ",samples[i,1]," with BBCAnalyzer","<br>",sep=""), add = TRUE)
                return()
            }
        }
        frequency<-list()
        for(i in 1:length(samples[,1])){
            if(file.exists(paste(input$output2,samples[i,1],".frequency.txt",sep=""))){
                test<-read.table(paste(input$output2,samples[i,1],".frequency.txt",sep=""),stringsAsFactors = F,header=T)
                if(length(test[1,])==22){
                    frequency[[i]]<-read.table(paste(input$output2,samples[i,1],".frequency.txt",sep=""),stringsAsFactors = F,header=T,
                                               colClasses = c("character","numeric","character",rep("numeric",19))) 
                }
                if(length(test[1,])==26){
                    frequency[[i]]<-read.table(paste(input$output2,samples[i,1],".frequency.txt",sep=""),stringsAsFactors = F,header=T,
                                               colClasses = c("character","numeric","character",rep("numeric",19),rep("character",4))) 
                }
            }
        }
        test_list<-list()
        for(i in 1:length(calling)){
            input_vector<-c()
            for(j in 1:length(calling[[i]][,1])){
                input_vector[j]<-paste(calling[[i]][j,1],calling[[i]][j,2],sep=";")
            }
            test_list[[i]]<-input_vector
        }
        target_common<-Reduce(intersect,test_list)
        if(length(target_common)==0){
            shinyjs::html("text", paste("Selected samples share no previously analyzed, intersecting target region","<br>",sep=""), add = TRUE)
            return()
        }
        if(length(target_common)>0){
            for(i in 1:length(test_list)){
                calling[[i]]<-calling[[i]][test_list[[i]]==target_common,]
                frequency[[i]]<-frequency[[i]][test_list[[i]]==target_common,]
            }
            if(input$relative_choice2=="Relative"){
                relative<-TRUE
            }
            if(input$relative_choice2=="Absolute"){
                relative<-FALSE
            }
            if(input$per_sample_choice2=="Sample"){
                per_sample<-TRUE
            }
            if(input$per_sample_choice2=="Position"){
                per_sample<-FALSE
            }
            qual_lower_bound<-input$qual_bounds2[1]
            qual_upper_bound<-input$qual_bounds2[2]
            
            if(input$vcf_input2=="No"){
                if(per_sample==TRUE&&relative==FALSE){
                    plotCountsPerSample(samples,input$output2,calling,input$known_file2,
                                        frequency,qual_lower_bound,qual_upper_bound,
                                        input$marks2)
                    output$plot <- renderImage({outfile <- tempfile(fileext='.png')
                    
                    png(outfile, width=800, height=400*length(samples[,1]))
                    par(mar=rep(0,4))
                    layout(matrix(1:length(samples[,1]),ncol=1,byrow=T))
                    for(i in 1:length(samples[,1])){
                        plot1<-readPNG(paste(input$output2,"/",samples[i,1],".png",sep=""))
                        plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
                        rasterImage(plot1,0,0,1,1)
                    }
                    dev.off()
                    
                    list(src = outfile,alt = "Live plotting failed")
                    }, deleteFile = TRUE)
                }
                if(per_sample==TRUE&&relative==TRUE){
                    plotCountsPerSampleRelative(samples,input$output2,calling,
                                                input$known_file2,frequency,qual_lower_bound,
                                                qual_upper_bound,input$marks2)
                    
                    output$plot <- renderImage({outfile <- tempfile(fileext='.png')
                    
                    png(outfile, width=800, height=400*length(samples[,1]))
                    par(mar=rep(0,4))
                    layout(matrix(1:length(samples[,1]),ncol=1,byrow=T))
                    for(i in 1:length(samples[,1])){
                        plot1<-readPNG(paste(input$output2,"/",samples[i,1],".png",sep=""))
                        plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
                        rasterImage(plot1,0,0,1,1)
                    }
                    dev.off()
                    
                    list(src = outfile,alt = "Live plotting failed")
                    }, deleteFile = TRUE)
                }
                if(per_sample==FALSE&&relative==FALSE){
                    plotCountsPerTarget(samples,input$output2,calling,input$known_file2,
                                        frequency,qual_lower_bound,qual_upper_bound,
                                        input$marks2)
                    output$plot <- renderImage({outfile <- tempfile(fileext='.png')
                    
                    png(outfile, width=800, height=300*length(target_common))
                    par(mar=rep(0,4))
                    layout(matrix(1:length(target_common),ncol=1,byrow=T))
                    for(i in 1:length(target_common)){
                        plot1<-readPNG(paste(input$output2,"/",target_common[i],".png",sep=""))
                        plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
                        rasterImage(plot1,0,0,1,1)
                    }
                    dev.off()
                    
                    list(src = outfile,alt = "Live plotting failed")
                    }, deleteFile = TRUE)
                }
                if(per_sample==FALSE&&relative==TRUE){
                    plotCountsPerTargetRelative(samples,input$output2,calling,
                                                input$known_file2,frequency,qual_lower_bound,
                                                qual_upper_bound,input$marks2)
                    output$plot <- renderImage({outfile <- tempfile(fileext='.png')
                    
                    png(outfile, width=800, height=300*length(target_common))
                    par(mar=rep(0,4))
                    layout(matrix(1:length(target_common),ncol=1,byrow=T))
                    for(i in 1:length(target_common)){
                        plot1<-readPNG(paste(input$output2,"/",target_common[i],".png",sep=""))
                        plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
                        rasterImage(plot1,0,0,1,1)
                    }
                    dev.off()
                    
                    list(src = outfile,alt = "Live plotting failed")
                    }, deleteFile = TRUE)
                }
            }
            if(input$vcf_input2=="Yes"){
                for(i in 1:length(calling)){
                    if(length(calling[[i]][1,])!=19){
                        shinyjs::html("text", paste("vcf file information not available for all samples","<br>",sep=""), add = TRUE)
                        return()
                    }
                    if(length(frequency[[i]][1,])!=26){
                        shinyjs::html("text", paste("vcf file information not available for all selected samples","<br>",sep=""), add = TRUE)
                        return()
                    }
                }
                if(per_sample==TRUE&&relative==FALSE){
                    plotCountsPerSampleVcf(samples,input$output2,calling,input$known_file2,
                                           frequency,qual_lower_bound,qual_upper_bound,
                                           input$marks2)
                    output$plot <- renderImage({outfile <- tempfile(fileext='.png')
                    
                    png(outfile, width=800, height=400*length(samples[,1]))
                    par(mar=rep(0,4))
                    layout(matrix(1:length(samples[,1]),ncol=1,byrow=T))
                    for(i in 1:length(samples[,1])){
                        plot1<-readPNG(paste(input$output2,"/",samples[i,1],".png",sep=""))
                        plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
                        rasterImage(plot1,0,0,1,1)
                    }
                    dev.off()
                    
                    list(src = outfile,alt = "Live plotting failed")
                    }, deleteFile = TRUE)
                }
                if(per_sample==TRUE&&relative==TRUE){
                    plotCountsPerSampleRelativeVcf(samples,input$output2,calling,
                                                   input$known_file2,frequency,
                                                   qual_lower_bound,qual_upper_bound,
                                                   input$marks2)
                    output$plot <- renderImage({outfile <- tempfile(fileext='.png')
                    
                    png(outfile, width=800, height=400*length(samples[,1]))
                    par(mar=rep(0,4))
                    layout(matrix(1:length(samples[,1]),ncol=1,byrow=T))
                    for(i in 1:length(samples[,1])){
                        plot1<-readPNG(paste(input$output2,"/",samples[i,1],".png",sep=""))
                        plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
                        rasterImage(plot1,0,0,1,1)
                    }
                    dev.off()
                    
                    list(src = outfile,alt = "Live plotting failed")
                    }, deleteFile = TRUE)
                }
                if(per_sample==FALSE&&relative==FALSE){
                    plotCountsPerTargetVcf(samples,input$output2,calling,input$known_file2,
                                           frequency,qual_lower_bound,qual_upper_bound,
                                           input$marks2)
                    output$plot <- renderImage({outfile <- tempfile(fileext='.png')
                    
                    png(outfile, width=800, height=300*length(target_common))
                    par(mar=rep(0,4))
                    layout(matrix(1:length(target_common),ncol=1,byrow=T))
                    for(i in 1:length(target_common)){
                        plot1<-readPNG(paste(input$output2,"/",target_common[i],".png",sep=""))
                        plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
                        rasterImage(plot1,0,0,1,1)
                    }
                    dev.off()
                    
                    list(src = outfile,alt = "Live plotting failed")
                    }, deleteFile = TRUE)
                }
                if(per_sample==FALSE&&relative==TRUE){
                    plotCountsPerTargetRelativeVcf(samples,input$output2,calling,
                                                   input$known_file2,frequency,
                                                   qual_lower_bound,qual_upper_bound,
                                                   input$marks2)
                    output$plot <- renderImage({outfile <- tempfile(fileext='.png')
                    
                    png(outfile, width=800, height=300*length(target_common))
                    par(mar=rep(0,4))
                    layout(matrix(1:length(target_common),ncol=1,byrow=T))
                    for(i in 1:length(target_common)){
                        plot1<-readPNG(paste(input$output2,"/",target_common[i],".png",sep=""))
                        plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
                        rasterImage(plot1,0,0,1,1)
                    }
                    dev.off()
                    
                    list(src = outfile,alt = "Live plotting failed")
                    }, deleteFile = TRUE)
                }
            }
            shinyjs::html("text", paste("Creating plots successful!",sep=""), add = TRUE)
            
        }
        
        
        
        
        
    })
}
)