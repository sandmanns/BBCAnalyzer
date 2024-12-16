#Determine single target bases from target regions:
determineTarget<-function(targetRegions){
    message("Determine Target")
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
            
            start<-targetRegions[targetLine,2]+1
            end<-targetRegions[targetLine,3]
            
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
    message("Analyze Reads")
    bases<-list()
    quality<-list()
    
    for(j in 1:length(samples[,1])){
        bases[[j]]<-data.frame()
        quality[[j]]<-data.frame()
    }
    counter<-0
    for(j in 1:length(samples[,1])){
        message("Sample ",j," out of ",length(samples[,1]))
        bases_insert<-data.frame()
        quality_insert<-data.frame()
        counter<-counter+1
        long_inserts<-0
        for(i in 1:length(target[,1])){
            message("Position ",i," out of ",length(target[,1]))
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
                for(k in 1:length(bam_info_rname)){ 
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
