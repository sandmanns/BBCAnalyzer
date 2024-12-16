#Variant calling:
#define threshold for frequency of variant
reportVariants<-function(samples,
                         directory2,
                         frequency_weight,
                         frequency_threshold,
                         target){
    message("Report Variants")
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
        message("Sample ",j)
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
    message("Report Variants")
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
        message("Sample ",j)
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
