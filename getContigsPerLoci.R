#!/usr/bin/env Rscript


args <- commandArgs(TRUE)

if (length(args)==0){
  stop("no argument !")
}

library(foreach)

library(doParallel)


########### arguments ###############
#####################################

#output directory
home<-args[1]
setwd(home)

#input table of Contigs with their infos.
mytable<-args[2]

registerDoParallel(cores=args[3])



########### function ###############
#this function will group Contigs following a shared ID for each type of locus (genic,antisense,intergenic)
#####################################

#argument "myContigs" : initial table 
#other arguments : columns names from the initial table, used during the processing
#it returns a data frame
#remark : if the name of the columns from the input data has changed, you just have to change the values here for each argument
#If you want to give this function directly an input data frame, put also in comment the line with read.delim
getContigsPerLoci<-function(myContigs,
                               is_mapped="is_mapped",
                               chromosome="chromosome",
                               start="start",
                               end="end",
                               strand="strand",
                               EnsemblGeneCol="gene",
                               HugoGeneCol="HUGO_ID_gene",
                               Ensembl_AS_GeneCol="as_gene",
                               Hugo_AS_GeneCol="HUGO_ID_as_gene",
                               gene_5p="gene_5p",
                               gene_3p="gene_3p",
                               pvalue="pvalue",
                               log2FC="log2FC"){
  
    myContigs<-read.delim(myContigs,sep="\t",header=T)
    
    myContigs<-myContigs[which(myContigs[is_mapped]==T),]
    
    #for each contig, create a new ID with chrs,strand,5'gene,3'gene, separated by "&"
    myContigs$new_ID<-paste(myContigs[,chromosome],myContigs[,strand],myContigs[,gene_5p],myContigs[,gene_3p],sep="&")
    
    #unique on the new IDs
    myIDs<-unique(as.character(myContigs[,"new_ID"]))
     
    #retrieve for each ID, all Contigs
    result<-foreach(i=1:length(myIDs), .combine=rbind) %dopar% {
    
        ContigsOfOneID<-myContigs[which(myContigs["new_ID"]==myIDs[i]),]
        
        NumberOfContigs<-nrow(ContigsOfOneID)
        
        min_start<-min(ContigsOfOneID[start])
        
        max_end<-max(ContigsOfOneID[end])
        
        best_P_value<-min(ContigsOfOneID[pvalue])
        
        min_log2FC<-min(ContigsOfOneID[log2FC])
        
        max_log2FC<-max(ContigsOfOneID[log2FC])
        
        if(abs(min_log2FC)>max_log2FC){best_log2FC<-min_log2FC}else{best_log2FC<-max_log2FC}
        
        #if the Contigs overlap a sense gene, locus_ID = sense gene ID
        if(as.character(ContigsOfOneID[1,EnsemblGeneCol])!="none"){
          
          locus_ID<-as.character(ContigsOfOneID[1,EnsemblGeneCol])
          gene_name<-as.character(ContigsOfOneID[1,HugoGeneCol])
          locusType<-"genic"
          
        #if the Contigs don't overlap a sense gene and overlap an antisense gene, locus_ID = antisense gene ID
        }else if(as.character(ContigsOfOneID[1,EnsemblGeneCol])=="none" & as.character(ContigsOfOneID[1,Ensembl_AS_GeneCol])!="none"){
          
          locus_ID<-as.character(ContigsOfOneID[1,Ensembl_AS_GeneCol])
          gene_name<-as.character(ContigsOfOneID[1,Hugo_AS_GeneCol])
          locusType<-"antisense"
         
        #if the Contigs don't overlap a sense gene and an antisense gene, they are intergenics, locus_ID = chrs & strand & 5'-gene & 3'-gene
        }else{
          
          #intergenic Contigs may have missing neighbors ; if both are missing, just keep the chromosome
          locus_ID<-gsub("&none&none$","",as.character(ContigsOfOneID$new_ID[1]))
          gene_name="none"
          locusType<-"intergenic"
          
        }
        
                #re-assign to the existing data frame a$line_in_SAMt the right line, the new values
                 data.frame(locus_ID=locus_ID,
                            gene_name=gene_name,
                            NumberOfContigs=NumberOfContigs,
                            best_log2FC=best_log2FC,
                            best_P_value=best_P_value,
                            chromosome=as.character(ContigsOfOneID[1,chromosome]),
                            start=min_start,
                            end=max_end,
                            strand=as.character(ContigsOfOneID[1,strand]),
                            locus_type=locusType,stringsAsFactors= F)
     }
    
    #plot the results
    # loci_frame<-data.frame(genic=length(result[which(result$locus_type=="genic"),1]),antisense=length(result[which(result$locus_type=="antisense"),1]),intergenic=length(result[which(result$locus_type=="intergenic"),1]))
    # pdf("ContigsPerLoci.pdf",width=15,height=10)
    # 
    # barplot(as.matrix(loci_frame),ylim=c(0,max(loci_frame)+max(loci_frame)/4),
    #         ylab="loci number",
    #         xlab="locus type",
    #         main=paste("Loci with differentially expressed Contigs (total = ",nrow(result),")",sep=""),
    #         width=1)
    # dev.off()
  
    return(result)
}

#example
myresult<-getContigsPerLoci(mytable)
write.table(myresult,file="ContigsPerLoci.tsv",sep="\t",row.names=F, col.names=T, quote=F)
