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

nb_cores<-args[3]



########### function ###############
#this function will clusterize Contigs following their locus (genic,antisense,intergenic)
#####################################

#argument "myContigs" : initial table 
#other arguments : columns names from the initial table, used during the processing
#nb_cores : number of cores to use during the processing
#remark : if the name of the columns from the input data has changed, you just have to change the values here for each argument
#it returns a data frame
clusterContigs<-function(myContigs,
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
                         log2FC="log2FC",
                         nb_merged_kmers="nb_merged_kmers",
                         nb_cores){
  
  myContigs<-read.delim(myContigs,sep="\t",header=T)
  
  myContigs<-myContigs[which(myContigs[is_mapped]==T),]
  
  #for each contig, create a new ID with chrs, strand, 5'gene, 3'gene, separated by "&" (for genic or antisense contigs, we just keep sense gene/antisense gene and strand)
  myContigs$new_ID<-paste(myContigs[,chromosome],myContigs[,strand],myContigs[,gene_5p],myContigs[,gene_3p],sep="&")
  
  for(line in 1:nrow(myContigs)){
    
    if(as.character(myContigs[line,EnsemblGeneCol])!="none"){
      
      myContigs[line,"new_ID"]<-paste(myContigs[line,EnsemblGeneCol],myContigs[line,strand],sep="&")
      
      
    }else if(as.character(myContigs[line,Ensembl_AS_GeneCol])!="none"){
      
      myContigs[line,"new_ID"]<-paste(myContigs[line,Ensembl_AS_GeneCol],myContigs[line,strand],sep="&")
      
    }
    
  }
  
  #unique on the new IDs
  myIDs<-unique(as.character(myContigs[,"new_ID"]))
  
  #set the number of cores to use
  if(nb_cores==""){nb_cores=getDoParWorkers()}
  
  print(paste("number of cores to use : ",nb_cores,sep=""))
  
  registerDoParallel(cores=nb_cores)
  
  #retrieve for each ID, all Contigs
  result<-foreach(i=1:length(myIDs), .combine=rbind) %dopar% {
    
    ContigsOfOneID<-myContigs[which(myContigs["new_ID"]==myIDs[i]),]
    
    #number of contigs in the locus
    nb_contigs<-nrow(ContigsOfOneID)
    
    min_start<-min(ContigsOfOneID[start])
    
    max_end<-max(ContigsOfOneID[end])
    
    ContigsOfOneID<-ContigsOfOneID[order(ContigsOfOneID[pvalue],-abs(ContigsOfOneID[log2FC])),]
    
    Best_P_value<-ContigsOfOneID[1,pvalue]
    
    log2FC_of_Best<-ContigsOfOneID[1,log2FC]
    
    meanCond1_of_Best<-ContigsOfOneID[1,grep("mean",value=T,names(ContigsOfOneID))][1]
    meanCond2_of_Best<-ContigsOfOneID[1,grep("mean",value=T,names(ContigsOfOneID))][2]
    
    locus_ID<-ContigsOfOneID[1,"new_ID"]
    
    #if the Contigs overlap a sense gene, locus_ID = sense gene ID
    if(as.character(ContigsOfOneID[1,EnsemblGeneCol])!="none"){
      
      gene_name<-paste(as.character(ContigsOfOneID[1,HugoGeneCol]))
      
      locusType<-"genic"
      
      #number of kmers in the locus
      nb_kmers<-sum(as.numeric(as.character(ContigsOfOneID[,nb_merged_kmers])))
      
      #if the Contigs don't overlap a sense gene and overlap an antisense gene, locus_ID = antisense gene ID
    }else if(as.character(ContigsOfOneID[1,EnsemblGeneCol])=="none" & as.character(ContigsOfOneID[1,Ensembl_AS_GeneCol])!="none"){
      
      gene_name<-paste(as.character(ContigsOfOneID[1,Hugo_AS_GeneCol]))
      
      locusType<-"antisense"
      
      #number of kmers in the locus
      nb_kmers<-sum(as.numeric(as.character(ContigsOfOneID[,nb_merged_kmers])))
      
      #if the Contigs don't overlap a sense gene and an antisense gene, they are intergenics, locus_ID = chrs & strand & 5'-gene & 3'-gene
    }else{
      
      #intergenic Contigs may have missing neighbors ; if both are missing, just keep the chromosome
      locus_ID<-gsub("&none&none$","",locus_ID)
      
      gene_name="none"
      
      locusType<-"intergenic"
      
      nb_kmers<-sum(as.numeric(as.character(ContigsOfOneID[,nb_merged_kmers])))
      
    }
    
    data.frame(locus_ID=locus_ID,
               gene_name=gene_name,
               nb_contigs=nb_contigs,
               nb_kmers=nb_kmers,
               Best_P_value=Best_P_value,
               log2FC_of_Best=log2FC_of_Best,
               meanCond1_of_Best,
               meanCond2_of_Best,
               chromosome=as.character(ContigsOfOneID[1,chromosome]),
               start=min_start,
               end=max_end,
               strand=as.character(ContigsOfOneID[1,strand]),
               locus_type=locusType,stringsAsFactors= F)
  }
  
  return(result)
}

#example
myresult<-clusterContigs(mytable)
write.table(myresult,file="ContigsPerLoci.tsv",sep="\t",row.names=F, col.names=T, quote=F)

