#!/usr/bin/env Rscript


args <- commandArgs(TRUE)

if (length(args)==0){
  stop("missing arguments !\nUsage : ./getContigsPerloci.R <output directory> <DiffContigsInfos.tsv> <nb cores>")
}

library(foreach)

library(doParallel)


########### arguments ###############
#####################################

#output directory
home<-args[1]
setwd(home)

#input table of Contigs with their infos.
initial_table<-args[2]

nb_cores<-args[3]



########### function ###############
#this function will clusterize Contigs following their locus (genic, antisense, intergenic, unmapped)
#####################################

#argument "initial_table" : DiffContigsInfos.tsv 
#other arguments : columns names from the initial table, used during the processing
#nb_cores : number of cores to use during the processing
#remark : if the name of the columns from the input data has changed, you just have to change the values here for each argument
#it returns a data frame
clusterContigs<-function(initial_table="",
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
                         nb_cores=1){
                                                
  initial_table<-read.delim(initial_table,sep="\t",header=T)
  
  mappedContigs<-initial_table[which(initial_table[is_mapped]==T),]
  
  unmappedContigs<-initial_table[which(initial_table[is_mapped]==F),]
  
  if(nrow(unmappedContigs)>0){unmappedContigs<-unmappedContigs[order(unmappedContigs[pvalue]),]}
  
  #for each contig, create a new ID with chrs, strand, 5'gene, 3'gene, separated by "&" (for genic/antisense contigs, we will just keep sense gene/antisense gene and strand)
  IDs_matrix<-as.matrix(data.frame(as.character(mappedContigs[,EnsemblGeneCol]),
                                   as.character(mappedContigs[,Ensembl_AS_GeneCol]),
                                   as.character(mappedContigs[,strand]),
                                   paste(mappedContigs[,chromosome],mappedContigs[,strand],mappedContigs[,gene_5p],mappedContigs[,gene_3p],sep="&")))
  
  getRefinedIDs<-function(EnsemblGeneCol,Ensembl_AS_GeneCol,strand,initial_id){

	    if(EnsemblGeneCol!="none"){
	      
	      return(paste(EnsemblGeneCol,strand,sep="&"))
	      
	      
	    }else if(Ensembl_AS_GeneCol!="none"){

	      return(paste(Ensembl_AS_GeneCol,strand,sep="&"))
	      
	    }else{

	      return(initial_id)

	    }

  }

  refined_IDs<-apply(IDs_matrix,1,function(IDs_matrix){

      getRefinedIDs(IDs_matrix[1],IDs_matrix[2],IDs_matrix[3],IDs_matrix[4])

  })

  rm(IDs_matrix);gc()

  mappedContigs$new_ID<-refined_IDs
  
  #unique on the new IDs
  myIDs<-unique(refined_IDs)
  
  print(paste("number of cores to use : ",nb_cores,sep=""))
  
  registerDoParallel(cores=nb_cores)
  
  #retrieve for each ID, all Contigs
  result<-foreach(i=1:length(myIDs), .combine=rbind) %dopar% {
    
	    ContigsOfOneID<-mappedContigs[which(mappedContigs["new_ID"]==myIDs[i]),]
	    
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
	    
	    #if the Contigs overlap a sense gene, gene name = sense gene ID
	    if(as.character(ContigsOfOneID[1,EnsemblGeneCol])!="none"){
	      
	      gene_name<-paste(as.character(ContigsOfOneID[1,HugoGeneCol]))
	      
	      locusType<-"genic"
	      
	      #number of kmers in the locus
	      nb_kmers<-sum(as.numeric(as.character(ContigsOfOneID[,nb_merged_kmers])))
	      
	      #if the Contigs don't overlap a sense gene and overlap an antisense gene, gene name = antisense gene ID
	    }else if(as.character(ContigsOfOneID[1,EnsemblGeneCol])=="none" & as.character(ContigsOfOneID[1,Ensembl_AS_GeneCol])!="none"){
	      
	      gene_name<-paste(as.character(ContigsOfOneID[1,Hugo_AS_GeneCol]))
	      
	      locusType<-"antisense"
	      
	      #number of kmers in the locus
	      nb_kmers<-sum(as.numeric(as.character(ContigsOfOneID[,nb_merged_kmers])))
	      
	      #if the Contigs don't overlap a sense gene and an antisense gene, they are intergenics, gene name = none
	    }else{
	      
	      #intergenic Contigs may have missing neighbors ; if both are missing, just keep the chromosome
	      locus_ID<-gsub("&none&none$","",locus_ID)
	      
	      gene_name<-"none"
	      
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
  
  #add the corresponding columns for the "locus" of unmapped contigs
  if(nrow(unmappedContigs)>0){
  
	  result<-rbind(result,data.frame(locus_ID="unmapped",
		                          gene_name="none",
		                          nb_contigs=nrow(unmappedContigs),
		                          nb_kmers=sum(as.numeric(as.character(unmappedContigs[,nb_merged_kmers]))),
		                          Best_P_value=unmappedContigs[1,pvalue],
		                          log2FC_of_Best=unmappedContigs[1,log2FC],
		                          meanCond1_of_Best=unmappedContigs[1,grep("mean",value=T,names(unmappedContigs))][1],
		                          meanCond2_of_Best=unmappedContigs[1,grep("mean",value=T,names(unmappedContigs))][2],
		                          chromosome="none",
		                          start=NA,
		                          end=NA,
		                          strand="none",
		                          locus_type="unmapped"))
  }
  
  return(result)
}

#example
myresult<-clusterContigs(initial_table=initial_table,nb_cores=nb_cores)
write.table(myresult,file="ContigsPerLoci.tsv",sep="\t",row.names=F, col.names=T, quote=F)
