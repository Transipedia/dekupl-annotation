#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

if (length(args)==0){
  stop("missing arguments !\nUsage : ./getSwitches.R <output directory> <DiffContigsInfos.tsv> <sample_conditions.tsv> <normalized_counts.tsv>")
}

library(DESeq2)

#### input data ####

home<-args[1]

#DiffContigsInfos from dekupl
all_contigs<-args[2]

#design file from dekupl
sample_conditions<-args[3]

#normalized counts (gene expression) from dekupl (Kallisto)
normalizedGeneCounts<-args[4]

all_contigs<-read.delim(all_contigs)

sample_conditions<-read.delim(sample_conditions)

normalizedGeneCounts<-read.delim(normalizedGeneCounts)

#### process data ####

#retrieve mapped contigs
mapped_DE_contigs<-all_contigs[which(all_contigs$gene%in%normalizedGeneCounts$id),]

#keep in normalized counts only rows in which we have mapped contigs
normalizedGeneCounts<-normalizedGeneCounts[which(normalizedGeneCounts$id%in%mapped_DE_contigs$gene),]

#contigs with their counts : keep only columns contig ID, gene ID and contig counts
#after the 37th column (log2FC), we have the counts
tab_counts_DEkupl<-mapped_DE_contigs[,c(1,22,38:(ncol(mapped_DE_contigs)))]

#intersect KALLISTO gene IDs & DEKUPL gene IDs
tab_counts_Kallisto<-merge(tab_counts_DEkupl,normalizedGeneCounts, by.x="gene", by.y="id", all.x=T, all.y=F)

#reorganize columns in order to have contig ID, gene ID & KALLISTO counts (same order as tab_counts_DEkupl)
tab_counts_Kallisto<-tab_counts_Kallisto[,c(2,1,(ncol(tab_counts_DEkupl)+1):(ncol(tab_counts_Kallisto)))]

#order both tables following the contig ID
tab_counts_Kallisto<-tab_counts_Kallisto[order(tab_counts_Kallisto$ID),] 

tab_counts_DEkupl<-tab_counts_DEkupl[order(tab_counts_DEkupl$ID),]

#keep the same header for both tables
names(tab_counts_Kallisto)[3:ncol(tab_counts_Kallisto)]<-names(tab_counts_DEkupl)[3:ncol(tab_counts_DEkupl)]

#prepare contigs with their counts for DESeq2 (row names = contig ID, and we keep only counts without any other columns)
rownames(tab_counts_DEkupl)<-tab_counts_DEkupl$ID 

tab_counts_DEkupl[,c(1,2)]<-NULL

#get conditions name
cond1<-as.character(sample_conditions[1,2])

cond2<-unique(as.character(sample_conditions[,2][sample_conditions[,2]!=cond1]))

#get number of samples for each condition
rep_cond1<-nrow(sample_conditions[which(sample_conditions[,2]==cond1),])

rep_cond2<-nrow(sample_conditions[which(sample_conditions[,2]==cond2),])

#set design
samples<-data.frame(row.names=names(tab_counts_DEkupl),condition=c(rep(cond1,rep_cond1),rep(cond2,rep_cond2)))

#create DESeqDataSet object from design & contigs DE
dds<-DESeqDataSetFromMatrix(countData=as.matrix(round(tab_counts_DEkupl)),colData=samples,design=~condition)

#compute normalization factor for each contig at each sample, thanks to their normalized gene counts from Kallisto
normFactors<-as.matrix((tab_counts_Kallisto[,3:ncol(tab_counts_Kallisto)]+1)/exp(rowMeans(log(tab_counts_Kallisto[,3:ncol(tab_counts_Kallisto)]+1))))

#allocation of normalization factors 
normalizationFactors(dds)<-normFactors

#estimate overdispersion parameters
#it's possible to have issues with estimateDispersions() if you have a low number of contigs ("dispersion trend not well captured")
#so, we use fitType="mean" instead of the default "parametric"
getDispersions<-function(my_object=""){
  
  dds<-try(estimateDispersions(my_object))
  
  if (class(dds)=="try-error"){
    
    cat("with fitType='parametric', the dispersion trend was not well captured by the function, we will use instead fitType='mean'")
    
    dds<-estimateDispersions(my_object,fitType="mean")
    
  }
  
  return(dds)
}

dds<-getDispersions(dds)

#binomiale negative test 
dds<-nbinomWaldTest(dds)

#results 
DESeq2Result<-results(dds)

#extract padj 
DESeq2Result<-DESeq2Result["padj"]

#make a custom table with contig ID,mean cond1, mean cond2, log2FC, padj, normalized counts for all libraries
new_result<-data.frame(id=row.names(DESeq2Result),DU_Pvalue=DESeq2Result$padj,row.names=NULL)

#merge the initial table of contigs with the result
all_contigs<-merge(all_contigs,new_result, by.x="ID", by.y="id", all.x=T, all.y=F)

#unmapped/intergenic contigs are given the value "NA"
all_contigs$DU_Pvalue[is.na(all_contigs$DU_Pvalue)]<-"NA"

write.table(all_contigs,file=paste(home,"DiffContigsInfos.tsv",sep=""),sep="\t",row.names=F, col.names=T, quote=F)

