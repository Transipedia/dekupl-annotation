#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

if (length(args)==0){
  stop("missing arguments !\nUsage : ./getIntronsByTranscripts.R <gff3 annotation>")
}

library(GenomicFeatures)

file<-args[1]

options(ucscChromosomeNames=F)

txdb<-makeTxDbFromGFF(file,format="gff")

IntronsByTranscripts<-unlist(intronsByTranscript(txdb,use.names=T))

transcriptsByGenes<-unlist(transcriptsBy(txdb,by="gene"))

transcriptsByGenes<-data.frame(transcript_id=transcriptsByGenes$tx_name,gene_id=names(transcriptsByGenes))

intronsGFF<-data.frame(transcript_id=names(IntronsByTranscripts),
                       chromosome=seqnames(IntronsByTranscripts),
                       source="annotation",
                       feature="intron",
                       start=start(IntronsByTranscripts),
                       end=end(IntronsByTranscripts),
                       empty1=".",
                       strand=strand(IntronsByTranscripts),
                       empty2=".",
                       ID=paste("ID=",names(IntronsByTranscripts),"_",start(IntronsByTranscripts),"_",end(IntronsByTranscripts),sep=""),
                       parent=paste(";Parent=",names(IntronsByTranscripts),sep=""),
                       transcript_id2=paste(";transcript_id=",names(IntronsByTranscripts),sep=""))

intronsGFF<-merge(intronsGFF,transcriptsByGenes, 
                  by.x="transcript_id",
                  by.y="transcript_id",
                  all.x=TRUE,
                  all.y=FALSE)

names(intronsGFF)[ncol(intronsGFF)]<-"gene_id"

intronsGFF$gene_id<-paste(";gene_id=",intronsGFF$gene_id,sep="")

intronsGFF<-intronsGFF[c("chromosome","source","feature","start","end","empty1","strand","empty2","ID","parent","gene_id","transcript_id2")]

intronsGFF$attributes<-paste(intronsGFF$ID,intronsGFF$parent,intronsGFF$gene_id,intronsGFF$transcript_id2,sep="")

intronsGFF<-subset(intronsGFF,select=-c(ID,parent,gene_id,transcript_id2))

write.table(intronsGFF,stdout(), sep='\t',row.names=F,col.names=F,quote=F)


