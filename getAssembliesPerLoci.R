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

#input table of assemblies with their infos.
mytable<-args[2]

registerDoParallel(cores=args[3])


#input table for some tests
#mytable<-"/home/marc/Documents/pipeline_test/DiffAssembliesInfos.txt"


########### function ###############
#this function will group assemblies following a shared ID for each type of locus (genic,antisense,intergenic)
#####################################

#argument "myAssemblies" : initial table 
#other arguments : columns names from the initial table, used during the processing
#it returns a data frame
#remark : if the name of the columns from the input data has changed, you just have to change the values here for each argument
#If you want to give this function directly an input data frame, put also in comment the line with read.delim
getAssembliesPerLoci<-function(myAssemblies,
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
  
    myAssemblies<-read.delim(myAssemblies,sep="\t",header=T)
    
    myAssemblies<-myAssemblies[which(myAssemblies[is_mapped]==T),]
    
    #for each assembly, create a new ID with chrs,strand,5'gene,3'gene, separated by "&"
    myAssemblies$new_ID<-paste(myAssemblies[,chromosome],myAssemblies[,strand],myAssemblies[,gene_5p],myAssemblies[,gene_3p],sep="&")
    
    #unique on the new IDs
    myIDs<-unique(as.character(myAssemblies[,"new_ID"]))
     
    #retrieve for each ID, all assemblies
    result<-foreach(i=1:length(myIDs), .combine=rbind) %dopar% {
    
        AssembliesOfOneID<-myAssemblies[which(myAssemblies["new_ID"]==myIDs[i]),]
        
        NumberOfAssemblies<-nrow(AssembliesOfOneID)
        
        min_start<-min(AssembliesOfOneID[start])
        
        max_end<-max(AssembliesOfOneID[end])
        
        best_P_value<-min(AssembliesOfOneID[pvalue])
        
        min_log2FC<-min(AssembliesOfOneID[log2FC])
        
        max_log2FC<-max(AssembliesOfOneID[log2FC])
        
        if(abs(min_log2FC)>max_log2FC){best_log2FC<-min_log2FC}else{best_log2FC<-max_log2FC}
        
        #if the assemblies overlap a sense gene, locus_ID = sense gene ID
        if(as.character(AssembliesOfOneID[1,EnsemblGeneCol])!="none"){
          
          locus_ID<-as.character(AssembliesOfOneID[1,EnsemblGeneCol])
          gene_name<-as.character(AssembliesOfOneID[1,HugoGeneCol])
          locusType<-"genic"
          
        #if the assemblies don't overlap a sense gene and overlap an antisense gene, locus_ID = antisense gene ID
        }else if(as.character(AssembliesOfOneID[1,EnsemblGeneCol])=="none" & as.character(AssembliesOfOneID[1,Ensembl_AS_GeneCol])!="none"){
          
          locus_ID<-as.character(AssembliesOfOneID[1,Ensembl_AS_GeneCol])
          gene_name<-as.character(AssembliesOfOneID[1,Hugo_AS_GeneCol])
          locusType<-"antisense"
         
        #if the assemblies don't overlap a sense gene and an antisense gene, they are intergenics, locus_ID = chrs & strand & 5'-gene & 3'-gene
        }else{
          
          #intergenic assemblies may have missing neighbors ; if both are missing, just keep the chromosome
          locus_ID<-gsub("&none&none$","",as.character(AssembliesOfOneID$new_ID[1]))
          gene_name="none"
          locusType<-"intergenic"
          
        }
        
                #re-assign to the existing data frame a$line_in_SAMt the right line, the new values
                 data.frame(locus_ID=locus_ID,
                            gene_name=gene_name,
                            NumberOfAssemblies=NumberOfAssemblies,
                            best_log2FC=best_log2FC,
                            best_P_value=best_P_value,
                            chromosome=as.character(AssembliesOfOneID[1,chromosome]),
                            start=min_start,
                            end=max_end,
                            strand=as.character(AssembliesOfOneID[1,strand]),
                            locus_type=locusType,stringsAsFactors= F)
     }
    
    #plot the results
    loci_frame<-data.frame(genic=length(result[which(result$locus_type=="genic"),1]),antisense=length(result[which(result$locus_type=="antisense"),1]),intergenic=length(result[which(result$locus_type=="intergenic"),1]))
    pdf("AssembliesPerLoci.pdf",width=15,height=10)

    barplot(as.matrix(loci_frame),ylim=c(0,max(loci_frame)+max(loci_frame)/4),
            ylab="loci number",
            xlab="locus type",
            main=paste("Loci with differentially expressed assemblies (total = ",nrow(result),")",sep=""),
            width=1)
    dev.off()
  
    return(result)
}

#example
myresult<-getAssembliesPerLoci(mytable)
write.table(myresult,file="AssembliesPerLoci.txt",sep="\t",row.names=F, col.names=T, quote=F)

