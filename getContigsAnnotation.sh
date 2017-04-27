#!/bin/bash

 usage() { echo -e "Usage: $0 -a <merged-diff-counts.tsv.gz (contigs from DEkupl)> -g <genome in fasta> -d <A_vs_B_DEGs.tsv (diff. genes from DEkupl)> -r <reference annotation (gff3 format)> -t <stranded data (choose between \"yes\" or \"no\")> -o <full path to output directory> -i <illumina adapters (you can use the file \"adapters.fa\" supplied with the program)> [options]\n\n\tOptions :\n
                  -b <path to bin/ of blast scripts (default : in \$PATH environment variable)>\n
                  -c <path to bedtools, preferentially 2.24 (default : in \$PATH environment variable)>\n
                  -j <GSNAP index name>\n
                  -k <path to directory of GSNAP index>\n
                  -m <path to bin/ of GSNAP (default : in \$PATH environment variable)>\n
                  -p <padj diff. gene threshold (default : 0.05)>\n
                  -s <path to samtools (default : in \$PATH environment variable)>\n
                  -n <thread number (default : 1)>\n\n\tResults :\n
                  - Table \"DiffContigsInfos.tsv\", summarizing for each contig, its location on the genome (if it's aligned), the neighborhood, the sequence alignment informations, and the differential expression informations\n
                  - BED \"file diff_contigs.bed\" for the visualization ; it contains useful informations from the summarization table.\n
                  
                  - Table \"contigsPerLoci.txt\" containing loci with differentially expressed contigs\n" 1>&2; exit 1;}

[[ $# -eq 0 ]] && usage

while getopts ":a:g:d:r:i:o:b:c:j:k:m:p:s:t:n:" opt; do
  case $opt in
  
      a)
      
	      DEkupl_result=$OPTARG
	      
	      echo -e "\ncontigs file from DEkupl is : $DEkupl_result\n" >&2
	      
              ;;
      
      g)
      
	      ref_fasta=$OPTARG
	      
	      echo -e "\nreference genome in fasta is : $ref_fasta\n" >&2
      
      	      ;;
      
      d)
      
              diff_genes=$OPTARG
 
              echo -e "\ndiff genes file from DEkupl is : $diff_genes\n" >&2
             
              ;;
      
      r)
             ref_annotation=$OPTARG
 
             echo -e "\nreference annotation in gff3 format is: $ref_annotation\n" >&2
      
      
             ;;
      
      i)
             illumina_adapters=$OPTARG
 
             echo -e "\nillumina adapters file is: $illumina_adapters\n" >&2
      
      
             ;;
      
      o)
      
             output_dir=$OPTARG
             
             echo -e "\noutput directory is : $output_dir\n" >&2
      
             ;;

      b)
      
             ncbi_blast_loc=$OPTARG

             ;;
      
      c)
      
             bedtools=$OPTARG
     
             ;;
      
      j)
      
             gsnap_index_name=$OPTARG
     
             ;;
      
      k)
      
             gsnap_index_dir=$OPTARG
      
             ;;
      
      m)
      
             GSNAP_loc=$OPTARG

             ;;
      
      
      p)
      
             padj_threshold=$OPTARG

             ;;
      
      s)
      
            samtools=$OPTARG
       
            ;;
      t)
      
            stranded=$OPTARG
       
            ;;            
            
      
      n)
      
            threads=$OPTARG

            ;;
      
      
      #invalid options (options not in the list)
      ######################
      
      
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    
  esac
done

##### scripts/functions that will be used by the main script
############################################################

getContigsPerLoci=$(dirname "$0")/getContigsPerLoci.R
blast=$(dirname "$0")/blast.sh
source $(dirname "$0")/getAnnotationFunctions.sh

##### check whether programs are supplied/installed
#we look for samtools, bedtools, blast directory, and gsnap directory
#####################################################################

#if samtools is supplied, print the version ; otherwise, we check in the $PATH environment variable ; if it's still missing, exit
if [ "$samtools" == "" ];then samtools=$(which samtools);if [ "$samtools" == "" ];then echo "samtools not found !";exit ;else echo "samtools version is : $($samtools --version) " ; fi ; else echo -e "samtools version is : $($samtools --version)\n" ;fi

#if bedtools is supplied, print the version ; otherwise, we check in the $PATH environment variable ; if it's still missing, exit
if [ "$bedtools" == "" ];then bedtools=$(which bedtools);if [ "$bedtools" == "" ];then echo "bedtools not found !";exit ;else echo "bedtools version is : $($bedtools -version) " ; fi ; else echo -e "bedtools version is : $($bedtools -version)\n" ;fi

#if blast scripts directory is supplied, print the location ; otherwise, we check in the $PATH environment variable ; if it's still missing, exit
if [ "$ncbi_blast_loc" == "" ];then ncbi_blast_loc=$(which blastn);if [ "$ncbi_blast_loc" == "" ];then echo "blast scripts directory not found !"; exit ; else ncbi_blast_loc=$(dirname $ncbi_blast_loc);echo -e "blast scripts directory is : ${ncbi_blast_loc}\n" ; fi; else echo -e "blast scripts directory is : $ncbi_blast_loc\n" ;fi
ncbi_blast_loc="${ncbi_blast_loc}/"
ncbi_blast_loc=$(echo "$ncbi_blast_loc" |sed 's/\/\//\//g')

#if gsnap directory is supplied, print the location ; otherwise, we check in the $PATH environment variable ; if it's still missing, exit
if [ "$GSNAP_loc" == "" ];then GSNAP_loc=$(which gsnap);if [ "$GSNAP_loc" == "" ];then echo "gsnap directory not found !"; exit ; else GSNAP_loc=$(dirname $GSNAP_loc);echo -e "gsnap directory is : ${GSNAP_loc}\n" ; fi; else echo -e "gsnap directory is : $GSNAP_loc\n" ;fi
GSNAP_loc="${GSNAP_loc}/"
GSNAP_loc=$(echo "$GSNAP_loc" |sed 's/\/\//\//g')


#if one of useful argument is missing, exit with the usage 
 
if [ "$DEkupl_result" == "" ] || [ "$ref_fasta" == "" ] || [ "$diff_genes" == "" ] || [ "$ref_annotation" == "" ] || [ "$output_dir" == "" ] || [ "$illumina_adapters" == "" ] || [ "$samtools" == "" ] || [ "$bedtools" == "" ] || [ "$ncbi_blast_loc" == "" ] || [ "$GSNAP_loc" == "" ] || ([ "$stranded" != "yes" ] && [ "$stranded" != "no" ]); then
      
echo -e "\none required argument is missing or is wrong (check also if all the programs are installed) !!\n"
      usage
      exit 1
fi


echo -e "stranded contigs : $stranded\n"


if [ "$threads" == "" ];then threads=1;fi
echo -e "threads number is : $threads\n"

if [ "$padj_threshold" == "" ];then padj_threshold=0.05;fi


#2.26 is buggy !! use 2.24
#bedtools="/home/marc/Downloads/bedtools2/bin/bedtools"

##### processing of the input variables
########################################

output_dir="${output_dir}/"

output_dir=$(echo "$output_dir" | sed 's/\/\//\//g')

if [ ! -d $output_dir ];then mkdir $output_dir ;fi

mapping_output="${output_dir}mapping_output/"

if [ ! -d $mapping_output ];then mkdir $mapping_output ;fi

#final tab file
FinalTable="${output_dir}DiffContigsInfos.tsv"
if [ -f $FinalTable ];then rm $FinalTable ;fi

diff_contigs_bed="${output_dir}diff_contigs.bed"
if [ -f $diff_contigs_bed ];then rm $diff_contigs_bed ;fi

#processing of the gff
grep -v "^#" $ref_annotation | sort -k1,1 -k4,4n >${output_dir}ref_annotation.tmp

ref_annotation=${output_dir}ref_annotation.tmp

zcat $DEkupl_result | awk 'NR==1{OFS="\t";print "ID",$0;exit}' >${output_dir}diffexFileHeader.txt

#convert the DEkupl contigs output in fasta, and at the same time put the tag of reference as an ID for joinings
zcat $DEkupl_result | awk 'NR>1{print}' | awk 'OFS="\t"{print $3,$0}' | tee ${output_dir}diffexFile.txt | awk 'OFS="\n"{print ">"$1,$3}' >${output_dir}OriginalFastaTags.fa 

#col 1 =ID ; col 2 = contig seq : will be used in joinings
cat ${output_dir}OriginalFastaTags.fa | paste - - | awk 'OFS="\t"{print $1,$2}' | LANG=en_EN sort -k1,1 >${output_dir}OriginalFastaTags.tmp


cat ${output_dir}diffexFileHeader.txt ${output_dir}diffexFile.txt >${output_dir}diffexFile.tmp && mv ${output_dir}diffexFile.tmp ${output_dir}diffexFile.txt && rm ${output_dir}diffexFileHeader.txt

##### blast of tags against adapters, and retain only those not matching
########################################################################

#if [ ! -f ${output_dir}nomatch_in_adapters.fa ];then 

echo -e "\n==== Filter out contigs matching in adapters ====\n"

echo -e "  blast on adapters !\n"

bash $blast -q ${output_dir}OriginalFastaTags.fa -s $illumina_adapters -o $output_dir -p $ncbi_blast_loc -d "illumina_adapters" -n $threads || { echo "blast on adapters failure!!" 1>&2; exit; }

#for checkings
awk '!seen[$1]++' "${output_dir}raw_blast.alignment_2.txt" | sed $'s/ /\t/g' >${output_dir}match_in_adapters.txt

grep "^>" ${output_dir}OriginalFastaTags.fa | sed 's/>//g' | sort >${output_dir}originaltagIDs.txt

comm -23 ${output_dir}originaltagIDs.txt <(cut -f 1 ${output_dir}raw_blast.alignment_2.txt | LANG=en_EN sort -u) | awk 'OFS="\t"{print ">"$1}' | LANG=en_EN sort -k1,1 >${output_dir}nomatch_in_adapters.txt && rm ${output_dir}originaltagIDs.txt

#reconstruction of the fasta with only tags not matching in adapters
LANG=en_EN join -t $'\t' -11 -21 ${output_dir}nomatch_in_adapters.txt ${output_dir}OriginalFastaTags.tmp | sed 's/\t/\n/g' >${output_dir}nomatch_in_adapters.tmp && rm ${output_dir}nomatch_in_adapters.txt && mv ${output_dir}nomatch_in_adapters.tmp ${output_dir}nomatch_in_adapters.fa
	

#fi

##### mapping of the contigs
###############################
 
#if no alignment file in BAM format for the contigs, build it 
#if [ ! -f ${mapping_output}contigs.bam ];then


echo -e "\n==== Mapping of contigs on the genome ====\n"

echo -e "number of contigs to align : $(($(wc -l ${output_dir}nomatch_in_adapters.fa|awk '{print $1}')/2))\n"

start_date=$(date)

#if no GSNAP index, build it 
if [ "$gsnap_index_name" == "" ] || [ "$gsnap_index_dir" == "" ];then

	echo -e "  no GASNAP index for the genome, we're going to make it !"

	gsnap_index_name="genome_index"

	gsnap_index_dir=$mapping_output 

	${GSNAP_loc}gmap_build -D $mapping_output -d $gsnap_index_name $ref_fasta  || { echo "gsnap db building failure!!" 1>&2; exit; }


fi

gsnap_index_dir="${gsnap_index_dir}/"
gsnap_index_dir=$(echo "$gsnap_index_dir" |sed 's/\/\//\//g')

${GSNAP_loc}gsnap -t $threads -A sam -N 1 -D $gsnap_index_dir -d $gsnap_index_name ${output_dir}nomatch_in_adapters.fa |$samtools view -bh >${mapping_output}contigs.bam || { echo "gsnap mapping failure!!" 1>&2; exit; }
end_date=$(date)

echo -e "\nstart mapping of contigs : $start_date\n"
echo -e "\nend mapping of contigs : $end_date\n"
	  
		  
		  
		  
		  
#fi
  
################## for each alignment line, reconstruction of the tag structure (splices, exons...)
###################################################################################################

echo -e "\n==== Reconstruction of the tag structure (splices, exons...) from the alignment ====\n"

#this process could be parallelized (split the file in many, run the while-loop on them, store the results in tmp files 1,2,3..., at the end merge the results)

start_date=$(date)

#from the alignment, get a bed and a table of the contigs with infos about the alignment
parseBam ${mapping_output}contigs.bam $FinalTable $diff_contigs_bed

$samtools view -f 4 ${mapping_output}contigs.bam | sort -k1,1 | sort -u -k1,1 | awk -v start_line=$(($(wc -l $diff_contigs_bed | awk '{print $1}')+1)) 'BEGIN{a=start_line}OFS="\t"{print a,$1;a=a+1}' >${output_dir}OriginalUnmappedTags.txt

end_date=$(date)

echo -e "\nstart parse BAM : $start_date\n"
echo -e "\nend parse BAM : $end_date\n"

############## blast of unmapped seq and add the good ones (full length alignment) in the final table/bed
#########################################################################################################

if [ -f ${output_dir}OriginalUnmappedTags.txt ];then

start_date=$(date)


 #if there are unmapped contigs, we'll try to align them with blast (actually we are looking for those with full length)
 if [[ $(wc -l ${output_dir}OriginalUnmappedTags.txt|awk '{print $1}') -gt 0 ]];then

	echo -e "\n==== blast of unmapped contigs : retrieving those mapping with their full length ====\n"

        #this function will rebuild a fasta file from the unmapped contigs
	getFastaFromUnmappedTags ${output_dir}OriginalUnmappedTags.txt ${output_dir}OriginalFastaTags.tmp ${output_dir}unmapped_fasta1.fa 

	bash $blast -q ${output_dir}unmapped_fasta1.fa -s $ref_fasta -o $output_dir -p $ncbi_blast_loc -d "human_genome" -n $threads

	#if alignment length = tag length, it's a good case missed by the mapper, we make a unique on the tag ID !
	#print tag ID,chr,start,end,strand,nb_mismatch,gaps , sort by mismaches, then gaps, then ID, then chromosome
	awk 'OFS="\t"{if($4==$13){if($14=="plus"){$14="+"}else{$14="-"};print $1,$2,$9,$10,$14,$5,$6}}' ${output_dir}raw_blast.alignment_2.txt | sort -k6,6n -k7,7n -k1,1 -k2,2 | awk '!seen[$1]++'>${output_dir}UnmappedEntirelyFoundByBlast.txt

	awk 'OFS="\t"{if($4==$13){if($14=="plus"){$14="+"}else{$14="-"};print $1,$2,$9,$10,$14,$5,$6}}' ${output_dir}raw_blast.alignment_2.txt >${output_dir}UnmappedEntirelyFoundByBlast.tmp1


	if [ -f ${output_dir}UnmappedEntirelyFoundByBlast.tmp2 ];then rm ${output_dir}UnmappedEntirelyFoundByBlast.tmp2 ; fi

        #re-insert the found alignments in the final table, after computing
	while read line ;do

	  mismaches=$(echo "$line" | cut -f 6)
	  gaps=$(echo "$line" | cut -f 7)
	  ID=$(echo "$line" | cut -f 1)
	  
	  #if there's an alignment as good as the best, it's a multihit
	  multihit=$(awk -v nb_mismatch=$mismaches -v gaps=$gaps -v ID=$ID 'BEGIN{a=0}{if($1==ID && $6==nb_mismatch && $7==gaps){a=a+1}}END{print a}' ${output_dir}UnmappedEntirelyFoundByBlast.tmp1)
	  
	  #1st block : print tag ID,lineInSAM,chr,start,end,strand,nb_mismatch,gaps,nb_hit
	  #2nd block : reorganize the columns in order to add them to the final table, SNV is added just after ins : line_in_SAM,tagID,is_mapped,mapper,chr,start,end,junction,nb_junction,exon_coord_list,other_split,strand,nb_hit,mis,del,ins,SNV,clipped_5p,clipped_3p
	  echo -e "$line\t$multihit" | awk 'OFS="\t"{sub("-","\t",$1);print}' | awk 'OFS="\t"{if($7>0 || $8>0){SNV="T"}else{SNV="F"};print $2,$1,"T","Blast",$3,$4,$5,"none",0,$3":"$4"-"$5,"F",$6,$9,$7,$8,0,SNV,0,0}' | tee -a ${output_dir}UnmappedEntirelyFoundByBlast.tmp2 | awk 'OFS="\t"{print}' >>$FinalTable

	done < ${output_dir}UnmappedEntirelyFoundByBlast.txt
	
	rm ${output_dir}UnmappedEntirelyFoundByBlast.tmp1 ${output_dir}UnmappedEntirelyFoundByBlast.txt

	if [ -f ${output_dir}UnmappedEntirelyFoundByBlast.tmp2 ];then

	  #put the result in the bed file 
	  awk 'OFS="\t"{if($7>$6){$6=$6-1;print $5,$6,$7,"LineInSam="$1";ID="$2";nb_hit="$13";nM="$14";del="$15";ins=0;clipped_5p=0;clipped_3p=0",1,$12,$6,$7,"255,0,0",1,$7-$6,0}else{$7=$7-1;print $5,$7,$6,"LineInSam="$1";ID="$2";nb_hit="$13";nM="$14";del="$15";ins=0;clipped_5p=0;clipped_3p=0",1,$12,$7,$6,"0,0,255",1,$6-$7,0}}' ${output_dir}UnmappedEntirelyFoundByBlast.tmp2 >>$diff_contigs_bed

          #contigs still unmapped after the previous method
	  comm -23 <(cut -f 1 ${output_dir}OriginalUnmappedTags.txt | LANG=en_EN sort) <(cut -f 1 ${output_dir}UnmappedEntirelyFoundByBlast.tmp2|LANG=en_EN sort ) | LANG=en_EN join -t $'\t' -11 -21 - <(LANG=en_EN sort -k1,1 ${output_dir}OriginalUnmappedTags.txt) >${output_dir}unmapped_tags_2.txt && rm ${output_dir}UnmappedEntirelyFoundByBlast.tmp2
	
	else

          #if the previous method doesn't give results, just re-use the initial unmapped contigs
          cat ${output_dir}OriginalUnmappedTags.txt >${output_dir}unmapped_tags_2.txt

        fi


  fi

end_date=$(date)

echo -e "\nstart blast of unmapped contigs on genome : $start_date\n"
echo -e "\nend blast of unmapped contigs on genome : $end_date\n"

fi

sort -k1,1 -k2,2n $diff_contigs_bed >${diff_contigs_bed}.tmp && mv ${diff_contigs_bed}.tmp $diff_contigs_bed

########## insert diffex infos in the bed file 
###############################################

LANG=en_EN join -t $'\t' -11 -21 <(awk 'OFS="\t"{initial_line=$0;ID_col=$4;print ID_col,initial_line}' $diff_contigs_bed | awk 'OFS="\t"{split ($1,x,";");a=x[2]; print a,$0}' | sed 's/^ID=//g' | awk 'OFS="\t"{$2="";print $0}'|LANG=en_EN sort -k1,1) <(awk 'NR>1{printf $1"\t";printf "%.3e\t",$5;printf "%.0f\t",$6;printf "%.0f\t",$7;printf "%.2f\n",$8}' ${output_dir}diffexFile.txt|LANG=en_EN sort -k1,1) | awk 'OFS="\t"{$5=$5";pval="$14";meanA="$15";meanB="$16";log2FC="$17;print}' | cut -f2-13 >${diff_contigs_bed}.tmp1 && mv ${diff_contigs_bed}.tmp1 ${diff_contigs_bed}


sort -k1,1 -k2,2n $diff_contigs_bed >${diff_contigs_bed}.tmp && mv ${diff_contigs_bed}.tmp $diff_contigs_bed


########## search max and min abs(log2FC) to scale the bed color (RGB) on them
#########################################################################

#max value for upregulated features
max_up_log2FC=$(zcat $DEkupl_result | awk 'NR>1{if($7 >0 || $7 <0){print $7}}'| awk '{if($1!~/-/){print}}' | sort -nr | head -n1)

#max value for downregulated features
max_down_log2FC=$(zcat $DEkupl_result | awk 'NR>1{if($7 >0 || $7 <0){print $7}}'| awk '{if($1~/-/){print}}' | sed 's/-//g'| sort -nr | head -n1)

#max absolute value
max_abs_log2FC=$(echo -e "${max_up_log2FC}\n${max_down_log2FC}"|sort -nr |head -n1)

#min value for upregulated features
min_up_log2FC=$(zcat $DEkupl_result | awk 'NR>1{if($7 >0 || $7 <0){print $7}}' | awk '{if($1!~/-/){print}}' | sort -n | head -n1)

#min value for downregulated features
min_down_log2FC=$(zcat $DEkupl_result |awk 'NR>1{if($7 >0 || $7 <0){print $7}}' | awk '{if($1~/-/){print}}' | sed 's/-//g'| sort -n | head -n1)

#min absolute value
min_abs_log2FC=$(echo -e "${min_up_log2FC}\n${min_down_log2FC}" | sort -n | head -n1)

#255,220,220 & 220,220,255 : light red and light blue
#255,0,0 & 0,0,255 : far red and far blue
#we will scale the color following the abs(log2FC)
#formula to range [min,max] to [a,b] :[((b-a)(x - min))/ (max - min)] +a

#join bed and diffex table (just column lineInSam and log2FC) by lineInSAM
LANG=en_EN join -t $'\t' -11 -21 <(awk 'OFS="\t"{initial_line=$0;ID_col=$4;print ID_col,initial_line}' $diff_contigs_bed | awk 'OFS="\t"{split ($1,x,";");a=x[2]; print a,$0}' | sed 's/^ID=//g' | awk 'OFS="\t"{$2="";print $0}'|LANG=en_EN sort -k1,1) <(awk 'NR>1{OFS="\t";print $1,$8}' ${output_dir}diffexFile.txt) | awk -v b=220 -v a=0 -v max=$max_abs_log2FC -v min=$min_abs_log2FC 'OFS="\t"{if($NF~/-/){abs_log2FC=-$14}else{abs_log2FC=$NF};RGB=220-((((b-a)*((0.5^(abs_log2FC))-(0.5^(min))))/((0.5^(max))-(0.5^(min))))+a);if($7=="+"){printf $0"\t";printf "255,";printf "%.0f",RGB;printf ",";printf "%.0f\n",RGB};if($7=="-"){printf $0"\t";printf "%.0f",RGB;printf ",";printf "%.0f",RGB;printf ",255\n"}}' | awk 'OFS="\t"{$10=$NF;print}' | cut -f 2-13 >${diff_contigs_bed}.tmp1 && mv ${diff_contigs_bed}.tmp1 ${diff_contigs_bed}

sort -k1,1 -k2,2n $diff_contigs_bed >${diff_contigs_bed}.tmp && mv ${diff_contigs_bed}.tmp $diff_contigs_bed


########## research of antisense & neighbours + location of the tags in the gene (UTR, exon...) ###########
###########################################################################################################

echo -e "\n==== research of the neighborhood of the contigs ====\n"

#genome size for bedtools 
$samtools view -H ${mapping_output}contigs.bam | grep -E -v "^@PG|^@HD" | awk 'OFS="\t"{print $2,$3}' | sed 's/SN\://;s/LN\://g' | sort -k1,1 >${mapping_output}genome_size.txt

start_date=$(date)

#searching antisense tags, retain only one target (sort by chr and start, then unique on attribute of the tag) (gives 2 cols in addition of the lineInSAM : Ensembl gene ID & HUGO ID)

if [ "$stranded" == "yes" ];then orientation_option="-s" ;dist_option="a";else orientation_option="";dist_option="ref"; fi

awk '{if($3=="gene"){print}}' $ref_annotation| $bedtools intersect -S -a $diff_contigs_bed -b - -loj -nonamecheck | sort -k1,1 -k2,2n | sort -u -k4,4 >${output_dir}antisense_tags.txt


if [ "$stranded" == "yes" ];then 

	paste -d'\t' <(awk 'OFS="\t"{print $4}' ${output_dir}antisense_tags.txt|awk -F';' '{print $1}' | awk '{sub("LineInSam=","",$1);print}') <(awk 'OFS="\t"{print $NF}' ${output_dir}antisense_tags.txt | awk -F';' 'OFS="\t"{print $1,$5}' | awk 'OFS="\t"{if($1=="."){$1="none";$2="none"};sub("ID=","",$1);sub("gene_name=","",$2);print}') >${output_dir}antisense_tags.tmp && mv ${output_dir}antisense_tags.tmp ${output_dir}antisense_tags.txt
	
else

	paste -d'\t' <(awk 'OFS="\t"{print $4}' ${output_dir}antisense_tags.txt | awk -F';' '{print $1}' | awk '{sub("LineInSam=","",$1);print}') <(awk 'OFS="\t"{print $NF}' ${output_dir}antisense_tags.txt | awk -F';' 'OFS="\t"{print "none","none"}') >${output_dir}antisense_tags.tmp && mv ${output_dir}antisense_tags.tmp ${output_dir}antisense_tags.txt

fi

LANG=en_EN join -t $'\t' -11 -21 <(LANG=en_EN sort -k1,1 $FinalTable) <(LANG=en_EN sort -k1,1 ${output_dir}antisense_tags.txt) >${FinalTable}.tmp && mv ${FinalTable}.tmp $FinalTable 


#searching sense tags (gives 2 cols in addition of the lineInSAM : Ensembl gene ID & HUGO ID)
awk '{if($3=="gene"){print}}' $ref_annotation | $bedtools intersect $orientation_option -a $diff_contigs_bed -b - -loj -nonamecheck | sort -k1,1 -k2,2n |sort -u -k4,4 >${output_dir}sense_tags.txt

paste -d'\t' <(awk 'OFS="\t"{print $4}' ${output_dir}sense_tags.txt|awk -F';' '{print $1}'|awk '{sub("LineInSam=","",$1);print}') <(awk 'OFS="\t"{print $NF}' ${output_dir}sense_tags.txt|awk -F';' 'OFS="\t"{print $1,$5}' | awk 'OFS="\t"{if($1=="."){$1="none";$2="none"};sub("ID=","",$1);sub("gene_name=","",$2);print}') >${output_dir}sense_tags.tmp && mv ${output_dir}sense_tags.tmp ${output_dir}sense_tags.txt

LANG=en_EN join -t $'\t' -11 -21 <(LANG=en_EN sort -k1,1 $FinalTable) <(LANG=en_EN sort -k1,1 ${output_dir}sense_tags.txt) >${FinalTable}.tmp && mv ${FinalTable}.tmp $FinalTable


#closest 5'end gene (gives 2 cols in addition of the lineInSAM : Ensembl gene ID & dist)
awk '{if($3=="gene"){print}}' $ref_annotation | $bedtools closest -nonamecheck -g ${mapping_output}genome_size.txt $orientation_option -io -fu -D $dist_option -t first -a $diff_contigs_bed -b - | sort -k1,1 -k2,2n | sort -u -k4,4 >${output_dir}closest_5end.txt

paste -d'\t' <(awk 'OFS="\t"{print $4}' ${output_dir}closest_5end.txt | awk -F';' '{print $1}'|awk '{sub("LineInSam=","",$1);print}') <(awk 'OFS="\t"{print $(NF-1)}' ${output_dir}closest_5end.txt | awk -F';' '{print $1}' | awk '{if($1=="."){$1="none"};sub("ID=","",$1);print}') <(awk 'OFS="\t"{print $(NF-1),$NF}' ${output_dir}closest_5end.txt | awk '{if($1=="."){$2="none"};print $2}') | awk 'OFS="\t"{if($3>0){$2="none";$3="none"}print}'>${output_dir}closest_5end.tmp && mv ${output_dir}closest_5end.tmp ${output_dir}closest_5end.txt

LANG=en_EN join -t $'\t' -11 -21 <(LANG=en_EN sort -k1,1 $FinalTable) <(LANG=en_EN sort -k1,1 ${output_dir}closest_5end.txt) >${FinalTable}.tmp && mv ${FinalTable}.tmp $FinalTable && rm ${output_dir}closest_5end.txt


#closest 3'end gene (gives 2 cols in addition of the lineInSAM : Ensembl gene ID & dist)
awk '{if($3=="gene"){print}}' $ref_annotation| $bedtools closest -nonamecheck -g ${mapping_output}genome_size.txt $orientation_option -io -fd -D $dist_option -t first -a $diff_contigs_bed -b - | sort -k1,1 -k2,2n | sort -u -k4,4 >${output_dir}closest_3end.txt

paste -d'\t' <(awk 'OFS="\t"{print $4}' ${output_dir}closest_3end.txt | awk -F';' '{print $1}' | awk '{sub("LineInSam=","",$1);print}') <(awk 'OFS="\t"{print $(NF-1)}' ${output_dir}closest_3end.txt | awk -F';' '{print $1}' | awk '{if($1=="."){$1="none"};sub("ID=","",$1);print}') <(awk 'OFS="\t"{print $(NF-1),$NF}' ${output_dir}closest_3end.txt |awk '{if($1=="."){$2="none"};print $2}') | awk 'OFS="\t"{if($3<0){$2="none";$3="none"}print}'>${output_dir}closest_3end.tmp && mv ${output_dir}closest_3end.tmp ${output_dir}closest_3end.txt

LANG=en_EN join -t $'\t' -11 -21 <(LANG=en_EN sort -k1,1 $FinalTable) <(LANG=en_EN sort -k1,1 ${output_dir}closest_3end.txt) >${FinalTable}.tmp && mv ${FinalTable}.tmp $FinalTable && rm ${output_dir}closest_3end.txt

#searching UTR tags (gives 1 col in addition of the lineInSAM : T or F)
awk '{if($3=="UTR"){print}}' $ref_annotation| $bedtools intersect $orientation_option -a $diff_contigs_bed -b - -loj -nonamecheck | sort -k1,1 -k2,2n | sort -u -k4,4 >${output_dir}UTR_tags.txt

paste -d'\t' <(awk 'OFS="\t"{print $4}' ${output_dir}UTR_tags.txt | awk -F';' '{print $1}' | awk '{sub("LineInSam=","",$1);print}') <(awk 'OFS="\t"{print $NF}' ${output_dir}UTR_tags.txt | awk -F';' '{if($1=="."){$1="F";print}else{print "T"}}') >${output_dir}UTR_tags.tmp && mv ${output_dir}UTR_tags.tmp ${output_dir}UTR_tags.txt

LANG=en_EN join -t $'\t' -11 -21 <(LANG=en_EN sort -k1,1 $FinalTable) <(LANG=en_EN sort -k1,1 ${output_dir}UTR_tags.txt) >${FinalTable}.tmp && mv ${FinalTable}.tmp $FinalTable 


#searching exonic tags ( gives 1 col in addition of the line_in_SAM : T or F)
awk '{if($3=="exon"){print}}' $ref_annotation | $bedtools intersect $orientation_option -a $diff_contigs_bed -b - -loj -nonamecheck | sort -k1,1 -k2,2n | sort -u -k4,4 >${output_dir}Exonic_tags.txt

paste -d'\t' <(awk 'OFS="\t"{print $4}' ${output_dir}Exonic_tags.txt | awk -F';' '{print $1}' | awk '{sub("LineInSam=","",$1);print}') <(awk 'OFS="\t"{print $NF}' ${output_dir}Exonic_tags.txt | awk -F';' '{if($1=="."){$1="F";print}else{print "T"}}') >${output_dir}Exonic_tags.tmp && mv ${output_dir}Exonic_tags.tmp ${output_dir}Exonic_tags.txt

LANG=en_EN join -t $'\t' -11 -21 <(LANG=en_EN sort -k1,1 $FinalTable) <(LANG=en_EN sort -k1,1 ${output_dir}Exonic_tags.txt) >${FinalTable}.tmp && mv ${FinalTable}.tmp $FinalTable

#searching intronic tags (it gives 1 col in addition of the line_in_SAM : T or F)
#join genic contigs, exonic contigs, and UTR contigs by line_in_SAM ; so we have 4 columns
#if the tag isn't genic, it cannot be intronic (intronic=F). In the case it is genic : if it's exonic or UTR, it cannot be intronic either(intronic=F). If it's not in the previous groups, it's intronic (intronic=T)
LANG=en_EN join -t $'\t' -11 -21 <(LANG=en_EN sort -k1,1 ${output_dir}sense_tags.txt) <(LANG=en_EN sort -k1,1 ${output_dir}Exonic_tags.txt) |LANG=en_EN join -t $'\t' -11 -21 - <(LANG=en_EN sort -k1,1 ${output_dir}UTR_tags.txt) | awk 'OFS="\t"{if($2=="none"){print $1,"F"}else{if($3=="T" || $4=="T"){print $1,"F"}else{print $1,"T"}}}'>${output_dir}Intronic_tags.txt

LANG=en_EN join -t $'\t' -11 -21 <(LANG=en_EN sort -k1,1 $FinalTable) <(LANG=en_EN sort -k1,1 ${output_dir}Intronic_tags.txt) >${FinalTable}.tmp && mv ${FinalTable}.tmp $FinalTable && rm ${output_dir}sense_tags.txt && rm ${output_dir}antisense_tags.txt && rm ${output_dir}Intronic_tags.txt && rm ${output_dir}Exonic_tags.txt && rm ${output_dir}UTR_tags.txt

#diff gene (gives 1 col)
awk -v padj_threshold=$padj_threshold 'NR>1{if($7!="NA" && $7<=padj_threshold){print $1}}' $diff_genes >${output_dir}diff_genes_ID.txt

# ! sensitive column ! : in final table we have to search sense gene (actually it's col 22); this give cols lineInSam & sense gene
awk 'NR>1{OFS="\t";if($22!="none"){print $1,$22}}' $FinalTable >${output_dir}genes_with_tags.txt

#join genes with contigs and diff_genes_ID by gene ID, and keep the lineInSam as first col ; this give cols lineInSam & gene_is_diff
LANG=en_EN join -t $'\t' -a1 -e'F' -12 -21 -o 1.1,2.1 <(LANG=en_EN sort -k2,2 ${output_dir}genes_with_tags.txt) <(LANG=en_EN sort -k2,2 ${output_dir}diff_genes_ID.txt) | awk 'OFS="\t"{if($2!="F"){$2="T"};print}' >${output_dir}diff_genes_with_tags.txt && rm ${output_dir}genes_with_tags.txt ${output_dir}diff_genes_ID.txt

#join the final table with diff_genes_ID
#! sensitive column ! : in final table we have to search gene_is_diff (actually it's col 31)
LANG=en_EN join -t $'\t' -a1 -e'F' -11 -21 <(LANG=en_EN sort -k1,1 $FinalTable) <(LANG=en_EN sort -k1,1 ${output_dir}diff_genes_with_tags.txt) | awk 'OFS="\t"{if($31==""){$31="F"}print}' >${FinalTable}.tmp && mv ${FinalTable}.tmp $FinalTable && rm ${output_dir}diff_genes_with_tags.txt

end_date=$(date)

echo -e "\nstart searching of contigs environment  : $start_date\n"
echo -e "\nend searching of contigs environment : $end_date\n"

#add the appropriates NAs (in addition of the existing lineInSam and tag ID column) in the final table for the unmatched contigs
#2 (LineInSam + ID) +29 NAs = 31 columns

if [ -f ${output_dir}unmapped_tags_2.txt ];then


	if [[ $(wc -l ${output_dir}unmapped_tags_2.txt |awk '{print $1}') -gt 0 ]];then

		while read line ;do 

		 empty=($(printf "%0.sNA\t" {1..29} |sed 's/\t$//g'))
		 
		 #is_mapped=F
		 empty[0]="F"
		 empty=$(echo ${empty[*]} |sed 's/ /\t/g')
		 
		 #line contain line in sam and ID
		 echo -e "${line}\t${empty}" >>$FinalTable
		   
		done < ${output_dir}unmapped_tags_2.txt
	
	fi

fi


#join the output table & the diff table by ID
#the ID in the final table is the 2nd column ; in diffex file, it's the 1st
#be carefull, there is a switch between col 1 and 2 after the join !! now the first column is the ID, and the second the lineInSam
LANG=en_EN join -t $'\t' -12 -21 <(LANG=en_EN sort -k2,2 $FinalTable) <(awk 'NR>1{print}' ${output_dir}diffexFile.txt | LANG=en_EN sort -k1,1) >${FinalTable}.tmp && mv ${FinalTable}.tmp $FinalTable

diff_ex_header=$(awk 'NR==1{print;exit}' ${output_dir}diffexFile.txt|cut -f 2-)
rm ${output_dir}diffexFile.txt

#print the complete header in a file
echo -e "ID\tLineInSam\tis_mapped\taligner\tchromosome\tstart\tend\tjunctions\tnb_junction\texon_coord_list\tother_split\tstrand\tnb_hit\tnb_mismatch\tnb_deletion\tnb_insertion\tSNV\tclipped_5p\tclipped_3p\tas_gene\tHUGO_ID_as_gene\tgene\tHUGO_ID_gene\tgene_5p\tgene_5p_dist\tgene_3p\tgene_3p_dist\tUTR\texonic\tintronic\tgene_is_diff\t${diff_ex_header}" >${output_dir}header.txt

#add the header in the final table
cat ${output_dir}header.txt $FinalTable >${FinalTable}.tmp && mv ${FinalTable}.tmp $FinalTable && rm ${output_dir}header.txt 


#remove the additional "tag" column
cut -f 1-33,35- $FinalTable >${FinalTable}.tmp && mv ${FinalTable}.tmp $FinalTable
sed 's/assembly/contig/g' $FinalTable >${FinalTable}.tmp && mv ${FinalTable}.tmp $FinalTable
echo -e "\n==== Clustering of the contigs by loci (genic/antisense/intergenic) ====\n"

#cluster contigs per loci
start_date=$(date)	
$getContigsPerLoci $output_dir $FinalTable $threads || { echo "R script failure !" 1>&2; exit; }
end_date=$(date)

echo -e "\nstart clustering of contigs : $start_date\n"
echo -e "\nend clustering of contigs : $end_date\n"

rm ${output_dir}ref_annotation.tmp ${output_dir}OriginalFastaTags.tmp
