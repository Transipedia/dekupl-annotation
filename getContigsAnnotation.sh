#!/bin/bash

usage() { echo -e "Usage: $0 <Required arguments> [Optional arguments]\n\n
 \tRequired arguments :\n
                  -a <merged-diff-counts.tsv.gz (contigs in \"{A}_vs_{B}_kmer_counts\" directory from Dekupl-run result)>\n
                  -g <genome in fasta (uncompressed). Use the same version as the annotation to avoid chromosome name issues>\n
                  -d <{A}vs{B}-DEGs.tsv (diff. genes in \"gene_expression\" directory from Dekupl-run result)>\n
                  -r <reference annotation (gff3 format, uncompressed). Use the same version as the genome file to avoid chromosome name issues>\n
                  -t <are the reads stranded ? (choose \"yes\" or \"no\")>\n
                  -e <normalized_counts.tsv (normalized gene counts in \"gene_expression\" directory from Dekupl-run result)>\n
                  -f <sample_conditions.tsv (design file in \"metadata\" directory from Dekupl-run result)>\n
                  -o <path to output directory>\n
                  -i <illumina adapters (you can use the file \"adapters.fa\" supplied with the program)>\n\n
\tOptional arguments :\n
                  -b <path to bin/ of blast scripts (default : in \$PATH environment variable)>\n
                  -c <path to bedtools, preferentially 2.24 (default : in \$PATH environment variable)>\n
                  -j <GSNAP genome index name (if the index of the genome has already been created, supply its name. Otherwise, the index is re-created with the name \"genome_index\")>\n
                  -k <path to directory of GSNAP genome index (if you have a former run, you can supply the full path of the \"temp_dir/mapping_output/\" inside your former output directory, in order to re-use the same genome index and save time)>\n
                  -m <path to bin/ of GSNAP (default : in \$PATH environment variable)>\n
                  -p <padj diff. gene threshold (default : 0.05)>\n
                  -q <contig color (choose 1 or 2 ; default : 1) > \n\n\t\t\t1 : contigs on forward strand are in red (contigs on reverse strand are in blue)\n\t\t\t2 : contigs on forward strand are in blue (contigs on reverse strand are in red)\n
                  -s <path to samtools (default : in \$PATH environment variable)>\n
                  -n <thread number (default : 1)>\n\n
\tResults :\n
                  - Table \"DiffContigsInfos.tsv\", summarizing for each contig, its location on the genome (if it's aligned), the neighborhood, the sequence alignment informations, and the differential expression informations\n
                  - BED file \"diff_contigs.bed\" for the visualization ; it contains useful informations from the summarization table.\n
                  
                  - Table \"contigsPerLoci.tsv\" (only for stranded data for now) containing loci with differentially expressed contigs \n" 1>&2; exit 1;}

[[ $# -eq 0 ]] && usage

while getopts ":a:g:d:r:i:o:b:c:j:k:m:p:s:t:n:e:f:q:" opt; do
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
            
      e)
      
            normalized_gene_counts=$OPTARG

            ;;            
      f)
      
            design=$OPTARG

            ;;
            
      q)
      
            contig_color=$OPTARG
          
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
chmod 755 $getContigsPerLoci

getSwitches=$(dirname "$0")/getSwitches.R
chmod 755 $getSwitches

getIntronsByTranscripts=$(dirname "$0")/getIntronsByTranscripts.R
chmod 755 $getIntronsByTranscripts

blast=$(dirname "$0")/blast.sh
chmod 755 $blast

chmod 755 $(dirname "$0")/getAnnotationFunctions.sh
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

#if one of the useful arguments is missing, exit with the usage 
if [ ! -f "$DEkupl_result" ] || [ ! -f "$ref_fasta" ] || [ ! -f "$diff_genes" ] || [ ! -f "$ref_annotation" ] || [ "$output_dir" == "" ] || [ ! -f "$illumina_adapters" ] || [ "$samtools" == "" ] || [ "$bedtools" == "" ] || [ "$ncbi_blast_loc" == "" ] || [ "$GSNAP_loc" == "" ] || ([ "$stranded" != "yes" ] && [ "$stranded" != "no" ]) || [ ! -f "$normalized_gene_counts" ] || [ ! -f "$design" ]; then
      
echo -e "\none required argument is missing or is wrong (check also if all the programs are installed) !!\n"
      usage
      exit 1
fi


echo -e "stranded contigs : $stranded\n"


if [ "$threads" == "" ];then threads=1;fi
echo -e "threads number is : $threads\n"

if [ "$padj_threshold" == "" ];then padj_threshold=0.05;fi

#set contig color
if [ "$contig_color" != "1" ] && [ "$contig_color" != "2" ];then contig_color=1;fi

if [ "$contig_color" == "1" ];then

  echo -e "\ncontig color is set to 1 : contigs on forward strand are in red (contigs on reverse strand are in blue)\n" >&2
  
  forward_contig_color="255,0,0"
  
  reverse_contig_color="0,0,255"
 
elif [ "$contig_color" == "2" ];then

  echo -e "\ncontig color is set to 2 : contigs on forward strand are in blue (contigs on reverse strand are in red)\n"
  
  forward_contig_color="0,0,255"
  
  reverse_contig_color="255,0,0"

fi


##### processing of the input variables
########################################

#output dir
output_dir="${output_dir}/"
output_dir=$(echo "$output_dir" | sed 's/\/\//\//g')
if [ ! -d $output_dir ];then mkdir $output_dir || { echo "cannot create $output_dir !!" 1>&2; exit; } ;fi

#temp dir (use to store temp files after the annotation, for checkings)
temp_dir="${output_dir}/temp_files/"
if [ ! -d $temp_dir ];then mkdir $temp_dir ;fi

mapping_output="${temp_dir}mapping_output/"
if [ ! -d $mapping_output ];then mkdir $mapping_output ;fi

#table summarizing for each contig, its location on the genome (if it's aligned), the neighborhood, the sequence alignment informations, and the differential expression informations
FinalTable="${output_dir}DiffContigsInfos.tsv"
if [ -f $FinalTable ];then rm $FinalTable ;fi

#contigs in BED12 format
diff_contigs_bed="${output_dir}diff_contigs.bed"
if [ -f $diff_contigs_bed ];then rm $diff_contigs_bed ;fi

#table containing loci (sense, antisense, intergenic, unmapped) with differentially expressed contigs
if [ -f ${output_dir}ContigsPerLoci.tsv ];then rm ${output_dir}ContigsPerLoci.tsv ;fi


#processing of the gff : construct the feature "intron" + sort
LANG=en_EN sort -k1,1 -k4,4n <(grep -v "^#" $ref_annotation) <($getIntronsByTranscripts $ref_annotation)>${temp_dir}ref_annotation.tmp

ref_annotation=${temp_dir}ref_annotation.tmp

zcat $DEkupl_result | awk 'NR==1{OFS="\t";print "ID",$0;exit}' >${temp_dir}diffexFileHeader.txt

#convert the DEkupl contigs output to fasta, and at the same time put the tag of reference as an ID for joinings
zcat $DEkupl_result | awk 'NR>1{print}' | awk 'OFS="\t"{print $3,$0}' | tee ${temp_dir}diffexFile.txt | awk 'OFS="\n"{print ">"$1,$3}' >${temp_dir}OriginalFastaContigs.fa 

#col 1 =ID ; col 2 = contig seq : will be used in joinings
cat ${temp_dir}OriginalFastaContigs.fa | paste - - | awk 'OFS="\t"{print $1,$2}' | LANG=en_EN sort -k1,1 >${temp_dir}OriginalFastaContigs.tmp

cat ${temp_dir}diffexFileHeader.txt ${temp_dir}diffexFile.txt >${temp_dir}diffexFile.tmp && mv ${temp_dir}diffexFile.tmp ${temp_dir}diffexFile.txt && rm ${temp_dir}diffexFileHeader.txt


##### blast of contigs against adapters, and retain only those not matching
########################################################################

echo -e "\n==== Filter out contigs matching in adapters ====\n"

echo -e "\nnumber of initial contigs $(($(wc -l ${temp_dir}OriginalFastaContigs.fa|awk '{print $1}')/2))\n"

start_date=$(date) 

bash $blast -q ${temp_dir}OriginalFastaContigs.fa -s $illumina_adapters -o $temp_dir -p $ncbi_blast_loc -d "illumina_adapters" -n $threads || { echo "blast on adapters failure!!" 1>&2; exit; }

#for checkings
awk '!seen[$1]++' "${temp_dir}raw_blast.alignment_2.txt" | sed $'s/ /\t/g' >${temp_dir}match_in_adapters.txt

grep "^>" ${temp_dir}OriginalFastaContigs.fa | sed 's/>//g' | LANG=en_EN sort >${temp_dir}originaltagIDs.txt

comm -23 ${temp_dir}originaltagIDs.txt <(cut -f 1 ${temp_dir}raw_blast.alignment_2.txt | LANG=en_EN sort -u) | awk 'OFS="\t"{print ">"$1}' | LANG=en_EN sort -k1,1 >${temp_dir}nomatch_in_adapters.txt && rm ${temp_dir}originaltagIDs.txt

#reconstruction of the fasta with only contigs not matching in adapters
LANG=en_EN join -t $'\t' -11 -21 ${temp_dir}nomatch_in_adapters.txt ${temp_dir}OriginalFastaContigs.tmp | sed 's/\t/\n/g' >${temp_dir}nomatch_in_adapters.tmp && rm ${temp_dir}nomatch_in_adapters.txt && mv ${temp_dir}nomatch_in_adapters.tmp ${temp_dir}nomatch_in_adapters.fa

end_date=$(date)

echo -e "\nstart blast of contigs on adapters : $start_date\n"
echo -e "\nend blast of contigs on adapters : $end_date\n"

#check if there are contigs to map
if [[ $(wc -l ${temp_dir}nomatch_in_adapters.fa|awk '{print $1}') -eq 0 ]];then

  echo -e "\nno contigs to map on the supplied genome, check your inputs !\n"
  
  exit
  
fi


##### mapping of the contigs
###############################
 
#remove pre-existing mapping file 
if [ -f ${mapping_output}contigs.bam ];then

  rm ${mapping_output}contigs.bam
  
fi

echo -e "\n==== Mapping of contigs on the genome ====\n"

echo -e "\nnumber of contigs to align : $(($(wc -l ${temp_dir}nomatch_in_adapters.fa|awk '{print $1}')/2))\n"

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

${GSNAP_loc}gsnap -t $threads -A sam -N 1 -D $gsnap_index_dir -d $gsnap_index_name -w 50000 ${temp_dir}nomatch_in_adapters.fa |$samtools view -bh >${mapping_output}contigs.bam || { echo "gsnap mapping failure!!" 1>&2; exit; }

end_date=$(date)

#check if there are no problems with the output bam file
if [[ $($samtools view ${mapping_output}contigs.bam |wc -l |awk '{print $1}') -eq 0 ]];then

  echo -e "\nno results for the mapping, check the supplied index directory !\n"
  
  exit


fi

echo -e "\nstart mapping of contigs : $start_date\n"
echo -e "\nend mapping of contigs : $end_date\n"


################## for each alignment line, reconstruction of the tag structure (splices, exons...)
###################################################################################################

echo -e "\n==== Reconstruction of the tag structure (splices, exons...) from the alignment ====\n"

start_date=$(date)

#from the alignment, get a bed and a table of the contigs with infos about the alignment
parseBam ${mapping_output}contigs.bam $FinalTable $diff_contigs_bed $forward_contig_color $reverse_contig_color

$samtools view -f 4 ${mapping_output}contigs.bam | LANG=en_EN sort -k1,1 | LANG=en_EN sort -u -k1,1 | awk -v start_line=$(($(wc -l $diff_contigs_bed | awk '{print $1}')+1)) 'BEGIN{a=start_line}OFS="\t"{print a,$1;a=a+1}' >${temp_dir}OriginalUnmappedContigs.txt

end_date=$(date)

echo -e "\nstart parse BAM : $start_date\n"
echo -e "\nend parse BAM : $end_date\n"


############## blast of unmapped seq and add the good ones (full length alignment) in the final table/bed
#########################################################################################################

#if there are unmapped contigs, try to align them on the supplied genome with another tool (blast here)
if [ -f ${temp_dir}OriginalUnmappedContigs.txt ];then

start_date=$(date)


 #if there are unmapped contigs, we'll try to align them with blast (actually we are looking for those with full length)
 if [[ $(wc -l ${temp_dir}OriginalUnmappedContigs.txt|awk '{print $1}') -gt 0 ]];then

	echo -e "\n==== blast of unmapped contigs : retrieving those mapping with their full length ====\n"

        #this function will rebuild a fasta file from the unmapped contigs
	getFastaFromUnmappedContigs ${temp_dir}OriginalUnmappedContigs.txt ${temp_dir}OriginalFastaContigs.tmp ${temp_dir}unmapped_fasta1.fa 

	bash $blast -q ${temp_dir}unmapped_fasta1.fa -s $ref_fasta -o $temp_dir -p $ncbi_blast_loc -d "human_genome" -n $threads

	#if alignment length = tag length, it's a good case missed by the mapper, we make a unique on the tag ID !
	#print tag ID,chr,start,end,strand,nb_mismatch,gaps , sort by mismaches, then gaps, then ID, then chromosome
	awk 'OFS="\t"{if($4==$13){if($14=="plus"){$14="+"}else{$14="-"};print $1,$2,$9,$10,$14,$5,$6}}' ${temp_dir}raw_blast.alignment_2.txt | LANG=en_EN sort -k6,6n -k7,7n -k1,1 -k2,2 | awk '!seen[$1]++'>${temp_dir}UnmappedEntirelyFoundByBlast.txt

	awk 'OFS="\t"{if($4==$13){if($14=="plus"){$14="+"}else{$14="-"};print $1,$2,$9,$10,$14,$5,$6}}' ${temp_dir}raw_blast.alignment_2.txt >${temp_dir}UnmappedEntirelyFoundByBlast.tmp1


	if [ -f ${temp_dir}UnmappedEntirelyFoundByBlast.tmp2 ];then rm ${temp_dir}UnmappedEntirelyFoundByBlast.tmp2 ; fi

        #re-insert the found alignments in the final table, after computing
	while read line ;do

	  mismaches=$(echo "$line" | cut -f 6)
	  gaps=$(echo "$line" | cut -f 7)
	  ID=$(echo "$line" | cut -f 1)
	  
	  #if there's an alignment as good as the best, it's a multihit
	  multihit=$(grep -P "$ID\t" ${temp_dir}UnmappedEntirelyFoundByBlast.tmp1 |awk -v nb_mismatch=$mismaches -v gaps=$gaps 'BEGIN{a=0}{if($6==nb_mismatch && $7==gaps){a=a+1}}END{print a}')
	  
	  #1st block : print tag ID,lineInSAM,chr,start,end,strand,nb_mismatch,gaps,nb_hit
	  #2nd block : reorganize the columns in order to add them to the final table, SNV is added just after ins : line_in_SAM,tagID,is_mapped,mapper,chr,start,end,junction,nb_junction,exon_coord,other_split,strand,nb_hit,mis,del,ins,SNV,clipped_5p,clipped_3p
	  echo -e "$line\t$multihit" | awk 'OFS="\t"{sub("-","\t",$1);print}' | awk 'OFS="\t"{if($7>0 || $8>0){SNV="T"}else{SNV="F"};print $2,$1,"T","Blast",$3,$4,$5,"none",0,$3":"$4"-"$5,"F",$6,$9,$7,$8,0,SNV,0,0}' | tee -a ${temp_dir}UnmappedEntirelyFoundByBlast.tmp2 | awk 'OFS="\t"{print}' >>$FinalTable

	done < ${temp_dir}UnmappedEntirelyFoundByBlast.txt
	
	rm ${temp_dir}UnmappedEntirelyFoundByBlast.tmp1 ${temp_dir}UnmappedEntirelyFoundByBlast.txt

	if [ -f ${temp_dir}UnmappedEntirelyFoundByBlast.tmp2 ];then

	  #put the result in the bed file 
	  awk -v forward_contig_color=$forward_contig_color -v reverse_contig_color=$reverse_contig_color 'OFS="\t"{if($7>$6){$6=$6-1;print $5,$6,$7,"LineInSam="$1";ID="$2";nb_hit="$13";nM="$14";del="$15";ins=0;clipped_5p=0;clipped_3p=0",1,$12,$6,$7,forward_contig_color,1,$7-$6,0}else{$7=$7-1;print $5,$7,$6,"LineInSam="$1";ID="$2";nb_hit="$13";nM="$14";del="$15";ins=0;clipped_5p=0;clipped_3p=0",1,$12,$7,$6,reverse_contig_color,1,$6-$7,0}}' ${temp_dir}UnmappedEntirelyFoundByBlast.tmp2 >>$diff_contigs_bed

          #contigs still unmapped after the previous method
	  comm -23 <(cut -f 1 ${temp_dir}OriginalUnmappedContigs.txt | LANG=en_EN sort) <(cut -f 1 ${temp_dir}UnmappedEntirelyFoundByBlast.tmp2|LANG=en_EN sort ) | LANG=en_EN join -t $'\t' -11 -21 - <(LANG=en_EN sort -k1,1 ${temp_dir}OriginalUnmappedContigs.txt) >${temp_dir}unmapped_contigs_2.txt && rm ${temp_dir}UnmappedEntirelyFoundByBlast.tmp2
	
	else

          #if the previous method didn't give results, just re-use the initial unmapped contigs
          cat ${temp_dir}OriginalUnmappedContigs.txt >${temp_dir}unmapped_contigs_2.txt

        fi


  fi

end_date=$(date)

echo -e "\nstart blast of unmapped contigs on genome : $start_date\n"
echo -e "\nend blast of unmapped contigs on genome : $end_date\n"

fi

LANG=en_EN sort -k1,1 -k2,2n $diff_contigs_bed >${diff_contigs_bed}.tmp && mv ${diff_contigs_bed}.tmp $diff_contigs_bed

#if there are contigs which are found in the supplied genome, process them
if [[ $(wc -l $diff_contigs_bed|awk '{print $1}') -gt 0 ]];then

        echo -e "\ncontigs are found in the supplied genome. We're going to process them\n"


	########## insert diffex infos in the bed file 
	###############################################

	LANG=en_EN join -t $'\t' -11 -21 <(awk 'OFS="\t"{initial_line=$0;ID_col=$4;print ID_col,initial_line}' $diff_contigs_bed | awk 'OFS="\t"{split ($1,x,";");a=x[2]; print a,$0}' | sed 's/^ID=//g' | awk 'OFS="\t"{$2="";print $0}'|LANG=en_EN sort -k1,1) <(awk 'NR>1{printf $1"\t";printf "%.3e\t",$5;printf "%.0f\t",$6;printf "%.0f\t",$7;printf "%.2f\n",$8}' ${temp_dir}diffexFile.txt|LANG=en_EN sort -k1,1) | awk 'OFS="\t"{$5=$5";pval="$14";meanA="$15";meanB="$16";log2FC="$17;print}' | cut -f2-13 >${diff_contigs_bed}.tmp1 && mv ${diff_contigs_bed}.tmp1 ${diff_contigs_bed}

	LANG=en_EN sort -k1,1 -k2,2n $diff_contigs_bed >${diff_contigs_bed}.tmp && mv ${diff_contigs_bed}.tmp $diff_contigs_bed


	########## search max and min abs(log2FC) to scale the bed color (RGB) on them
	#########################################################################

	#max value for upregulated features
	max_up_log2FC=$(zcat $DEkupl_result | awk 'NR>1{if($7 >0 || $7 <0){print $7}}'| awk '{if($1!~/-/){print}}' | LANG=en_EN sort -nr | head -n1)

	#max value for downregulated features
	max_down_log2FC=$(zcat $DEkupl_result | awk 'NR>1{if($7 >0 || $7 <0){print $7}}'| awk '{if($1~/-/){print}}' | sed 's/-//g'| LANG=en_EN sort -nr | head -n1)

	#max absolute value
	max_abs_log2FC=$(echo -e "${max_up_log2FC}\n${max_down_log2FC}"| LANG=en_EN sort -nr |head -n1)

	#min value for upregulated features
	min_up_log2FC=$(zcat $DEkupl_result | awk 'NR>1{if($7 >0 || $7 <0){print $7}}' | awk '{if($1!~/-/){print}}' | LANG=en_EN sort -n | head -n1)

	#min value for downregulated features
	min_down_log2FC=$(zcat $DEkupl_result |awk 'NR>1{if($7 >0 || $7 <0){print $7}}' | awk '{if($1~/-/){print}}' | sed 's/-//g'| LANG=en_EN sort -n | head -n1)

	#min absolute value
	min_abs_log2FC=$(echo -e "${min_up_log2FC}\n${min_down_log2FC}" | LANG=en_EN sort -n | head -n1)

	#255,220,220 & 220,220,255 : light red and light blue
	#255,0,0 & 0,0,255 : far red and far blue
	#we will scale the color following the abs(log2FC)
	#formula to range [min,max] to [a,b] :[((b-a)(x - min))/ (max - min)] +a

	#join bed and diffex table (just column ID and log2FC) by ID, then scale the color of each line on the abs(log2FC)
	LANG=en_EN join -t $'\t' -11 -21 <(awk 'OFS="\t"{initial_line=$0;ID_col=$4;print ID_col,initial_line}' $diff_contigs_bed | awk 'OFS="\t"{split ($1,x,";");a=x[2]; print a,$0}' | sed 's/^ID=//g' | awk 'OFS="\t"{$2="";print $0}'|LANG=en_EN sort -k1,1) <(awk 'NR>1{OFS="\t";print $1,$8}' ${temp_dir}diffexFile.txt) | \
	 awk -v contig_color=$contig_color -v b=220 -v a=0 -v max=$max_abs_log2FC -v min=$min_abs_log2FC 'OFS="\t"{if($NF~/-/){abs_log2FC=-$14}else{abs_log2FC=$NF};RGB=220-((((b-a)*((0.5^(abs_log2FC))-(0.5^(min))))/((0.5^(max))-(0.5^(min))))+a);if($7=="+"){if(contig_color==1){printf $0"\t";printf "255,";printf "%.0f",RGB;printf ",";printf "%.0f\n",RGB}else{printf $0"\t";printf "%.0f",RGB;printf ",";printf "%.0f",RGB;printf ",255\n"}};if($7=="-"){if(contig_color==1){printf $0"\t";printf "%.0f",RGB;printf ",";printf "%.0f",RGB;printf ",255\n"}else{printf $0"\t";printf "255,";printf "%.0f",RGB;printf ",";printf "%.0f\n",RGB}}}' | awk 'OFS="\t"{$10=$NF;print}' | cut -f 2-13 >${diff_contigs_bed}.tmp1 && mv ${diff_contigs_bed}.tmp1 ${diff_contigs_bed}

	LANG=en_EN sort -k1,1 -k2,2n $diff_contigs_bed >${diff_contigs_bed}.tmp && mv ${diff_contigs_bed}.tmp $diff_contigs_bed


	########## research of antisense & neighbours + location of the contigs in the gene (UTR, exon...) ###########
	###########################################################################################################

	echo -e "\n==== research of the neighborhood of the contigs ====\n"

	#genome size for bedtools 
	$samtools view -H ${mapping_output}contigs.bam | grep -E -v "^@PG|^@HD" | awk 'OFS="\t"{print $2,$3}' | sed 's/SN\://;s/LN\://g' | LANG=en_EN sort -k1,1 >${mapping_output}genome_size.txt

	start_date=$(date)

	#searching antisense gene, retain only one target (sort by chr and start, then unique on attribute of the tag) (gives 2 cols in addition of the lineInSAM : Ensembl gene ID & HUGO ID)
	if [ "$stranded" == "yes" ];then orientation_option="-s" ;dist_option="a";else orientation_option="";dist_option="ref"; fi

	awk '{if($3=="gene"){print}}' $ref_annotation| $bedtools intersect -S -a $diff_contigs_bed -b - -loj -nonamecheck | LANG=en_EN sort -k1,1 -k2,2n | LANG=en_EN sort -u -k4,4 >${temp_dir}antisense_contigs.txt || { echo "searching of antisense gene failure (bedtools instersect -S )!!" 1>&2; exit; }


	if [ "$stranded" == "yes" ];then 

		#remark : for the last blocks of commands (in the 2nd argument of the "paste"), we use the table a[] to store the gff attributes, and we keep the result only if we have the regex "gene_name" (case insensitive)
		paste -d'\t' <(awk 'OFS="\t"{print $4}' ${temp_dir}antisense_contigs.txt|awk -F';' '{print $1}' | awk '{sub("LineInSam=","",$1);print}') \
		<(awk 'OFS="\t"{print $NF}' ${temp_dir}antisense_contigs.txt | awk -F';' 'OFS="\t"{print $1,$0}' | awk 'OFS="\t"{if($1=="."){$1="none";$2="none"};print}' | awk 'OFS="\t"{if($1=="none"){print}else{split($2,a,";");IGNORECASE=1;for(i=1;i<=length(a);i++){if(a[i]~/^gene_name/||a[i]~/^Name/){b=a[i]}};print $1,b}}' | awk 'BEGIN{IGNORECASE=1}OFS="\t"{sub("ID=","");sub("gene_name=","");print}') >${temp_dir}antisense_contigs.tmp && mv ${temp_dir}antisense_contigs.tmp ${temp_dir}antisense_contigs.txt

	else

		paste -d'\t' <(awk 'OFS="\t"{print $4}' ${temp_dir}antisense_contigs.txt | awk -F';' '{print $1}' | awk '{sub("LineInSam=","",$1);print}') <(awk 'OFS="\t"{print $NF}' ${temp_dir}antisense_contigs.txt | awk -F';' 'OFS="\t"{print "none","none"}') >${temp_dir}antisense_contigs.tmp && mv ${temp_dir}antisense_contigs.tmp ${temp_dir}antisense_contigs.txt

	fi

	#normally, if the bed file isn't empty, you shouldn't have an empty file (because it's a "left outer join") ; if it's the case, probably there's no concordance between the genome & the annotation !
	if [[ $(wc -l ${temp_dir}antisense_contigs.txt|awk '{print $1}') -eq 0 ]];then echo -e "\ncheck if your genome & annotation files have concordant chromosomes !\n";exit;fi

	LANG=en_EN join -t $'\t' -11 -21 <(LANG=en_EN sort -k1,1 $FinalTable) <(LANG=en_EN sort -k1,1 ${temp_dir}antisense_contigs.txt) >${FinalTable}.tmp && mv ${FinalTable}.tmp $FinalTable 


	#searching sense gene (gives 2 cols in addition of the lineInSAM : Ensembl gene ID & HUGO ID)
	awk '{if($3=="gene"){print}}' $ref_annotation | $bedtools intersect $orientation_option -a $diff_contigs_bed -b - -loj -nonamecheck | LANG=en_EN sort -k1,1 -k2,2n | LANG=en_EN sort -u -k4,4 >${temp_dir}sense_contigs.txt || { echo "searching of sense gene failure (bedtools instersect -s )!!" 1>&2; exit; }

	#remark : for the last blocks of commands (in the 2nd argument of the "paste"), we use the table a[] to store the gff attributes, and we keep the result only if we have the regex "gene_name" (case insensitive)
	paste -d'\t' <(awk 'OFS="\t"{print $4}' ${temp_dir}sense_contigs.txt|awk -F';' '{print $1}'|awk '{sub("LineInSam=","",$1);print}') \
	<(awk 'OFS="\t"{print $NF}' ${temp_dir}sense_contigs.txt|awk -F';' 'OFS="\t"{print $1,$0}' | awk 'OFS="\t"{if($1=="."){$1="none";$2="none"};print}' | awk 'OFS="\t"{if($1=="none"){print}else{split($2,a,";");IGNORECASE=1;for(i=1;i<=length(a);i++){if(a[i]~/^gene_name/||a[i]~/^Name/){b=a[i]}};print $1,b}}' | awk 'BEGIN{IGNORECASE=1}OFS="\t"{sub("ID=","");sub("gene_name=","");print}') >${temp_dir}sense_contigs.tmp && mv ${temp_dir}sense_contigs.tmp ${temp_dir}sense_contigs.txt

	LANG=en_EN join -t $'\t' -11 -21 <(LANG=en_EN sort -k1,1 $FinalTable) <(LANG=en_EN sort -k1,1 ${temp_dir}sense_contigs.txt) >${FinalTable}.tmp && mv ${FinalTable}.tmp $FinalTable


	#closest 5'end gene (gives 2 cols in addition of the lineInSAM : Ensembl gene ID & dist)
	awk '{if($3=="gene"){print}}' $ref_annotation | $bedtools closest -nonamecheck -g ${mapping_output}genome_size.txt $orientation_option -io -fu -D $dist_option -t first -a $diff_contigs_bed -b - | LANG=en_EN sort -k1,1 -k2,2n | LANG=en_EN sort -u -k4,4 >${temp_dir}closest_5end.txt || { echo "searching of closest 5' gene failure (bedtools closest -io -fu )!!" 1>&2; exit; }

	paste -d'\t' <(awk 'OFS="\t"{print $4}' ${temp_dir}closest_5end.txt | awk -F';' '{print $1}'|awk '{sub("LineInSam=","",$1);print}') <(awk 'OFS="\t"{print $(NF-1)}' ${temp_dir}closest_5end.txt | awk -F';' '{print $1}' | awk '{if($1=="."){$1="none"};sub("ID=","",$1);print}') <(awk 'OFS="\t"{print $(NF-1),$NF}' ${temp_dir}closest_5end.txt | awk '{if($1=="."){$2="none"};print $2}') | awk 'OFS="\t"{if($3>0){$2="none";$3="none"}print}'>${temp_dir}closest_5end.tmp && mv ${temp_dir}closest_5end.tmp ${temp_dir}closest_5end.txt

	LANG=en_EN join -t $'\t' -11 -21 <(LANG=en_EN sort -k1,1 $FinalTable) <(LANG=en_EN sort -k1,1 ${temp_dir}closest_5end.txt) >${FinalTable}.tmp && mv ${FinalTable}.tmp $FinalTable && rm ${temp_dir}closest_5end.txt


	#closest 3'end gene (gives 2 cols in addition of the lineInSAM : Ensembl gene ID & dist)
	awk '{if($3=="gene"){print}}' $ref_annotation| $bedtools closest -nonamecheck -g ${mapping_output}genome_size.txt $orientation_option -io -fd -D $dist_option -t first -a $diff_contigs_bed -b - | LANG=en_EN sort -k1,1 -k2,2n | LANG=en_EN sort -u -k4,4 >${temp_dir}closest_3end.txt || { echo "searching of closest 3' gene failure (bedtools closest -io -fd )!!" 1>&2; exit; }

	paste -d'\t' <(awk 'OFS="\t"{print $4}' ${temp_dir}closest_3end.txt | awk -F';' '{print $1}' | awk '{sub("LineInSam=","",$1);print}') <(awk 'OFS="\t"{print $(NF-1)}' ${temp_dir}closest_3end.txt | awk -F';' '{print $1}' | awk '{if($1=="."){$1="none"};sub("ID=","",$1);print}') <(awk 'OFS="\t"{print $(NF-1),$NF}' ${temp_dir}closest_3end.txt |awk '{if($1=="."){$2="none"};print $2}') | awk 'OFS="\t"{if($3<0){$2="none";$3="none"}print}'>${temp_dir}closest_3end.tmp && mv ${temp_dir}closest_3end.tmp ${temp_dir}closest_3end.txt

	LANG=en_EN join -t $'\t' -11 -21 <(LANG=en_EN sort -k1,1 $FinalTable) <(LANG=en_EN sort -k1,1 ${temp_dir}closest_3end.txt) >${FinalTable}.tmp && mv ${FinalTable}.tmp $FinalTable && rm ${temp_dir}closest_3end.txt

	#searching UTR contigs (gives 1 col in addition of the lineInSAM : T or F)
	awk '$3 ~ "UTR"{print}' $ref_annotation| $bedtools intersect $orientation_option -a $diff_contigs_bed -b - -loj -nonamecheck | LANG=en_EN sort -k1,1 -k2,2n | LANG=en_EN sort -u -k4,4 >${temp_dir}UTR_contigs.txt || { echo "searching of contigs in UTRs failure (bedtools intersect )!!" 1>&2; exit; }

	paste -d'\t' <(awk 'OFS="\t"{print $4}' ${temp_dir}UTR_contigs.txt | awk -F';' '{print $1}' | awk '{sub("LineInSam=","",$1);print}') <(awk 'OFS="\t"{print $NF}' ${temp_dir}UTR_contigs.txt | awk -F';' '{if($1=="."){$1="F";print}else{print "T"}}') >${temp_dir}UTR_contigs.tmp && mv ${temp_dir}UTR_contigs.tmp ${temp_dir}UTR_contigs.txt

	LANG=en_EN join -t $'\t' -11 -21 <(LANG=en_EN sort -k1,1 $FinalTable) <(LANG=en_EN sort -k1,1 ${temp_dir}UTR_contigs.txt) >${FinalTable}.tmp && mv ${FinalTable}.tmp $FinalTable 


	#searching exonic contigs ( gives 1 col in addition of the line_in_SAM : T or F)
	awk '{if($3=="exon"){print}}' $ref_annotation | $bedtools intersect $orientation_option -a $diff_contigs_bed -b - -loj -nonamecheck | LANG=en_EN sort -k1,1 -k2,2n | LANG=en_EN sort -u -k4,4 >${temp_dir}Exonic_contigs.txt || { echo "searching of contigs in exons failure (bedtools intersect )!!" 1>&2; exit; }

	paste -d'\t' <(awk 'OFS="\t"{print $4}' ${temp_dir}Exonic_contigs.txt | awk -F';' '{print $1}' | awk '{sub("LineInSam=","",$1);print}') <(awk 'OFS="\t"{print $NF}' ${temp_dir}Exonic_contigs.txt | awk -F';' '{if($1=="."){$1="F";print}else{print "T"}}') >${temp_dir}Exonic_contigs.tmp && mv ${temp_dir}Exonic_contigs.tmp ${temp_dir}Exonic_contigs.txt

	LANG=en_EN join -t $'\t' -11 -21 <(LANG=en_EN sort -k1,1 $FinalTable) <(LANG=en_EN sort -k1,1 ${temp_dir}Exonic_contigs.txt) >${FinalTable}.tmp && mv ${FinalTable}.tmp $FinalTable
	
	
	#searching intronic contigs ( gives 1 col in addition of the line_in_SAM : T or F)
	awk '{if($3=="intron"){print}}' $ref_annotation | $bedtools intersect $orientation_option -a $diff_contigs_bed -b - -loj -nonamecheck | LANG=en_EN sort -k1,1 -k2,2n | LANG=en_EN sort -u -k4,4 >${temp_dir}Intronic_contigs.txt || { echo "searching of contigs in introns failure (bedtools intersect )!!" 1>&2; exit; }

	paste -d'\t' <(awk 'OFS="\t"{print $4}' ${temp_dir}Intronic_contigs.txt | awk -F';' '{print $1}' | awk '{sub("LineInSam=","",$1);print}') <(awk 'OFS="\t"{print $NF}' ${temp_dir}Intronic_contigs.txt | awk -F';' '{if($1=="."){$1="F";print}else{print "T"}}') >${temp_dir}Intronic_contigs.tmp && mv ${temp_dir}Intronic_contigs.tmp ${temp_dir}Intronic_contigs.txt

	LANG=en_EN join -t $'\t' -11 -21 <(LANG=en_EN sort -k1,1 $FinalTable) <(LANG=en_EN sort -k1,1 ${temp_dir}Intronic_contigs.txt) >${FinalTable}.tmp && mv ${FinalTable}.tmp $FinalTable


	#diff gene (gives 1 col)
	awk -v padj_threshold=$padj_threshold 'NR>1{if($7!="NA" && $7<=padj_threshold){print $1}}' $diff_genes >${temp_dir}diff_genes_ID.txt

	# ! sensitive column ! : in final table we have to search sense gene (actually it's col 22); this give cols lineInSam & sense gene
	awk 'NR>1{OFS="\t";if($22!="none"){print $1,$22}}' $FinalTable >${temp_dir}genes_with_contigs.txt

	#join genes with contigs and diff_genes_ID by gene ID, and keep the lineInSam as first col ; this give cols lineInSam & gene_is_diff
	LANG=en_EN join -t $'\t' -a1 -e'F' -12 -21 -o 1.1,2.1 <(LANG=en_EN sort -k2,2 ${temp_dir}genes_with_contigs.txt) <(LANG=en_EN sort -k2,2 ${temp_dir}diff_genes_ID.txt) | awk 'OFS="\t"{if($2!="F"){$2="T"};print}' >${temp_dir}diff_genes_with_contigs.txt && rm ${temp_dir}genes_with_contigs.txt ${temp_dir}diff_genes_ID.txt

	#join the final table with diff_genes_ID
	#! sensitive column ! : in final table we have to search gene_is_diff (actually it's col 31)
	LANG=en_EN join -t $'\t' -a1 -e'F' -11 -21 <(LANG=en_EN sort -k1,1 $FinalTable) <(LANG=en_EN sort -k1,1 ${temp_dir}diff_genes_with_contigs.txt) | awk 'OFS="\t"{if($31==""){$31="F"}print}' >${FinalTable}.tmp && mv ${FinalTable}.tmp $FinalTable && rm ${temp_dir}diff_genes_with_contigs.txt

	end_date=$(date)

	echo -e "\nstart searching of contigs environment  : $start_date\n"
	echo -e "\nend searching of contigs environment : $end_date\n"

fi

#add the appropriates NAs (in addition of the existing lineInSam and tag ID column) in the final table for the unmatched contigs
#2 (LineInSam + ID) +29 NAs = 31 columns

if [ -f ${temp_dir}unmapped_contigs_2.txt ];then

	if [[ $(wc -l ${temp_dir}unmapped_contigs_2.txt |awk '{print $1}') -gt 0 ]];then

		while read line ;do 

		 empty=($(printf "%0.sNA\t" {1..29} |sed 's/\t$//g'))
		 
		 #is_mapped=F
		 empty[0]="F"
		 empty=$(echo ${empty[*]} |sed 's/ /\t/g')
		 
		 #line contain line in sam and ID
		 echo -e "${line}\t${empty}" >>$FinalTable
		   
		done < ${temp_dir}unmapped_contigs_2.txt
	
	fi

fi

#join the output table & the diff table by ID
#the ID in the final table is the 2nd column ; in diffex file, it's the 1st
#be carefull, there is a switch between col 1 and 2 after the join !! now the first column is the ID, and the second the lineInSam
LANG=en_EN join -t $'\t' -12 -21 <(LANG=en_EN sort -k2,2 $FinalTable) <(awk 'NR>1{print}' ${temp_dir}diffexFile.txt | LANG=en_EN sort -k1,1) >${FinalTable}.tmp && mv ${FinalTable}.tmp $FinalTable

diff_ex_header=$(awk 'NR==1{print;exit}' ${temp_dir}diffexFile.txt|cut -f 2-)
rm ${temp_dir}diffexFile.txt

#print the complete header in a file
echo -e "ID\tLineInSam\tis_mapped\taligner\tchromosome\tstart\tend\tjunctions\tnb_junction\texon_coord\tother_split\tstrand\tnb_hit\tnb_mismatch\tnb_deletion\tnb_insertion\tSNV\tclipped_5p\tclipped_3p\tas_gene\tHUGO_ID_as_gene\tgene\tHUGO_ID_gene\tgene_5p\tgene_5p_dist\tgene_3p\tgene_3p_dist\tUTR\texonic\tintronic\tgene_is_diff\t${diff_ex_header}" >${temp_dir}header.txt

#add the header in the final table
cat ${temp_dir}header.txt $FinalTable >${FinalTable}.tmp && mv ${FinalTable}.tmp $FinalTable && rm ${temp_dir}header.txt 

#remove the additional "tag" column
cut -f 1-33,35- $FinalTable >${FinalTable}.tmp && mv ${FinalTable}.tmp $FinalTable
sed 's/assembly/contig/g' $FinalTable >${FinalTable}.tmp && mv ${FinalTable}.tmp $FinalTable


############## compute the differential usage
#############################################

echo -e "\n==== Computing of the differential usage ====\n"

start_date=$(date)

$getSwitches $output_dir $FinalTable $design $normalized_gene_counts || { echo "getSwitches.R script failure !" 1>&2; exit; }

end_date=$(date)

echo -e "\nstart computing of differential usage : $start_date\n"
echo -e "\nend computing of differential usage : $end_date\n"


############## Clustering of the contigs by loci (genic/antisense/intergenic/unmapped)
######################################################################################

if [ "$stranded" == "yes" ];then
       
        echo -e "\n==== Clustering of the contigs by loci (genic/antisense/intergenic) ====\n"
	
	#cluster contigs per loci
	start_date=$(date)	
	$getContigsPerLoci $output_dir $FinalTable $threads || { echo "getContigsPerLoci.R script failure !" 1>&2; exit; }
	end_date=$(date)

	echo -e "\nstart clustering of contigs : $start_date\n"
	echo -e "\nend clustering of contigs : $end_date\n"

fi

#delete some temp files
rm ${temp_dir}ref_annotation.tmp ${temp_dir}OriginalFastaContigs.tmp ${temp_dir}sense_contigs.txt ${temp_dir}antisense_contigs.txt ${temp_dir}Intronic_contigs.txt ${temp_dir}Exonic_contigs.txt ${temp_dir}UTR_contigs.txt

echo -e "\n**** Annotation done ****\n"

