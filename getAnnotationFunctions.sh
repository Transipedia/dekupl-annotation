#!/bin/bash

################# functions ###########
#######################################


#this function will parse the modified bed file from parseBam() to give informations on the alignments
#it creates a bed file and a table with the informations
#modifed bed : 12 classical bed columns + CIGAR + MD tag + NH tag + other_split + line number
parseModifiedBed(){

        modifiedBed=$1
        
        mapper_name=$2
        
        forward_contig_color=$3
        
        reverse_contig_color=$4
        
        file_name=$(basename $modifiedBed)
        
        temp_dir="$(dirname $modifiedBed)/"
        
        file_number=$(echo ${file_name//_/ }|awk '{print $1}')
        
        if [ -f ${temp_dir}${file_number}_subTable ];then rm ${temp_dir}${file_number}_subTable ;fi
        if [ -f ${temp_dir}${file_number}_subBed ];then rm ${temp_dir}${file_number}_subBed ;fi
    

        j=1
 	while read line;do
 	
	
	   LineInSam=$(echo "$line"|cut -f 17)
	
	   #echo -e "line in SAM is : $LineInSam"

	   ID=$(echo "$line"|cut -f 4)

	   #echo -e "ID is : $ID"

	   chromosome=$(echo "$line"|cut -f 1)
	   
	   start=$(echo "$line"|cut -f 2)
	   
	   end=$(echo "$line"|cut -f 3)
	   
	   strand=$(echo "$line"|cut -f 6)
	   
	   if [ "$strand" == "+" ];then color=$forward_contig_color ; else color=$reverse_contig_color ;fi
	   
	   other_split=$(echo "$line"|cut -f 16)
	   
	   nb_mismatch=$(echo "$line"|cut -f 14)
	   
	   nb_hit=$(echo "$line" |cut -f 15)
	   
	   CIGAR=$(echo "$line"|cut -f 13)
	   
	   
	   split_nb_insertion=($(echo $CIGAR|grep -o -E "[0-9]+I"|sed 's/I//g'))
	 
	   nb_insertion=0
	   if [[ ${#split_nb_insertion[*]} -ne 0 ]];then
	 
		 for i in $(seq 0 $((${#split_nb_insertion[*]}-1)));do
		 
		    nb_insertion=$((nb_insertion+${split_nb_insertion[i]}))
		    
		 done
	 
	   fi
	   
	   #echo -e "insertions are : ${nb_insertion}"
	 
	   split_nb_deletion=($(echo $CIGAR|grep -o -E "[0-9]+D"|sed 's/D//g' ))
	 
	   nb_deletion=0
	   if [[ ${#split_nb_deletion[*]} -ne 0 ]];then
	 
		 for i in $(seq 0 $((${#split_nb_deletion[*]}-1)));do
		 
		    nb_deletion=$((nb_deletion+${split_nb_deletion[i]}))
		    
		 done
	 
	   fi
	   
	   #echo -e "deletions are : ${nb_deletion}"
	   
	   #flag contigs with SNV
	   if [[ $nb_mismatch -gt 0 ]] || [[ $nb_deletion -gt 0 ]] || [[ $nb_insertion -gt 0 ]];then SNV="T" ; else SNV="F";fi

	   #looking for 5'/3' soft/hard clipped for strand +
	   if [ "$strand" == "+" ];then
	 
		 nb_5SoftClipped=$(echo $CIGAR|grep -o -E "^[0-9]+S"|sed 's/S//g')
		 
		 #echo -e "strand + ; split 5' soft clipped is : $nb_5SoftClipped"
		 
		 
		 nb_5HardClipped=$(echo $CIGAR|grep -o -E "^[0-9]+H"|sed 's/H//g')
		 
		 #echo -e "strand + ; split 5' hard clipped is : $nb_5HardClipped"
		 
		 
		 #we merge 5' soft and 5' hard clipped
		 nb_5SoftClipped=$((nb_5SoftClipped+nb_5HardClipped+0))

		 #echo -e "strand + ; total 5' clipped is : $nb_5SoftClipped"



		 nb_3SoftClipped=$(echo $CIGAR|grep -o -E "[0-9]+S$"|sed 's/S//g')
		 
		 #echo -e "strand + ; split 3' soft clipped is : $nb_3SoftClipped"
		 
		 
		 nb_3HardClipped=$(echo $CIGAR|grep -o -E "[0-9]+H$"|sed 's/H//g')
		 
		 #echo -e "strand + ; split 3' hard clipped is : $nb_3HardClipped"

		 nb_3SoftClipped=$((nb_3SoftClipped+nb_3HardClipped+0))

		 #echo -e "strand + ; total 3' clipped is : $nb_3SoftClipped"


	   #looking for 5'/3' soft/hard clipped for strand -
	   else

		 nb_5SoftClipped=$(echo $CIGAR|grep -o -E "[0-9]+S$"|sed 's/S//g')
		 
		 #echo -e "strand - ; split 5' soft clipped is : $nb_5SoftClipped"
		 
		 
		 nb_5HardClipped=$(echo $CIGAR|grep -o -E "[0-9]+H$"|sed 's/H//g')
		 
		 #echo -e "strand - ; split 5' hard clipped is : $nb_5HardClipped"
		 
		 #we merge 5' soft and hard clipped
		 nb_5SoftClipped=$((nb_5SoftClipped+nb_5HardClipped+0))

		 #echo -e "strand - ; total 5' clipped is : $nb_5SoftClipped"


		 nb_3SoftClipped=$(echo $CIGAR|grep -o -E "^[0-9]+S"|sed 's/S//g')
		 
		 #echo -e "strand - ; split 3' soft clipped is : $nb_3SoftClipped"
		 
		 
		 nb_3HardClipped=$(echo $CIGAR|grep -o -E "^[0-9]+H"|sed 's/H//g')

		 
		 #echo -e "strand - ; split 3' hard clipped is : $nb_3HardClipped"

		 #we merge 3' soft and hard clipped
		 nb_3SoftClipped=$((nb_3SoftClipped+nb_3HardClipped+0))

		 #echo -e "strand - ; total 5' clipped is : $nb_3SoftClipped"

	      
	   fi
	   
	   #all exons lengths (array)
	   split_exons_length=($(echo "$line"|cut -f 11|sed 's/,/ /g'))
	   
	   #all exons relative starts (array)
	   split_exons_RelStart=($(echo "$line"|cut -f 12|sed 's/,/ /g'))
	   
	   #initialization of coordinates of all exons
	   exon_coord_list=()
	   
	   #storage of all exons absolute starts and end in order to use them to infer coordinates of junctions
	   all_exons_AbsStart=()
	   all_exons_AbsEnd=()
	  
	   for i in $(seq 0 $((${#split_exons_length[*]}-1)));do
	   
	      one_exon_AbsStart=$((start+${split_exons_RelStart[i]}))
	      all_exons_AbsStart+=($one_exon_AbsStart)
	      
	      one_exon_AbsEnd=$((one_exon_AbsStart+${split_exons_length[i]}))
	      all_exons_AbsEnd+=($one_exon_AbsEnd)
	   
	      exon_coord_list+=("${chromosome}:${one_exon_AbsStart}-${one_exon_AbsEnd}")
	   
	   done
	   
	   exon_coord_list=$(echo ${exon_coord_list[*]}|sed 's/ /+/g')
	   
	   #initialization of coordinates of all junctions
	   junction_coord=()
	   
	   if [[ ${#all_exons_AbsStart[*]} -gt 1 ]];then
	   
	      nb_junction=$((${#all_exons_AbsEnd[*]}-1))
	   
	      for i in $(seq 0 $((nb_junction-1)));do
	      
	      
		  one_junc_start=$((${all_exons_AbsEnd[i]}+1))
		
		  one_junc_end=$((${all_exons_AbsStart[i+1]}-1))
		  
		  junction_coord+=("${chromosome}:${one_junc_start}-${one_junc_end}")
	      
	      done
	      
	   else
	   
	      nb_junction=0
	   
	      junction_coord="none"
	   
	   fi
	   
	   junction_coord=$(echo ${junction_coord[*]}|sed 's/ /+/g')
	   
	   #echo -e "exons coordinates are : ${exon_coord_list}"
	   #echo -e "junctions coordinates are : ${junction_coord}"
	   
	   #echo -e "\n --- \n"
	   
	   #all exons lengths (sep by "," in order to re-use them in the bed)
	   split_exons_length=$(echo ${split_exons_length[*]}|sed 's/ /,/g')
	   
	   #all relative starts (sep by "," in order to re-use them in the bed)
	   split_exons_RelStart=$(echo ${split_exons_RelStart[*]}|sed 's/ /,/g')
	   
	   #we add +1 at the start for the final table because it was in bed format (0-based)
	   echo -e "${LineInSam}\t${ID}\tT\t${mapper_name}\t${chromosome}\t$((${start}+1))\t${end}\t${junction_coord}\t${nb_junction}\t${exon_coord_list}\t$other_split\t${strand}\t${nb_hit}\t${nb_mismatch}\t${nb_deletion}\t${nb_insertion}\t$SNV\t${nb_5SoftClipped}\t${nb_3SoftClipped}" >>${temp_dir}${file_number}_subTable
	   
	   
	   echo -e "${chromosome}\t${start}\t${end}\tLineInSam=${LineInSam};ID=${ID};nb_hit=${nb_hit};nM=${nb_mismatch};del=${nb_deletion};ins=${nb_insertion};clipped_5p=${nb_5SoftClipped};clipped_3p=${nb_3SoftClipped}\t1\t${strand}\t${start}\t${end}\t${color}\t$((${nb_junction}+1))\t${split_exons_length}\t${split_exons_RelStart}" >>${temp_dir}${file_number}_subBed
	   
	   #echo -e "${chromosome}\t${start}\t${end}\tLineInSam=${LineInSam};ID=${ID};nb_hit=${nb_hit};nM=${nb_mismatch};del=${nb_deletion};ins=${nb_insertion};clipped_5p=${nb_5SoftClipped};clipped_3p=${nb_3SoftClipped}\t1\t${strand}\t${start}\t${end}\t${color}\t$((${nb_junction}+1))\t${split_exons_length}\t${split_exons_RelStart}"
	   
	j=$((j+1))  
	
	#if [[ j -eq 76 ]];then exit ;fi
	
	
	done < $modifiedBed
	
	
}


#this function will build, from a bam file, a modified bed file (12 columns + column of CIGAR + column of MD tag + column flag + multihit column + line number)
#it will call parseModifiedBed() in parallel on each sub modified bed to parse them, and return a complete bed and table
#1st argument : bam file
#2nd argument : output table
#3rd argument : output bed 
#4th argument color of contigs on forward strand
#5th argument color of contigs on reverse strand
parseBam(){
      
        bam_file=$1
        
        temp_dir="$(dirname $bam_file)/bam_file_processing/"
        
        if [ -d $temp_dir ];then 
        
           rm -rf $temp_dir
           
        else
        
          mkdir $temp_dir
        
        fi
        
        output_table=$2
        
        output_bed=$3
        
        forward_contig_color=$4
        
        reverse_contig_color=$5
     
        if [ -f $output_table ];then rm $output_table ; fi
        if [ -f $output_bed ];then rm $output_bed ; fi
        
        $samtools view -H $bam_file >${temp_dir}sam_header.txt
        
        #we convert primary alignment in col ID+bed12, and we keep the longer seq for each ID (useful if we have chimeric alignments).
        #remark : we cannot have the extra SAM tags in the bed12 file, so we store the samtools result before the conversion in ${temp_dir}TAGS.tmp
        #in the last pipe, we remove /1 or /2 that could be added by bedtools for read1 or read2
	$samtools view -F 4 -F 0x100 $bam_file |awk 'OFS="\t"{$(NF+1)=length($10);print}' |LANG=en_EN sort -k 1,1 -k10,10nr|LANG=en_EN sort -u -k 1,1|awk 'OFS="\t"{$NF="";print}' |LANG=en_EN sort -k1,1 -k3,3|tee ${temp_dir}TAGS.tmp |cat ${temp_dir}sam_header.txt - |$samtools view -bh|$bedtools bamtobed -bed12 -i stdin|awk 'OFS="\t"{print $4,$0}' |awk 'OFS="\t"{gsub(/\/2$/,"",$1);gsub(/\/1$/,"",$1);gsub(/\/2$/,"",$5);gsub(/\/1$/,"",$5);print}' >${temp_dir}primary_aligment.txt

        #get col ID, CIGAR, MD & NH tags (if these last tags are missing, the value "NA" is given)
        awk -F '\t' '{OFS="\t";ID=$1;CIGAR=$6;MD="";NH="";for(i=12;i<=NF;i++){if($i~/^MD\:Z\:.*/){gsub(/MD\:Z\:/,"",$i);if($i~/[0-9]/ || $i~/[A-Z]/){gsub(/[0-9]/,"",$i);gsub(/\^/,"",$i);if($i==""){$i=0}else{$i=length($i)}};MD=$i};if($i~/NH\:i\:[0-9]+/){gsub(/NH\:i\:/,"",$i);NH=$i}};if(MD==""){MD="NA"};if(NH==""){NH="NA"};print ID,CIGAR,MD,NH}' ${temp_dir}TAGS.tmp >${temp_dir}TAGS.txt && rm ${temp_dir}TAGS.tmp 
	
	#join col ID+bed12 with the CIGAR, NH tag, & MD tag (no need to sort)
	LANG=en_EN join -t $'\t' -11 -21 ${temp_dir}primary_aligment.txt ${temp_dir}TAGS.txt >${temp_dir}modifiedBed12.tmp
	
	#looking for chimeric alignments with the flag
	$samtools view -F 4 -F 0x100 -f 0x800 $bam_file >${temp_dir}chimeric_split1.tmp 
	
	#if there's no such flag, keep the ID, and put in the second column False for all the contigs 
	if [[ $(wc -l ${temp_dir}chimeric_split1.tmp|awk '{print $1}' ) -eq 0 ]];then 
	
	   awk 'OFS="\t"{print $1,"F"}' ${temp_dir}primary_aligment.txt >${temp_dir}chimeric_split1.txt && rm ${temp_dir}chimeric_split1.tmp
	   
	#otherwise, keep the ID, and put in the second column True   
	else
	
	   awk 'OFS="\t"{print $1,T}' ${temp_dir}chimeric_split1.tmp >${temp_dir}chimeric_split1.txt && rm ${temp_dir}chimeric_split1.tmp
	
	fi

        #looking for chimeric alignments with the number of seq (with GSNAP, assembly is split in reads 1 & 2 on 2 chromosomes, and these don't have the flag 0x800)
	IDs=($($samtools view -F 4 -f 0x1 $bam_file|cut -f 1|LANG=en_EN sort -u))

        #keep only ID & chromosome
	$samtools view -f 0x1 $bam_file|cut -f 1,3 >${temp_dir}chimeric_split2.tmp
	
	#if we have splits, we will check if they are chimeric 
	if [[ $(wc -l  ${temp_dir}chimeric_split2.tmp|awk '{print $1}') -gt 0 ]];then 

		if [ -f ${temp_dir}chimeric_split2.txt ];then rm ${temp_dir}chimeric_split2.txt;fi

		for i in $(seq 0 $((${#IDs[*]}-1)));do

		   One_ID=${IDs[$i]}
		  
		   all_chr=$(grep "$One_ID" ${temp_dir}chimeric_split2.tmp |cut -f 2|LANG=en_EN sort -u|wc -l)
		   
		   
		   #if we have more than 1 chromosome for the assembly, it's chimeric
		   if [[ $all_chr -gt 1 ]];then
		  
		   
		      echo  -e "${One_ID}\tT" >>${temp_dir}chimeric_split2.txt
		    
		   else
		    
		      echo  -e "${One_ID}\tF" >>${temp_dir}chimeric_split2.txt
		    
		   fi

		   
		done
		
	
	else
	
		awk 'OFS="\t"{print $1,"F"}' ${temp_dir}primary_aligment.txt >${temp_dir}chimeric_split2.txt
	
	fi
	
	rm ${temp_dir}chimeric_split2.tmp ${temp_dir}primary_aligment.txt
	
	#concatenate both types of chimeric
	cat ${temp_dir}chimeric_split1.txt ${temp_dir}chimeric_split2.txt >${temp_dir}other_split.txt
	
	#unique on the ID : keep preferentially IDs with the value T
	LANG=en_EN sort -k 1,1 -k 2,2r ${temp_dir}other_split.txt |LANG=en_EN sort -u -k 1,1 >${temp_dir}other_split.tmp && mv ${temp_dir}other_split.tmp ${temp_dir}other_split.txt
	
	#reconstruct the modified bed (12 classical bed columns + CIGAR + MD tag + NH tag + other_split + line number)
	LANG=en_EN join -t $'\t' -11 -21 <(LANG=en_EN sort -k1,1 ${temp_dir}modifiedBed12.tmp) <(LANG=en_EN sort -k1,1 ${temp_dir}other_split.txt)| cut -f 2-|awk 'OFS="\t"{print $0,NR}' >${temp_dir}modifiedBed12.txt

        mapper_name=$(grep "@PG" ${temp_dir}sam_header.txt|cut -f 2|sed 's/ID://g')

        total_line=$(tac ${temp_dir}modifiedBed12.txt|head -n 1|awk '{print $NF}')
       
        echo -e "number of aligned contigs to parse is : $total_line"
     
        split_lines=$(echo $((total_line/threads))|awk '{print int($1)}')
        
        if [[ $split_lines -eq 0 ]];then split_lines=$threads ; fi
     
        echo -e "the file will be split, and each sub-file will contain at most $split_lines lines to process in parallel\n"
       
        #removing pre-existing input sub-files 
        find $temp_dir -name "*_subfile.txt" -type f -delete
       
        #split the modified bed in many chunks, then run the function parseModifiedBed() on them
        awk -v split_lines=$split_lines -v temp_dir=$temp_dir 'NR%split_lines==1{OFS="\t";x=++i"_subfile.txt"}{OFS="";print > temp_dir x}' ${temp_dir}modifiedBed12.txt
       
        #removing pre-existing sub-results
        find $temp_dir -name "*_subBed" -type f -delete
        find $temp_dir -name "*_subTable" -type f -delete
       
        export -f parseModifiedBed
       
        find $temp_dir -name "*_subfile.txt" |grep -E "[0-9]+_subfile.txt" |xargs -I {} -P $threads bash -c "parseModifiedBed {} "$mapper_name" "$forward_contig_color" "$reverse_contig_color""
      
        cat $(find $temp_dir -name "*_subBed") |LANG=en_EN sort -k 1,1 -k 2,2n >$output_bed && find $temp_dir -name "*_subBed" -type f -delete
        cat $(find $temp_dir -name "*_subTable") >$output_table && find $temp_dir -name "*_subTable" -type f -delete
       
        rm -rf $temp_dir 
       
}

#this function will rebuild a fasta file from the unmapped contigs
#1st argument : unmapped contigs (1st col = line in sam ; 2nd col = tag ID)
#2nd argument : initial fasta (1st column = >tag ID ; 2nd column = seq, it's generated by this script)
#3rd argument : output name
getFastaFromUnmappedContigs(){

  unmapped_seq=$1
  initial_fasta=$2
  output_fasta=$3
  
  
  #put tag ID before line in SAM, in a temp file (trick to conserve the line in SAM with the tag ID, in order to re-insert the result in the final table 
  awk 'OFS="\t"{print $2,$1}' $unmapped_seq|LANG=en_EN sort -k1,1 >${unmapped_seq}.tmp
  
  #1st block : sort by tag ID unmapped contigs
  #2nd block : sort by tag ID original fasta
  #3rd block : join both by tag ID
  #add the line in sam just after the tag ID (it will be used laparseModifiedBedter for the joining)

  LANG=en_EN join -t $'\t' -11 -21 <(awk '{print ">"$2}' $unmapped_seq |LANG=en_EN sort -k1,1) <(LANG=en_EN sort -k1,1 $initial_fasta)|sed 's/>//g'| LANG=en_EN join -t $'\t' -11 -21 ${unmapped_seq}.tmp - |awk 'OFS="\t"{print ">"$1"-"$2,$3}'|sed 's/\t/\n/g' >$output_fasta

}
