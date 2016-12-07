#!/bin/bash
#http://blog.nextgenetics.net/?e=8
 
 
 usage() { echo -e "Usage: $0 -q <query in fasta> -s <target in fasta> -o <output_dir> -d <db name> [options]\n\n\tOptions :\n
                  -t <\"n\" or \"p\"> : nucleotide or protein blast (default : \"n\").\n
                  -b <\"blastn\" or \"tblasx\" or \"blast\" or \"blastx\" or \"tblastn\"> : type of comparison (default : blastn).\n
                  -c <\"yes\" or \"no\" (default : \"yes\")> : filtering for low-complexity, repeats, etc. If \"no\", \"-dust no\" is set.\n
                  -e <evalue threshold (default : 1e-4)>\n
                  -m <max target seq (default : 500)>\n
                  -n <thread number (default : 1)>\n\n\tResults :\n
                  -p <location of blast scripts> : path to blast scripts (default : search in path environment variable).\n
                  -w <word size (default : 11)>\n
                  
                  - Blast table format 6 (raw_blast.alignment_2.txt), with the following columns :\n
                  \t\"qseqid\", \"sseqid\", \"pident\", \"length\", \"mismatch\", \"gaps\", \"qstart\", \"qend\", \"sstart\", \"send\", \"evalue\", \"bitscore\", \"qlen\", \"sstrand\".\n
                  - File with raw alignment (raw_blast.alignment_1.txt)\n" 1>&2; exit 1;}


[[ $# -eq 0 ]] && usage


while getopts ":q:s:o:t:b:e:p:w:d:m:n:c:" opt; do
  case $opt in
  
      q)
      
	      query=$OPTARG
	      echo "query file is: $query" >&2
	      
              ;;
      
      s)
      
	      subject=$OPTARG
	      echo "target file is: $subject" >&2
      
              ;;
      
      o)
      
             output=$OPTARG
 
             echo "output directory is: $output" >&2
             
             ;;
      
      t)
      
      	     type=$OPTARG
      
             ;;
      
      b)
      
             blast=$OPTARG
      
             ;;

      e)
      
             evalue=$OPTARG

             ;;
      
      p)
      
             program_path=$OPTARG
     
             ;;
      
      w)
      
             word_size=$OPTARG

             ;;
      
      d)
      
            db_name=$OPTARG
     
            echo "db name is: $db_name" >&2
       
            ;;
      
      m)
      
            max_target_seqs=$OPTARG

            ;;
      
      c)
      
            filtering=$OPTARG

            ;;
      
      n)
      
            num_threads=$OPTARG

            ;;
      
      
      #invalid options (options not in the list)
      ######################
      
      
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    
  esac
done

#if path to blast scripts isn't set, look in the path environment variable
if [ "$program_path" == "" ];then

	program_path=$(which blastn)

        if [ "$program_path" != "" ];then program_path=$(dirname $program_path);fi

fi
echo "location of blast scripts is: $program_path"

 #if one of the arguments is missing, exit with the usage (can also use "-z")
 
if [ "$query" == "" ] || [ "$subject" == "" ] || [ "$output" == "" ] || [ "$db_name" == "" ] || [ "$program_path" == "" ]; then
      
echo -e "\none required argument is missing or is wrong, check again !!\n"
      usage
      exit 1
fi

#create output directory         
output="${output}/"
output=$(echo $output |sed 's/\/\//\//g')
if [ ! -d $output ]; then mkdir $output ||{ echo "impossible to create output directory ! " 1>&2; exit 1; } ; fi 


if [ "$type" == "" ];then type="n";fi
if [ "$type" != "n" ] ;then

  if [ "$type" != "p" ];then
   echo "wrong parameter for type of sequence"
   echo $type
   echo "it should be n or p"
   exit 1
  fi
fi

if [ "$type" == "p" ];then type="prot";fi
if [ "$type" == "n" ];then type="nucl";fi
echo "type of sequence is: $type"

#-dust no : no filter
if [ "$filtering" == "no" ];then 
  
    filtering="-dust no"
    
    echo "no filtering ($filtering) !"

else

    echo "filtering is set !"
    filtering="" 
    
fi

if [ "$max_target_seqs" == "" ];then max_target_seqs=500;fi
echo "max_target_seqs is : $max_target_seqs" 

if [ "$evalue" == "" ];then evalue=1e-4;fi
echo "e-value threshold is: $evalue" 

if [ "$num_threads" == "" ];then num_threads=1;fi
echo "num_threads is : $num_threads"

if [ "$word_size" == "" ];then word_size=11;fi
echo "word size is: $word_size"

if [ "$blast" == "" ];then blast="blastn";fi
echo "type of blast is: $blast"


program_path="${program_path}/"
program_path=$(echo $program_path |sed 's/\/\//\//g')

if ! [[ $(ls ${output} |grep "${db_name}\..*"    ) ]];then 

  echo -e "\nbuilding blast database\n" 
   
  ${program_path}makeblastdb -in $subject -dbtype $type -title $output -out ${output}$db_name ||{ echo "making db 1 blast failure ! " 1>&2;rm ${output}/$db_name ; exit 1; }


fi

echo -e "\nmaking blast now\n"
    
#${program_path}${blast} -query $query -db ${output}$db_name -evalue $evalue $filtering -word_size $word_size -max_target_seqs $max_target_seqs -task ${blast} -outfmt "0" -out ${output}raw_blast.alignment_1.txt -num_threads $num_threads 2> /dev/null  || { echo "blast failure 1" 1>&2;exit 1;}



#cat $query|paste - - |awk -v split_lines=$split_lines -v output_dir=${output} 'NR%split_lines==1{OFS="\t";x=++i"_blastsubfile.tmp"}{OFS="";print > output_dir x}'

#for i in $(find ${output} -name "*_blastsubfile.txt");do
#  
#   sed 's/\t/\n/g' $i >

#done

#if it takes too long, split the files, and process them in parallel
${program_path}${blast} -query $query -db ${output}$db_name -evalue $evalue $filtering -word_size $word_size -max_target_seqs $max_target_seqs -task ${blast} -outfmt "6 qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore qlen sstrand" -out ${output}raw_blast.alignment_2.txt -num_threads $num_threads 2> /dev/null || { echo "blast failure 2" 1>&2;exit 1;}

echo -e "alignment done !\n"  
    
    
    
    

