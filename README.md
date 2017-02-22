

Required dependencies :

	- bash
	- awk
	- R libraries foreach & doParallel
	- bedtools >= 2.24 (preferentially 2.24)
	- Blast 2.5.0+
	- GSNAP

========================================================


Usage: ./getContigsAnnotation.sh -a <contigs from DEkupl> -g <genome in fasta> -d <diff genes from DEkupl> -r <reference annotation (gff3 format)> -o <output_dir> -i <path to illumina adapters> [options]

	Options :

                  -b        path to bin/ of blast scripts (default : in $PATH environment variable)

                  -c        path to bedtools, preferentially 2.24 (default : in $PATH environment variable)

                  -j        GSNAP index name

                  -k        path to directory of GSNAP index

                  -m        path to bin/ of GSNAP (default : in $PATH environment variable)

                  -p        padj diff. gene threshold (default : 0.05)

                  -s        path to samtools (default : in $PATH environment variable)

                  -n        thread number (default : 1)

	Results :

                  - Table "DiffContigsInfos.tsv", summarizing for each assembly, its location on the genome (if it's aligned), the neighborhood, the sequence alignment informations, and the differential expression informations

                  - BED "file diff_contigs.bed" for the visualization ; it contains useful informations from the summarization table.

                  
                  - Table "ContigsPerLoci.tsv" containing loci with differentially expressed contigs
          
========================================================                  
                  
Steps of the annotation : 

	1- Contigs from merged-diff-counts.tsv are converted in fasta.
	
	2- Contigs which align on the adapters are filtered out.
	
	3- The clean contigs are mapped on the genome with GSNAP.
	
	4- The BAM file is parsed to obtain a primary BED12 file (diff_contigs.bed), and a table summarizing the alignment (DiffContigsInfos.tsv). This step is parallelized.
	
	5- Unmapped contigs with GSNAP are aligned on the genome with Blast. Those with full length alignment are pushed in the BED12 and the  summarizing table.
	
	6- The BED12 file is enriched with infos from merged-diff-counts.tsv.
	
	7- Unaligned contigs are pushed in the summarizing table, with the appropriates values "NA", "F" or "none".
	
    8- The summarizing table is enriched with neighborhood of the contigs (antisense, closest genes, etc).
	
	9- From the summarizing table, a table of the contigs grouped by loci (genic, antisense, intergenic), and a pdf are built (DiffContigsInfos.tsv & ContigsPerLoci.pdf).                  
	
========================================================

Description of the outputs :

 	1- diff_contigs.bed : BED12 file of aligned contigs with GSNAP/Blast, with strand-specific color (red : strand + ; blue : strand -), with intensity scaled on the abs(log2FC).
 	
 	2- DiffContigsInfos.tsv : summary of all supplied contigs, no matter that they are aligned or not. The file "Contigs_Ontology_24_11_2016.xlsx" supplied with this README gives the ontology used. The column "TERM" in this xlsx file gives all the columns added by the annotation in the summarizing table.
 	
 	3- ContigsPerLoci.tsv : contigs grouped by loci (genic, antisense, intergenic). The locus ID for a genic/antisense locus is the HUGO ID. For an intergenic locus, we have the concatenation of : chromosome, strand, 5'-gene,3'-gene, separated by "&".
 	
 	5- ContigsPerLoci.pdf : Result of "ContigsPerLoci.tsv" in barplots
 	
========================================================

Miscellaneous :

  

  	- If the file "nomatch_in_adapters.fa" containing unmapped contigs on the adapters is created, the step of blast on adapters will be avoided. If you want to re-run it anyway, put in comment the if-statement.
  
  	- If the alignemnt file of contigs is already created, this step will be avoided. If you want to re-run it anyway, put in comment the if-statement.
  	
  	- If Blast data bases for the human genome and the adapters is already created, this step in the script blast.sh will be avoided.
  	
  	- The Blast of unmapped contigs on the genome (after the mapping with GSNAP) can take long time even though it is parallelized ; in script blast.sh, you can split the input file in many, and run blast on each sub-file (with xargs) to speed up.
