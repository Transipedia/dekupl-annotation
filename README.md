# DE-kupl annotation

DE-kupl annotation is part of the DE-kupl package, and performs annotations of DE contigs identified by DE-kupl.

## Required dependencies :

* bash (version >= 4.3.46)
* awk (version >= 4.1.3)
* R (version >= version 3.2.3) with libraries foreach, doParallel, GenomicFeatures & DESeq2
* bedtools (versions 2.24, 2.26 ; don't use version 2.25 !!)
* Blast (version >= 2.5.0+)
* GSNAP (version >= 2016-11-07)
* samtools (version >= 1.3)

Note: needs sufficient memory to index the genome (as much as gmap_build needs, a few tens of GB's).

## Usage : 

Usage: ./getContigsAnnotation.sh < Required arguments > [Optional arguments]


 	Required arguments :

                  -a   <merged-diff-counts.tsv.gz (contigs in "{A}_vs_{B}_kmer_counts" directory from Dekupl-run result)>

                  -g   <genome in fasta (uncompressed). Use the same version as the annotation to avoid chromosome name issues>

                  -d   <{A}vs{B}-DEGs.tsv (diff. genes in "gene_expression" directory from Dekupl-run result)>

                  -r   <reference annotation (gff3 format, uncompressed). Use the same version as the genome file to avoid chromosome name issues>

                  -t   <are the reads stranded ? (choose "yes" or "no")>

                  -e   <normalized_counts.tsv (normalized gene counts in "gene_expression" directory from Dekupl-run result)>

                  -f   <sample_conditions.tsv (design file in "metadata" directory from Dekupl-run result)>

                  -o   <path to output directory>

                  -i   <illumina adapters (you can use the file "adapters.fa" supplied with the program)>


	Optional arguments :

                  -b   <path to bin/ of blast scripts (default : in $PATH environment variable)>

                  -c   <path to bedtools, preferentially 2.24 (default : in $PATH environment variable)>

                  -j   <GSNAP genome index name (if the index of the genome has already been created, supply its name. Otherwise, the index is re-created with the name "genome_index")>

                  -k   <path to directory of GSNAP genome index (if you have a former run, you can supply the full path of the \"mapping_output/\" inside your former output directory, in order to re-use the same genome index and save time)>

                  -m   <path to bin/ of GSNAP (default : in $PATH environment variable)>

                  -p   <padj diff. gene threshold (default : 0.05)>

                  -s   <path to samtools (default : in $PATH environment variable)>

                  -n   <thread number (default : 1)>

## Warnings :

- After the downloading of `dekupl-annotation`, make `chmod 755 *.{sh,R}` into the directory to avoid rights issues.
- Use the same genome & annotation files version as the one of the file `transcripts.fa` used by Kallisto in Dekupl-run (e.g by default it's gencode 24, if it has changed, you should take this into account here, to avoid incoherences in the results because of non-maching IDs).

## Results :

- Table `DiffContigsInfos.tsv`, summarizing for each assembly, its location on the genome (if it's aligned), the neighborhood, the sequence alignment informations, and the differential expression informations
- BED `file diff_contigs.bed` for the visualization ; it contains useful informations from the summarization table.
- Table `ContigsPerLoci.tsv` (only for stranded data for now) containing loci with differentially expressed contigs
          
## Steps of the annotation : 

1. Contigs from `merged-diff-counts.tsv.gz` are converted in fasta.	
2. Contigs which align on the adapters are filtered out.
3. The clean contigs are mapped on the genome with GSNAP.
4. The BAM file is parsed to obtain a primary BED12 file (`diff_contigs.bed`), and a table summarizing the alignment (`DiffContigsInfos.tsv`). This step is parallelized.
5. Unmapped contigs with GSNAP are aligned on the genome with Blast. Those with full length alignment are pushed in the BED12 and the  summarizing table.
6. The BED12 file is enriched with infos from merged-diff-counts.tsv.
7. Unaligned contigs are pushed in the summarizing table, with the appropriates values "NA", "F" or "none".
8. The summarizing table is enriched with neighborhood of the contigs (antisense, closest genes, etc).
9. From the summarizing table, a table of the contigs grouped by loci (genic, antisense, intergenic, unmapped) is built (`DiffContigsInfos.tsv`).                  

## Description of the outputs :

1. **`diff_contigs.bed`** : BED12 file of aligned contigs with GSNAP/Blast, with strand-specific color (red : strand + ; blue : strand -), with intensity scaled on the abs(log2FC).
2. **`DiffContigsInfos.tsv`** : summary of all supplied contigs, no matter that they are aligned or not. The file `Contigs_Ontology_24_11_2016.xlsx` supplied with this README gives the ontology used. The column "TERM" in this xlsx file gives all the columns added by the annotation in the summarizing table.
3. **`ContigsPerLoci.tsv`** : contigs grouped by loci (genic, antisense, intergenic, unmapped). The locus ID for a genic/antisense locus is the Ensembl ID followed by the strand (separated by "&"). For an intergenic locus, we have the concatenation of : chromosome, strand, 5'-gene, 3'-gene (separated by "&").
 	
## Miscellaneous :

- If Blast data bases for the human genome and the adapters are already created, this step in the script blast.sh will be avoided.
- The Blast of unmapped contigs on the genome (after the mapping with GSNAP) can take long time even though it is parallelized ; in script blast.sh, you can split the input file in many, and run blast on each sub-file (with xargs) to speed up.
