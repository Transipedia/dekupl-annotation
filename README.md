![dekupl-annot logo](dekupl-annot-logo.png)

[![Build Status](https://travis-ci.org/Transipedia/dekupl-annotation.svg?branch=new-version)](https://travis-ci.org/Transipedia/dekupl-annotation) [![docker pull](https://img.shields.io/docker/pulls/transipedia/dekupl-annotation.svg)](https://hub.docker.com/r/transipedia/dekupl-annotation/) [![conda install](https://anaconda.org/transipedia/dekupl-annotation/badges/downloads.svg)](https://anaconda.org/Transipedia/dekupl-annotation)

DE-kupl annotation is part of the DE-kupl package, and performs annotations of DE contigs identified by DE-kupl.

- [Usage](#usage)
    - [Creating the index](#creating-the-index)
    - [Annotating contigs](#annotating-contigs)
- [Tutorial & toys](#tutorial--toys)
- [Installation](#installation)
    - [Option 1: Use dekupl-annotation with conda](#option-1-use-dekupl-annotation-with-conda)
    - [Option 2: Use dekupl-annotation with Docker](#option-2-use-dekupl-annotation-with-docker)
    - [Option 3: Use dekupl-annotation with singularity](#option-3-use-dekupl-annotation-with-singularity)
    - [Option 4: Install from the sources (not recommended)](#option-4-install-from-the-sources-not-recommended)
- [Output files](#output-files)
- [Ontology](#ontology)
- [Dev environnement](#dev-environnement)

## Usage

### Creating the index

To use DEkupl Annotation, you first need to build an index directory, that will prepare all
the file for annotation.

Required file are a genome as a multi-FASTA and annotations in a GFF3 format. Both
of this file can be downloaded from [Ensembl](https://www.ensembl.org/info/data/ftp/index.html).

Input files can be gzipped.

```
Usage:
    dkpl index -g gff_file -f genome.fasta -i index_dir

  Mandatory Arguments:
      -a,--annotations FILE   GFF annotation file
      -g,--genome FILE        Genome if FASTA format
      -i,--index DIR          Output index directory.

  Optional Arguments:
      -h,--help           show this help message and exit
      -t,--threads INT    Number of threads
      --star              Index the genome with STAR (for discovery of chimeric RNA)
                          It needs ~30gb of RAM for Human genome
```

You can use the `--star` option to create a index of the genome with STAR that will
be used during the annotation for searching chimeric RNA (fusions, circRNA, etc).
STAR index is loaded into memory and requires up to 30gb for the human genome.

### Annotating contigs

Once an index has been constructed you can annotate contigs from a DEkupl-run.
The minimal input to run DEkupl-annotation is the contigs file
(merged-diff-counts.tsv.gz) produced by DEkupl-run. If the input files are stranded,
you need to add the `--stranded` option, in order to get the annotations on the right strand.

```
Usage:
    dkpl annot -i index_dir/ merged-diff-counts.tsv.gz

Options:
  Requiered Arguments:
      -i,--index DIR      path to the index directory (created with dekupl-annotation index)

  Input/Output:
      -o,--output DIR     path to the output directory (default: "DEkupl_annotation/")
      -d,--deg FILE       (Optional) {A}vs{B}-DEGs.tsv (diff. genes in "gene_expression" directory from Dekupl-run result)
      --norm-gene-counts FILE
                          (Optional) Normalized gene counts
      --sample-conditions FILE
                          (Optional) Sample conditions. First column is sample name,
                          second column is sample condition)
      --repeat            FILE
                    			(Optional) Repeat sequences reference. A human repeat reference is available in the toy directory.

  Optional Arguments:
      -t,--threads INT    Number of threads (for GSNAP)
      -s,--stranded       RNA-Seq is strand-specific.
      -p,--deg-padj       padj diff. gene threshold (default : 0.05)
      --max-splice-length
                          Splice with greater length are considered as chimeric junctions (default 1000000)
      --contig-color INT  Contig color mode (default 1):
                            1 : contigs on forward strand are in red (contigs on reverse strand are in blue)
		                        2 : contigs on forward strand are in blue (contigs on reverse strand are in red)
      --kamrat
                          Allow input tweaked by Kamrat analysis (with an additional score column)
      -h,--help           show this help message and exit
      -v,--verbose        print additional debug messages
```

Output files will be placed under the `DEkupl_annotation` directory unless you specify
another output directory with `-o` option.

Extra file can be supplied to complete the annotation process (see ontolgy table) :

* differentially expressed genes (`--deg`) for DEG analyzer
* nomarlized gene counts (`--norm-gene-counts`) and sample conditions (`--sample-conditions`) for Switches analyzer.

If the index was created with the `--star` option, then dkpl-annot will look for chimeric junctions.


## Tutorial & toys

Toy files are available with this repository to test dkpl-annot.

```
dkpl index -g toy/references/GRCh38-chr22.fa.gz -a toy/references/GRCh38-chr22.gff.gz -i test_index
dkpl annot -i test_index --deg toy/dkpl-run/DEGs.tsv.gz --norm-gene-counts toy/dkpl-run/normalized_counts.tsv --sample-conditions toy/dkpl-run/sample_conditions_full.tsv toy/dkpl-run/merged-diff-counts.tsv.gz
```
## Installation

We recommand tu use [singularity](https://singularity.lbl.gov/) to install dekupl-annotation, but you can also use Docker, and manual installation.

### Option 1: Use dekupl-annotation with singularity

One can create a singularity container from the docker image.

- **Step 1: Build Singularity image**
    ```
    singularity build dekupl-annotation.simg docker://transipedia/dekupl-annotation:1.0.3
    ```

It's advised to mount some volumes (input/output directories). To mount the "/store" volume you should use "--bind /store:/mnt". That way, you can access the /store directory by using /mnt. For example, if your input files for the index are located at /store/references/{files} :

- **Step 2a: Run dekupl-annotation index with mounted volumes**
    ```
    singularity run --bind /store:/mnt ./dekupl-annotation.simg index -g /mnt/references/GRCh38-chr22.fa.gz \
    -a /mnt/references/GRCh38-chr22.gff.gz -i mnt/index    
    ```
Then, if your DEkupl_run files are located at /store/DEkupl_result/{files} :

- **Step 3a: Run dekupl-annotation annot with mounted volumes**
    ```
    singularity run --bind /store:/mnt ./dekupl-annotation.simg annot -i /mnt/index /mnt/DEkupl_result/case_vs_control_kmer_counts/merged-diff-counts.tsv.gz -o ./ --sample-conditions /mnt/DEkupl_result/metadata/sample_conditions_full.tsv --deg /mnt/DEkupl_result/gene_expression/normalvstumor-DEGs.tsv --norm-gene-counts /mnt/DEkupl_result/gene_expression/normalized_counts.tsv --repeat /mnt/references/homo_sapiens/RMBlastLib.fasta    
    ```


Alternatively (Not advised !), you can move every input file in the directory where dekupl-annotation is running (subdirectories allowed).

- **Step 2b: Run dekupl-annotation with every input file in dekupl-annotation.simg directory**
    ```
    singularity run ./dekupl-annotation.simg index -g references/GRCh38-chr22.fa.gz \
    -a toy/references/GRCh38-chr22.gff.gz -i index
    singularity run ./dekupl-annotation.simg annot -i index DEkupl_result/case_vs_control_kmer_counts/merged-diff-counts.tsv.gz -o ./ --sample-conditions DEkupl_result/metadata/sample_conditions_full.tsv --deg DEkupl_result/gene_expression/normalvstumor-DEGs.tsv --norm-gene-counts DEkupl_result/gene_expression/normalized_counts.tsv --repeat references/homo_sapiens/RMBlastLib.fasta    
    ```


### Option 2: Use dekupl-annotation with Docker

    - **Step 1: Retrieve the docker image.**
        ```
        docker pull transipedia/dekupl-annotation:1.0.3
        ```
    - **Step 2: Run dekupl-annotation**.
        You will need to mount some volumes (input and output directories). If your input files for the index are located at /store/references/{files}
        ```
        docker run --rm -v /store/:/data/ \
        transipedia/dekupl-annotation index -g /data/references/GRCh38-chr22.fa.gz \
        -a /data/references/GRCh38-chr22.gff.gz -i /data/index
        ```


### Option 3: Install from the sources (not recommended)

- **Step 1: Install dependancies**. Before using Dekupl-annotation, install these dependencies:
    - **Required**: bash (version >= 4.3.46), R (version >= version 3.2.3) with libraries DESeq2, GSNAP (version >= 2016-11-07), samtools (version >= 1.3) & blast (version >= 2.5.0+)
    - **Optional** : STAR (version >= 2.5.3) for chimeric RNA
    - Installing dependancies on Linux Debian
    ```
    apt-get install cpanminus libdist-zilla-perl gmap samtools ncbi-blast+ rna-star
    Rscript install_r_packages.R # Install DESeq2 from bioconductor
    ```
- **Step 2: Install dekupl-annot**
    - **Global install**: The following command, will clone the repository and install dkpl-annot globaly with dzil and cpanm.
    ```
    git clone -branch ALtweaks https://github.com/Transipedia/dekupl-annotation.git && cd dekupl-annotation
    dzil install --install-command 'cpanm .'
    ```
    - **Local install**: For local install you need to use the `-l LOCAL_DIR` parameter of cpanm. Then you need to make sure that the Perl library that have been installed locally are available to the path using the `PERL5LIB` environnement variable.
    ```
    dzil install --install-command 'cpanm -l $HOME/.local .'
    export PERL5LIB=$HOME/.local/lib/perl5:$PERL5LIB
    ```

## Output files

- Table `DiffContigsInfos.tsv`, summarizing for each contig, its location on the genome (if it's aligned), the neighborhood, the sequence alignment informations, and the differential expression informations.

- BED file `diff_contigs.bed` for the visualization ; it contains useful informations from the summarization table. BED12 file of aligned contigs with GSNAP/Blast, with strand-specific color (red : strand + ; blue : strand - ; grey : unstranded). Color intensity is scaled on the log2FC value.

- Table `ContigsPerLoci.tsv`, contigs grouped by loci (genic, antisense, intergenic, unmapped). The locus ID for a genic/antisense locus is the Ensembl ID followed by the strand (separated by "&"). For an intergenic locus, we have the concatenation of : chromosome, strand, 5'-gene, 3'-gene (separated by "&").

## Ontology

| Term                    | Type  | Analyzer    | File Source                    | Description                                                                                        |
| ----------------------- | ----- | ----------- | ------------------------------ | -------------------------------------------------------------------------------------------------- |
| tag                     | Str   | Contigs     | contigs                        | kmer of reference                                                                                  |
| nb_merged_kmers         | Int   | Contigs     | contigs                        | Number of k-mer merged to produce the contig                                                       |
| contig                  | Str   | Contigs     | contigs                        | Sequence of the contig                                                                             |
| contig_size             | Int   | Contigs     | contigs                        | Size (in nucleotides) of the contig                                                                |
| pvalue                  | Float | Contigs     | contigs                        | Pvalue statistics from dekupl-run                                                                  |
| meanA                   | Float | Contigs     | contigs                        | Average counts for condition A                                                                     |
| meanB                   | Float | Contigs     | contigs                        | Average counts for condition B                                                                     |
| log2FC                  | Float | Contigs     | contigs                        | Log2 Fold-Change between samples in condition A vs B                                               |
| is_mapped               | Bool  | Bam         | bam                            | Is the contig mapped to at least one location                                                      |
| line_in_sam             | Int   | Bam         | bam                            | Line number in the SAM/BAM                                                                         |
| chromosome              | Str   | Bam         | bam                            | Chromosome                                                                                         |
| start                   | Int   | Bam         | bam                            | Begining of the alignment on the reference                                                         |
| end                     | Int   | Bam         | bam                            | End of the alignment on the reference                                                              |
| strand                  | Char  | Bam         | bam                            | Strand of the alignment (+/-). set to NA is unstranded data.                                       |
| cigar                   | Str   | Bam         | bam                            | CIGAR string from the SAM alignment.                                                               |
| nb_insertion            | Int   | Bam         | bam                            | Number of insertions in the alignment (from cigar)                                                 |
| nb_deletion             | Int   | Bam         | bam                            | Number of deletions in the alignment (from cigar)                                                  |
| nb_splice               | Int   | Bam         | bam                            | Number of splices in the alignment (from cigar)                                                    |
| nb_snv                  | Bool  | Bam         | bam                            | Number of SNV in the contigs (computed as the number of mismatches minus indels)                   |
| clipped_3p              | Bool  | Bam         | bam                            | Number of clipped bases (soft/hard) from 3prim contig                                              |
| clipped_5p              | Bool  | Bam         | bam                            | Number of clipped bases (soft/hard) from 5prim contig                                              |
| is_clipped_3p           | Bool  | Bam         | bam                            | True if contig's 3prim is soft/hard clipped (from cigar)                                           |
| is_clipped_5p           | Bool  | Bam         | bam                            | True if contig's 5prim is soft/hard clipped (from cigar)                                           |
| query_cover             | Float | Bam         | bam                            | Fraction of the query that have been aligned to the reference                                      |
| alignment_identity      | Float | Bam         | bam                            | Fraction of exact match over the query alignment length (splices do not count)                     |
| nb_hit                  | Int   | Bam         | bam                            | Number of alignment given for the contig (NH field)                                                |
| nb_mismatches           | Int   | Bam         | bam                            | Number of mismatches in the alignment (NM field)                                                   |
| is_chimeric             | Bool  | ChimericRNA | Chimeric.out.junction (STAR)   | The contig contains a chimeric junctions                                                           |
| chimeric_junctions      | Str   | ChimericRNA | Chimeric.out.junction (STAR)   | List of chimeric junctions (format: `chr1:pos1:strand1|chr2:pos2:trand`)                           |
| gene_id                 | Str   | Annotations | gff                            | Overlapping gene ID (from GFF ID field). On the same strand if --stranded option.                  |
| gene_symbol             | Str   | Annotations | gff                            | Overlapping gene symbol (from GFF Name field). On the same strand if --stranded option.            |
| gene_strand             | Char  | Annotations | gff                            | Overlapping gene strand (+/-). On the same strand if --stranded option.                            |
| gene_biotype            | Str   | Annotations | gff                            | Overlapping gene biotype (from GFF biotype field). On the same strand if --stranded option.        |
| as_gene_id              | Str   | Annotations | gff                            | Overlapping antisense gene ID (from GFF ID field). Only defined if --stranded option.              |
| as_gene_symbol          | Str   | Annotations | gff                            | Overlapping antisense gene symbol (from GFF Name field). Only defined if --stranded option.        |
| as_gene_strand          | Char  | Annotations | gff                            | Overlapping antisense gene strand (+/-). Only defined if --stranded option.                        |
| as_gene_biotype         | Str   | Annotations | gff                            | Overlapping antisense gene biotype (from GFF biotype field). Only defined if --stranded option.    |
| upstream_gene_id        | Str   | Annotations | gff                            | Nearest downstream gene ID. Same strand if --stranded option, downstream otherwise.                |
| upstream_gene_symbol    | Str   | Annotations | gff                            | Nearest downstream gene symbol. Same strand if --stranded option, downstream otherwise.            |
| upstream_gene_strand    | Str   | Annotations | gff                            | Nearest downstream gene strand (+/-).                                                              |
| upstream_gene_biotype   | Str   | Annotations | gff                            | Nearest downstream gene biotype.                                                                   |
| downstream_gene_id      | Str   | Annotations | gff                            | Nearest upstream gene ID. Same strand if --stranded option, upstream otherwise.                    |
| downstream_gene_symbol  | Str   | Annotations | gff                            | Nearest upstream gene symbol. Same strand if --stranded option, upstream otherwise.                |
| downstream_gene_strand  | Str   | Annotations | gff                            | Nearest upstream gene strand (+/-).                                                                |
| downstream_gene_biotype | Str   | Annotations | gff                            | Nearest upstream gene biotype.                                                                     |
| exonic                  | Bool  | Annotations | gff                            | Overlap between an exon and the contig. Same strand if --stranded option, both strand otherwise.   |
| intronic                | Bool  | Annotations | gff                            | Overlap between an intron and the contig. Same strand if --stranded option, both strand otherwise. |
| gene_is_diff            | Bool  | DEG         | DEGs                           | The main annotated gene (gene_id ontology) is differantially expressed                             |
| du_pvalue               | Float | Switches    | DEGs                           | P-value for contig differential usage                                                              |
| du_stats                | Float | Switches    | gene_counts, sample_conditions | Differential usage statistic                                                                       |

## Dev environnement

For developpement, `git clone` this repository and enter the project folder.

Then, add the local dir to the PERL5LIB env var to use the modules locally.

```
export PERL5LIB=$PWD/lib:$PERL5LIB
```
