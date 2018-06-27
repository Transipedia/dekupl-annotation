# DE-kupl Annotation

DE-kupl annotation is part of the DE-kupl package, and performs annotations of DE contigs identified by DE-kupl.

## Installation

Dependencies are cpan-minus (aka cpanm) and Dist::Zilla :

```
apt-get install cpanminus libdist-zilla-perl
```

The install rankvar with dzil and cpanm.

```
git clone https://gitlab.seq.one/workset/rankvar2.git && cd rankvar2
dzil install --install-command 'cpanm .'
```

## Ontology
| Term               | Type  | Analyzer    | File Source                    | Description                                                                                        |
| ------------------ | ----- | ----------- | ------------------------------ | -------------------------------------------------------------------------------------------------- |
| tag                | Str   | Contigs     | contigs                        | kmer of reference                                                                                  |
| nb_merged_kmers    | Int   | Contigs     | contigs                        | Number of k-mer merged to produce the contig                                                       |
| contig             | Str   | Contigs     | contigs                        | Sequence of the contig                                                                             |
| contig_size        | Int   | Contigs     | contigs                        | Size (in nucleotides) of the contig                                                                |
| pvalue             | Float | Contigs     | contigs                        | Pvalue statistics from dekupl-run                                                                  |
| meanA              | Float | Contigs     | contigs                        | Average counts for condition A                                                                     |
| meanB              | Float | Contigs     | contigs                        | Average counts for condition B                                                                     |
| log2FC             | Float | Contigs     | contigs                        | Log2 Fold-Change between samples in condition A vs B                                               |
| is_mapped          | Bool  | Bam         | bam                            | Is the contig mapped to at least one location                                                      |
| line_in_sam        | Int   | Bam         | bam                            | Line number in the SAM/BAM                                                                         |
| chromosome         | Str   | Bam         | bam                            | Chromosome                                                                                         |
| start              | Int   | Bam         | bam                            | Begining of the alignment on the reference                                                         |
| end                | Int   | Bam         | bam                            | End of the alignment on the reference                                                              |
| strand             | Char  | Bam         | bam                            | Strand of the alignment (+/-). set to NA is unstranded data.                                       |
| nb_insertion       | Int   | Bam         | bam                            | Number of insertions in the alignment (from cigar)                                                 |
| nb_deletion        | Int   | Bam         | bam                            | Number of deletions in the alignment (from cigar)                                                  |
| nb_splice          | Int   | Bam         | bam                            | Number of splices in the alignment (from cigar)                                                    |
| clipped_3p         | Bool  | Bam         | bam                            | True if contig's 3prim is soft/hard clipped (from cigar)                                           |
| clipped_5p         | Bool  | Bam         | bam                            | True if contig's 5prim is soft/hard clipped (from cigar)                                           |
| query_cover        | Float | Bam         | bam                            | Fraction of the query that have been aligned to the reference                                      |
| alignment_identity | Float | Bam         | bam                            | Fraction of exact match over the query alignment length (splices do not count)                     |
| nb_hit             | Int   | Bam         | bam                            | Number of alignment given for the contig (NH field)                                                |
| nb_mismatches      | Int   | Bam         | bam                            | Number of mismatches in the alignment (NM field)                                                   |
| gene_id            | Str   | Annotations | gff                            | Overlapping gene ID (from GFF ID field). On the same strand if --stranded option.                  |
| gene_symbol        | Str   | Annotations | gff                            | Overlapping gene symbol (from GFF Name field). On the same strand if --stranded option.            |
| gene_strand        | Char  | Annotations | gff                            | Overlapping gene strand (+/-). On the same strand if --stranded option.                            |
| gene_biotype       | Str   | Annotations | gff                            | Overlapping gene biotype (from GFF biotype field). On the same strand if --stranded option.        |
| as_gene_id         | Str   | Annotations | gff                            | Overlapping antisense gene ID (from GFF ID field). Only defined if --stranded option.              |
| as_gene_symbol     | Str   | Annotations | gff                            | Overlapping antisense gene symbol (from GFF Name field). Only defined if --stranded option.        |
| as_gene_strand     | Char  | Annotations | gff                            | Overlapping antisense gene strand (+/-). Only defined if --stranded option.                        |
| as_gene_biotype    | Str   | Annotations | gff                            | Overlapping antisense gene biotype (from GFF biotype field). Only defined if --stranded option.    |
| 5p_gene_id         | Str   | Annotations | gff                            | Nearest 5’ gene ID. Same strand if --stranded option, downstream otherwise.                       |
| 5p_gene_symbol     | Str   | Annotations | gff                            | Nearest 5’ gene symbol. Same strand if --stranded option, downstream otherwise.                   |
| 5p_gene_strand     | Str   | Annotations | gff                            | Nearest 5’ gene gene strand (+/-).                                                                |
| 3p_gene_id         | Str   | Annotations | gff                            | Nearest 3’ gene ID. Same strand if --stranded option, upstream otherwise.                         |
| 3p_gene_symbol     | Str   | Annotations | gff                            | Nearest 3’ gene symbol. Same strand if --stranded option, upstream otherwise.                     |
| 3p_gene_strand     | Str   | Annotations | gff                            | Nearest 3’ gene gene strand (+/-).                                                                |
| exonic             | Bool  | Annotations | gff                            | Overlap between an exon and the contig. Same strand if --stranded option, both strand otherwise.   |
| intronic           | Bool  | Annotations | gff                            | Overlap between an intron and the contig. Same strand if --stranded option, both strand otherwise. |
| gene_is_diff       | Bool  | DEG         | DEGs                           | The main annotated gene (gene_id ontology) is differantially expressed                             |
| du_pvalue          | Float | Switches    | DEGs                           | P-value for contig differential usage                                                              |
| du_stats           | Float | Switches    | gene_counts, sample_conditions | Differential usage statistic                                                                       |

## Dev environnement

For developpement, `git clone` this repository and enter the project folder.

Then, add the local dir to the PERL5LIB env var to use the modules locally.

```
export PERL5LIB=$PWD/lib:$PERL5LIB
```

You are ready to code!