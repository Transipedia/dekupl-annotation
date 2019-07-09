use strict;
use warnings;

use Test::More tests => 15;
use DEkupl::IntervalQuery;
use DEkupl::GenomicInterval;
use DEkupl::Annotations;
use DEkupl::Annotations::Exon;
use DEkupl::Annotations::Gene;
use DEkupl::Analyzer::Annotations;
use Inline::Files 0.68;
use File::Temp;
use Test::Exception;

my $query = DEkupl::GenomicInterval->new(
    'chr'   => '12',
    'start' => 9,
    'end'   => 22,
);

my $query2 = DEkupl::GenomicInterval->new(
    'chr'   => '12',
    'start' => 12,
    'end'   => 18,
);

my $geneA = DEkupl::Annotations::Gene->new(
  'chr'   => '12',
  'id' => 'geneA',
);

my $geneA_exon1 = DEkupl::Annotations::Exon->new(
  'chr' => '12',
  'start' => 1,
  'end' => 4,
  'gene' => $geneA,
);

my $geneA_exon2 = DEkupl::Annotations::Exon->new(
  'chr' => '12',
  'start' => 42,
  'end' => 58,
);

my $geneB = DEkupl::Annotations::Gene->new(
  'chr'   => '12',
  'id' => 'geneB',
);

my $geneB_exon1 = DEkupl::Annotations::Exon->new(
  'chr' => '12',
  'start' => 7,
  'end' => 25,
  'gene' => $geneB,
);

my $geneC = DEkupl::Annotations::Gene->new(
  'chr'   => '12',
  'id' => 'geneC',
);

my $geneC_exon1 = DEkupl::Annotations::Exon->new(
  'chr' => '12',
  'start' => 11,
  'end' => 32,
  'gene' => $geneC,
);

{
  my @results = ($geneA);

  my $candidate = DEkupl::Analyzer::Annotations::_selectBestCandidate(\@results,$query);

  is($candidate->{gene}->id, 'geneA');
  is($candidate->{is_exonic}, 0);
  is($candidate->{is_intronic}, 1);
}

# Rule 1 : exonic over non-exonic genes
{
  my @results = ($geneA,$geneB,$geneB_exon1);

  my $candidate = DEkupl::Analyzer::Annotations::_selectBestCandidate(\@results,$query);

  is($candidate->{gene}->id, 'geneB');
  is($candidate->{is_exonic}, 1);
  is($candidate->{is_intronic}, 0);
}

# Rule 2 : longer exonic overlap 
{
  my @results = ($geneA,$geneC,$geneC_exon1,$geneB,$geneB_exon1);

  my $candidate = DEkupl::Analyzer::Annotations::_selectBestCandidate(\@results,$query);

  is($candidate->{gene}->id, 'geneB');
  is($candidate->{is_exonic}, 1);
  is($candidate->{is_intronic}, 0);
}

# Rule 3 : longer gene if overlapping interval 
{
  my @results = ($geneA,$geneC,$geneC_exon1,$geneB,$geneB_exon1);

  my $candidate = DEkupl::Analyzer::Annotations::_selectBestCandidate(\@results,$query2);

  is($candidate->{gene}->id, 'geneC');
  is($candidate->{is_exonic}, 1);
  is($candidate->{is_intronic}, 0);
}

# Test on real GFF file
{
  my ($fh, $gff_file) = File::Temp::tempfile( SUFFIX => '.gff', UNLINK => 0);
  while(<GFF>) {print $fh $_;}
  $fh->close;

  my $annotations = DEkupl::Annotations->new();
  $annotations->loadFromGFF($gff_file);

  my $interval_query = DEkupl::IntervalQuery->new();
  $interval_query->loadAnnotations($annotations);

  # chr11:5,225,504-5,225,560
  my $query = DEkupl::GenomicInterval->new(
      'chr'   => 'chr11',
      'start' => 5225504,
      'end'   => 5225560,
      'strand' => '-',
  );

  my $results   = $interval_query->fetchByRegion($query);
  my $candidate = DEkupl::Analyzer::Annotations::_selectBestCandidate($results, $query);

  is($candidate->{gene}->id, 'ENSG00000244734');
  is($candidate->{is_exonic}, 1);
  is($candidate->{is_intronic}, 0);
}

__GFF__
##gff-version 3
#description: evidence-based annotation of the human genome (GRCh38), version 31 (Ensembl 97)
#provider: GENCODE
#contact: gencode-help@ebi.ac.uk
#format: gff3
#date: 2019-06-27
##sequence-region chr1 1 248956422
chr11	HAVANA	gene	5225464	5229395	.	-	.	ID=ENSG00000244734.4;gene_id=ENSG00000244734.4;gene_type=protein_coding;gene_name=HBB;level=2;hgnc_id=HGNC:4827;havana_gene=OTTHUMG00000066678.8
chr11	HAVANA	transcript	5225464	5227071	.	-	.	ID=ENST00000335295.4;Parent=ENSG00000244734.4;gene_id=ENSG00000244734.4;transcript_id=ENST00000335295.4;gene_type=protein_coding;gene_name=HBB;transcript_type=protein_coding;transcript_name=HBB-201;level=2;protein_id=ENSP00000333994.3;transcript_support_level=1;hgnc_id=HGNC:4827;tag=CAGE_supported_TSS,basic,MANE_Select,appris_principal_1,CCDS;ccdsid=CCDS7753.1;havana_gene=OTTHUMG00000066678.8;havana_transcript=OTTHUMT00000495006.2
chr11	HAVANA	exon	5226930	5227071	.	-	.	ID=exon:ENST00000335295.4:1;Parent=ENST00000335295.4;gene_id=ENSG00000244734.4;transcript_id=ENST00000335295.4;gene_type=protein_coding;gene_name=HBB;transcript_type=protein_coding;transcript_name=HBB-201;exon_number=1;exon_id=ENSE00001829867.2;level=2;protein_id=ENSP00000333994.3;transcript_support_level=1;hgnc_id=HGNC:4827;tag=CAGE_supported_TSS,basic,MANE_Select,appris_principal_1,CCDS;ccdsid=CCDS7753.1;havana_gene=OTTHUMG00000066678.8;havana_transcript=OTTHUMT00000495006.2
chr11	HAVANA	CDS	5226930	5227021	.	-	0	ID=CDS:ENST00000335295.4;Parent=ENST00000335295.4;gene_id=ENSG00000244734.4;transcript_id=ENST00000335295.4;gene_type=protein_coding;gene_name=HBB;transcript_type=protein_coding;transcript_name=HBB-201;exon_number=1;exon_id=ENSE00001829867.2;level=2;protein_id=ENSP00000333994.3;transcript_support_level=1;hgnc_id=HGNC:4827;tag=CAGE_supported_TSS,basic,MANE_Select,appris_principal_1,CCDS;ccdsid=CCDS7753.1;havana_gene=OTTHUMG00000066678.8;havana_transcript=OTTHUMT00000495006.2
chr11	HAVANA	start_codon	5227019	5227021	.	-	0	ID=start_codon:ENST00000335295.4;Parent=ENST00000335295.4;gene_id=ENSG00000244734.4;transcript_id=ENST00000335295.4;gene_type=protein_coding;gene_name=HBB;transcript_type=protein_coding;transcript_name=HBB-201;exon_number=1;exon_id=ENSE00001829867.2;level=2;protein_id=ENSP00000333994.3;transcript_support_level=1;hgnc_id=HGNC:4827;tag=CAGE_supported_TSS,basic,MANE_Select,appris_principal_1,CCDS;ccdsid=CCDS7753.1;havana_gene=OTTHUMG00000066678.8;havana_transcript=OTTHUMT00000495006.2
chr11	HAVANA	exon	5226577	5226799	.	-	.	ID=exon:ENST00000335295.4:2;Parent=ENST00000335295.4;gene_id=ENSG00000244734.4;transcript_id=ENST00000335295.4;gene_type=protein_coding;gene_name=HBB;transcript_type=protein_coding;transcript_name=HBB-201;exon_number=2;exon_id=ENSE00001057381.1;level=2;protein_id=ENSP00000333994.3;transcript_support_level=1;hgnc_id=HGNC:4827;tag=CAGE_supported_TSS,basic,MANE_Select,appris_principal_1,CCDS;ccdsid=CCDS7753.1;havana_gene=OTTHUMG00000066678.8;havana_transcript=OTTHUMT00000495006.2
chr11	HAVANA	CDS	5226577	5226799	.	-	1	ID=CDS:ENST00000335295.4;Parent=ENST00000335295.4;gene_id=ENSG00000244734.4;transcript_id=ENST00000335295.4;gene_type=protein_coding;gene_name=HBB;transcript_type=protein_coding;transcript_name=HBB-201;exon_number=2;exon_id=ENSE00001057381.1;level=2;protein_id=ENSP00000333994.3;transcript_support_level=1;hgnc_id=HGNC:4827;tag=CAGE_supported_TSS,basic,MANE_Select,appris_principal_1,CCDS;ccdsid=CCDS7753.1;havana_gene=OTTHUMG00000066678.8;havana_transcript=OTTHUMT00000495006.2
chr11	HAVANA	exon	5225464	5225726	.	-	.	ID=exon:ENST00000335295.4:3;Parent=ENST00000335295.4;gene_id=ENSG00000244734.4;transcript_id=ENST00000335295.4;gene_type=protein_coding;gene_name=HBB;transcript_type=protein_coding;transcript_name=HBB-201;exon_number=3;exon_id=ENSE00001600613.2;level=2;protein_id=ENSP00000333994.3;transcript_support_level=1;hgnc_id=HGNC:4827;tag=CAGE_supported_TSS,basic,MANE_Select,appris_principal_1,CCDS;ccdsid=CCDS7753.1;havana_gene=OTTHUMG00000066678.8;havana_transcript=OTTHUMT00000495006.2
chr11	HAVANA	CDS	5225598	5225726	.	-	0	ID=CDS:ENST00000335295.4;Parent=ENST00000335295.4;gene_id=ENSG00000244734.4;transcript_id=ENST00000335295.4;gene_type=protein_coding;gene_name=HBB;transcript_type=protein_coding;transcript_name=HBB-201;exon_number=3;exon_id=ENSE00001600613.2;level=2;protein_id=ENSP00000333994.3;transcript_support_level=1;hgnc_id=HGNC:4827;tag=CAGE_supported_TSS,basic,MANE_Select,appris_principal_1,CCDS;ccdsid=CCDS7753.1;havana_gene=OTTHUMG00000066678.8;havana_transcript=OTTHUMT00000495006.2
chr11	HAVANA	stop_codon	5225598	5225600	.	-	0	ID=stop_codon:ENST00000335295.4;Parent=ENST00000335295.4;gene_id=ENSG00000244734.4;transcript_id=ENST00000335295.4;gene_type=protein_coding;gene_name=HBB;transcript_type=protein_coding;transcript_name=HBB-201;exon_number=3;exon_id=ENSE00001600613.2;level=2;protein_id=ENSP00000333994.3;transcript_support_level=1;hgnc_id=HGNC:4827;tag=CAGE_supported_TSS,basic,MANE_Select,appris_principal_1,CCDS;ccdsid=CCDS7753.1;havana_gene=OTTHUMG00000066678.8;havana_transcript=OTTHUMT00000495006.2
chr11	HAVANA	five_prime_UTR	5227022	5227071	.	-	.	ID=UTR5:ENST00000335295.4;Parent=ENST00000335295.4;gene_id=ENSG00000244734.4;transcript_id=ENST00000335295.4;gene_type=protein_coding;gene_name=HBB;transcript_type=protein_coding;transcript_name=HBB-201;exon_number=1;exon_id=ENSE00001829867.2;level=2;protein_id=ENSP00000333994.3;transcript_support_level=1;hgnc_id=HGNC:4827;tag=CAGE_supported_TSS,basic,MANE_Select,appris_principal_1,CCDS;ccdsid=CCDS7753.1;havana_gene=OTTHUMG00000066678.8;havana_transcript=OTTHUMT00000495006.2
chr11	HAVANA	three_prime_UTR	5225464	5225597	.	-	.	ID=UTR3:ENST00000335295.4;Parent=ENST00000335295.4;gene_id=ENSG00000244734.4;transcript_id=ENST00000335295.4;gene_type=protein_coding;gene_name=HBB;transcript_type=protein_coding;transcript_name=HBB-201;exon_number=3;exon_id=ENSE00001600613.2;level=2;protein_id=ENSP00000333994.3;transcript_support_level=1;hgnc_id=HGNC:4827;tag=CAGE_supported_TSS,basic,MANE_Select,appris_principal_1,CCDS;ccdsid=CCDS7753.1;havana_gene=OTTHUMG00000066678.8;havana_transcript=OTTHUMT00000495006.2