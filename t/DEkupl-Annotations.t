use strict;
use warnings;

use Test::More tests => 13;
use DEkupl::Annotations;
use DEkupl::Annotations::Exon;
use DEkupl::Annotations::Gene;
use Inline::Files 0.68;
use File::Temp;

{
  my $gene = DEkupl::Annotations::Gene->new(chr => "1", strand => '+', id => "TOTO"); 

  is($gene->chr,"1");
  is($gene->strand,'+');
  is($gene->id,"TOTO");

  my $exon = DEkupl::Annotations::Exon->new(start => 12, end => 26, gene => $gene); 

  is($exon->start,12);
  is($exon->end,26);
  is($gene->allExons,1);

  is($gene->start,12);
  is($gene->end,26);
}

{
  # Create a temp GTF annotation file
  my $gtf_file = new File::Temp( SUFFIX => '.gtf', UNLINK => 1);
  while(<ANNOTATIONS>) {print $gtf_file $_;}
  close $gtf_file;

  # Create a annotation object and load the GTF
  my $annotations = DEkupl::Annotations->new();
  $annotations->loadFromFile($gtf_file->filename,'gtf');
  is(scalar $annotations->allGenes,2);
  my $geneA = $annotations->getGene('geneA');
  is($geneA->start,9);
  is($geneA->end,24);
  is(scalar $geneA->allExons, 2);
  $annotations->removeGene('geneB');
  is(scalar $annotations->allGenes,1);
}

__ANNOTATIONS__
1	StringTie	transcript	10	25	1000	+	.	gene_id "STRG.1"; transcript_id "STRG.1.1"; cov "105.907524"; FPKM "128.387177";
1	StringTie	exon	10	15	1000	+	.	gene_id "geneA"; transcript_id "geneA.1"; exon_number "1"; cov "42.406654";
1	StringTie	exon	20	25	1000	+	.	gene_id "geneA"; transcript_id "geneA.1"; exon_number "2"; cov "161.272552";
2	StringTie	exon	5	10	1000	-	.	gene_id "geneB"; transcript_id "geneB.1"; exon_number "2"; cov "161.272552";
2	StringTie	exon	15	20	1000	-	.	gene_id "geneB"; transcript_id "geneB.1"; exon_number "2"; cov "161.272552";
