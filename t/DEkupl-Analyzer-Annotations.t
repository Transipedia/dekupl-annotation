use strict;
use warnings;

use Test::More tests => 12;

use DEkupl::GenomicInterval;
use DEkupl::Annotations::Exon;
use DEkupl::Annotations::Gene;
use DEkupl::Analyzer::Annotations;

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
