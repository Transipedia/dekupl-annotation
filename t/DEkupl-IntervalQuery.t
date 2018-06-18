use strict;
use warnings;

use Test::More tests => 1;

use DEkupl::Annotations::Exon;
use DEkupl::Annotations::Gene;

use DEkupl::IntervalQuery;

{
  my $gene = DEkupl::Annotations::Gene->new(chr => "1", strand => '+', id => "TOTO"); 

  my $exon = DEkupl::Annotations::Exon->new(start => 12, end => 26, gene => $gene); 

  my $interval_query = DEkupl::IntervalQuery->new();

  $interval_query->addInterval($gene);
  $interval_query->addInterval($exon);

  my $search = DEkupl::GenomicInterval->new(chr => "1", start => 7, end => 23, strand => '+');

  my $results = $interval_query->fetchByRegion($search);

  is(scalar @{$results}, 2);
}