use strict;
use warnings;

use Test::More tests => 5;

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

  my ($three_prim_res,$dist_three) = $interval_query->fetchNearest3prim(
    DEkupl::GenomicInterval->new(
      chr => "1",
      start => 1,
      end => 3,
      strand => '+'
    )
  );

  is($three_prim_res->id, 'TOTO');
  is($dist_three,9);

  my ($five_prim_res,$dist_five) = $interval_query->fetchNearest3prim(
    DEkupl::GenomicInterval->new(
      chr => "1",
      start => 100,
      end => 200,
      strand => '+'
    )
  );

  is($three_prim_res->id, 'TOTO');

  my ($three_prim_res_rev,$three_prim_res_rev_dist) = $interval_query->fetchNearest3prim(
    DEkupl::GenomicInterval->new(
      chr => "1",
      start => 1,
      end => 3,
      strand => '-'
    )
  );

  is($three_prim_res_rev, undef);
  
}