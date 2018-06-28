use strict;
use warnings;

use Test::More tests => 10;

use DEkupl::Utils;
use DEkupl::Analyzer::BAM;

my $sam_line = "AACATTCAGGTTGTTTTTTTTTTTTTGTTTC	0	22	20127177	40	8S18M1I30M	*	0	0	AACTCAACATTCAGGTTGTTTTTTTTTTTTTGTTTCTAAGTTTTTGCCC	*	MD:Z:48	NH:i:1	HI:i:1	NM:i:1	SM:i:40	XQ:i:40	X2:i:0	XO:Z:UU	XG:Z:M";
my $parsed_sam_line = DEkupl::Utils::parseSAMLine($sam_line);

is(DEkupl::Analyzer::BAM::_getStrandFromFlag(22,1),'-');
is(DEkupl::Analyzer::BAM::_getStrandFromFlag(22),undef);

my $cig_stats = DEkupl::Analyzer::BAM::_computeStatsFromCigar($parsed_sam_line->{cigar});

# my %cig_stats = (
#     ref_aln_length    => 0, # Length of the alignment on the reference (from the first based aligned to the last, including deletion and splice)
#     query_aln_length  => 0, # Length of the alignment on the query (including deletion)
#     nb_match          => 0,
#     nb_insertion      => 0, # Number of insertion in the query (I cigar element)
#     nb_deletion       => 0, # Number of deletiong in the query (D cigar element)
#     nb_splice         => 0, # Number of splice in the query (N cigar element),
#     is_clipped_5p     => 0,
#     is_clipped_3p     => 0,
# );

is($cig_stats->{ref_aln_length}, 48);
is($cig_stats->{query_aln_length}, 49);
is($cig_stats->{nb_match}, 48);
is($cig_stats->{nb_insertion}, 1);
is($cig_stats->{nb_deletion}, 0);
is($cig_stats->{nb_splice}, 0);
is($cig_stats->{is_clipped_5p}, 1);
is($cig_stats->{is_clipped_3p}, 0);
