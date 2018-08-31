use strict;
use warnings;

use Test::More tests => 4;
use DEkupl::Analyzer::Switches;
use Inline::Files 0.68;
use File::Temp;
use Test::Exception;

{
  my $valid_sample_cond_file = new File::Temp( SUFFIX => '.tsv', UNLINK => 1);
  while(<VALIDSAMPLECOND>) {print $valid_sample_cond_file $_;}
  close $valid_sample_cond_file;
  my @sample_names = qw(SRR2966453 SRR2966454 SRR2966455 SRR2966474 SRR2966475 SRR2966476);
  my $ret = DEkupl::Analyzer::Switches::verifySampleConditionsFile($valid_sample_cond_file,\@sample_names);
  is($ret, 0);
  @sample_names = qw(SRR2966453 SRR2966454 SRR2966455 SRR2966474 SRR2966475);
  
  dies_ok {DEkupl::Analyzer::Switches::verifySampleConditionsFile($valid_sample_cond_file,\@sample_names)} 'Unknown sample name';
}

{
  my $sample_cond_file = new File::Temp( SUFFIX => '.tsv', UNLINK => 1);
  while(<SPACESEPSAMPLECOND>) {print $sample_cond_file $_;}
  close $sample_cond_file;

  my @sample_names = qw(SRR2966453 SRR2966454 SRR2966455 SRR2966474 SRR2966475 SRR2966476);

  dies_ok { DEkupl::Analyzer::Switches::verifySampleConditionsFile($sample_cond_file,\@sample_names) } 'Malformed file';
}


{
  my $sample_cond_file = new File::Temp( SUFFIX => '.tsv', UNLINK => 1);
  while(<MISSINGSAMPEPSAMPLECOND>) {print $sample_cond_file $_;}
  close $sample_cond_file;

  my @sample_names = qw(SRR2966453 SRR2966454 SRR2966455 SRR2966474 SRR2966475 SRR2966476);

  dies_ok { DEkupl::Analyzer::Switches::verifySampleConditionsFile($sample_cond_file,\@sample_names) } 'Missing sample';
}

__VALIDSAMPLECOND__
sample	condition	normalization_factor
SRR2966453	D0D1	1.26283604091483
SRR2966454	D0D1	1.09777764671068
SRR2966455	D0D1	1.04686452518516
SRR2966474	D6D7	1.00455068135207
SRR2966475	D6D7	1.03033580675732
SRR2966476	D6D7	1.08107691202423
__SPACESEPSAMPLECOND__
sample condition normalization_factor
SRR2966453 D0D1 1.26283604091483
SRR2966454 D0D1 1.09777764671068
SRR2966455 D0D1 1.04686452518516
SRR2966474 D6D7 1.00455068135207
SRR2966475 D6D7 1.03033580675732
SRR2966476 D6D7 1.08107691202423
__MISSINGSAMPEPSAMPLECOND__
sample condition normalization_factor
SRR2966453 D0D1 1.26283604091483
SRR2966454 D0D1 1.09777764671068
SRR2966474 D6D7 1.00455068135207
SRR2966475 D6D7 1.03033580675732
SRR2966476 D6D7 1.08107691202423
