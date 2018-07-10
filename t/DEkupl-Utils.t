use strict;
use warnings;

use Test::More tests => 4;
use DEkupl::Utils;

{
  my $id = DEkupl::Utils::parseEnsemblID('gene:ENSG00000186092');
  is($id,'ENSG00000186092');
  my ($type,$id2) = DEkupl::Utils::parseEnsemblID('gene:ENSG00000186092');
  is($type,'gene')
}

{
  my $id = DEkupl::Utils::getAtomicGeneID('ENSG00000186092.12');
  is($id,'ENSG00000186092');
  $id = DEkupl::Utils::getAtomicGeneID('ENSG00000186092');
  is($id,'ENSG00000186092');
}

# check dependencies versions
{
  DEkupl::Utils::checkSamtoolsVersion();
  DEkupl::Utils::checkGSNAPVersion();
}