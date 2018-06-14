use strict;
use warnings;

use Test::More tests => 1;
use DEkupl::ContigsDB;
use File::Temp qw/ tempdir /;

my $tempdir = tempdir(CLEANUP => 1);
my $contigs_db = DEkupl::ContigsDB->new(db_folder => $tempdir);


# Save a contig into the DB
{
  my $contig = {
  tag => 'AAA',
  info => 'TOTO'
  };
  $contigs_db->saveContig($contig);
}

# Load the contig
{
  my $contig = $contigs_db->loadContig('AAA');
  is($contig->{info}, "TOTO");
}

