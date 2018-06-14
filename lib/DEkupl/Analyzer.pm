package DEkupl::Analyzer;
# ABSTRACT: Base class for dekupl analyzers

use Moose::Role;

has 'contigs_db' => (
  is => 'ro',
  isa => 'DEkupl::ContigsDB'
);

# return header columns to be printed in the output file
requires 'getHeaders';

# return values for a given contig passed in argument
requires 'getValues';

1;