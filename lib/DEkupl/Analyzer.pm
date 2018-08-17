package DEkupl::Analyzer;
# ABSTRACT: Base class for dekupl analyzers

use Moose::Role;

with 'DEkupl::Base';

has 'contigs_db' => (
  is => 'ro',
  isa => 'DEkupl::ContigsDB'
);

has 'is_stranded' => (
  is => 'ro',
  isa => 'Bool',
  required => 1,
);

# return header columns to be printed in the output file
requires 'getHeaders';

# return values for a given contig passed in argument
requires 'getValues';

1;