package DEkupl::Genome;
# ABSTRACT: A genome (loaded from FASTA)

use Moose;
use DEkupl::Utils;

has 'references_length' => (
  traits => ['Hash'],
  is => 'ro',
  isa => 'HashRef[Int]',
  default => sub { {} },
  handles => {
    getReferenceLength  => 'get',
    allReferencesLength => 'values',
    allReferences       => 'keys',
    hasReference        => 'exists',
  },
);

has 'fasta_file' => (
  is => 'ro',
  isa => 'Str',
  required => 1,
  trigger => sub {
    my $self = shift;
    my $seq_it = DEkupl::Utils::fastaFileIterator($self->fasta_file);
    while(my $entry = $seq_it->()) {
      my ($ref_name) = $entry->{name} =~ /^(\S+)/;
      $self->references_length->{$ref_name} = length $entry->{seq};
    }
  }
);

sub genomeLength {
  my $self = shift;
  my $genome_length = 0;
  foreach my $l ($self->allReferencesLength) {
    $genome_length += $l;
  }
  return $genome_length;
}

no Moose;
__PACKAGE__->meta->make_immutable;