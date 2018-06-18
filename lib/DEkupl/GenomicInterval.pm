package DEkupl::GenomicInterval;
# ABSTRACT: Base class for all genomic based object

use Moose;
use DEkupl::Utils;

has chr => (
  isa => 'Str',
  is  => 'rw',
  required => 1,
);

has start => (
  isa => 'Int',
  is  => 'rw',
  required  => 1,
);

has end => (
  isa => 'Int',
  is  => 'rw',
  required  => 1,
  # Verify that 'end' is greater than start
  trigger   => sub {
    my $self = shift;
    if($self->end < $self->start) {
      die "End position (".$self->end.") must be greater (or equal) to start (".$self->start.")";
    }
  },
);

has strand => (
  isa => 'Strand',
  is  => 'rw',
  required => 1,
  default  => '+',
);

sub length {
  my $self = shift;
  return $self->end - $self->start + 1;
}

no Moose;
__PACKAGE__->meta->make_immutable;