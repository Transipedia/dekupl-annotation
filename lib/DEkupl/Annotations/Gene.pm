package DEkupl::Annotations::Gene;
# ABSTRACT: An Gene

use Moose;
use Moose::Util::TypeConstraints;

extends 'DEkupl::GenomicInterval';

use DEkupl::Utils;
use DEkupl::Annotations::Exon;

has exons_hash => (
  traits => ['Hash'],
  is     => 'ro',
  isa    => 'HashRef[DEkupl::Annotations::Exon]',
  default => sub { {} },
  handles => {
    allExons  => 'values',
    getExon   => 'get',
    nbExons   => 'count',
  },
);

has id  => (is => 'ro', isa => 'Str', required => 1);
has symbol => (is => 'ro', isa => 'Str');

# Compute start on the fly
has '+start' => (
  lazy  => 1,
  default => sub {
    my $self = shift;
    my $start = undef;
    foreach my $exon ($self->allExons) {
      if(!defined $start || $exon->start < $start) {
        $start = $exon->start;
      }
    }
    return $start;
  },
);

has '+end' => (
  lazy    => 1,
  default => sub {
    my $self = shift;
    my $end = undef;
    foreach my $exon ($self->allExons) {
      if(!defined $end || $exon->end > $end) {
        $end = $exon->end;
      }
    }
    return $end;
  },
);

sub addExon {
  my $self = shift;
  my $new_exon = shift;
  my $exon_key = _getExonKey($new_exon);
  my $exon     = $self->getExon($exon_key);
  # If there is already an exon with this key
  # we merge their transcript ids
  if(defined $exon) {
    foreach my $transcript ($new_exon->allTranscripts) {
      $exon->addTranscript($transcript);
    }
    return $exon;
  } else {
    $self->exons_hash->{$exon_key} = $new_exon;
    return $new_exon;
  }
}

sub sortedExons {
  my $self = shift;
  return sort { $a->start <=> $b->start } $self->allExons;
}


sub _getExonKey {
  my $exon = shift;
  return $exon->start.",".$exon->end;
}

no Moose;
__PACKAGE__->meta->make_immutable;

__END__
=head1 ACCESSORS

=head2 chr => 'String'

Getter for the gene's chr

=head2 strand => ['+'|'-']

Getter for the gene's strand

=head2 start => 'Int'

Getter for the gene's start position (retrieved 'on the fly' from the exons')

=head2 end => 'Int'

Getter for the gene's end position (retrieved 'on the fly' from the exons')

=head1 METHODS

=head2 new

  Arg [chr]     : 'Str'   - gene's chromosome
  Arg [strand]  : '(+|-)' - gene's strand
  Arg [id]      : 'Str'   - gene's id

Create a new 'DEkupl::Annotations::Gene' object

=head2 addExon('DEkupl::Annotations::Exon')

Add a new exon for this gene

=head2 nbExons => 'Int'

Return the number of exons

=head2 allExons() => Array['DEkupl::Annotations::Exon']

Return all gene's exons

=head2 sortedExons => Array['DEkupl::Annotations::Exon']

Retturn all gene's exons sorted by 'start' position