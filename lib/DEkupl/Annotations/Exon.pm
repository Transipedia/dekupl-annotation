package DEkupl::Annotations::Exon;
# ABSTRACT: An Exon

use Moose;

extends 'DEkupl::GenomicInterval';

use DEkupl::Annotations::Gene;

has gene  => (
  is        => 'ro',
  isa       => 'DEkupl::Annotations::Gene',
  required  => 1,
  trigger   => sub {
    my $self = shift;
    $self->gene->addExon($self);
  }
);

# Use gene value for chr
has '+chr' => (
  lazy    => 1,
  default => sub {
    my $self = shift;
    return $self->gene->chr;
  },
);

# Use gene value for strand
has '+strand' => (
  lazy    => 1,
  default => sub {
    my $self = shift;
    return $self->gene->strand;
  },
);

has transcripts => (
  traits  => ['Array'],
  is      => 'ro',
  isa     => 'ArrayRef[Str]',
  default => sub { [] },
  handles => {
    addTranscript   => 'push',
    allTranscripts  => 'elements',
  },
);

# Remove a transcript given its id
sub removeTranscript($) {
  my $self = shift;
  my $transcript  = shift;
  my @transcripts = grep { $_ ne $transcript } @{$self->transcripts};
  $self->transcripts(\@transcripts);
}

no Moose;
__PACKAGE__->meta->make_immutable;

__END__

=head1 ACCESSORS

=head2 chr => 'String'

Getter for the exon chr (retrieved from the gene)

=head2 strand => ['+'|'-']

Getter for the exon strand (retrieved from the gene)

=head2 start => 'Int'

Getter/setter for the exon's start

=head2 end => 'Int'

Getter/setter for the exon's end. End should be greater or equal to start

=head1 METHODS

=head2 new

  Arg [start] : 'Int' - exon's start position
  Arg [end]   : 'Int' - exon's end position
  Arg [gene]  : 'DEkupl::Annotations::Gene' - exon's gene

Create a new 'DEkupl::Annotations::Exon' object

=head2 addTranscript('transcript_id')

Add a new transcript id for this exon

=head2 removeTranscript('transcript_id')

Remove a transcript_id given its name

=head2 allTranscripts() => Array['transcript_id']

Return all transcript ids