package DEkupl::Analyzer::ChimericRNA;
# ABSTRACT: Append Chimeric RNA information from STAR output

use Moose;

use DEkupl::Utils;

with 'DEkupl::Analyzer';

# Chimeric.out.junction file from STAR output
has 'chimeric_file' => (
  is => 'ro',
  isa => 'Str',
  required => 1,
);

my @columns = ('is_chimeric', 'chimeric_junctions');

sub BUILD {
  my $self = shift;
  my $fh = DEkupl::Utils::getReadingFileHandle($self->chimeric_file);

  my %gene_is_diff;

  # Keep counts for verbose prints
  my $nb_chimeric_junctions = 0;

  while(<$fh>) {
    chomp;
    # column 1: chromosome of the donor
    # column 2: first base of the intron of the donor (1-based)
    # column 3: strand of the donor
    # column 4: chromosome of the acceptor
    # column 5: first base of the intron of the acceptor (1-based)
    # column 6: strand of the acceptor
    # column 7: junction type: -1=encompassing junction (between the mates), 1=GT/AG, 2=CT/AC
    # column 8: repeat length to the left of the junction
    # column 9: repeat length to the right of the junction
    # Columns 10-14 describe the alignments of the two chimeric segments, it is SAM like. Alignments are given with respect to the (+) strand
    # column 10: read name
    # column 11: first base of the first segment (on the + strand) column 12: CIGAR of the first segment
    # column 13: first base of the second segment
    # column 14: CIGAR of the second segment
    my ($chr1, $pos1, $strand1, $chr2, $pos2, $strand2, $type, $repeat_length_left, $repeat_length_right, $tag) = split "\t", $_;

    # Load contig from the DB
    my $contig = $self->contigs_db->loadContig($tag);
    # The contig might have been deleted by adapter removal
    next if !defined $contig;

    $nb_chimeric_junctions++;
    
    # Set values
    $contig->{is_chimeric} = DEkupl::Utils::booleanEncoding(1);
    push @{$contig->{chimeric_junctions}}, "$chr1:$pos1:$strand1|$chr2:$pos2:$strand2";

    # Save contig
    $self->contigs_db->saveContig($contig);
  }

  $self->verboseLog("$nb_chimeric_junctions chimeric junctions found in STAR output");
}

sub getHeaders {
  my $self = shift;
  return @columns;
}

sub getValues {
  my $self = shift;
  my $contig = shift;
  my @values;
  $contig->{is_chimeric} = DEkupl::Utils::booleanEncoding(0) if !defined $contig->{is_chimeric};
  my $chimeric_junctions = defined $contig->{chimeric_junctions}? join(',', @{$contig->{chimeric_junctions}}) : $DEkupl::Utils::NA_value;
  return ($contig->{is_chimeric}, $chimeric_junctions);
}


no Moose;
__PACKAGE__->meta->make_immutable;