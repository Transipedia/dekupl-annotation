package DEkupl::IntervalQuery;
# ABSTRACT: Store and query genomics intervals.

use Moose;
use Set::IntervalTree;

has 'interval_trees' => (
  is => 'ro',
  isa => 'HashRef[Set::IntervalTree]',
  default => sub { {} },
  traits => ['Hash'],
  handles => {
    _getIntervalTree => 'get',
    _setIntervalTree => 'set'
  },
);

sub addInterval {
  my $self = shift;
  my $genomic_interval = shift;

  my $interval_tree = $self->getIntervalTree($genomic_interval->chr,$genomic_interval->strand);

  # If there is no already existing IntervalTree for this ("chr","strand") pair
  if(!defined $interval_tree) {
    # We create a new one
    $interval_tree = Set::IntervalTree->new;

    # We add this new interval tree with the others
    $self->setIntervalTree($genomic_interval->chr,$genomic_interval->strand,$interval_tree);
  }

  # We insert the given interval in the IntervalTree
  # pos_end +1 because Interval tree use [a,b) intervals
  $interval_tree->insert($genomic_interval, $genomic_interval->start, $genomic_interval->end + 1);
}

=head2 loadAnnotations

  Arg [1] : Annotations     - DEkupl::Annotations object
  Description: Load annotations into the intervalQuery

=cut

sub loadAnnotations {
  my $self = shift;
  my $annotations = shift;

  # TODO Verify that annotations are indeed DEkupl::Anotations

  foreach my $gene ($annotations->allGenes) {
    # Skip empty genes
    next if $gene->nbExons == 0;

    # Add gene to the interval query
    $self->addInterval($gene);
    
    # Add exons to the interval query
    foreach my $exon ($gene->allExons) {
      $self->addInterval($exon);
    }
  }
}

=head2 fetchByRegion

  Arg [1] : GenomicInterval     - Genomic Interval
  Arg [2] : (Optional) Boolean  - Windowed query, only return intervals which
                                  are completely contained in the queried region.

=cut

sub fetchByRegion {
  my ($self,$genomic_interval,$windowed) = @_;

  my $interval_tree = $self->getIntervalTree($genomic_interval->chr,$genomic_interval->strand);
  
  if(defined $interval_tree) {
    if(defined $windowed && $windowed) {
      # pos_end +1 because Interval tree use [a,b) intervals
      return $interval_tree->fetch_window($genomic_interval->start,$genomic_interval->end+1);
    } else {
      # pos_end +1 because Interval tree use [a,b) intervals
      return $interval_tree->fetch($genomic_interval->start,$genomic_interval->end+1);
    }
  }

  return [];
}

=head1 PRIVATE METHODS

=head2 getIntervalTree 

  Arg [1] : String             - Chromosome
  Arg [2] : (Optional)         - Strand
  Description : Return the Set::IntervalTree reference for the chromosome and strand (Default : 1)
  ReturnType  : Set::IntervalTree

=cut

sub getIntervalTree {
  my ($self,$chr,$strand) = @_;
  my $interval_key = _getIntervalTreeKey($chr,$strand);
  return $self->_getIntervalTree($interval_key);
}

=head2 setIntervalTree
  Arg [1] : String             - Chromosome
  Arg [2] : (Optional) Integer - Strand
  Arg [3] : Set::IntervalTree  - Interval tree
  Description : Add an Set::IntervalTree object for a specific ("chr","strand") pair.
                Strand is set to 1 if none (or undef) is provided
=cut

sub setIntervalTree {
  my ($self,$chr,$strand,$interval_tree) = @_;
  my $interval_key = _getIntervalTreeKey($chr,$strand);
  $self->_setIntervalTree($interval_key => $interval_tree);
}

=head2 _getIntervalTreeKey

  Arg [1] : String             - Chromosome
  Arg [2] : (Optional) Integer - Strand
  Description : Static method that return and unique key for the ("chr","strand") pair passed in arguements.
                Strand is set to 1 if none (or undef) is provided
  ReturnType  : String

=cut

sub _getIntervalTreeKey {
  my ($chr,$strand) = @_;
  return "$chr"."@"."$strand";
}

no Moose;
__PACKAGE__->meta->make_immutable;