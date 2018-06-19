package DEkupl::Analyzer::Annotations;
# ABSTRACT: Append annotation informations to contigs

use Moose;

use DEkupl::Utils;

with 'DEkupl::Analyzer';

# This is a object used to query annotations loaded into memory
has 'interval_query' => (
  is => 'ro',
  isa => 'DEkupl::IntervalQuery',
  required => 1,
);

# TODO This should be a hash to auto generate the documentation!
my @columns = (
  'gene_id',
  'gene_symbol'
);

sub BUILD {
  my $self = shift;

  my $contigs_it = $self->contigs_db->contigsIterator();

  while(my $contig = $contigs_it->()) {

    if($contig->{is_mapped}) {

      # TODO we should add a non strand-specific mode!
      # Query annotations
      my $query = DEkupl::GenomicInterval->new(
        chr => $contig->{chromosome},
        start => $contig->{start},
        end => $contig->{end},
        strand => $contig->{strand},
      );

      my $results = $self->interval_query->fetchByRegion($query);

      foreach my $res (@{$results}) {
        
        my $res_type = ref($res);

        # TODO we should do a special treatment when there is multiple genes overlapping
        # the position. Usually we should choose the one that is 'protein_coding' over
        # a non_conding gene!
        if($res_type eq 'DEkupl::Annotations::Gene') {
          $contig->{gene_id} = $res->id;
          $contig->{gene_symbol} = $res->symbol;
        }
      }

      # Save contig
      $self->contigs_db->saveContig($contig);
    }
  }
}


sub getHeaders {
  my $self = shift;
  return @columns;
}

sub getValues {
  my $self = shift;
  my $contig = shift;
  my @values = map { defined $contig->{$_}? $contig->{$_} : $DEkupl::Utils::NA_value } @columns;
  return @values;
}

no Moose;
__PACKAGE__->meta->make_immutable;