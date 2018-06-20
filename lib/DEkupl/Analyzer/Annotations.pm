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
  'gene_symbol',
  'gene_biotype',
  'as_gene_id',
  'as_gene_symbol',
  'as_gene_biotype',
  'exonic',
  'intronic',
);

sub BUILD {
  my $self = shift;

  my $contigs_it = $self->contigs_db->contigsIterator();

  while(my $contig = $contigs_it->()) {

    if($contig->{is_mapped}) {

      my ($fwd_results,$rv_results);

      my $query = DEkupl::GenomicInterval->new(
          chr => $contig->{chromosome},
          start => $contig->{start},
          end => $contig->{end},
      );

      if($self->is_stranded) {
        # TODO we should add a non strand-specific mode!
        # Query annotations
        $query->strand($contig->{strand});
        $fwd_results = $self->interval_query->fetchByRegion($query);
        $query->strand(DEkupl::Utils::reverseStrand($contig->{strand}));
        $rv_results = $self->interval_query->fetchByRegion($query);
      } else {
        $query->strand('+');
        $fwd_results = $self->interval_query->fetchByRegion($query);
        $query->strand('-');
        push @{$fwd_results}, @{$self->interval_query->fetchByRegion($query)};
        $rv_results = [];
      }

      my $exonic = 0;
      my $intronic;
      foreach my $strand (('fwd','rv')) {
        my @results = $strand eq 'fwd'? @{$fwd_results} : @{$rv_results};
        foreach my $res (@results) {
          my $res_type = ref($res);

          # TODO we should do a special treatment when there is multiple genes overlapping
          # the position. Usually we should choose the one that is 'protein_coding' over
          # a non_conding gene!
          if($res_type eq 'DEkupl::Annotations::Gene') {
            if($strand eq 'fwd') {
              $contig->{gene_id} = $res->id;
              $contig->{gene_symbol} = $res->symbol;
              $contig->{gene_biotype} = $res->biotype;
            } else {
              $contig->{as_gene_id} = $res->id;
              $contig->{as_gene_symbol} = $res->symbol;
              $contig->{as_gene_biotype} = $res->biotype;
            }
          } elsif($res_type eq 'DEkupl::Annotations::Exon') {
            $exonic = 1;
            # The contig overlap the exon and the intron
            if ($query->start < $res->start || $query->end > $res->end) {
              $intronic = 1;
            } elsif(!defined $intronic) {
              $intronic = 0;
            }
          }
        }
      }
      
      $intronic = 1 if !defined $intronic; # We have found no exons, therefor we are fully intronic.

      $contig->{exonic} = DEkupl::Utils::booleanEncoding($exonic);
      $contig->{intronic} = DEkupl::Utils::booleanEncoding($intronic);

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