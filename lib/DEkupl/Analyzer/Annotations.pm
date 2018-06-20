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
  '5p_gene_id',
  '5p_gene_symbol',
  '5p_gene_dist',
  '3p_gene_id',
  '3p_gene_symbol',
  '3p_gene_dist',
  'exonic',
  'intronic',
);

sub BUILD {
  my $self = shift;

  my $contigs_it = $self->contigs_db->contigsIterator();

  while(my $contig = $contigs_it->()) {

    if($contig->{is_mapped}) {

      my ($fwd_results,$rv_results);
      my ($five_prim_result, $five_prim_dist, $three_prim_result, $three_prim_dist);

      my $query = DEkupl::GenomicInterval->new(
          chr => $contig->{chromosome},
          start => $contig->{start},
          end => $contig->{end},
      );

      if($self->is_stranded) {
        # Query annotations from contig strand
        $query->strand($contig->{strand});
        $fwd_results        = $self->interval_query->fetchByRegion($query);

        # If we have no results on the contig regions, we try to find the nearset neighbors.
        if(scalar @{$fwd_results} == 0) {
          ($five_prim_result,$five_prim_dist)   = $self->interval_query->fetchNearest5prim($query);
          ($three_prim_result,$three_prim_dist) = $self->interval_query->fetchNearest3prim($query);
        }
        
        # Set RV results with contig reverse strand
        $query->strand(DEkupl::Utils::reverseStrand($contig->{strand}));
        $rv_results = $self->interval_query->fetchByRegion($query);

      # else : unstranded case
      } else {
        # Add annotations on FWD strand
        $query->strand('+');
        $fwd_results = $self->interval_query->fetchByRegion($query);

        # Append annotations from RV strand
        $query->strand('-');
        push @{$fwd_results}, @{$self->interval_query->fetchByRegion($query)};

        # If we have no results on the contig regions, we try to find the nearset neighbors.
        if(scalar @{$fwd_results} == 0) {
          # Query forward strand
          $query->strand('+');
          my ($five_prim_result_fwd, $five_prim_dist_fwd) = $self->interval_query->_fetchNearestDown($query);
          my ($three_prim_result_fwd, $three_prim_dist_fwd) = $self->interval_query->_fetchNearestUp($query);

          # Query reverse strand
          $query->strand('-');
          my ($five_prim_result_rv, $five_prim_dist_rv) = $self->interval_query->_fetchNearestDown($query);
          my ($three_prim_result_rv, $three_prim_dist_rv) = $self->interval_query->_fetchNearestUp($query);

          # Select the closest 5prim gene between the two strand
          # If both are defined we choose the closest
          # Otherwise we chose the one that is defined
          if(defined $five_prim_result_fwd && defined $five_prim_result_rv) { 
            if($five_prim_dist_fwd < $five_prim_dist_rv) {
              ($five_prim_result,$five_prim_dist) = ($five_prim_result_fwd, $five_prim_dist_fwd);
            } else {
              ($five_prim_result,$five_prim_dist) = ($five_prim_result_rv, $five_prim_dist_rv);
            }
          } elsif(defined $five_prim_result_fwd) {
            ($five_prim_result,$five_prim_dist) = ($five_prim_result_fwd, $five_prim_dist_fwd);
          } elsif(defined $five_prim_result_rv) {
            ($five_prim_result,$five_prim_dist) = ($five_prim_result_rv, $five_prim_dist_rv);
          }
          
          # Select the closest 3prim gene between the two strand
          # If both are defined we choose the closest
          # Otherwise we chose the one that is defined
          if(defined $three_prim_result_fwd && defined $three_prim_result_rv) {
            if($three_prim_dist_fwd < $three_prim_dist_rv) {
              ($three_prim_result,$three_prim_dist) = ($three_prim_result_fwd, $three_prim_dist_fwd);
            } else {
              ($three_prim_result,$three_prim_dist) = ($three_prim_result_rv, $three_prim_dist_rv);
            }
          } elsif(defined $three_prim_result_fwd) {
            ($three_prim_result,$three_prim_dist) = ($three_prim_result_fwd, $three_prim_dist_fwd);
          } elsif(defined $three_prim_result_rv) {
            ($three_prim_result,$three_prim_dist) = ($three_prim_result_rv, $three_prim_dist_rv);
          }
        }

        # Set empty RV results
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
            } elsif($strand eq '5prim') {
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

      $contig->{exonic}   = DEkupl::Utils::booleanEncoding($exonic);
      $contig->{intronic} = DEkupl::Utils::booleanEncoding($intronic);

      # Set 5prim annotations
      if(defined $five_prim_result) {
        # If we have match an exon, we get the gene object
        if(ref($five_prim_result) eq 'DEkupl::Annotations::Exon') {
          $five_prim_result = $five_prim_result->gene;
        }
        $contig->{'5p_gene_id'}     = $five_prim_result->id;
        $contig->{'5p_gene_symbol'} = $five_prim_result->symbol;
        $contig->{'5p_gene_dist'}   = $five_prim_dist;
      }

      # Set 3prim annotations
      if(defined $three_prim_result) {
        # If we have match an exon, we get the gene object
        if(ref($three_prim_result) eq 'DEkupl::Annotations::Exon') {
          $three_prim_result = $three_prim_result->gene;
        }
        $contig->{'3p_gene_id'}     = $three_prim_result->id;
        $contig->{'3p_gene_symbol'} = $three_prim_result->symbol;
        $contig->{'3p_gene_dist'}   = $three_prim_dist;
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