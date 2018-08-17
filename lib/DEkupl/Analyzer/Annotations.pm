package DEkupl::Analyzer::Annotations;
# ABSTRACT: Append annotation informations to contigs

use Moose;

use List::Util qw(min max);

use DEkupl::Utils;

with 'DEkupl::Analyzer';

# This is a object used to query annotations loaded into memory
has 'interval_query' => (
  is => 'ro',
  isa => 'DEkupl::IntervalQuery',
  required => 1,
);

has 'loci_file' => (
  is => 'ro',
  isa => 'Str',
);

# TODO This should be a hash to auto generate the documentation!
my @columns = (
  'gene_id',
  'gene_symbol',
  'gene_strand',
  'gene_biotype',
  'as_gene_id',
  'as_gene_symbol',
  'as_gene_strand',
  'as_gene_biotype',
  'upstream_gene_id',
  'upstream_gene_strand',
  'upstream_gene_symbol',
  'upstream_gene_dist',
  'downstream_gene_id',
  'downstream_gene_strand',
  'downstream_gene_symbol',
  'downstream_gene_dist',
  'exonic',
  'intronic',
);

sub BUILD {
  my $self = shift;

  my $contigs_it = $self->contigs_db->contigsIterator();

  my %loci;

  my $nb_annot            = 0;
  my $nb_annot_as         = 0;
  my $nb_annot_upstream   = 0;
  my $nb_annot_downstream = 0;

  while(my $contig = $contigs_it->()) {

    if($contig->{is_mapped}) {

      my ($fwd_results,$rv_results);
      my ($upstream_result, $upstream_dist, $downstream_result, $downstream_dist);

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
          ($upstream_result,$upstream_dist)   = $self->interval_query->fetchNearest5prim($query);
          ($downstream_result,$downstream_dist) = $self->interval_query->fetchNearest3prim($query);
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
          my ($upstream_result_fwd, $upstream_dist_fwd)     = $self->interval_query->_fetchNearestDown($query);
          my ($downstream_result_fwd, $downstream_dist_fwd) = $self->interval_query->_fetchNearestUp($query);

          # Query reverse strand
          $query->strand('-');
          my ($upstream_result_rv, $upstream_dist_rv)     = $self->interval_query->_fetchNearestDown($query);
          my ($downstream_result_rv, $downstream_dist_rv) = $self->interval_query->_fetchNearestUp($query);

          # Select the closest 5prim gene between the two strand
          # If both are defined we choose the closest
          # Otherwise we chose the one that is defined
          if(defined $upstream_result_fwd && defined $upstream_result_rv) { 
            if($upstream_dist_fwd < $upstream_dist_rv) {
              ($upstream_result,$upstream_dist) = ($upstream_result_fwd, $upstream_dist_fwd);
            } else {
              ($upstream_result,$upstream_dist) = ($upstream_result_rv, $upstream_dist_rv);
            }
          } elsif(defined $upstream_result_fwd) {
            ($upstream_result,$upstream_dist) = ($upstream_result_fwd, $upstream_dist_fwd);
          } elsif(defined $upstream_result_rv) {
            ($upstream_result,$upstream_dist) = ($upstream_result_rv, $upstream_dist_rv);
          }
          
          # Select the closest 3prim gene between the two strand
          # If both are defined we choose the closest
          # Otherwise we chose the one that is defined
          if(defined $downstream_result_fwd && defined $downstream_result_rv) {
            if($downstream_dist_fwd < $downstream_dist_rv) {
              ($downstream_result,$downstream_dist) = ($downstream_result_fwd, $downstream_dist_fwd);
            } else {
              ($downstream_result,$downstream_dist) = ($downstream_result_rv, $downstream_dist_rv);
            }
          } elsif(defined $downstream_result_fwd) {
            ($downstream_result,$downstream_dist) = ($downstream_result_fwd, $downstream_dist_fwd);
          } elsif(defined $downstream_result_rv) {
            ($downstream_result,$downstream_dist) = ($downstream_result_rv, $downstream_dist_rv);
          }
        }

        # Set empty RV results
        $rv_results = [];
      }

      # Select the best gene candidate for each strand
      foreach my $strand (('fwd','rv')) {
        # Get the right results array
        my @results = $strand eq 'fwd'? @{$fwd_results} : @{$rv_results};

        my $candidate = _selectBestCandidate(\@results, $query);

        # We have found a candidate gene!
        if(defined $candidate->{gene}) {
          if($strand eq 'fwd') {
            $nb_annot++;
            _setContigGeneInfo($contig,$candidate->{gene});
            $contig->{exonic}   = DEkupl::Utils::booleanEncoding($candidate->{is_exonic});
            $contig->{intronic} = DEkupl::Utils::booleanEncoding($candidate->{is_intronic});
          } elsif($strand eq 'rv') {
            $nb_annot_as++;
            _setContigGeneInfo($contig,$candidate->{gene},'as');
          }
        }
      }

      # Set 5prim annotations
      if(defined $upstream_result) {
        # If we have match an exon, we get the gene object
        if(ref($upstream_result) eq 'DEkupl::Annotations::Exon') {
          $upstream_result = $upstream_result->gene;
        }
        $nb_annot_upstream++;
        _setContigGeneInfo($contig,$upstream_result,'upstream');
        $contig->{'upstream_gene_dist'}   = $upstream_dist;
      }

      # Set 3prim annotations
      if(defined $downstream_result) {
        # If we have match an exon, we get the gene object
        if(ref($downstream_result) eq 'DEkupl::Annotations::Exon') {
          $downstream_result = $downstream_result->gene;
        }
        $nb_annot_downstream++;
        _setContigGeneInfo($contig,$downstream_result,'downstream');
        $contig->{'downstream_gene_dist'}   = $downstream_dist;
      }

      # Adding locus inforamations to the loci hash
      my ($locus_id, $locus_type) = _getLocusIdFromContig($contig);
      my $locus = $loci{$locus_id};
      # If the locus is already know we update its values
      if(defined $locus) {
        $locus->{best_pvalue} = $contig->{pvalue} if $contig->{pvalue} < $locus->{best_pvalue};
        $locus->{best_log2FC} = $contig->{log2FC} if abs($contig->{log2FC}) > abs($locus->{best_log2FC});
        $locus->{start}       = $contig->{start}  if $contig->{start} < $locus->{start};
        $locus->{end}         = $contig->{end}    if $contig->{end} > $locus->{end};
        $locus->{number_of_contigs}++;
      # This is a new locus
      } else {
        $loci{$locus_id} = {
          locus_id          => $locus_id,
          gene_symbol       => $contig->{gene_symbol},
          number_of_contigs => 1,
          best_log2FC       => $contig->{log2FC},
          best_pvalue       => $contig->{pvalue},
          chromosome        => $contig->{chromosome},
          start             => $contig->{start},
          end               => $contig->{end},
          strand            => $contig->{strand},
          locus_type        => $locus_type,
        };
      }

      # Save contig
      $self->contigs_db->saveContig($contig);
    }
  }

  $self->verboseLog("$nb_annot contigs annotated with a gene");
  $self->verboseLog("$nb_annot_as contigs annotated with an antisens gene");
  $self->verboseLog("$nb_annot_upstream contigs annotated with an upstream gene");
  $self->verboseLog("$nb_annot_downstream contigs annotated with an downstream gene");

  # Print locus to loci files
  my @loci_headers = qw(locus_id gene_symbol number_of_contigs best_log2FC best_pvalue chromosome start end strand locus_type);
  if(defined $self->loci_file) {
    my $loci_fh = DEkupl::Utils::getWritingFileHandle($self->loci_file) if defined $self->loci_file;
    print $loci_fh join("\t", @loci_headers), "\n";
    foreach my $locus (values %loci) {
      print $loci_fh join("\t", map { defined $locus->{$_}? $locus->{$_} : $DEkupl::Utils::NA_value } @loci_headers), "\n";
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

sub _getLocusIdFromContig {
  my $contig = shift;
  my $locus_id;
  my $locus_type;

  if(defined $contig->{gene_id}) {
    $locus_id   = $contig->{gene_id};
    $locus_type = 'genic';
  } elsif(defined $contig->{as_gene_id}) {
    $locus_id   = $contig->{as_gene_id};
    $locus_type = 'antisense';
  } else {
    # chr9&-&ENSG00000177910.7&ENSG00000271923.1
    $locus_id = join('&',
      $contig->{chromosome},
      defined $contig->{strand}? $contig->{strand} : $DEkupl::Utils::NA_value,
      defined $contig->{downstream_gene_id}?  $contig->{downstream_gene_id} : $DEkupl::Utils::NA_value,
      defined $contig->{upstream_gene_id}?    $contig->{upstream_gene_id}   : $DEkupl::Utils::NA_value,
    );
    $locus_type = 'intergenic';
  }

  return ($locus_id, $locus_type);
}

sub _setContigGeneInfo {
  my $contig  = shift;
  my $gene    = shift;
  my $type    = shift;
  my $prefix = "gene_";
  if(defined $type) {
    $prefix = $type.'_'.$prefix;
  }
  $prefix = "" unless defined $prefix;
  $contig->{$prefix."id"}      = $gene->id;
  $contig->{$prefix."strand"}  = $gene->strand;
  $contig->{$prefix."symbol"}  = $gene->symbol;
  $contig->{$prefix."biotype"} = $gene->biotype;
}

sub _selectBestCandidate {
  my $results = shift;
  my $query   = shift;
  my $exonic = 0;
  my $intronic;
  my $exonic_overlap_length;
  my $current_gene;

  foreach my $res (@{$results}) {
    my $res_type = ref($res);
    # TODO we should do a special treatment when there is multiple genes overlapping
    # the position. Usually we should choose the one that is 'protein_coding' over
    # a non_conding gene!
    if($res_type eq 'DEkupl::Annotations::Gene') {
      # We do not override possible exonic overlapping genes
      if(!$exonic) {
        $current_gene = $res;
      }
    } elsif($res_type eq 'DEkupl::Annotations::Exon') {
      # Genes with exons overlap take over non-exonic gene annotations
      # If multiple exons are overlapping, we take the one with the largest
      # overlapping length
      my $overlap = min($res->end,$query->end) - max($res->start,$query->start) + 1;
      if(!defined $exonic_overlap_length || $overlap > $exonic_overlap_length) {
        $exonic_overlap_length = $overlap;
        # The contig overlap the exon and the intron (only for fwd annotation)
        $intronic = ($query->start < $res->start || $query->end > $res->end)? 1 : 0;
        $current_gene = $res->gene;
      # Both candidates have the same overlapping length, we select the one with the longer gene length
      } elsif(defined $exonic_overlap_length && $overlap == $exonic_overlap_length) {
        # Selelect the longest gene
        if($res->gene->length > $current_gene->length) {
          # The contig overlap the exon and the intron (only for fwd annotation)
          $intronic = ($query->start < $res->start || $query->end > $res->end)? 1 : 0;
          $current_gene = $res->gene;
        }
      }
      $exonic = 1;
    }
  }

  $intronic = 1 if !defined $intronic; # We have found no exons, therefor we are fully intronic.

  return {
    'is_exonic' => $exonic,
    'is_intronic' => $intronic,
    'gene' => $current_gene,
  };
}

no Moose;
__PACKAGE__->meta->make_immutable;