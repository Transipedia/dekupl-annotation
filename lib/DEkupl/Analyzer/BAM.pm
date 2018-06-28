package DEkupl::Analyzer::BAM;
# ABSTRACT: Append Alignment informations to contigs

use Moose;

use DEkupl::Utils;

with 'DEkupl::Analyzer';

has 'bam_file' => (
  is => 'ro',
  isa => 'Str',
  required => 1,
);

my @columns = (
  'is_mapped',
  'line_in_sam',
  'chromosome',
  'start',
  'end',
  'strand',
  'cigar',
  'nb_insertion',
  'nb_deletion',
  'nb_splice',
  'clipped_3p',
  'clipped_5p',
  'query_cover',
  'alignment_identity',
  'nb_hit',
  'nb_mismatch',
);

my %flags = (
  MULTIPLE_SEGMENTS => 1,
  PROPERLY_ALIGNED => 2,
  UNMAPPED => 4,
  NEXT_UNMAPPED => 8,
  REVERSE_COMPLEMENTED => 16,
  NEXT_REVERSE_COMPLEMENTED => 32,
  FIRST_SEGMENT => 64,
  LAST_SEGMENT => 128,
  SECONDARY_ALIGNMENT => 256,
  QUALITY_CONTROLS_FAILED => 512,
  PCR_DUPLICATED => 1024,
  SUPPLEMENTARY_ALIGNMENT => 2048,
);

# TODO Split this into function to be more testable!
sub BUILD {
  my $self = shift;
  my $bam_it = DEkupl::Utils::bamFileIterator($self->bam_file);
  my $i = 0;
  while(my $sam_line = $bam_it->()) {

    # Slect only primary alignments
    next if $sam_line->{flag} & $flags{SECONDARY_ALIGNMENT} || $sam_line->{flag} & $flags{SUPPLEMENTARY_ALIGNMENT};
 
    # The tag is the query name
    my $tag = $sam_line->{qname};

    # Load contig from the DB
    my $contig = $self->contigs_db->loadContig($tag);

    $contig->{is_mapped} = DEkupl::Utils::booleanEncoding(!($sam_line->{flag} & $flags{UNMAPPED}));
    $contig->{line_in_sam} = $i;

    # We skip the rest of annotation of the contig is unmapped
    next if $sam_line->{flag} & $flags{UNMAPPED};

    $contig->{chromosome} = $sam_line->{rname};
    $contig->{start} = $sam_line->{pos};
    $contig->{cigar} = $sam_line->{original_cigar};

    # Set strand if we are in 'strand-specific' mode
    $contig->{strand} = _getStrandFromFlag($sam_line->{flag}, $self->is_stranded);

    # Compute alignment length from cigar
    # Extracting information from CIGAR element
    my %cig_stats = %{_computeStatsFromCigar($sam_line->{cigar})};

    $contig->{nb_insertion} = $cig_stats{nb_insertion};
    $contig->{nb_deletion}  = $cig_stats{nb_deletion}; 
    $contig->{nb_splice}    = $cig_stats{nb_splice};
    $contig->{end}          = $sam_line->{pos} + $cig_stats{ref_aln_length} - 1;
    $contig->{clipped_3p}   = DEkupl::Utils::booleanEncoding($cig_stats{is_clipped_3p});
    $contig->{clipped_5p}   = DEkupl::Utils::booleanEncoding($cig_stats{is_clipped_5p});

    # Compute the query cover (fraction of based aligned)
    # As in BLAST
    $contig->{query_cover}  = $cig_stats{query_aln_length} / length($sam_line->{seq});

    # Extracting information from extented SAM fields
    my $nb_hit      = $sam_line->{extended_fields}->{NH};
    my $nb_mismatch = $sam_line->{extended_fields}->{NM};
    
    $contig->{nb_hit}       = $nb_hit if defined $nb_hit;
    $contig->{nb_mismatch}  = $nb_mismatch if defined $nb_mismatch;

    # Compute alignemnt identity
    $contig->{alignment_identity} = ($cig_stats{nb_match} - $nb_mismatch) / $cig_stats{query_aln_length};
    
    # Save contig
    $self->contigs_db->saveContig($contig);
    $i++;
  }
}

sub getHeaders {
  my $self = shift;
  return @columns;
}

sub getValues {
  my $self = shift;
  my $contig = shift;
  my @values = map { defined $contig->{$_}? $contig->{$_} : 'NA' } @columns;
  return @values;
}

sub _getStrandFromFlag {
  my $flag        = shift;
  my $is_stranded = shift;
  my $strand;
  if($is_stranded) {
    if($flag & $flags{REVERSE_COMPLEMENTED}) {
      $strand = '-';
    } else {
      $strand = '+';
    }
  }
  return $strand;
}

sub _computeStatsFromCigar {
  my $cigar = shift;
  my %cig_stats = (
    ref_aln_length    => 0, # Length of the alignment on the reference (from the first based aligned to the last, including deletion and splice)
    query_aln_length  => 0, # Length of the alignment on the query (including deletion)
    nb_match          => 0,
    nb_insertion      => 0, # Number of insertion in the query (I cigar element)
    nb_deletion       => 0, # Number of deletiong in the query (D cigar element)
    nb_splice         => 0, # Number of splice in the query (N cigar element),
    is_clipped_5p     => 0,
    is_clipped_3p     => 0,
  );

  if(@{$cigar} > 1) {
    $cig_stats{is_clipped_5p} = 1 if $cigar->[0]->{op} =~ /[SH]/;
    $cig_stats{is_clipped_3p} = 1 if $cigar->[$#{$cigar}]->{op} =~ /[SH]/; 
  }
     
  foreach my $cigel (@{$cigar}) {
    if($cigel->{op} =~ /[MX=]/ ) {
      $cig_stats{ref_aln_length}    += $cigel->{nb};
      $cig_stats{query_aln_length}  += $cigel->{nb};
      $cig_stats{nb_match}          += $cigel->{nb};
    } elsif($cigel->{op} eq 'I') {
      $cig_stats{nb_insertion}      += 1;
      $cig_stats{query_aln_length}  += $cigel->{nb};
    } elsif($cigel->{op} eq 'N') {
      $cig_stats{nb_splice}         += 1;
      $cig_stats{ref_aln_length}    += $cigel->{nb};
    } elsif($cigel->{op} eq 'D') {
      $cig_stats{nb_deletion}       += 1;
      $cig_stats{ref_aln_length}    += $cigel->{nb};
    }
  }
  
  return \%cig_stats;
}

no Moose;
__PACKAGE__->meta->make_immutable;