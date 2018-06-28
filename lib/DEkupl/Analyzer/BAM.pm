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
    if($self->is_stranded) {
      my $strand;
      if($sam_line->{flag} & $flags{REVERSE_COMPLEMENTED}) {
        $strand = '-';
      } else {
        $strand = '+';
      }
      $contig->{strand} = $strand;
    }

    # Compute alignment length from cigar
    # Extracting information from CIGAR element
    my $ref_aln_length = 0;   # Length of the alignment on the reference (from the first based aligned to the last, including deletion and splice)
    my $query_aln_length = 0; # Length of the alignment on the query (including deletion)
    my $nb_match = 0;
    my $nb_insertion = 0;     # Number of insertion in the query (I cigar element)
    my $nb_deletion = 0;      # Number of deletiong in the query (D cigar element)
    my $nb_splice = 0;        # Number of splice in the query (N cigar element)
    foreach my $cigel (@{$sam_line->{cigar}}) {
      if($cigel->{op} =~ /[MX=]/ ) {
        $ref_aln_length += $cigel->{nb};
        $query_aln_length += $cigel->{nb};
        $nb_match += $cigel->{nb};
      } elsif($cigel->{op} eq 'I') {
        $nb_insertion++;
        $query_aln_length += $cigel->{nb};
      } elsif($cigel->{op} eq 'N') {
        $nb_splice++;
        $ref_aln_length += $cigel->{nb};
      } elsif($cigel->{op} eq 'D') {
        $nb_deletion++;
        $ref_aln_length += $cigel->{nb};
      }
    }

    $contig->{nb_insertion} = $nb_insertion;
    $contig->{nb_deletion} = $nb_deletion; 
    $contig->{nb_splice} = $nb_splice;
    $contig->{end} = $sam_line->{pos} + $ref_aln_length - 1;

    # Check if alignment is clipped (Hard or Soft)
    my $clipped_5p = 0;
    my $clipped_3p = 0;
    if(@{$sam_line->{cigar}} > 1) {
      $clipped_5p = 1 if $sam_line->{cigar}->[0]->{op} =~ /[SH]/;
      $clipped_3p = 1 if $sam_line->{cigar}->[$#{$sam_line->{cigar}}]->{op} =~ /[SH]/; 
    }
    $contig->{clipped_3p} = DEkupl::Utils::booleanEncoding($clipped_3p);
    $contig->{clipped_5p} = DEkupl::Utils::booleanEncoding($clipped_5p);

    # Compute the query cover (fraction of based aligned)
    # As in BLAST
    my $query_cover = $query_aln_length / length($sam_line->{seq});
    $contig->{query_cover} = $query_cover;
    

    # Extracting information from extented SAM fields
    my $nb_hit = $sam_line->{extended_fields}->{NH};
    my $nb_mismatch = $sam_line->{extended_fields}->{NM};
    
    $contig->{nb_hit} = $nb_hit if defined $nb_hit;
    $contig->{nb_mismatch} = $nb_mismatch if defined $nb_mismatch;

    # Compute alignemnt identity
    my $alignment_identity = ($nb_match - $nb_mismatch) / $query_aln_length;
    $contig->{alignment_identity} = $alignment_identity;
    
    # Save contig
    $self->contigs_db->saveContig($contig);
    $i++;
  }
  # Load contigs?
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

# sub _isFlagged {
#   my $value = shift;
#   my $flag = shift;
#   return $value & $flag;
# }


no Moose;
__PACKAGE__->meta->make_immutable;