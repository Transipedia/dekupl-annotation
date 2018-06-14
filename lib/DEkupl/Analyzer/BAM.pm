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

sub BUILD {
  my $self = shift;
  my $bam_it = $self->bamFileIterator();
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
    $contig->{chromosome} = $sam_line->{rname};
    $contig->{start} = $sam_line->{pos};

    # Compute alignment length from cigar
    my $alignement_length = 0;
    map { $alignement_length += $_->{nb} if $_->{op} =~ /[MDNX=]/ } @{$sam_line->{cigar}};

    $contig->{end} = $sam_line->{pos} + $alignement_length - 1;

    my $nb_hit = $sam_line->{extended_fields}->{NH};
    $contig->{nb_hit} = $nb_hit if defined $nb_hit;

    my $nb_mismatch = $sam_line->{extended_fields}->{NM};
    $contig->{nb_mismatch} = $nb_mismatch if defined $nb_mismatch;


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

sub bamFileIterator {
  my $self = shift;
  my $region = shift;
  $region = "" if !defined $region;
  my $file = $self->bam_file;
  open(my $fh, "-|", "samtools view $file $region" )or die "Cannot open $file, check if samtools are installed.";

  return sub {
    my $line = <$fh>;
    if($line) {
      return _parseSAMLine($line);
    }
    return $line;
  }

}

sub _parseSAMLine {
  my $line = shift;
  my ($qname,$flag,$rname,$pos,$mapq,$cigar,$rnext,$pnext,$tlen,$seq,$qual,@others) = split("\t",$line);
  my @cigar_hash = map { { op => substr($_,-1), nb => substr($_,0,length($_)-1)} } $cigar =~ /(\d+\D)/g;
  my %extended_fields = map { my ($id, $t, $v) = split ':', $_; $id => $v; } @others;
  return {
    qname => $qname,
    flag => $flag,
    rname => $rname,
    pos => $pos,
    mapq => $mapq,
    cigar => \@cigar_hash,
    original_cigar => $cigar,
    rnext => $rnext,
    pnext => $pnext,
    tlen => $tlen,
    seq => $seq,
    qual => $qual,
    extended_fields => \%extended_fields,
  };
}

# sub _isFlagged {
#   my $value = shift;
#   my $flag = shift;
#   return $value & $flag;
# }


no Moose;
__PACKAGE__->meta->make_immutable;