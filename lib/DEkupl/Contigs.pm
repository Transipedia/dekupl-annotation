package DEkupl::Contigs;
# ABSTRACT: Manipulate contigs from contigs_file

use Moose;

use DEkupl::Utils;

with 'DEkupl::Analyzer';

has 'contigs_file' => (
  is => 'ro',
  isa => 'Str',
  required => 1,
);

has 'sample_names' => (
  is => 'rw',
  isa => 'ArrayRef[Str]',
  traits => ['Array'],
  handles => {
    all_samples => 'elements',
  },
);

has 'min_abs_log2FC' => (
  is => 'rw',
  isa => 'Num'
);

has 'max_abs_log2FC' => (
  is => 'rw',
  isa => 'Num'
);

my @columns = (
  'tag',
  'nb_merged_kmers',
  'contig',
  'contig_size',
  'pvalue',
  'meanA',
  'meanB',
  'log2FC'
);

sub BUILD {
  my $self = shift;
  my $fh = DEkupl::Utils::getReadingFileHandle($self->contigs_file);
  my $header = <$fh>;

  chomp $header;
  my $info =  _splitLine($header);

  # Load sample names
  my @samples;
  foreach my $s (@{$info->{counts}}) {
    push @samples, $s;
  }
  $self->sample_names(\@samples);

  close($fh);
}

sub generateFasta {
  my $self = shift;
  my $fasta_file = shift;
  my $fhi = DEkupl::Utils::getReadingFileHandle($self->contigs_file);
  my $fho = DEkupl::Utils::getWritingFileHandle($fasta_file);
  my $header = <$fhi>;
  while(<$fhi>) {
    chomp;
    # nb_merged_kmers   contig  tag   pvalue   meanA   meanB   log2FC   sample1   sample2 ...
    my $contig = _splitLine($_);
    print $fho ">".$contig->{tag}."\n".$contig->{contig},"\n";
  }
  close($fho);
  close($fhi);
}

sub loadContigsDB {
  my $self = shift;
  my $contigs_db = shift;
  my $fh = DEkupl::Utils::getReadingFileHandle($self->contigs_file);
  my $header = <$fh>;

  my $i = 0;
  my $min_abs_log2FC = 0;
  my $max_abs_log2FC = 0;

  while(<$fh>) {
    chomp;
    my $contig = _splitLine($_);
    $contig->{contig_size} = length $contig->{contig};
    $contigs_db->saveContig($contig);
    
    # Setting min-max absoluted log2FC (used for coloring contigs)
    my $abs_log2FC = abs($contig->{log2FC});
    if(!defined $min_abs_log2FC) {
      $min_abs_log2FC = $abs_log2FC;
      $max_abs_log2FC = $abs_log2FC;
    } else {
      $min_abs_log2FC = $abs_log2FC if $abs_log2FC < $min_abs_log2FC;
      $max_abs_log2FC = $abs_log2FC if $abs_log2FC > $max_abs_log2FC;
    }
    
    # print STDERR "$i contigs loaded\n" if $i % 100 == 0;
    $i++;
  }

  $self->min_abs_log2FC($min_abs_log2FC);
  $self->max_abs_log2FC($max_abs_log2FC);
  
  close($fh);
}

sub getHeaders {
  my $self = shift;
  my $fh = DEkupl::Utils::getReadingFileHandle($self->contigs_file);
  my $header = <$fh>;
  chomp $header;
  my $contig_headers = _splitLine($header);
  return (@columns, @{$contig_headers->{counts}});
}

sub getValues {
  my $self = shift;
  my $contig = shift;
  my @values = map { $contig->{$_} } @columns;
  push @values, @{$contig->{counts}};
  return @values;
}

sub _splitLine {
  my $line = shift;

  # nb_merged_kmers   contig  tag   pvalue   meanA   meanB   log2FC   sample1   sample2 ...
  my ($nb_merged_kmers, $contig, $tag, $pvalue, $meanA, $meanB, $log2FC, @counts) = split "\t", $line;

  return {
      'tag' => $tag,
      'nb_merged_kmers' => $nb_merged_kmers,
      'contig' => $contig,
      'pvalue' => $pvalue,
      'meanA' => $meanA,
      'meanB' => $meanB,
      'log2FC' => $log2FC,
      'counts' => \@counts,
  };
}

no Moose;
__PACKAGE__->meta->make_immutable;