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

my @columns = (
  'tag',
  'nb_merged_kmers',
  'assembly',
  'pvalue',
  'meanA',
  'meanB',
  'log2FC'
);

# sub BUILD {
#   my $self = shift;
#   # Load contigs?
# }

sub generateFasta {
  my $self = shift;
  my $fasta_file = shift;
  my $fhi = DEkupl::Utils::getReadingFileHandle($self->contigs_file);
  my $fho = DEkupl::Utils::getWritingFileHandle($fasta_file);
  my $header = <$fhi>;
  while(<$fhi>) {
    chomp;
    # nb_merged_kmers   assembly  tag   pvalue   meanA   meanB   log2FC   sample1   sample2 ...
    my $contig = _splitLine($_);
    print $fho ">".$contig->{tag}."\n".$contig->{assembly},"\n";
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
  while(<$fh>) {
    chomp;
    my $contig = _splitLine($_);
    $contigs_db->saveContig($contig);
    # print STDERR "$i contigs loaded\n" if $i % 100 == 0;
    $i++;
  }
  
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

  # nb_merged_kmers   assembly  tag   pvalue   meanA   meanB   log2FC   sample1   sample2 ...
  my ($nb_merged_kmers, $assembly, $tag, $pvalue, $meanA, $meanB, $log2FC, @counts) = split "\t", $line;

  return {
      'tag' => $tag,
      'nb_merged_kmers' => $nb_merged_kmers,
      'assembly' => $assembly,
      'pvalue' => $pvalue,
      'meanA' => $meanA,
      'meanB' => $meanB,
      'log2FC' => $log2FC,
      'counts' => \@counts,
  };
}

no Moose;
__PACKAGE__->meta->make_immutable;