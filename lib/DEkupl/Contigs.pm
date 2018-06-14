package DEkupl::Contigs;
# ABSTRACT: Manipulate contigs from contigs_file

use Moose;

use DEkupl::Utils;

has 'contigs_file' => (
  is => 'ro',
  isa => 'Str',
  required => 1,
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
    my @f = split "\t", $_;
    print $fho ">".$f[2]."\n".$f[1],"\n";
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
    # nb_merged_kmers   assembly  tag   pvalue   meanA   meanB   log2FC   sample1   sample2 ...
    my ($nb_merged_kmers, $assembly, $tag, $pvalue, $meanA, $meanB, $log2FC, @counts) = split "\t", $_;
    my $contig = {
      'tag' => $tag,
      'assembly' => $assembly,
      'pvalue' => $pvalue,
      'meanA' => $meanA,
      'meanB' => $meanB,
      'log2FC' => $log2FC,
      'counts' => \@counts,
    };
    $contigs_db->saveContig($contig);
    # print STDERR "$i contigs loaded\n" if $i % 100 == 0;
    $i++;
  }
  
  close($fh);
}

no Moose;
__PACKAGE__->meta->make_immutable;