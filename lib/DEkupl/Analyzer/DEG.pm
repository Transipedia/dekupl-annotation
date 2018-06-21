package DEkupl::Analyzer::DEG;
# ABSTRACT: Append DEG (Differential Expressed Gene) informations to contigs

use Moose;

use DEkupl::Utils;

with 'DEkupl::Analyzer';

has 'deg_file' => (
  is => 'ro',
  isa => 'Str',
  required => 1,
);

has 'max_padj' => (
  is => 'ro',
  isa => 'Num',
  default => 0.05,
);

has 'min_log2fc' => (
  is => 'ro',
  isa => 'Num',
  default => 0,
);

my @columns = ('gene_is_diff');

sub BUILD {
  my $self = shift;
  my $fh = DEkupl::Utils::getReadingFileHandle($self->deg_file);

  my %gene_is_diff;

  my $headers = <$fh>;
  while(<$fh>) {
    chomp;
    #gene_id baseMean        log2FoldChange  lfcSE   stat    pvalue  padj
    my ($gene_id, $base_mean, $log2FC, $lfcSE, $stat, $pvalue, $padj) = split "\t", $_;

    # Remove .* extentensions
    $gene_id = DEkupl::Utils::getAtomicGeneID($gene_id);

    next if $padj eq 'NA' || $log2FC eq 'NA';
    
    my $is_diff = $padj <= $self->max_padj && abs($log2FC) >= $self->min_log2fc? 1 : 0;
    
    $gene_is_diff{$gene_id} = $is_diff;
  }

  my $contigs_it = $self->contigs_db->contigsIterator();
  while(my $contig = $contigs_it->()) {
    if(defined $contig->{gene_id}) {
      my $is_diff = $gene_is_diff{$contig->{gene_id}};
      $contig->{gene_is_diff} = DEkupl::Utils::booleanEncoding($is_diff) if defined $is_diff;

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
  my @values = map { defined $contig->{$_}? $contig->{$_} : 'NA' } @columns;
  return @values;
}


no Moose;
__PACKAGE__->meta->make_immutable;