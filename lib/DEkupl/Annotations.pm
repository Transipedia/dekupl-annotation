package DEkupl::Annotations;
# ABSTRACT: A collection of annotations

use Moose;
use Carp;

use DEkupl::Utils;
use DEkupl::Annotations::Exon;
use DEkupl::Annotations::Gene;

=head1 SYNOPSIS

A collection of annotations (gene => [exon => transcript])
Annotations are 0-based.

=cut

# We store genes in a hash where the id is the uniq gene_id
has genes_hash => (
  traits => ['Hash'],
  is  => 'ro',
  isa => 'HashRef[DEkupl::Annotations::Gene]',
  default => sub { {} },
  handles => {
    allGenes    => 'values',
    allGeneIds  => 'keys',
    getGene     => 'get',
    nbGenes     => 'count',
    removeGene  => 'delete',
  },
);

my %gene_features = (
  'gene' => 'coding gene',
  'pseudogene' => 'pseudogene',
  'ncRNA_gene' => 'non_coding gene',
);

sub sortedGenes {
  my $self = shift;
  return sort { $a->chr cmp $b->chr || $a->start <=> $b->start } $self->allGenes;
}

sub addGene($) {
  my $self = shift;
  my $gene = shift;
  $self->genes_hash->{$gene->id} = $gene;
}

=head2 loadFromGFF
  $annotations->loadFromGFF($gtf_file);

=cut

sub loadFromGFF {
  my ($self,$gtf_file) = @_;

  # Now we read the annotations and load them into memory
  my $gff_it = DEkupl::Utils::gffFileIterator($gtf_file,'gff3');

  # We supposed that GFF are corectly sorted
  # First genes, then mRNA, then exons
  my $gene_id;
  while(my $annot = $gff_it->()) {

    print STDERR "Loading annotations from chr ".$annot->{chr}."...\n" if($annot->{feature} eq 'chromosome');

    my $id = $annot->{attributes}->{ID};
    # next if !defined $id;

    # We only consider exon annotations
    if(defined $gene_features{$annot->{feature}}) {
      # print STDERR ".JFCYKUIMLCJFVKBILOML.K\n";
      $gene_id = $id;
      my $gene_symbol = $annot->{attributes}->{Name};

      my $gene = $self->getGene($gene_id);
      if(!defined $gene) {
        my %gene_info = (
          chr     => $annot->{chr},
          strand  => $annot->{strand},
          id      => $gene_id,
          symbol  => $gene_symbol,
        );
        my $biotype = $annot->{attributes}->{biotype};
        $gene_info{biotype} = $biotype if defined $biotype;
        # If this gene is not defined yet, we create a new entry
        $gene = DEkupl::Annotations::Gene->new(%gene_info);
        $self->addGene($gene);
      }
    }

    # print STDERR "TATA $gene_id for exon\n";

    # print STDERR "TOTO ".$annot->{feature}."\n";

    if($annot->{feature} eq 'exon') {
      # print STDERR "Get gene $gene_id for exon\n";

      my $gene = $self->getGene($gene_id);

      # print STDERR "ADDing exon to gene ".$gene->id."\n";

      # The exon adds itself to the gene
      my $exon = DEkupl::Annotations::Exon->new(
        start     => $annot->{start} - 1,
        end       => $annot->{end} - 1,
        gene      => $gene,
        # transcripts => [$annot->{attributes}->{Parent}],
        transcripts => [],
      );
    }
  }
}

=head2 loadFromGTF
  $annotations->loadFromGTF($gtf_file);

=cut

sub loadFromGTF {
  my ($self,$gtf_file) = @_;

  # Now we read the annotations and load them into memory
  my $gtf_it = DEkupl::Utils::gffFileIterator($gtf_file,'gtf');
  while(my $annot = $gtf_it->()) {

    # We only consider exon annotations
    next if $annot->{feature} ne 'exon';

    my $gene_id = $annot->{attributes}->{gene_id};
    next if !defined $gene_id;

    my $gene = $self->getGene($gene_id);

    # If this gene is not defined yet, we create a new entry
    if(!defined $gene) {
      $gene = DEkupl::Annotations::Gene->new(
        chr     => $annot->{chr},
        strand  => $annot->{strand},
        id      => $gene_id,
      );
      $self->addGene($gene);
    }

    # The exon adds itself to the gene
    my $exon = DEkupl::Annotations::Exon->new(
      start     => $annot->{start} - 1,
      end       => $annot->{end} - 1,
      gene      => $gene,
      transcripts => [$annot->{attributes}->{transcript_id}],
    );
  }
}

no Moose;
__PACKAGE__->meta->make_immutable;

__END__
=head1 METHODS

=head2 new

Create a new 'DEkupl::Annotations' object

=head2 addGene('DEkupl::Annotations::Gene')

Add a new gene to the annotations

=head2 getGene('gene_id') => 'DEkupl::Annotations::Gene'

Given a gene_id return the gene

=head2 removeGene('gene_id')

Given a gene_id remove this gene from the collection (if it exists)

=head2 nbGenes => 'Int'

Return the number of genes

=head2 allGenes => Array['DEkupl::Annotations::Gene']

Return all genes contained in the annotations

=head2 allGeneIds => Array['gene_id']

Return all gene ids contained in the annotations

=head2 sortedGenes => Array['DEkupl::Annotations::Gene']

Return all genes sorted by 'chr' and 'start' position

=head2 loadGTF('gtf_file')

Load annotations contained if the gtf_file passed in argument

=head2 saveGTF('gtf_file')

Save the annotations in the gtf_file passed in argument