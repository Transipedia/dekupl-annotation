package DEkupl::App::Index;
# ABSTRACT: dkpl-index binary

use Moose;
use Getopt::Long;
use File::Temp qw/ tempdir /;

use DEkupl;
use DEkupl::Genome;


sub BUILD {
  my $self = shift;

  my ($help);
  my ($gff_file, $fasta_file, $index_dir);

  # TODO create the fasta transcriptome based on the GFF file
  # TODO index the fasta transcriptome with kallisto

  GetOptions(
      "help"             => \$help,
      # Input / Output
      "a|annotations=s"    => \$gff_file,
      "g|genome=s"         => \$fasta_file,
      "i|index=s"           => \$index_dir
  ) or usage();

  usage() if ($help);

  if(!defined $gff_file) {
    print STDERR "Missing gff file (-a,--annotations)\n";
    usage() && exit(1);
  }

  if(!defined $fasta_file) {
    print STDERR "Missing genome file (-g,--genome)\n";
    usage() && exit(1);
  }

  if(!defined $index_dir) {
    print STDERR "Missing index directory (-i,--index)\n";
    usage() && exit(1);
  }

  # Create output dir
  if(-e $index_dir) {
    print STDERR "Output directory already exists\n";
    exit(1);
  } else {
    mkdir $index_dir;
  }

  my $logs_dir    = "$index_dir/logs";
  mkdir $logs_dir;
  my $index_fasta = "$index_dir/genome.fa.gz";
  my $index_gff   = "$index_dir/annotations.gff3.gz";
  
  my $index_gff_fh = DEkupl::Utils::getWritingFileHandle($index_gff);

  print STDERR "Parsing FASTA file to load reference names\n";
  my $genome = DEkupl::Genome->new(fasta_file => $fasta_file);

  
  print STDERR "Saving gzipped FASTA file to index dir\n";
  if($fasta_file =~ /.gz$/) {
    system("cp $fasta_file $index_fasta");
  } else  {
    system("cat $fasta_file | gzip -c > $index_fasta");
  }

  foreach my $ref ($genome->allReferences) {
    print STDERR "reference found: $ref\n";
  }
  
  my $gff_it = DEkupl::Utils::gffFileIterator($gff_file,'gff3');

  # print headers
  my $gff_fh = DEkupl::Utils::getReadingFileHandle($gff_file);
  while(<$gff_fh>) {
    if($_ =~ /^#/) {
      print $index_gff_fh $_;
    } else {
      last;
    }
  }

  print STDERR "Parsing annotations to verify GFF format\n";
  my $nb_genes = 0;
  my $nb_exons = 0;
  my $nb_matching_ref = 0;
  while(my $annot = $gff_it->()) {
    # TODO Check that the sorting is okay
    # Sort by gene!!!
    if($genome->hasReference($annot->{chr})) {
      $nb_matching_ref++;
      if($annot->{feature} eq 'gene') {
        $nb_genes++;
      } elsif($annot->{feature} eq 'exon') {
        $nb_exons++;
      }
      print $index_gff_fh $annot->{_original_line},"\n";
    } else {
      print STDERR "Missing reference sequence ". $annot->{chr} . " in the genome\n";
    }
  }
  close($index_gff_fh);
  print STDERR "Found $nb_genes genes, and $nb_exons exons.\n";

  # TODO print an info json file with version and stuff!!!
  my $gsnap_log = "$logs_dir/gsnap.log";
  print STDERR "Creating GSNAP index (logs in $gnsap_log)\n";
  system("gmap_build -D $index_dir -d gsnap --gunzip $index_fasta 2> $gsnap_log");
}

sub usage {

  my $usage =<<END;
Usage:
    dkpl indx -g gff_file -f genome.fasta [-t transcript.fa] -i index_dir

  Options:
      -a,--annotations FILE   GFF annotation file
      -g,--genome FILE        Genome if FASTA format
      -i,--index DIR          Output index directory.

  Optional Arguments:
      -h,--help           show this help message and exit
END

  print $usage;
}

no Moose;
__PACKAGE__->meta->make_immutable;