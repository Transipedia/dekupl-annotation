package DEkupl::App::Index;
# ABSTRACT: dkpl-index binary

use Moose;
use Getopt::Long;
use File::Temp qw/ tempdir /;

use DEkupl;
use DEkupl::Annotations;
use DEkupl::Genome;

with 'DEkupl::Base';

sub BUILD {
  my $self = shift;

  my ($help);
  my ($gff_file, $fasta_file, $index_dir, $star);
  my $nb_threads = 1;

  # TODO create the fasta transcriptome based on the GFF file
  # TODO index the fasta transcriptome with kallisto

  GetOptions(
      "help"             => \$help,
      # Input / Output
      "a|annotations=s" => \$gff_file,
      "g|genome=s"      => \$fasta_file,
      "i|index=s"       => \$index_dir,
      "star"            => \$star,
      "t|threads=i"      => \$nb_threads,
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

  my $step = 0;
  my $logs_dir    = "$index_dir/logs";
  my $index_fasta = "$index_dir/genome.fa.gz";
  my $index_gff   = "$index_dir/annotations.gff3.gz";

  # Create output dir
  if(-e $index_dir) {
    $self->log("Output directory already exists");
  } else {
    mkdir $index_dir;
  }

  mkdir $logs_dir if ! -e $logs_dir;

  DEkupl::Utils::printStep(\$step, "Parsing FASTA file to load reference names");
  my $genome = DEkupl::Genome->new(fasta_file => $fasta_file);
  foreach my $ref ($genome->allReferences) {
    print STDERR "reference found: $ref\n";
  }
  
  if(!-e $index_fasta) {
    DEkupl::Utils::printStep(\$step, "Saving gzipped FASTA file to index dir");
    if($fasta_file =~ /.gz$/) {
      system("cp $fasta_file $index_fasta");
    } else {
      system("cat $fasta_file | gzip -c > $index_fasta");
    }
  } else {
    DEkupl::Utils::printStep(\$step, "Skip FASTA file (already exists)");
  }
  
  if(-e $index_gff) {
    DEkupl::Utils::printStep(\$step, "Skip GFF file (already exists)");
  } else {
    DEkupl::Utils::printStep(\$step, "Saving GFF file (already exists)");
    my $index_gff_fh = DEkupl::Utils::getWritingFileHandle($index_gff);
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
    my $nb_gene_names = 0;
    my $nb_gene_biotype = 0;
    my %missing_ref;
    my %ids; # in order to verify that each gene ID is uniq!
    while(my $annot = $gff_it->()) {

      # TODO Check that the sorting is okay
      # Sort by gene!!!
      if($genome->hasReference($annot->{chr})) {
        $nb_matching_ref++;
        if(defined $DEkupl::Annotations::gene_features{$annot->{feature}}) {
          my $id = $annot->{attributes}->{ID};
          if(defined $id) {
            if(defined $ids{$id}) {
              die "ERROR: Duplicated gene ID ($id) in GFF annotations";
            }
            $ids{$id} = 1;
          }
          $nb_genes++;
          $nb_gene_names++ if defined $annot->{attributes}->{Name};
          $nb_gene_biotype++ if defined $annot->{attributes}->{biotype};
        } elsif($annot->{feature} eq 'exon') {
          $nb_exons++;
        }
        # Print any line that match the FASTA references set
        print $index_gff_fh $annot->{_original_line},"\n";
      } else {
        if(!defined $missing_ref{$annot->{chr}}) {
          print STDERR "Missing reference sequence ". $annot->{chr} . " in the genome\n";
        }
        $missing_ref{$annot->{chr}} = 1;
      }
    }
    close($index_gff_fh);

    die "ERROR: No gene feature found in GFF file" if $nb_genes == 0;
    die "ERROR: Reference sequences between genome (FASTA) and annotations (GFF) are not matching" if $nb_matching_ref == 0;

    print STDERR "Found $nb_genes genes, and $nb_exons exons.\n";
    print STDERR "No gene symbol found ('Name' attribute in GFF)\n" if $nb_gene_names == 0;
    print STDERR "No gene biotype found ('biotype' attribute in GFF)\n" if $nb_gene_biotype == 0;
  }
  
  if(-e "$index_dir/gsnap") {
    DEkupl::Utils::printStep(\$step, "Skipping GSNAP index");
  } else {
    # TODO print an info json file with version and stuff!!!
    my $gsnap_log = "$logs_dir/gsnap.log";
    DEkupl::Utils::printStep(\$step, "Creating GSNAP index (logs in $gsnap_log). This can take a while...");
    system(join(' ',
      "gmap_build",
      "-D $index_dir",
      "-d gsnap",
      "--gunzip",
      $index_fasta,
      ">> $gsnap_log 2>&1"
    )) == 0 or die "Error in GSNAP index building (see logs in $gsnap_log)";
  }

  if($star) {
    if(-e "$index_dir/star") {
      DEkupl::Utils::printStep(\$step, "Skipping STAR index)");
    } else {
      mkdir "$index_dir/star";
      my $star_log = "$logs_dir/star.log";
      DEkupl::Utils::printStep(\$step, "Creating STAR index. This can take a while");
      my $tmp_fasta = File::Temp->new(UNLINK => 1, SUFFIX => '.fa');
      system("gunzip -c $index_fasta > $tmp_fasta");

      # TODO We should use the Annotations when construction STAR index!!!
      my $command = join(' ',
        "STAR",
        "--runThreadN $nb_threads",
        "--runMode genomeGenerate",
        "--genomeDir $index_dir/star",
        "--genomeFastaFiles $tmp_fasta",
        "2> $star_log"
      );
      system($command) == 0 or die "Error in STAR index building (see logs in $star_log)";
    }
  }
}

sub usage {

  my $usage =<<END;
Usage:
    dkpl indx -g gff_file -f genome.fasta -i index_dir

  Mandatory Arguments:
      -a,--annotations FILE   GFF annotation file
      -g,--genome FILE        Genome if FASTA format
      -i,--index DIR          Output index directory.

  Optional Arguments:
      -h,--help           show this help message and exit
      -t,--threads INT    Number of threads
      --star              Index the genome with STAR (for discovery of chimeric RNA)
                          It needs ~30gb of RAM for Human genome
END

  print $usage;
}

no Moose;
__PACKAGE__->meta->make_immutable;