package DEkupl::App::Annot;
# ABSTRACT: dkpl-annot binary

use Moose;
use Getopt::Long;
use File::Temp qw/ tempdir /;

use DEkupl;
use DEkupl::GSNAP;
use DEkupl::Contigs;
use DEkupl::ContigsDB;
use DEkupl::Annotations;
use DEkupl::IntervalQuery;

use DEkupl::Analyzer::BAM;
use DEkupl::Analyzer::Annotations;
use DEkupl::Analyzer::DEG;
use DEkupl::Analyzer::Switches;

sub BUILD {
  my $self = shift;

  my ($help);
  my ($contigs_file, $index_dir, $deg_file, $version);
  my $normalized_gene_counts_file;
  my $sample_conditions_file;
  my $output_dir = "DEkupl_annotation";
  my $index_name = "GRCh38-chr22";
  my $nb_threads = 1;
  my $tmp_dir = "/tmp";
  my $debug = 1;
  my $is_stranded = 0;
  my $deg_padj_threshold = 0.05;

  GetOptions(
      "help"             => \$help,
      # Input / Output
      #"c|contigs=s"      => \$contigs_file,
      "o|output=s"       => \$output_dir,
      "i|index=s"        => \$index_dir,
      "d|deg=s"          => \$deg_file,
      "norm-gene-counts=s" => \$normalized_gene_counts_file,
      "sample-conditions=s" => \$sample_conditions_file,

      # Options
      "t|threads=i"      => \$nb_threads,
      "s|stranded"       => \$is_stranded,
      "p|deg-padj"       => \$deg_padj_threshold,
      "version"          => \$version,
  ) or pod2usage(-verbose => 1);

  usage() if ($help);

  if($version) {
    print STDERR "rankvar2 ", $DEkupl::VERSION, "\n";
    exit 0;
  }

  $contigs_file = shift @ARGV;

  if(!defined $contigs_file) {
    print STDERR "Missing contig file (-c,--contigs)\n";
    usage() && exit(1);
  }

  if(!defined $index_dir) {
    print STDERR "Missing index directory (-i,--index)\n";
    usage() && exit(1);
  }

  # Create output dir
  mkdir $output_dir;

  # Find files from index
  # TODO we should have a class for that
  # that load the config file and check that everything is ok
  my $gff_file = "$index_dir/annotations.gff3.gz";
  my @analyzers;

  # Create contigs object
  my $contigs = DEkupl::Contigs->new(contigs_file => $contigs_file, is_stranded => $is_stranded);
  push @analyzers, $contigs;

  my $gsnap = DEkupl::GSNAP->new(
    index_dir   => $index_dir,
    index_name  => "gsnap",
    nb_threads  => $nb_threads,
  );

  # Create FASTA file
  my $fasta_file = "$output_dir/contigs.fa.gz";
  if(-e $fasta_file && $debug) {
    print STDERR "Skipping FASTA file\n";
  } else {
    print STDERR "Generating FASTA file\n";
  $contigs->generateFasta($fasta_file);
  }

  # Generate BAM file
  my $bam_file = "$output_dir/contigs.bam";
  if(-e $bam_file && $debug) {
    print STDERR "Skipping GSNAP mapping\n";
  } else {
    print STDERR "Running GSNAP\n";
    $gsnap->generateBam($fasta_file,$bam_file);
  }

  # Create contigs DB
  #my $tempdir = tempdir(CLEANUP => 1);
  my $tempdir = tempdir("$tmp_dir/dkplannot_tmp.XXXXX",CLEANUP => 1);
  my $contigs_db = DEkupl::ContigsDB->new(db_folder => $tempdir);
  print STDERR "Loading contigs DB into $tempdir\n";
  $contigs->loadContigsDB($contigs_db);

  print STDERR "Parsing BAM file\n";
  my $bed_file = "$output_dir/diff_contigs.bed.gz";
  my $bam_analyzer = DEkupl::Analyzer::BAM->new(
    contigs     => $contigs,
    contigs_db  => $contigs_db,
    bam_file    => $bam_file,
    is_stranded => $is_stranded,
    bed_file    => $bed_file,
  );
  push @analyzers, $bam_analyzer;

  # TODO: We should do that in separate space, in order to unload annotations from
  # memory as soon as the annotation is done remove them after annotation.
  print STDERR "Loading annotations into memory\n";
  my $annotations = DEkupl::Annotations->new();
  $annotations->loadFromGFF($gff_file,'gff3');

  print STDERR "Loading annotations to the interval tree\n";
  my $interval_query = DEkupl::IntervalQuery->new();
  $interval_query->loadAnnotations($annotations);

  print STDERR "Annotating contigs\n";
  my $loci_file = "$output_dir/ContigsPerLoci.tsv.gz";
  my $annot_analyzer = DEkupl::Analyzer::Annotations->new(
    contigs_db      => $contigs_db,
    interval_query  => $interval_query,
    is_stranded     => $is_stranded,
    loci_file       => $loci_file,
  );
  push @analyzers, $annot_analyzer;

  # If we have DEG (Differentially expressed genes, we add appened information to the contigs)
  if(defined $deg_file) {
    print STDERR "Adding DEG informations\n";
    my $deg_analyzer = DEkupl::Analyzer::DEG->new(
      contigs_db  => $contigs_db,
      deg_file    => $deg_file,
      is_stranded => $is_stranded,
      max_padj    => $deg_padj_threshold,
    );
    push @analyzers, $deg_analyzer;
  }

  if(defined $normalized_gene_counts_file && defined $sample_conditions_file) {
    print STDERR "Computing switches\n";
    my @sample_names = $contigs->all_samples;
    my $switches_analyzer = DEkupl::Analyzer::Switches->new(
      contigs_db                  => $contigs_db,
      is_stranded                 => $is_stranded,
      sample_names                => \@sample_names,
      normalized_gene_counts_file => $normalized_gene_counts_file,
      sample_conditions_file      => $sample_conditions_file,
    );
    push @analyzers, $switches_analyzer;
  }


  # TODO we should print automatic headers with descriptions of each fields as
  # in the VCF format!!!
  my $contigs_info_file = "$output_dir/DiffContigsInfos.tsv";
  my $contigs_info_fh   = DEkupl::Utils::getWritingFileHandle($contigs_info_file);

  # Print headers
  print $contigs_info_fh join("\t",
    map { $_->getHeaders() } @analyzers
  ),"\n";

  # Print values
  my $contigs_it = $contigs_db->contigsIterator;
  while (my $contig = $contigs_it->()) {
    print $contigs_info_fh join("\t",
      map { $_->getValues($contig) } @analyzers
    ),"\n";
  }

  close($contigs_info_fh);

}

sub usage {

  my $usage =<<END;
Usage:
    dkpl annot -i index_dir/ merged-diff-counts.tsv.gz

Options:
  Requiered Arguments:
      -i,--index DIR      path to the index directory (created with dkplannot index)

  Input/Output:
      -o,--output DIR     path to the output directory (default: "DEkupl_annotation/")
      -d,--deg FILE       (Optional) {A}vs{B}-DEGs.tsv (diff. genes in "gene_expression" directory from Dekupl-run result)
      --norm-gene-counts FILE 
                          (Optional) Normalized gene counts
      --sample-conditions FILE 
                          (Optional) Sample conditions. First column is sample name,
                          second column is sample condition)

  Optional Arguments:
      -t,--threads INT    Number of threads (for GSNAP)
      -s,--stranded       RNA-Seq is strand-specific.
      -p,--deg-padj       padj diff. gene threshold (default : 0.05)

      -h,--help           show this help message and exit
END

  print $usage;
}

no Moose;
__PACKAGE__->meta->make_immutable;