package DEkupl::App::Annot;
# ABSTRACT: dkpl-annot binary

use Moose;
use Getopt::Long;
use File::Temp qw/ tempdir /;

use DEkupl;
use DEkupl::Utils;
use DEkupl::RemoveAdapters;
use DEkupl::Contigs;
use DEkupl::ContigsDB;
use DEkupl::Annotations;
use DEkupl::IntervalQuery;

use DEkupl::Aligner::GSNAP;
use DEkupl::Aligner::STAR;

use DEkupl::Analyzer::BAM;
use DEkupl::Analyzer::Annotations;
use DEkupl::Analyzer::DEG;
use DEkupl::Analyzer::Switches;
use DEkupl::Analyzer::ChimericRNA;

sub BUILD {
  my $self = shift;

  my ($help,$verbose);
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
  my $contig_color_mode = 1;
  my $max_splice_length = 100000;
  my $adapters_fasta;

  GetOptions(
      "help"             => \$help,
      "v|verbose"        => \$verbose,
      # Input / Output
      #"c|contigs=s"      => \$contigs_file,
      "o|output=s"       => \$output_dir,
      "i|index=s"        => \$index_dir,
      "d|deg=s"          => \$deg_file,
      "norm-gene-counts=s" => \$normalized_gene_counts_file,
      "sample-conditions=s" => \$sample_conditions_file,

      # Options
      "t|threads=i"         => \$nb_threads,
      "s|stranded"          => \$is_stranded,
      "p|deg-padj=f"        => \$deg_padj_threshold,
      "max-splice-length=i" => \$max_splice_length,
      "contig-color=i"      => \$contig_color_mode,
      "adapters=s"          => \$adapters_fasta,
      "version"             => \$version,
  ) or pod2usage(-verbose => 1);

  usage() if ($help);

  if($version) {
    print STDERR "rankvar2 ", $DEkupl::VERSION, "\n";
    exit 0;
  }

  $contigs_file = shift @ARGV;

  if(!defined $contigs_file) {
    print STDERR "Missing contig file (aka. merged-diff-counts.tsv.gz)\n";
    usage() && exit(1);
  }

  if(!defined $index_dir) {
    print STDERR "Missing index directory (-i,--index)\n";
    usage() && exit(1);
  }

  if($contig_color_mode != 1 && $contig_color_mode != 2) {
    die "Invalid value for --contig-color option (must be 1 or 2)\n";
  }

  # Create output dir
  unless (-e $output_dir) {
     mkdir $output_dir  or die "Cannot create the output directory";
  }

  my $step = 0;

  # Find files from index
  # TODO we should have a class for that
  # that load the config file and check that everything is ok
  my $gff_file = "$index_dir/annotations.gff3.gz";
  my @analyzers;

  # Create contigs object
  my $contigs = DEkupl::Contigs->new(
    verbose       => $verbose,
    contigs_file  => $contigs_file, 
    is_stranded   => $is_stranded
  );
  push @analyzers, $contigs;

  # Create contigs DB
  my $tempdir = tempdir("$tmp_dir/dkplannot_tmp.XXXXX", CLEANUP => 1);
  my $contigs_db = DEkupl::ContigsDB->new(db_folder => $tempdir);
  DEkupl::Utils::printStep(\$step,"Loading contigs DB into $tempdir");
  $contigs->loadContigsDB($contigs_db);

  # Remove contigs matching adapters
  my $remove_adapters = DEkupl::RemoveAdapters->new(verbose => $verbose, fasta => $adapters_fasta);
  DEkupl::Utils::printStep(\$step,"Remove contigs matching adapters");
  $remove_adapters->removeAdapterContigs($contigs_db);

  # Create FASTA file
  my $fasta_file = "$output_dir/contigs.fa.gz";
  if(-e $fasta_file && $debug) {
    DEkupl::Utils::printStep(\$step,"Skipping FASTA file");
  } else {
    DEkupl::Utils::printStep(\$step,"Generating FASTA file");
    $contigs->generateFasta($fasta_file);
  }

  # Generate BAM with GSNAP
  my $bam_file = "$output_dir/contigs.bam";
  if(-e $bam_file && $debug) {
    DEkupl::Utils::printStep(\$step,"Skipping GSNAP mapping");
  } else {
    DEkupl::Utils::printStep(\$step,"Running GSNAP");
    my $gsnap = DEkupl::Aligner::GSNAP->new(
      index_dir   => $index_dir,
      index_name  => "gsnap",
      nb_threads  => $nb_threads,
    );
    $gsnap->generateBam($fasta_file, $bam_file);
  }

  # Annotating contigs with BAM from GSNAP
  DEkupl::Utils::printStep(\$step,"Parsing BAM file");
  my $bed_file = "$output_dir/diff_contigs.bed.gz";
  my $bam_analyzer = DEkupl::Analyzer::BAM->new(
    verbose           => $verbose,
    contigs           => $contigs,
    contigs_db        => $contigs_db,
    bam_file          => $bam_file,
    is_stranded       => $is_stranded,
    bed_file          => $bed_file,
    contig_color_mode => $contig_color_mode,
  );
  push @analyzers, $bam_analyzer;

  # Generate Chimeric junctions with STAR (if index is available)
  if(-e "$index_dir/star") {
    my $chimeric_file = "$output_dir/STAR/Chimeric.out.junction";
    if(-e $chimeric_file) {
      DEkupl::Utils::printStep(\$step,"Skipping STAR mapping");
    } else {
      DEkupl::Utils::printStep(\$step,"Running STAR");
      my $star = DEkupl::Aligner::STAR->new(
        index_dir   => "$index_dir/star",
        nb_threads  => $nb_threads,
        verbose     => $verbose,
      );
      mkdir "$output_dir/STAR" if !-e "$output_dir/STAR/";
      $star->generateChimericJunctions($fasta_file,"$output_dir/STAR");
    }
    my $chimeric_analyzer = DEkupl::Analyzer::ChimericRNA->new(
      verbose           => $verbose,
      contigs_db        => $contigs_db,
      is_stranded       => $is_stranded,
      chimeric_file     => $chimeric_file,
      max_splice_length => $max_splice_length,
    );
    push @analyzers, $chimeric_analyzer;
  }

  # TODO: We should do that in separate space, in order to unload annotations from
  # memory as soon as the annotation is done remove them after annotation.
  DEkupl::Utils::printStep(\$step,"Loading annotations into memory");
  my $annotations = DEkupl::Annotations->new(verbose => $verbose);
  $annotations->loadFromGFF($gff_file,'gff3');

  DEkupl::Utils::printStep(\$step,"Loading annotations to the interval tree");
  my $interval_query = DEkupl::IntervalQuery->new();
  $interval_query->loadAnnotations($annotations);

  DEkupl::Utils::printStep(\$step,"Annotating contigs");
  my $loci_file = "$output_dir/ContigsPerLoci.tsv.gz";
  my $annot_analyzer = DEkupl::Analyzer::Annotations->new(
    verbose         => $verbose,
    contigs_db      => $contigs_db,
    interval_query  => $interval_query,
    is_stranded     => $is_stranded,
    loci_file       => $loci_file,
  );
  push @analyzers, $annot_analyzer;

  # If we have DEG (Differentially expressed genes, we add appened information to the contigs)
  if(defined $deg_file) {
    DEkupl::Utils::printStep(\$step,"Adding DEG informations");
    my $deg_analyzer = DEkupl::Analyzer::DEG->new(
      verbose     => $verbose,
      contigs_db  => $contigs_db,
      deg_file    => $deg_file,
      is_stranded => $is_stranded,
      max_padj    => $deg_padj_threshold,
    );
    push @analyzers, $deg_analyzer;
  }

  if(defined $normalized_gene_counts_file && defined $sample_conditions_file) {
    DEkupl::Utils::printStep(\$step,"Computing switches");
    my @sample_names = $contigs->all_samples;
    my $switches_analyzer = DEkupl::Analyzer::Switches->new(
      verbose                     => $verbose,
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

  DEkupl::Utils::printStep(\$step,"Printing final output in $contigs_info_file");
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
      --max-splice-length 
                          Splice with greater length are considered as chimeric junctions (default 100000)
      --contig-color INT  Contig color mode (default 1):
                            1 : contigs on forward strand are in red (contigs on reverse strand are in blue)
		                        2 : contigs on forward strand are in blue (contigs on reverse strand are in red)
      --adapters STR      FASTA file with sequence to filter out (default : illumina adapters)
      -h,--help           show this help message and exit
      -v,--verbose        print additionnal debug messages
END

  print $usage;
}

no Moose;
__PACKAGE__->meta->make_immutable;