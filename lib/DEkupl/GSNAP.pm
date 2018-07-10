package DEkupl::GSNAP;
# ABSTRACT: Run GSNAP Aligner

use Moose;
use File::Temp;

has 'index_dir' => (
  is => 'ro',
  isa => 'Str',
  required => 1,
);

has 'index_name' => (
  is => 'ro',
  isa => 'Str',
  required => 1,
);

has 'nb_threads' => (
  is => 'rw',
  isa => 'Int',
  default => 1
);

# TODO We should verify GSNAP index directory
sub BUILD {
  my $self = shift;
  DEkupl::Utils::checkSamtoolsVersion();
  DEkupl::Utils::checkGSNAPVersion();
}

sub generateBam {
  my $self = shift;
  my $fata_input_file = shift;
  my $output_file = shift;

  my $logs = File::Temp->new(UNLINK => 0, SUFFIX => '.log');
  print STDERR "GSNAP logs are located in $logs\n";

  my $command = join(" ",
    "gsnap -t ", $self->nb_threads,
    "-A sam", # -A, --format=STRING            Another format type, other than default.
    "-N 1", # -N, --novelsplicing=INT              Look for novel splicing (0=no (default), 1=yes)
    "-D", $self->index_dir, #  -D, --dir=directory            Genome directory.
    "-d", $self->index_name, # -d, --db=STRING                Genome database
    #"-w 50000", # -w, --localsplicedist=INT            Definition of local novel splicing event (default 200000)
    "--gunzip", # Uncompress gzipped input files
    $fata_input_file,
    "2> $logs",
    "| samtools view -bh >", # Convert output SAM to BAM format with samtools
    $output_file
  );

  print STDERR "Executing GSNAP\nCommand: $command\n",
  
  # run GSNAP verify that the execution of GSNAP ended well
  system($command) == 0 or die("GSNAP alignement failed, see logs in $logs");
  unlink $logs;
}

no Moose;
__PACKAGE__->meta->make_immutable;