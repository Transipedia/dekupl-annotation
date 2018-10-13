package DEkupl::Aligner::STAR;
# ABSTRACT: Run STAR Aligner

use Moose;
use File::Temp;

with 'DEkupl::Base';

has 'index_dir' => (
  is => 'ro',
  isa => 'Str',
  required => 1,
);

has 'nb_threads' => (
  is => 'rw',
  isa => 'Int',
  default => 1,
);

has 'max_splice_length' => (
  is => 'ro',
  isa => 'Int',
  default => 100000,
);

# TODO We should verify GSNAP index directory
sub BUILD {
  my $self = shift;
  DEkupl::Utils::checkSTARVersion();
}

sub generateChimericJunctions {
  my $self = shift;
  my $fasta_input_file = shift;
  my $output_directory = shift;

  my $logs = File::Temp->new(UNLINK => 0, SUFFIX => '.log');
  $self->verboseLog("STAR logs are located in $logs");

  # Create a temp fasta that cut contigs longer than 650bp that causes
  # an error in STAR
  # FIXME This shoud be better handled !!!
  my $fasta_tmp = File::Temp->new(UNLINK => 0, SUFFIX => '.fa.gz');
  close $fasta_tmp;
  $self->verboseLog("Tmp fasta files with truncated contigs in $fasta_tmp");
  my $fasta_fh = DEkupl::Utils::getWritingFileHandle($fasta_tmp);
  my $seq_it = DEkupl::Utils::fastaFileIterator($fasta_input_file);
  while(my $entry = $seq_it->()) {
    if(length($entry->{seq}) > 650) {
      $entry->{seq} = substr $entry->{seq}, 0, 650;
    # We found bug in STAR for short contigs, causing segmentation fault
    # Also having contigs > 50 will result in more accurate chimeric detection
    } elsif(length($entry->{seq}) < 50) {
      next;
    }
    print $fasta_fh ">".$entry->{name}."\n".$entry->{seq}."\n";
  }
  close($fasta_fh);

  my $command = join(" ",
    "STAR",
    "--runThreadN", $self->nb_threads,
    "--genomeDir", $self->index_dir,
    "--readFilesIn", $fasta_tmp,
    "--readFilesCommand gunzip -c",
    "--outFileNamePrefix $output_directory/",
    "--twopassMode Basic",
    "--outSAMtype None",
    "--outReadsUnmapped None",
    #"--chimSegmentMin 30", # stringent-mode
    #"--chimJunctionOverhangMin 12", # stringent-mode

    # These parameters are taken from STAR-Fusion
    # See https://github.com/STAR-Fusion/STAR-Fusion/wiki
    #"--alignMatesGapMax 100000",
    "--alignIntronMax", $self->max_splice_length,
    "--chimSegmentMin 12",
    "--chimJunctionOverhangMin 12",
    "--alignSJDBoverhangMin 10",
    "--chimSegmentReadGapMax 3",
    "--alignSJstitchMismatchNmax 5 -1 5 5",
    "2> $logs",
    ">> $logs"
  );

  $self->verboseLog("Executing STAR ...");
  $self->verboseLog("Command: $command");
  
  # run GSNAP verify that the execution of GSNAP ended well
  system($command) == 0 or die("STAR alignement failed, see logs in $logs");
  unlink $logs;

  return "$output_directory/Chimeric.out.junction";
}

no Moose; 
__PACKAGE__->meta->make_immutable;