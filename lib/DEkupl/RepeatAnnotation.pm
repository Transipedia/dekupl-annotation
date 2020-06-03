package DEkupl::RepeatAnnotation;
# ABSTRACT: Run Blast Aligner & Enthropy Script to find Repeat Sequences

use Moose;
use File::Temp;

with 'DEkupl::Base';

has 'nb_threads' => (
  is => 'rw',
  isa => 'Int',
  default => 1
);

has 'fasta' => (
  is => 'ro',
);

my @columns = ('rep_type');

sub BUILD {
  my $self = shift;
  DEkupl::Utils::checkBlastnVersion();
  DEkupl::Utils::checkMakeblastdbVersion();
}

sub _createBlastIndex {
  my $self = shift;
  my $fasta = $self->fasta;

  my $blast_index_basename = File::Temp->new(UNLINK => 1, SUFFIX => '');
  my $logs = File::Temp->new(UNLINK => 0, SUFFIX => '.log');
  $self->verboseLog("Create BLAST Index, basename $blast_index_basename");
  my $command = join(" ",
    "makeblastdb",
    "-in $fasta",
    "-dbtype nucl",
    "-out $blast_index_basename",
    "2> $logs",
    ">> $logs",
  );
  system($command) == 0 or die("Blast db creation failed, see logs in $logs");
  unlink $logs;
  return $blast_index_basename;
}

sub alignRepeats {
  my $self = shift;
  my $contigs_db = shift;
  my $results = shift;
  # Create blast index
  my $blast_db = $self->_createBlastIndex();

  # Create temporary files
  my $contigs_fasta = File::Temp->new(UNLINK => 0, SUFFIX => '.fa');
  my $logs          = File::Temp->new(UNLINK => 0, SUFFIX => '.log');

  # Generate FASTA file with contigs from contigs db
  $contigs_db->generateFasta($contigs_fasta);

  # Run Blastn
  my $command = join(" ",
    "blastn",
    "-query $contigs_fasta",
    "-db $blast_db",
    "-word_size 11",
    "-max_hsps 1",
    "-max_target_seqs 1",
    "-dust no",
    "-soft_masking false",
    "-task blastn",
    "-evalue 1e-3",
    '-outfmt "6 qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore qlen sstrand"',
    "-out $results/resultblast.txt",
    "-num_threads " . $self->nb_threads,
    "2> $logs"
  );
  system($command) == 0 or die("Blast of contigs against repeats failed, see logs in $logs");
  unlink $logs;
  unlink $contigs_fasta;

  #Selecting contigs aligned to the Repeat database
  my %aligned_contigs;
  my %contig_dict;
  {
    open(my $fh, '<', "$results/resultblast.txt") or die("Cannot open $results");
    while(<$fh>) {
      chomp;
      # ACACTCTTTCCCTACACGACGCTGTTCCATCT        SINE1   100.000 32      0       0       1       32      1       32      1.77e-13        59.0    32      plus
      my ($contig_id, $repeat_id, $p_identity) = split "\t", $_;
      $aligned_contigs{$contig_id}++;
      $contig_dict{$contig_id} = $repeat_id;
    }
  }

  my $contigs_it = $contigs_db->contigsIterator;
  while(my $contig = $contigs_it->()) {
    if(exists($contig_dict{$contig->{tag}})) {
      my $contig_load = $contigs_db->loadContig($contig->{tag});
      $contig_load->{rep_type} = $contig_dict{$contig->{tag}};
      $contigs_db->saveContig($contig_load);
    } else {
      my $contig_load = $contigs_db->loadContig($contig->{tag});
      $contig_load->{rep_type} = "NA";
      $contigs_db->saveContig($contig_load);
    }
  }

}

sub _log2 {
  my $self = shift;
  return log($self) / log(2);
}

sub _stringEntropy {
  my $self = shift;
  my @chars = split("",$self);
  my $nchar = length($self);
  my %counts;
  $counts{$_}++ foreach @chars;
  my @unic_count;
  while (my ($key, $value) = each(%counts)) {
        $value = $value/$nchar;
        push(@unic_count,$value);
        };

  my $entropy = 0;
  foreach (@unic_count) {
    if ($_ != 0 ) {
      my $log2 = _log2($_);
      $entropy = $entropy + $_*$log2;
    };
  };
  $entropy = -$entropy;
  return $entropy;
}

sub contigsEntropy {
  my $self = shift;
  my $results = shift;
  my $single_repeat_file = "$results/ContigsEntropy.csv";
  my $fh = DEkupl::Utils::getWritingFileHandle($single_repeat_file);
  my $contigs_it = $self->contigsIterator;
  while(my $contig = $contigs_it->()) {
    my $contig_entropy = _stringEntropy($contig->{contig});
    my $simple_repeat="";
    my $contig_load = $self->loadContig($contig->{tag});
    if ($contig_entropy <= 1.5) {
      $simple_repeat="Simple Repeat";
      $contig_load->{rep_type} = "Simple Repeat";
    }
    $self->saveContig($contig_load);
    my $CE_line = join("\t",
      $contig->{tag},
      $contig->{contig},
      $contig_entropy,
      $simple_repeat
    );
    print $fh $CE_line,"\n";
  }
  close($fh);
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
