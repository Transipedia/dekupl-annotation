package DEkupl::Utils;

use strict;
use warnings;

use Carp;
use Term::ANSIColor;

use Moose::Util::TypeConstraints;

subtype 'Strand',
  as 'Str',
  where { $_ =~ /^(\+|-)$/ },
  message { "Strand must be either '+' or '-'" };

our $NA_value = 'NA';

# scale a a value x defined on a range [min,max] to the range [a,b]
#
#        (b-a)(x - min)
# f(x) = --------------  + a
#          max - min
sub scaleValue {
  my ($a,$b,$min,$max,$x) = @_;
  return (($b-$a)*($x-$min))/($max-$min)+$a;
}

sub reverseStrand {
  my $strand = shift;
  if($strand eq '+') {
    return '-';
  }
  return '+';
}

sub getReadingFileHandle {
    my $file = shift;
    my $gzip = shift;
    $gzip = 0 if !defined $gzip;
    my $fh;
    if($gzip || $file =~ /\.gz$/) {
        open($fh, "gunzip -c $file |") or die "Cannot open $file";
    } else {
        open($fh, $file) or die "Cannot open $file";
    }
    return $fh;
}

sub getWritingFileHandle {
    my $file = shift;
    my $gzip = shift;
    $gzip = 0 if !defined $gzip;
    my $fh;
    if($gzip || $file =~ /\.gz$/) {
        open($fh, "| gzip -c > $file") or die "Cannot open $file";
    } else {
        open($fh, '>', $file) or die "Cannot open $file";
    }
    return $fh;
}

sub slurpFile {
    my $file = shift;
    my $data;
    {
        my $fh = getReadingFileHandle($file);
        local $/ = undef;
        $data = <$fh>;
        close $fh;
    }
    return $data;
}

sub saveToFile {
    my $file = shift;
    my $data = shift;
    my $fh = getWritingFileHandle($file);
    print $fh $data, "\n";
    close($fh);
}

sub booleanEncoding {
    my $v = shift;
    my $encoded_v = $v? 'T' : 'F';
}

sub bamFileIterator {
  my $file = shift;
  my $region = shift;
  $region = "" if !defined $region;
  
  open(my $fh, "-|", "samtools view $file $region" )or die "Cannot open $file, check if samtools are installed.";

  return sub {
    my $line = <$fh>;
    if($line) {
      return DEkupl::Utils::parseSAMLine($line);
    }
    return $line;
  }

}

=head2 gffFileIterator

manage GFF3 and GTF2 file format

Example:

  my $it = gffFileIterator($file,'type');
  while (my $annot = $it->()) {
    print "chr    : $annot->{chr}
           start  : $annot->{start}
           end    : $annot->{end}";
  }

Return a hashref with the annotation parsed:

  { chr         => 'field_1',
    source      => 'field_2',
    feature     => 'field_3',
    start       => 'field_4',
    end         => 'field_5',
    score       => 'field_6',
    strand      => 'field_7',
    frame       => 'field_8'
    attributes  => { 'attribute_id' => 'attribute_value', ...},
    seek_pos    => 'Seek position of this line in the file',
  }

gffFileIterator is B<5x faster than Bio-Perl> Bio::Tools::

=cut

sub gffFileIterator {
  my $file = shift;
  my $type = shift;
  return getFileIterator(file => $file, type => $type);
}

=head2 getFileIterator

Generic method to parse files.

=cut

sub getFileIterator {
  my %args = @_;
  my $file = $args{file};
  my $type = $args{type};
  my $skip = $args{skip};
  my $header_regex = $args{header_regex};
  my $parsing_method = $args{parsing_method};
  my @parsing_arguments = ();


  croak "Missing arguments in getFileIterator" if !defined $file;

  if(!defined $parsing_method && defined $type) {
    if($type =~ /gff3/i || $type eq 'gtf' || $type eq 'gff2') {
      $header_regex = '^#';
      $parsing_method = sub { return parseGFFLine(@_,$type) };
      push (@parsing_arguments,$type);
    } else {
      croak "Undefined format type";
    }
  }

  # set defaults
  #$header_regex = '^#' if !defined $header_regex;
  $parsing_method = sub{ return { line => shift } } if !defined $parsing_method;

  # We get a filehandle on the file
  my $fh = getReadingFileHandle($file);

  # $curr_pos will hold the SEEK_POS of the current_line
  my $curr_pos = tell($fh);
  # $line will hold the content of the current_line
  my $line = <$fh>;

  # if we need to skip lines
  if(defined $skip) {
    for(my $i = 0; $i < $skip; $i++) {
      $curr_pos = tell($fh);
      $line = <$fh>;
    }
  }

  # Skip line that match a specific regex
  if(defined $header_regex) {
    while($line && $line =~ /$header_regex/) {
      $curr_pos = tell($fh);
      $line = <$fh>;
    }
  }

  # The iterator itself (an unnamed subroutine
  return sub {
    # Skip lines starting with # even after the header was removed
    while(defined $line && $line =~ /^#/) {
      $line = <$fh>;
    }
    if (defined $line) { 
      chomp $line;
      # We parse the line with the appropriate methd
      my $parsed_line = $parsing_method->($line,@parsing_arguments);
      $parsed_line->{_original_line} = $line;

      $curr_pos = tell($fh);
      $line = <$fh>; # Get next line
      return $parsed_line;
    } else {
      return undef;
    }
  };
}

sub parseGFFLine {
  my $line = shift;
  my $type = shift;
  my $attribute_split;
  if($type =~ /gff3/i) {
    $attribute_split = sub {my $attr = shift; return $attr =~ /(\S+)=(.*)/;};
  } elsif ($type eq 'gtf' || $type eq 'gff2') {
    $attribute_split = sub {my $attr = shift; return $attr  =~ /(\S+)\s+"(.*)"/;};
  } else {
    die "Undefined GFF format (must be either gff3,gtf, or gff2)";
  }
  my($chr,$source,$feature,$start,$end,$score,$strand,$frame,$attributes) = split("\t",$line);

  # Load attributes into a hash
  my %attributes_hash;
  if(defined $attributes) {
    my @attributes_tab = split(";",$attributes);
    foreach my $attr (@attributes_tab) {
      my ($k,$v) = $attribute_split->($attr);
      $attributes_hash{$k} = $v;
    }
  }

  # This is a hack for GENCODE GFF files that use difference attribute field than ENSEMBL for gene names and biotypes
  my %conv_attr = (
    'gene_name' => 'Name',
    'gene_type' => 'biotype',
  );

  foreach my $attr (keys %conv_attr) {
    if(defined $attributes_hash{$attr} && !defined $attributes_hash{$conv_attr{$attr}}) {
      $attributes_hash{$conv_attr{$attr}} = $attributes_hash{$attr};
    }
  }
  
  return { 
    chr        => $chr,
    source     => $source,
    feature    => $feature,
    start      => $start,
    end        => $end,
    score      => $score,
    strand     => $strand,
    frame      => $frame,
    attributes => \%attributes_hash,
  };
}

sub parseSAMLine {
  my $line = shift;
  my ($qname,$flag,$rname,$pos,$mapq,$cigar,$rnext,$pnext,$tlen,$seq,$qual,@others) = split("\t",$line);
  my @cigar_hash = map { { op => substr($_,-1), nb => substr($_,0,length($_)-1)} } $cigar =~ /(\d+\D)/g;
  my %extended_fields = map { my ($id, $t, $v) = split ':', $_; $id => $v; } @others;
  return {
    qname => $qname,
    flag => $flag,
    rname => $rname,
    pos => $pos,
    mapq => $mapq,
    cigar => \@cigar_hash,
    original_cigar => $cigar,
    rnext => $rnext,
    pnext => $pnext,
    tlen => $tlen,
    seq => $seq,
    qual => $qual,
    extended_fields => \%extended_fields,
  };
}

sub getAtomicGeneID {
  my $gene_id = shift;
  my ($atomic_gene_id, $version) = $gene_id =~ /^(\S+)\.(\d+)$/;
  if(defined $atomic_gene_id && defined $version) {
    return $atomic_gene_id;
  }
  return $gene_id;
}

sub parseEnsemblID {
  my $ensembl_id = shift;
  my ($type,$id) = $ensembl_id =~ /(\S+):(\S+)/;
  if(defined $id) {
    return ($type,$id);
  }
  return $ensembl_id;
}

sub fastaFileIterator {
  my $file = shift;

  # Get file handle for $file
  my $fh = getReadingFileHandle($file);

  # Read prev line for FASTA because we dont know the number
  # of line used for the sequence
  my $prev_line = <$fh>;
  chomp $prev_line;
  return sub {
    my ($name,$seq,$qual);
    if(defined $prev_line) {
      ($name) = $prev_line =~ />(.*)$/;
      $prev_line = <$fh>;
      # Until we find a new sequence identifier ">", we
      # concatenate the lines corresponding to the sequence
      while(defined $prev_line && $prev_line !~ /^>/) {
        chomp $prev_line;
        $seq .= $prev_line;
        $prev_line = <$fh>;
      }
      return {name => $name, seq => $seq, qual => $qual};
    } else {
      return undef;
    }
  };
}

=head2 isVersionGreaterOrEqual($v1,$v2)

Return true is version number v1 is greater than v2

=cut

sub isVersionGreaterOrEqual($$) {
  my ($v1,$v2) = @_;
  my @v1_nums = split(/\./,$v1);
  my @v2_nums = split(/\./,$v2);
  for(my $i = 0; $i < @v1_nums; $i++) {
    if($v1_nums[$i] >= $v2_nums[$i]) {
      return 1;
    } else {
      return 0;
    }
  }
  if(scalar @v2_nums > @v1_nums) {
    return 0;
  } else {
    return 1;
  }
}

sub checkSamtoolsVersion {
  my $min_version = "1.3";
  if(system("samtools --version > /dev/null") == 0) {
    my ($version) = `samtools --version` =~ /samtools\s(\S+)/;
    die "samtools version ($version) is outdated (version >= $min_version)" if !isVersionGreaterOrEqual($version,$min_version);
    print STDERR "Using SAMtools version $version\n";
  } else {
    die "samtools is not acessible in the \$PATH";
  }
}

sub checkGSNAPVersion {
  my $min_version = "2016.11.07";
  if(system("gsnap --version 2> /dev/null > /dev/null") == 0) {
    my ($version) = `gsnap --version 2> /dev/null` =~ /version\s(\S+)/;
    $version =~ s/-/\./;
    die "gsnap version ($version) is outdated (version >= $min_version)" if !isVersionGreaterOrEqual($version,$min_version);
    print STDERR "Using GSNAP version $version\n";
  } else {
    die "samtools is not acessible in the \$PATH";
  }
}

sub checkBlastnVersion {
  my $min_version = "2.5.0";
  if(system("blastn -version > /dev/null") == 0) {
    my ($version) = `blastn -version` =~ /blastn:\s(\S+)\+/;
    die "blastn version ($version) is outdated (version >= $min_version)" if !isVersionGreaterOrEqual($version,$min_version);
    print STDERR "Using blastn version $version\n";
  } else {
    die "blastn is not acessible in the \$PATH";
  }
}

sub checkMakeblastdbVersion {
  my $min_version = "2.5.0";
  if(system("makeblastdb -version > /dev/null") == 0) {
    # makeblastdb: 2.7.1+
    # Package: blast 2.7.1, build Feb 16 2018 14:27:59
    my ($version) = `makeblastdb -version` =~ /makeblastdb:\s(\S+)\+/;
    die "makeblastdb version ($version) is outdated (version >= $min_version)" if !isVersionGreaterOrEqual($version,$min_version);
    print STDERR "Using makeblastdb version $version\n";
  } else {
    die "makeblastdb is not acessible in the \$PATH";
  }
}

sub checkSTARVersion {
  my $min_version = "2.5.3";
  if(system("STAR --version > /dev/null") == 0) {
    # STAR_2.5.3a
    my ($version) = `STAR --version` =~ /STAR_(\S+)[a-z]$/;
    die "STAR version ($version) is outdated (version >= $min_version)" if !isVersionGreaterOrEqual($version,$min_version);
    print STDERR "Using STAR version $version\n";
  } else {
    die "STAR is not acessible in the \$PATH";
  }
}

sub printStep {
  my ($step, $message) = @_;
  print STDERR color('bold blue');
  print STDERR "[Step ".++$$step."]";
  print STDERR color('reset');
  print STDERR " $message\n";
}

1;