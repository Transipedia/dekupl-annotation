package DEkupl::Utils;

use strict;
use warnings;

use Carp;

use Moose::Util::TypeConstraints;

subtype 'Strand',
  as 'Str',
  where { $_ =~ /^(\+|-)$/ },
  message { "Strand must be either '+' or '-'" };

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

1;