package DEkupl::Utils;

use strict;
use warnings;

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

1;