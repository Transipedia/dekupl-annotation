package DEkupl::ContigsDB;
# ABSTRACT: Store informations for contigs

use Moose;
use JSON::XS;

use DEkupl::Utils;

has 'db_folder' => (
  is => 'ro',
  isa => 'Str',
  required => 1,
);

sub saveContig {
  my $self = shift;
  my $contig = shift;
  my $tag = $self->_getTag($contig);
  my $file = $self->_getContigFile($tag);
  my $json = encode_json $contig;
  DEkupl::Utils::saveToFile($file, $json);
}

sub loadContig {
  my $self = shift;
  my $tag = shift;
  my $file = $self->_getContigFile($tag);
  my $contig;
  if(-e $file) {
    my $json = DEkupl::Utils::slurpFile($file);
    $contig = decode_json $json;
  }
  return $contig;
}

sub _getTag {
  my $self = shift;
  my $contig = shift;
  my $tag = $contig->{tag};
  die("Missing tag entry in contig") unless defined $tag;
  return $tag;
}

sub _getContigFile {
  my $self = shift;
  my $tag = shift;
  #return $self->db_folder . "/" . $tag . ".gz";
  return $self->db_folder . "/" . $tag;
}

no Moose;
__PACKAGE__->meta->make_immutable;