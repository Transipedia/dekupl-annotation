package DEkupl::Base;
# ABSTRACT: Base class for dekupl modules

use Moose::Role;

has 'verbose' => (
  is => 'ro',
  isa => 'Bool',
  default => 0,
);

sub verboseLog {
  my $self = shift;
  my $message = shift;
  print STDERR "[".$self->meta->name."] $message\n" if $self->verbose;
}

1;