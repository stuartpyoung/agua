#!/usr/bin/perl -w

my $input = <STDIN>;
print "input not supplied. Exiting\n" if not defined $input;

use JSON;
my $jsonParser = JSON->new->allow_nonref;
$input =~ s/([^"])true([^"])/$1"true"$2/g;
$input =~ s/([^"])false([^"])/$1"false"$2/g;

my $object = $jsonParser->decode($input);
print to_json($object, {pretty => 1, indent => 2})
