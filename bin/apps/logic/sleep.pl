#!/usr/bin/perl -w

use strict;

my $sleep = $ARGV[1];
$sleep = 2 if not defined $sleep or not $sleep;
print "Running $0\n";
print "Sleeping $sleep seconds\n";
sleep($sleep);
print "Completed $0\n";
print "\n------------------------------------------------------------\n";
print "---[completed sleep: failed]---";
print "\n------------------------------------------------------------\n";
