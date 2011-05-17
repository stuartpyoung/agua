#!/usr/bin/perl -w

# runner.pl runs the runned.pl

my $runned = "runned.pl";

print "Executing $runned...\n";
my $stdoutfile = "stdout.txt";
open(STDOUT, ">$stdoutfile") or die "Can't open stdout file: $stdoutfile\n";
exec ("perl $runned");
#print "completed $0\n";