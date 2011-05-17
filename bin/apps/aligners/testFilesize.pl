#!/usr/bin/perl -w

use strict;

$| = 1;

my $logfile = "/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland/1/chr22/check/eland-1-check/check.log";


#### AVOID MAX FILE SIZE OF 49152 ON LINUX
#use POSIX;
my $fh;
#	sysopen $fh, $logfile, O_RDONLY;
#open(LOG, "<$logfile") or die "Can't open logfile: $logfile\n";

open($fh, "<$logfile") or die "Can't open logfile: $logfile\n";
my $filesize = -s $logfile;
print "filesize: $filesize\n";
my $contents;
read($fh, $contents, $filesize);
print "length contents: ", length($contents), "\n";

#
#my $counter = 0;
#foreach my $line ( @lines )
#{
#    $counter++;
#    my ($filesize, $file) = $line =~ /^\S+\t(\S+)\t[^\t]+\t[^\t]+\t([^\t]+)/;
#                and exit if not defined $file;
#                and exit if not defined $filesize;
#}
