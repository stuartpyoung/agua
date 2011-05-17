#!/usr/bin/perl -w
use strict;

# save a ton of memory by not storing any of your records in memory
# http://www.perlmonks.org/?node_id=53920

die "Usage: $0 col[,col[...]] [file[,...]]\n" unless @ARGV;
my @cols= map { $_-1 } split/,/,shift;
my @lines= <>;
my @sort= map { my $x=join"\0"x5,(split/\s*,\s*/)[@cols];
        $x =~ s/(^|[^\d.])(\d+)/$1.pack("N",$2)/eg; $x } @lines;
print @lines[  sort { $sort[$a] cmp $sort[$b] } 0..$#sort  ];