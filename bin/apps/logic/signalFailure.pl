#!/usr/bin/perl -w

my $label = $ARGV[0];
print "Please supply label argument\n" and exit if not defined $label;
my $labels = "subjob1,subjob2";
print qq{	#### SEND JOB COMPLETION SIGNAL
------------------------------------------------------------
---[completed $label: incomplete $labels]---
------------------------------------------------------------
};