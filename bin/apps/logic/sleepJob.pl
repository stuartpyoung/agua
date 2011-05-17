#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../lib";
use lib "$Bin/../../../lib/external";

use TestCluster;

print "Running $0\n";
my $test = TestCluster->new(
    {
        maxjobs     =>  10,
        cluster     =>  "LSF",
        queue       =>  "small",
        walltime    =>  24,
        submit      =>  1
    }
);
$test->sleepJob();

print "Doing sleeps\n";
my $sleep = 3;
my $counter = $sleep;
while ( $counter > 0 )
{
    print "$counter\n";
    $counter--;
    sleep(1);
}
print "Completed $0\n";
