#!/usr/bin/perl -w

my $sleep = $ARGV[0];
$sleep = 5 if not defined $sleep;

while ( 1 )
{
    qsum();
    sleep($sleep);
}

sub qsum {
    print "Doing qsum\n";
    my $command = "qstat -f";
    my $output = `$command`;


    my ($nodes, $jobs) = $output =~ /^(.+)\n#{10,}\n.+?PENDING JOBS.+?\n#{10,}\n(.+)$/ms;
    print "nodes: $nodes\n";
    print "jobs: $jobs\n";

    my @lines = split "\n", $jobs;
    my $jobcounts = {};
    foreach my $line ( @lines )
    {
        my ($jobname) = $line =~ /^\s*\S+\s+\S+\s+(\S+)/;
        if ( exists $jobcounts->{$jobname} )
        {
            $jobcounts->{$jobname}++;
        }
        else
        {
            $jobcounts->{$jobname} = 1;
        }
    }

    print "=================================================================================\n";
    print "$nodes\n";
    foreach my $key ( sort keys %$jobcounts )
    {
        print "$key: $jobcounts->{$key}\n";
    }
}