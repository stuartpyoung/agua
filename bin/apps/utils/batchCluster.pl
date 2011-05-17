#!/usr/bin/perl -w

#### DEBUG
$DEBUG = 1;

#### TIME
my $time = time();

=head2

	APPLICATION     batchCluster

    PURPOSE

        SPLIT LINES INTO A SEPARATE FILE PER UNIQUE VALUE IN THE USER-DEFINED COLUMN

    INPUT

        1. LOCATION OF A TAB-SEPARATED BATCH FILE CONTAINING COLUMN VALUES

        2. COMMA-SEPARATED LIST OF ARGUMENT (COLUMN) NAMES

        3. LOCATION OF APPLICATION TO BE RUN IN BATCH MODE

		4. OUTPUT DIR FOR WRITING SHELL SCRIPTS

		5. 

    OUTPUT

        1. RUN APPLICATION IN BATCH MODE IN SERIES USING THE ARGUMENTS

			PROVIDED IN EACH LINE OF THE BATCH FILE. EACH COLUMN IN THE

			BATCH FILE REPRESENTS AN ARGUMENT VALUE AND THE COLUMN NAMES

			ARE THE ARGUMENT AS IT SHOULD BE PRESENTED TO THE APPLICATION

    USAGE

    ./batch.pl  <--batchfile String> <--application String> <--columns String> [--separator String] [-h]

        --batchfile		:   Location of .tsv (tab-separated values) batch file
        --columns       :   Names of columns in batch file
        --application	:   Location of application to be run in batch mode
		--outputdir		:	Location to print STDOUT files
		--dot			:	Print the entry number every 'dot' batch entries
		--qsub			: 	Full path to qsub/msub
		--queue			:	Name of queue, e.g., "-q gsmall"
        --help			:   Print help info


    EXAMPLES

/nethome/syoung/base/bin/utils/batchCluster.pl \
--batchfile /nethome/syoung/base/pipeline/jbrowse/ucsc/conf/chromosome-gff-biodb-to-json-batchfile.txt \
--outputdir /nethome/syoung/base/pipeline/jbrowse/ucsc/runs/biodb \
--columns "--executable,--inputdir,--outputdir,--jsonfile" \
--application /nethome/syoung/base/apps/agua/0.4/bin/apps/run-biodb-to-json.pl \
&> /nethome/syoung/base/pipeline/jbrowse/ucsc/runs/biodb/batch-biodb.out




EXAMPLE BATCH FILE CONTENTS:

/nethome/syoung/base/pipeline/jbrowse1/ucsc/genome-gtf/affy-exon-probes-gtf	1	/nethome/syoung/base/pipeline/jbrowse1/ucsc/chromosome-gtf/	affy-exon-probes-	gff
/nethome/syoung/base/pipeline/jbrowse1/ucsc/genome-gtf/broad-histone-gtf	1	/nethome/syoung/base/pipeline/jbrowse1/ucsc/chromosome-gtf/	broad-histone-	gff
...


=cut

use strict;

#### USE LIBRARY
use FindBin qw($Bin);
use lib "$Bin/../../../lib";

#### FLUSH BUFFER
$| = 1;

#### INTERNAL MODULES
use Sampler;
use Timer;
use Util;
use Conf::Agua;

#### EXTERNAL MODULES
use Term::ANSIColor qw(:constants);
use Getopt::Long;
use File::Path;

#### SET DEFAULT SEPARATOR
my $DEFAULT_SEPARATOR = "\\s+";

#### SAVE ARGUMENTS
my @arguments = @ARGV;

#### GET OPTIONS
my $batchfile;
my $application;
my $columns;
my $outputdir;
my $help;

#### CLUSTER OPTIONS
my $qsub = "/usr/local/bin/msub";
my $jobs = 30;
my $cpus = 1;
my $sleep = 5;
my $queue = '';
my $dot = 1;

if ( not GetOptions (
    'batchfile=s' => \$batchfile,
    'application=s' => \$application,
    'columns=s' => \$columns,
    'outputdir=s' => \$outputdir,
    'dot=i' => \$dot,
    'qsub=s' => \$qsub,
    'cpus=s' => \$cpus,
    'queue=s' => \$queue,
    'help' => \$help
) )
{ print "Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### SET VALUE TO ZERO IF batchfile NOT DEFINED
#### OTHERWISE, SET VALUE TO THE batchfile (ORDINAL) NUMBER
#### WHERE 1 IS THE LAST batchfile
die "batchfile not defined (Use --help for usage)\n" if not defined $batchfile;
die "application not defined (Use --help for usage)\n" if not defined $application;
die "columns not defined (Use --help for usage)\n" if not defined $columns;
die "outputdir not defined (Use --help for usage)\n" if not defined $outputdir;


#### GET OUTPUT SUBDIR
my $output_subdir = next_subdir($outputdir, "batchrun");
print "output_subdir: $output_subdir\n";

#### CREATE OUTPUT SUBDIR
if ( not -d $output_subdir )
{
	File::Path::mkpath($output_subdir) or die "Can't create output subdir: $output_subdir\n";
}

#### GENERATE ARRAY OF JOBS
my @commands;
open(FILE, $batchfile) or die "Can't open batch file: $batchfile\n";
my $batch_counter = 0;
while ( <FILE> )
{
	next if $_ =~ /^\s*$/ or $_ =~ /^#/; 
	$batch_counter++;

	my @values = split "\t", $_;
	my @arguments = split ",", $columns;
	if ( $#values != $#arguments )
	{
		print "No. of arguments (" , $#arguments + 1, ") is not the same as no. of values (" , $#values + 1, ")\n";
		print "Batch file line: $_\n";
		next;
	}

	my $command = "$application \\\n";
	for (my $i = 0; $i < $#values + 1; $i++ )
	{
		$values[$i] =~ s/\n$// if $ == $#values;
		$command .= qq{$arguments[$i] "$values[$i]" \\\n};
	}

	my $appstub = $application;
	$appstub =~ s/\.[^\.]{1,5}$//;
	$command .= "&> $output_subdir/$appstub.$batch_counter" if defined $stdout;
	$command =~ s/\\\n$// if not defined $stdout;

	push @commands, $command;
}

#### SUBMIT ARRAY OF JOBS TO CLUSTER
my $qsub_counter = 0;
my $pids = [];
my $scriptfiles = [];

foreach my $command ( @commands )
{	
    $qsub_counter++;
    print "$qsub_counter\n" if $qsub_counter % $dot == 0;

	#### CREATE AN ADDITIONAL SUBDIR FOR EACH COMMAND RUN
	my $final_subdir = "$output_subdir/$qsub_counter";
	if ( not -d $final_subdir )
	{
		File::Path::mkpath($final_subdir) or die "Can't create final subdir: $final_subdir\n";
	}

	#### MOVE TO OUTPUT SUBDIR
	chdir("$final_subdir") or die "Can't move to output subdir: $final_subdir\n";

	#### PRINT SHELL SCRIPT	
	my $scriptfile = "$final_subdir/batchscript-$qsub_counter.sh";
	my $usagefile = "$final_subdir/usage-$qsub_counter.txt";
	my $shellscript = qq{#!/bin/sh

## Request N CPUs (instead of nodes). CPUs are always on one machine.
#PBS -l ncpus=$cpus
#PBS -j oe

HOST=`hostname`
echo \$HOST

$command

echo "PBS_ID: " \$PBS_ID

/usr/local/bin/qstat -f \$PBS_ID > $usagefile`;

echo "Printed usage file: $usagefile

};

	open(SHFILE, ">$scriptfile") or die "Can't open script file: $scriptfile\n";
	print SHFILE $shellscript;
	close(SHFILE);
	chmod(0777, $scriptfile) or die "Can't chmod 0777 script file: $scriptfile\n";

	print "scriptfile printed: $scriptfile\n";
	print `cat $scriptfile`;
	#exit;

	### RUN SHELL SCRIPT ON CLUSTER
	my $qsub_command = "$qsub $queue $scriptfile";
	print "\n$qsub_command***\n";

	#### RUN QSUB AND STORE PID
	my $pid = `$qsub_command`;
	($pid) = $pid =~ /^\s+(\d+)\s+/ms;

	push @$pids, $pid;
	print "Currently ", scalar(@$pids), "in list: @$pids\n";
	print `date`;
	while ( scalar(@$pids) >= $jobs )
	{
		sleep($sleep);
		$pids = Sampler::monitor_jobs($pids);   
	}

	#### CLEAN UP
	#`rm -fr $scriptfile`;
	push @$scriptfiles, $scriptfile;
}

#### WAIT TIL ALL JOBS ARE FINISHED
print "Waiting for remaining jobs to finish (", scalar(@$pids), ")\n";
print "sleep: $sleep seconds\n";
my $sleep_counter = 0;
while ( defined $pids and scalar(@$pids) > 0 )
{
	sleep($sleep);
	$pids = Sampler::monitor_jobs($pids);

	#### PRINT PROGRESS 100 DOTS PER LINE
	$sleep_counter++;
	print "." if $sleep_counter %100 != 0;
	print "\n" if $sleep_counter %100 == 0;
}
print "All jobs completed\n";


#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "\nRun time: $runtime\n";
print "Completed $0\n";
print Util::datetime(), "\n";
print "****************************************\n\n\n";
exit;

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#									SUBROUTINES
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


sub next_subdir
{

    my $directory   =   shift;
    my $subdir    =   shift;

	return if not defined $directory;
	return if not defined $subdir;

	$directory =~ s/\/$//;

	### INCREMENT THE COUNTER UNTIL SUBDIR NOT FOUND
	my $i = 0;
	$i++ while ( -e "$directory/$subdir.$i" );


	return "$directory/$subdir.$i";
}


sub usage
{
	print `perldoc $0`;
    exit;
}

