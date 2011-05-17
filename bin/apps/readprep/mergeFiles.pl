#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();

=head2

	APPLICATION     mergeFiles

    PURPOSE

        1. SAMPLE mode FROM A LIBRARY 

    INPUT

		1. FASTA OR FASTA.GZ FILE

    OUTPUT

		OPTION: total

			1. TOTAL mode IN ALL READ FILES IN INPUT DIRECTORY

		OPTION: sample

			1. size NUMBER OF FASTA FILES, RECORDS RANDOMLY SELECTED FROM

				FILES IN THE INPUT DIRECTORY

    USAGE

    ./mergeFiles.pl <--inputdir String> <--mode String> [-h]

    --inputdir		:   /full/path/to/inputdir
    --outputdir     :   /full/path/to/outputdir
    --mode          :   Functional mode (total, samples)
    --indices		:   Numbers of merge files to crate
	--cluster		:   Type of cluster, e.g., LSF, PBS
	--maxjobs		:   Maximum number of concurrent maxjobs (DEFAULT = 30)
	--cleanup		:   Remove all script files after completion
	--cleanup		:   Period to sleep between polls to determine run completion
	--dot			:   Print counter every 'dot' number of records
    --help			:   print help info


    EXAMPLES

/nethome/bioinfo/apps/agua/0.5/bin/apps/readprep/mergeFiles.pl \
--inputdir /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/SRX000600/fixedwidth,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/SRX000601/fixedwidth,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/SRX000602/fixedwidth,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/SRX000603/fixedwidth,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/SRX001539/fixedwidth,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/SRX001540/fixedwidth \
--label min3length26 \
--outputdir  /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/samples/100M \
--mode fastq \
--indices 1-29 \
--paired \
--cluster LSF \
--maxjobs 600 \
--queue small \
&> /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/samples/100M/mergeFiles.out


=cut

use strict;

#### FLUSH BUFFER
$| = 1;

#### EXTERNAL MODULES
use Data::Dumper;
use File::Path;
use FindBin qw($Bin);
use Getopt::Long;
use Term::ANSIColor qw(:constants);

#### USE LIBRARY
use lib "$Bin/../../../lib";	
use lib "$Bin/../../../lib/external";	

#### INTERNAL MODULES
use Monitor;
use ReadPrep;
use Sampler;
use Timer;
use Util;
use Conf::Agua;

#### GET modeS
use Getopt::Long;
my $inputdir;
my $label;
my $tempdir;
my $outputdir;
my $indices;
my $max_index;
my $paired;
my $cluster;
my $maxjobs = 30;
my $sleep = 5;
my $queue;
my $cleanup;
my $help;
GetOptions (
    'inputdir=s'	=> \$inputdir,
    'label=s'		=> \$label,
    'outputdir=s'	=> \$outputdir,
    'indices=s'		=> \$indices,
    'max_index=i'	=> \$max_index,
    'paired'		=> \$paired,
    'cluster=s' 	=> \$cluster,
    'maxjobs=i' 	=> \$maxjobs,
    'queue=s' 		=> \$queue,
    'cleanup' 		=> \$cleanup,
    'sleep=i' 		=> \$sleep,
    'help'			=> \$help             
) or die "No arguments provided. Use --help for usage\n";

usage() if defined $help;

#### CHECK INPUTS
die "inputdir not defined (use --help for usage)\n" if not defined $inputdir;
die "label not defined (use --help for usage)\n" if not defined $label;
die "outputdir not defined (use --help for usage)\n" if not defined $outputdir;
die "cluster not defined (use --help for usage)\n" if not defined $cluster;
die "maxjobs not defined (use --help for usage)\n" if not defined $maxjobs;
#die "queue not defined (use --help for usage)\n" if not defined $queue;
#die "cleanup not defined (use --help for usage)\n" if not defined $cleanup;
#die "sleep not defined (use --help for usage)\n" if not defined $sleep;

#### CREATE OUTPUT DIR IF NOT PRESENT
File::Path::mkpath($outputdir) if not -d $outputdir;

#### GET FILES TO BE MERGED AND TARGET FILES
my $mergefiles = Sampler::mergeFiles(
	{
		'inputdir'	=>	$inputdir,
		'label'		=>	$label,
		'outputdir'	=>	$outputdir,
		'max_index'	=>	$max_index,
		'indices'	=>	$indices,
		'tempdir'	=>	$tempdir,
		'paired'	=>	$paired
	}
);

#### INITIALISE CLUSTER OBJECT
my $clusterObject = ReadPrep->new(
	{
		cluster 	=> $cluster,
		queue 		=> $queue,
		cleanup 	=> $cleanup,
		maxjobs     => $maxjobs,
		sleep       => $sleep
	}
);

my $jobs = [];

my @outputfiles = keys ( %$mergefiles );
@outputfiles = sort by_last_number @outputfiles;

my $counter = 0;
foreach my $outputfile ( @outputfiles  )
{
	print "mergeFiles.pl    outputfile: $outputfile\n";	

	my $inputfiles = $mergefiles->{$outputfile};

	#### SET LABEL
	my ($filename) = $outputfile =~ /([^\/]+)$/;
	$filename =~ s/\.gz$//;
	$filename =~ s/\.zip$//;
	$filename =~ s/\.fastq$//;

	my $label = "mergeFiles-$counter-$filename";

	#### GET CUMULATIVE MERGE COMMANDS
	my $commands = [];
	push @$commands, "rm -fr $outputfile";
	foreach my $inputfile ( @$inputfiles )
	{
		push @$commands, "cat $inputfile >> $outputfile";
	}
#exit;

	my $job = $clusterObject->setJob($commands, $label, $outputdir);
	push @$jobs, $job;

	$counter++;
}

#### RUN MERGE JOBS
print "mergeFiles.pl    Running ", scalar(@$jobs), " jobs\n";
$clusterObject->runJobs($jobs, $label );

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "\nRun time: $runtime\n";
print "Completed $0\n";
print Timer::datetime(), "\n";
print "****************************************\n\n\n";
exit;

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#									SUBROUTINES
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


sub by_last_number
{
	my ($aa) = $a =~ /(\d+)\.reads/;
	my ($bb) = $b =~ /(\d+)\.reads/;

	$aa <=> $bb;
}

sub usage
{
	print GREEN;
	print `perldoc $0`;
	print RESET;
	exit;
}


