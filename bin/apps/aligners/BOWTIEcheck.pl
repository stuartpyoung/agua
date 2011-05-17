#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     BOWTIE

    PURPOSE

        WRAPPER SCRIPT FOR RUNNING BOWTIE:

            1. RUN ONE LANE AT A TIME (MATE PAIRED)

            2. OUTPUT TO A USER-SPECIFIED DIRECTORY

	VERSION		0.02

	HISTORY

		LATER: ADD SPLIT FILES, CONVERSION OF SAM TO BAM
			MERGING OF BAM, AND RECONVERSION TO SAM

		0.02 ADDED CHUNK-BY-CHROMOSOME IF referencedir SPECIFIED
			AND RUN CUFFLINKS COMMAND

		0.01 BASIC VERSION

    INPUT

        1. LIST OF 'INPUTFILES' (s_*_1_sequence.txt)

        2. LIST OF 'MATEFILES' (s_*_2_sequence.txt)

        3. OUTPUT DIRECTORY CONTAINING THE SEQUENCES

    OUTPUT

        1. 

    NOTES

        1. IN THE CASE OF MULTIPLE INPUT FILES, THE READS IN ALL

            FILES MUST BE THE SAME LENGTH BECAUSE A SINGLE s_*_pair.xml

            FILE WILL BE USED FOR THE READS MERGED INTO ONE FILE

    USAGE

    ./BOWTIEcheck.pl <--inputfiles String> [--matesfiles String] <--distance String> <--params String> <--outputdir String> <--referencedir String> <--label Integer> \
    [--queue String] [--maxjobs Integer] [--cpus Integer ] [--cluster String] [--walltime Integer] [--sleep Integer] [--tempdir String] [--cleanup] [--keep] [--verbose] [--help]

    --inputfiles    :   Comma-separated filenames (e.g., 's_1_1_sequence.txt,s_2_1_sequence.txt')
    --paired     :   Comma-separated '*sequence.txt' mate file names
    --params        :   Additional parameters to be passed to bowtie
    --outputdir     :   Create this directory and write output files to it
    --referencedir     :   Location of indexed referencedir file stub (e.g., /some/path/hg19)
    --species       :   Name of the referencedir species (e.g., 'human', 'mouse')
    --label         :   Submit jobs to cluster using this label
    --keep          :   Keep intermediate files
    --maxjobs       :   Max. number of concurrent cluster maxjobs (DEFAULT: 30)
    --cpus          :   Max. number of cpus per job (DEFAULT: 4)
    --cluster       :   Type of cluster job scheduler (e.g., 'PBS', 'SGE', 'LSF')
    --queue         :   Cluster queue options
	--walltime		:	Terminate job if run time exceeds this number of hours
    --cleanup       :   Delete *stdout.txt and *sh files once runs finished
    --sleep         :   Pause for this number of seconds waiting for *stdout.txt files
	--verbose       :	Print debug information
	--tempdir       :	Use this temporary output directory on execution host
    --help          :   print help info

    EXAMPLES

screen -S analysis1
mkdir -p /scratch/syoung/base/pipeline/bixby/run1/ln/bowtie-multi/analysis1/

cd /scratch/syoung/base/pipeline/bixby/run1/ln/bowtie-multi/analysis1/

perl /nethome/bioinfo/apps/agua/0.4/bin/apps/BOWTIEcheck.pl \
--inputfiles /nethome/bioinfo/data/sequence/labs/bixby/run1/s_1_sequence.txt,/nethome/bioinfo/data/sequence/labs/bixby/run1/s_2_sequence.txt \
--params "--solexa1.3-quals" \
--distance 200 \
--queue "-q ngs" \
--gtf /nethome/bioinfo/data/sequence/chromosomes/mouse/mm9/cufflinks/mm9-knownGene.gtf \
--referencedir /nethome/bioinfo/data/sequence/chromosomes/mouse/mm9/bowtie \
--outputdir /scratch/syoung/base/pipeline/bixby/run1/ln/bowtie-multi/analysis1 \
--label bixln \
--keep \
--cluster LSF \
--cpus 4 \
--maxjobs 1000 \
--species mouse \
--clean




TEST WITH SHORT FILE:

cd /nethome/syoung/base/pipeline/bixby/run1/ln/bowtie-multi/

perl /nethome/bioinfo/apps/agua/0.5/bin/apps/BOWTIEcheck.pl \
--inputfiles /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/samples/100M/1/min3length26-1.reads_1.sequence.txt \
--paired /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/samples/100M/1/min3length26-1.reads_2.sequence.txt \
--distance 200 \
--queue large \
--referencedir /nethome/bioinfo/data/sequence/chromosomes/human/hg19/bowtie/chr22 \
--outputdir /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/samples/100M/1/bowtie \
--keep \
--cluster LSF \
--cpus 4 \
--maxjobs 1000 \
--species human \
--label sample1-bowtie \
--reads 1000000 \
&> /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/samples/100M/1/bowtie/chr22/sample1-bowtie.out 



=cut

use strict;

#### EXTERNAL MODULES
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use Getopt::Long;
use FindBin qw($Bin);
use File::Path;
use File::Copy;

#### USE LIBRARY
use lib "$Bin/../../../lib";	
use lib "$Bin/../../../lib/external";	

#### INTERNAL MODULES
use BOWTIE;
use FileTools;
use Timer;
use Util;
use Conf::Agua;

##### STORE ARGUMENTS FOR PRINTING TO USAGE FILE LATER
my @arguments = @ARGV;
unshift @arguments, $0;

#### FLUSH BUFFER
$| =1;

#### GET CONF 
my $conf = Conf::Agua->new(inputfile=>"$Bin/../../../conf/default.conf");
my $bowtie = $conf->getKeyValue("applications", 'BOWTIE');
my $cluster = $conf->getKeyValue("agua", 'CLUSTERTYPE');
my $qstat = $conf->getKeyValue("cluster", 'QSTAT');
my $qsub = $conf->getKeyValue("cluster", 'QSUB'); #### /usr/local/bin/qsub

#### SET msub TO qsub FOR ARRAY JOB SUBMISSION
$qsub = "/usr/local/bin/qsub";
print "BOWTIEcheck.pl    Using bowtie version: $bowtie\n";

#### GET OPTIONS
my $min;
my $max;
my $paired;
my $replicates;
my $outputdir;
my $referencedir;
my $distance;
my $stdout;
my $params;
my $label;

# CLUSTER
my $maxjobs = 30;
my $sleep = 5;
my $verbose;
my $tempdir;
my $cpus = 1;
my $queue = "small";
my $walltime = 24; #### WALL TIME IN HOURS (INTEGER)
my $dot = 1;
my $cleanup;

my $help;
if ( not GetOptions (

	#### GENERAL
    'min=s' 		=> \$min,
    'max=s' 		=> \$max,
    'paired'   		=> \$paired,
    'replicates=s' 	=> \$replicates,
    'outputdir=s'   => \$outputdir,
    'referencedir=s'=> \$referencedir,
    'stdout=s' 		=> \$stdout,
    'distance=i' 	=> \$distance,
    'params=s'      => \$params,
    'label=s'       => \$label,

	#### CLUSTER
    'maxjobs=i'     => \$maxjobs,
    'cpus=i'        => \$cpus,
    'cluster=s' 	=> \$cluster,
    'queue=s'       => \$queue,
    'walltime=i'    => \$walltime,
    'sleep=s' 		=> \$sleep,
    'cleanup'       => \$cleanup,
    'verbose' 		=> \$verbose,
    'tempdir=s' 	=> \$tempdir,
    'help'          => \$help
) )
{ print "Use option --help for usage instructions.\n";  exit;    };


#### MAKE OUTPUT DIR IF NOT EXISTS
print "BOWTIEcheck.pl    Output directory is a file: $outputdir\n" and exit if -f $outputdir;
File::Path::mkpath($outputdir) if not -d $outputdir;
print "BOWTIEcheck.pl    Can't create output directory: $outputdir\n" if not -d $outputdir;

#### PRINT TO STDOUT IF DEFINED stdout
if ( defined $stdout )
{
	print "BOWTIEcheck.pl    Printing STDOUT to stdoutfile:\n\n$stdout\n\n";

	my ($stdout_path) = $stdout =~ /^(.+)\/[^\/]+$/;
	File::Path::mkpath($stdout_path) if not -d $stdout_path;
	open(STDOUT, ">$stdout") or die "Can't open STDOUT file: $stdout\n" if defined $stdout;
}

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "min not defined (Use --help for usage)\n" if not defined $min;
die "max not defined (Use --help for usage)\n" if not defined $max;
die "distance not defined (option --help for usage)\n" if not defined $distance;
die "outputdir not defined (Use --help for usage)\n" if not defined $outputdir;
die "referencedir not defined (Use --help for usage)\n" if not defined $referencedir;
die "label not defined (Use --help for usage)\n" if not defined $label;
die "replicates not defined (Use --help for usage)\n" if not defined $replicates;

#### DEBUG
print "BOWTIEcheck.pl    paired: $paired\n";
print "BOWTIEcheck.pl    outputdir: $outputdir\n";
print "BOWTIEcheck.pl    referencedir: $referencedir\n";
print "BOWTIEcheck.pl    cluster: $cluster\n";
print "BOWTIEcheck.pl    label: $label\n";
print "BOWTIEcheck.pl    params: $params\n" if defined $params;

#### INSTANTIATE BOWTIE
my $runBowtie = BOWTIE->new(
	{
		#### INPUTS (FROM USER)
		min 		=> $min,
		max 		=> $max,
		paired   	=> $paired,
		replicates 	=> $replicates,
		outputdir   => $outputdir,
		referencedir=> $referencedir,
		distance    => $distance,
		params      => $params,
		label       => $label,

		#### EXECUTABLES (FROM CONF)
		bowtie      => $bowtie,

		#### CLUSTER (PRESETS AND USER)
		cluster 	=> $cluster,
		queue 		=> $queue,
		walltime 	=> $walltime,
		qstat 		=> $qstat,
		maxjobs     => $maxjobs,
		cpus        => $cpus,
		sleep       => $sleep,
		cleanup     => $cleanup,
		verbose 	=> $verbose,
		tempdir 	=> $tempdir,
		dot         => $dot,
		qsub        => $qsub,

		command 	=>	\@arguments
	}
);

my ($completed, $sublabels, $missingfiles, $dubiousfiles);
($completed, $label, $sublabels, $missingfiles, $dubiousfiles) = $runBowtie->check();	

#### SEND JOB COMPLETION SIGNAL
print "\n------------------------------------------------------------\n";
print "---[completed $label: $completed $sublabels]---";
if ( scalar(@$missingfiles) > 0 )
{
	print "\n";
	print scalar(@$missingfiles);
	print " missing file:\n" if scalar(@$missingfiles) == 1;
	print " missing files:\n" if scalar(@$missingfiles) != 1;
	print @$missingfiles if scalar(@$missingfiles) > 0;	
}
if ( scalar(@$dubiousfiles) > 0 )
{
	print "\n";
	print scalar(@$dubiousfiles);
	print " dubious file:\n" if scalar(@$dubiousfiles) == 1;
	print " dubious files:\n" if scalar(@$dubiousfiles) != 1;
	print "lowerbound\tupperbound\tsizemultiple\taveragesize\tfilesize\tlocation\n";
	print @$dubiousfiles if scalar(@$dubiousfiles) > 0;	
}
print "\n------------------------------------------------------------\n";

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "BOWTIEcheck.pl    Run time: $runtime\n";
print "BOWTIEcheck.pl    Completed $0 @arguments\n";
print "BOWTIEcheck.pl    ";
print Timer::datetime(), "\n";
print "BOWTIEcheck.pl    ****************************************\n\n\n";
exit;


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#									SUBROUTINES
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


sub usage
{
	print GREEN;
	print `perldoc $0`;
	print RESET;
	exit;
}


