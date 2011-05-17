#!/usr/bin/perl -w

#### DEBUG
$DEBUG = 1;

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     NOVOALIGN

    PURPOSE

        WRAPPER SCRIPT FOR RUNNING NOVOALIGN:

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

    ./NOVOALIGN.pl <--inputfiles String> [--matesfiles String] <--distance String> <--params String> <--outputdir String> <--referencedir String> <--label Integer>
    [--keep] [--queue String] [--maxjobs Integer] [--cpus Integer ][--help]

    --inputfiles          :   Comma-separated '*sequence.txt' file names (e.g., 's_1_1_sequence.txt,s_2_1_sequence.txt')
    --matefiles           :   Comma-separated '*sequence.txt' mate file names
    --params              :   Additional parameters to be passed to novoalign
    --outputdir           :   Create this directory and write output files to it
    --referencedir           :   Location of indexed referencedir file
                              (NB: Use file stub, e.g., /some/path/hg19)
    --species             :   Name of the referencedir species (e.g., 'human', 'mouse')
    --label               :   Name to used to submit maxjobs to cluster
    --keep                :   Keep intermediate files
    --queue               :   Cluster queue options
    --maxjobs                :   Max. number of concurrent cluster maxjobs (DEFAULT: 30)
    --cpus                :   Max. number of cpus per job (DEFAULT: 4)
    --help                :   print help info

    EXAMPLES

/nethome/bioinfo/apps/agua/0.5/bin/apps/aligners/NOVOALIGN.pl \
--inputfiles /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/2000bp/100M/simpleheader/yoruba1-1.reads_1.sequence.txt \
--matefiles /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/2000bp/100M/simpleheader/yoruba1-1.reads_2.sequence.txt \
--referencedir /nethome/bioinfo/data/sequence/chromosomes/human/hg19/novoalign/chr22 \
--outputdir /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/2000bp/autochr22/novoalign \
--distance 2000 \
--deviation 300 \
--species human \
--label novoalign-autochr22 \
--maxjobs 1000 \
--queue large \
--cluster LSF \
--walltime 24 \
--cpus 1 \
--reads 1000000 \
--stdout /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/2000bp/autochr22/novoalign/autochr22-1.out 


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
use NOVOALIGN;
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
my $novoalign = $conf->getKeyValue("applications", 'NOVOALIGN');
my $samtools = $conf->getKeyValue("applications", 'SAMTOOLS');
my $cluster = $conf->getKeyValue("agua", 'CLUSTERTYPE');
my $qstat = $conf->getKeyValue("cluster", 'QSTAT');
my $qsub = $conf->getKeyValue("cluster", 'QSUB'); #### /usr/local/bin/qsub

#### SET msub TO qsub FOR ARRAY JOB SUBMISSION
$qsub = "/usr/local/bin/qsub";

#### GET OPTIONS
my $min;
my $max;
my $paired;
my $replicates;
my $outputdir;
my $referencedir;
my $stdout;
my $distance;
my $deviation;
my $params;
my $label;

# CLUSTER
my $maxjobs = 30;
my $sleep = 5;
my $verbose;
my $tempdir;
my $cpus = 1;
my $queue = "-q gsmall";
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
    'deviation=i' 	=> \$deviation,
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
print "NOVOALIGN.pl    Output directory is a file: $outputdir\n" and exit if -f $outputdir;
File::Path::mkpath($outputdir) if not -d $outputdir;
print "NOVOALIGN.pl    Can't create output directory: $outputdir\n" if not -d $outputdir;

#### PRINT TO STDOUT IF DEFINED stdout
if ( defined $stdout )
{
	print "NOVOALIGN.pl    Printing STDOUT to stdoutfile:\n\n$stdout\n\n";

	my ($stdout_path) = $stdout =~ /^(.+)\/[^\/]+$/;
	File::Path::mkpath($stdout_path) if not -d $stdout_path;
	open(STDOUT, ">$stdout") or die "Can't open STDOUT file: $stdout\n" if defined $stdout;
}

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
print "outputdir not defined (option --help for usage)\n" and exit  if not defined $outputdir;
print "referencedir not defined (option --help for usage)\n" and exit  if not defined $referencedir;
print "label not defined (option --help for usage)\n" and exit  if not defined $label;
print "deviation not defined for paired end insert (option --help for usage)\n" and exit if $paired and not defined $deviation;
print "distance not defined for paired end insert (option --help for usage)\n" and exit if $paired and not defined $distance;

#### DEBUG
print "NOVOALIGN.pl    outputdir: $outputdir\n";
print "NOVOALIGN.pl    referencedir: $referencedir\n";
print "NOVOALIGN.pl    cluster: $cluster\n";
print "NOVOALIGN.pl    params: $params\n" if defined $params;

#### INSTANTIATE NOVOALIGN
my $runNovoalign = NOVOALIGN->new(
	{
		#### INPUTS (FROM USER)
		min 		=> $min,
		max 		=> $max,
		paired   	=> $paired,
		replicates 	=> $replicates,
		distance   	=> $distance,
		deviation   => $deviation,
		params      => $params,
		outputdir   => $outputdir,
		referencedir=> $referencedir,
		label       => $label,

		#### EXECUTABLES (PRESETS FROM CONF)
		novoalign   => $novoalign,

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
($completed, $label, $sublabels, $missingfiles, $dubiousfiles) = $runNovoalign->check();	

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
print "NOVOALIGN.pl    Run time: $runtime\n";
print "NOVOALIGN.pl    Completed $0\n";
print "NOVOALIGN.pl    ";
print Timer::datetime(), "\n";
print "NOVOALIGN.pl    ****************************************\n\n\n";
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


