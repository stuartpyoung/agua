#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     binBam

    PURPOSE

		SPLIT UP *bam FILES INTO EQUALLY-SIZED CHROMOSOMAL RANGE BINS:

		-	USE UCSC BINNING SYSTEM 
			http://genomewiki.ucsc.edu/index.php/Bin_indexing_system
			Kent & al. The Human Genome Browser at UCSC", doi: 10.1101/gr.229102 . Genome Res. 2002. 12:996-1006

		-   DIVIDE UP *.bam FILE ENTRIES BY CHROMOSOMAL POSITION

		-	WIDTH OF THE BINS IS DEFINED BY THE SPECIFIED BIN LEVEL

		-   DEFAULT: BIN LEVEL 1 (8 x 64MB BINS)

		-   ADD binlevel, binnumber AND totalbins TO BINNED *bam FILE NAME


	VERSION		0.01

	HISTORY

		0.01 BASIC VERSION FOR INITIAL IMPLEMENTATION OF UCSC BINNING
			(Used when chromEnd is less than or equal to 536,870,912 BASES)

    INPUTS

        1. LOCATION OF *bam FILE

        2. BIN LEVEL

    OUTPUTS

        1. SERIES OF *bam FILES hit.binlevel1.bin1-8.bam

		2. ONE FOR EACH CHROMOSOME

    USAGE

    ./binBam.pl \
	<--inputdirs String>  <--outputdir String> <--binlevel String> \
	[--queue String] [--maxjobs Integer] [--cpus Integer ] [--cluster String] [--walltime Integer] [--sleep Integer] [--tempdir String] [--cleanup] [--keep] [--verbose] [--help]

    --inputdirs	:   Comma-separated list of directories containing
	                         chr* subdirs
    --outputdir :   Create this directory and write output files to it
    --binlevel  :   UCSC bin level (1-4 in initial implementation, chr < 536 Mbp)
    --queue     :   Cluster queue options
    --maxjobs   :   Max. number of concurrent cluster maxjobs (DEFAULT: 30)
    --cpus      :   Max. number of cpus per job (DEFAULT: 4)
    --help      :   print help info

    EXAMPLES


/nethome/bioinfo/apps/agua/0.5/bin/apps/snp/binBam.pl \
--inputdirs /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland2/1/chr22/hit.bam \
--binlevel 2 \
--reference chr22 \
--maxjobs 1000 \
--walltime 48 \
--queue small \
--outputdir /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland2/1/chr22/bins \
--cluster LSF 

\
--stdout /ngs/bioinfo/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/autochr22/eland/cumulative-1/binBam-stdout.txt


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
use UCSCBin;
use Timer;
use Util;
use Conf::Agua;

##### STORE ARGUMENTS TO PRINT TO FILE LATER
my @arguments = @ARGV;

#### FLUSH BUFFER
$| =1;

#### GET CONF 
my $conf = Conf::Agua->new(inputfile=>"$Bin/../../../conf/default.conf");
my $samtools = $conf->getKeyValue("applications", 'SAMTOOLS');
my $cluster = $conf->getKeyValue("agua", 'CLUSTERTYPE');
my $qstat = $conf->getKeyValue("cluster", 'QSTAT');
my $qsub = $conf->getKeyValue("cluster", 'QSUB'); #### /usr/local/bin/qsub

#### SET msub TO qsub FOR ARRAY JOB SUBMISSION
$qsub = "/usr/local/bin/qsub";

#### GET OPTIONS
# GENERAL
my $inputdirs;
my $outputdir;
my $reference;
my $binlevel;
my $stdout;

# CLUSTER
my $maxjobs = 30;
my $sleep = 5;
my $verbose;
my $tempdir;
my $cpus = 1;
my $queue;
my $walltime = 24; #### WALL TIME IN HOURS (INTEGER)
my $dot = 1;
my $cleanup;

my $help;
if ( not GetOptions (

	#### GENERAL
    'inputdirs=s'  => \$inputdirs,
    'outputdir=s'  => \$outputdir,
    'reference=s'   => \$reference,
    'binlevel=i'   	=> \$binlevel,
    'stdout=s' 		=> \$stdout,

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

#### PRINT TO STDOUT IF DEFINED stdout
if ( defined $stdout )
{
	print "binBam.pl    Printing STDOUT to stdoutfile:\n\n$stdout\n\n";

	my ($stdout_path) = $stdout =~ /^(.+)\/[^\/]+$/;
	File::Path::mkpath($stdout_path) if not -d $stdout_path;
	open(STDOUT, ">$stdout") or die "Can't open STDOUT file: $stdout\n" if defined $stdout;
}

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "queue not defined (option --help for usage)\n" if not defined $queue;
die "inputdirs not defined (option --help for usage)\n" if not defined $inputdirs;
#die "outputdir not defined (option --help for usage)\n" if not defined $outputdir;

#### CHECK INPUT DIRS
my @infiles = split ",", $inputdirs;
foreach my $infile ( @infiles )
{   
    print "binBam.pl    Can't find inputdir: $infile\n" and exit if not -d $infile;
}

#### SET hitRange EXECUTABLE
my $executable = "$Bin/hitRange.pl";

#### DEBUG

#### RETRIEVE COMMAND
my $command = "$0 @arguments";

#### INSTANTIATE UCSC BINS
my $binner = UCSCBin->new(
	{
		#### INPUTS (FROM USER)
		inputdirs  => \@infiles,
		outputdir  =>  $outputdir,
		reference   => $reference,
		samtools	=> $samtools,
		binlevel	=> $binlevel,
		executable	=> $executable,

		#### CLUSTER (PRESETS AND USER)
		command		=> $command,
		cluster 	=> $cluster,
		queue 		=> $queue,
		walltime 	=> $walltime,
		cpus        => $cpus,
		qstat 		=> $qstat,
		qsub        => $qsub,
		maxjobs     => $maxjobs,
		sleep       => $sleep,
		cleanup     => $cleanup,
		verbose 	=> $verbose,
		tempdir 	=> $tempdir,
		dot         => $dot
	}
);

#### RUN RNA PAIRED TRANSCRIPTOME ANALYSIS
print "binBam.pl    Doing binner->binBam()\n";
$binner->binBam();

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "binBam.pl    Run time: $runtime\n";
print "binBam.pl    Completed $0\n";
print "binBam.pl    ";
print Timer::datetime(), "\n";
print "binBam.pl    ****************************************\n\n\n";
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


