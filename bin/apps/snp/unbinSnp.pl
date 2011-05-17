#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     cumulativeSNP

    PURPOSE

		WE ARE PROVIDED A LIST OF DIRECTORIES CONTAINING chr* CHROMOSOME

		SUBDIRECTORIES, EACH CONTAINING A CHROMOSOME *.sam FILE. FOR EACH

		CHROMOSOME, STARTING WITH THE FIRST DIRECTORY AND GOING IN ORDER,

		FOR EACH DIRECTORY:

			1. CUMULATIVELY CONCAT THE FILES INTO A NEW merged.sam FILE

			2. COPY THE merged.sam FILE AND CONVERT IT INTO A merged.sam.gz FILE

			3. CALL SNPS USING THE merged.sam.gz FILE

			4. RETAIN THE merged.sam FILE FOR USE IN THE NEXT ITERATION

	VERSION		0.01

	HISTORY

		0.01 BASIC VERSION MERGING SAM FILES

    INPUTS

        1. 'N' NUMBER OF COMMA SEPARATED LIST OF 'INPUTFILES' (E.G., '/Full/path/to/hit.sam')

        2. SPECIES NAME TO DETERMINE CORRECT SAMTOOLS INDEXES

    OUTPUTS

        1. 'N' NUMBER OF SNP FILES IN INCREASING ORDER OF READ HIT QUANTITY

		2. ONE FOR EACH CHROMOSOME

    USAGE

    ./unbinSnp.pl \
	<--inputdirs String>  <--outputdir String> \
	<--reference String> <--label Integer> \
    [--keep] [--queue String] [--maxjobs Integer] [--cpus Integer ][--help]

    --inputdirs	:   Comma-separated list of directories containing chr* subdirs
    --outputdir :   Directory containing the binned *snp files
    --binlevel  :   UCSC bin level (1-4 in initial implementation, chr < 536 Mbp)
    --queue     :   Cluster queue options
    --maxjobs   :   Max. number of concurrent cluster maxjobs (DEFAULT: 30)
    --cpus      :   Max. number of cpus per job (DEFAULT: 4)
    --help      :   print help info

    EXAMPLES


perl /nethome/bioinfo/apps/agua/0.5/bin/apps/snp/unbinSnp.pl \
--inputdirs /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/autochr22/eland/1,/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/autochr22/eland/2 \
--filename hit.bam \
--outputdir /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland2/cumulative \
--walltime 48 \
--cluster LSF \
--queue small \
--cpus 1 \
--maxjobs 1000 \
--binlevel 2


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
use SNP;
use FileTools;
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
my $filename;
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
    'inputdirs=s'  	=> \$inputdirs,
    'outputdir=s'   => \$outputdir,
    'filename=s' 	=> \$filename,
    'binlevel=s' 	=> \$binlevel,
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
	print "unbinSnp.pl    Printing STDOUT to stdoutfile:\n\n$stdout\n\n";

	my ($stdout_path) = $stdout =~ /^(.+)\/[^\/]+$/;
	File::Path::mkpath($stdout_path) if not -d $stdout_path;
	open(STDOUT, ">$stdout") or die "Can't open STDOUT file: $stdout\n" if defined $stdout;
}

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "queue not defined (option --help for usage)\n" if not defined $queue;
die "inputdirs not defined (option --help for usage)\n" if not defined $inputdirs;
die "outputdir not defined (Use --help for usage)\n" if not defined $outputdir;
die "filename not defined (Use --help for usage)\n" if not defined $filename;

#### MAKE OUTPUT DIR IF NOT EXISTS
print "unbinSnp.pl    outputdir is a file: $outputdir\n" and exit if -f $outputdir;
File::Path::mkpath($outputdir) if not -d $outputdir;
print "unbinSnp.pl    Can't create output directory: $outputdir\n" and exit if not -d $outputdir;


#### CHECK INPUT DIRS
my @indirs = split ",", $inputdirs;
foreach my $indir ( @indirs )
{   
    print "unbinSnp.pl    Can't find inputdir: $indir\n" and exit if not -d $indir;
}

#### DEBUG
print "unbinSnp.pl    inputdirs: $inputdirs\n";
print "unbinSnp.pl    outputdir: $outputdir\n";
print "unbinSnp.pl    cluster: $cluster\n";

#### RETRIEVE COMMAND
my $command = "$0 @arguments";

#### INSTANTIATE cumulativeSNP
my $snp = SNP->new(
	{
		#### INPUTS (FROM USER)
		inputdirs  	=> \@indirs,
		outputdir   => $outputdir,
		filename   	=> $filename,
		binlevel   	=> $binlevel,
		samtools	=> $samtools,

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
print "unbinSnp.pl    Doing cumulativeSnp()\n";
$snp->unbinSnp();

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "unbinSnp.pl    Run time: $runtime\n";
print "unbinSnp.pl    Completed $0\n";
print "unbinSnp.pl    ";
print Timer::datetime(), "\n";
print "unbinSnp.pl    ****************************************\n\n\n";
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


