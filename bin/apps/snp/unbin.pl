#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     cumulativeUCSCBin

    PURPOSE

		1. UNBIN FILES THAT HAVE BEEN BINNED ACCORDING TO THE HIT RANGE

			OF READ ALIGNMENTS IN *bam FILES

		2. CREATE ONE CONCATENATED FILE PER CHROMOSOME

	VERSION		0.01

	HISTORY

		0.01 BASIC FILE CONCATENATION

    INPUTS

        1. COMMA SEPARATED LIST OF INPUT DIRECTORIES

        2. SUFFIX FOR CONCATENATED FILE

		3. OUTPUT DIRECTORY FOR CONCATENATED FILES

		4. BIN LEVEL ORIGINALLY USED TO BIN FILES

    OUTPUTS

        1. ONE CONCATENATED FILE FOR EACH CHROMOSOME PRINTED TO OUTPUT DIRECTORY

    USAGE

    ./unbin.pl \
	<--inputdirs String> <--outputdir String> <--binlevel String> <--suffix Integer> \
    [--cluster String] [--maxjobs Integer] [--queue Integer] [--cpus Integer ] \
	[--help]

    --inputdirs	:   Comma-separated list of directories containing chr* subdirs
    --outputdir :   Directory containing the binned *snp files
    --outputdir :   Directory containing the binned *snp files
    --binlevel  :   UCSC bin level (1-4 in initial implementation, chr < 536 Mbp)
    --queue     :   Cluster queue options
    --maxjobs   :   Max. number of concurrent cluster maxjobs (DEFAULT: 30)
    --cpus      :   Max. number of cpus per job (DEFAULT: 4)
    --help      :   print help info

    EXAMPLES



perl /nethome/bioinfo/apps/agua/0.5/bin/apps/snp/unbin.pl \
--filename hit.bam \
--suffix snp \
--outputdir /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland2/cumulative3 \
--inputdirs /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/autochr22/eland/1,\
/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/autochr22/eland/2,\
/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland2/3,\
/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland2/4,\
/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland2/5,\
/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland2/6,\
/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland2/7,\
/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland2/8,\
/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland2/9,\
/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland2/10,\
/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland2/11,\
/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland2/12,\
/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland2/13,\
/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland2/14,\
/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland2/15,\
/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland2/16,\
/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland2/17,\
/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland2/18,\
/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland2/19,\
/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland2/20,\
/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland2/21,\
/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland2/22,\
/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland2/23,\
/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland2/24,\
/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland2/25,\
/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland2/26,\
/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland2/27,\
/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland2/28,\
/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland2/29,\
/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland2/30,\
/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland2/31,\
/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland2/32,\
/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland2/33 \
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
use UCSCBin;
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
my $suffix;
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
	'suffix=s' 		=> \$suffix,
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
	print "unbin.pl    Printing STDOUT to stdoutfile:\n\n$stdout\n\n";

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
die "binlevel not defined (Use --help for usage)\n" if not defined $binlevel;
die "suffix not defined (Use --help for usage)\n" if not defined $suffix;
die "filename not defined (Use --help for usage)\n" if not defined $filename;

#### MAKE OUTPUT DIR IF NOT EXISTS
print "unbin.pl    outputdir is a file: $outputdir\n" and exit if -f $outputdir;
File::Path::mkpath($outputdir) if not -d $outputdir;
print "unbin.pl    Can't create output directory: $outputdir\n" and exit if not -d $outputdir;


#### CHECK INPUT DIRS
my @indirs = split ",", $inputdirs;
foreach my $indir ( @indirs )
{   
    print "unbin.pl    Can't find inputdir: $indir\n" and exit if not -d $indir;
}

#### DEBUG

#### RETRIEVE COMMAND
my $command = "$0 @arguments";

#### INSTANTIATE BINNER
my $binner = UCSCBin->new(
	{
		#### INPUTS (FROM USER)
		inputdirs  	=> \@indirs,
		outputdir   => $outputdir,
		filename   	=> $filename,
		binlevel   	=> $binlevel,
		suffix   	=> $suffix,
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

#### DO UNBIN
$binner->unbin();

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "unbin.pl    Run time: $runtime\n";
print "unbin.pl    Completed $0\n";
print "unbin.pl    ";
print Timer::datetime(), "\n";
print "unbin.pl    ****************************************\n\n\n";
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


