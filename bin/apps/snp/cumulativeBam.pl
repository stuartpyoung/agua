#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     cumulativeBam

    PURPOSE

		WE ARE PROVIDED A LIST OF DIRECTORIES CONTAINING chr* CHROMOSOME

		SUBDIRECTORIES, EACH CONTAINING A CHROMOSOME *.sam FILE. FOR EACH

		CHROMOSOME, STARTING WITH THE FIRST DIRECTORY AND GOING IN ORDER,

		FOR EACH DIRECTORY:

			1. CONCAT NEW *bam FILE WITH EXISTING merged-N.bam FILE

			2. REPEAT UNTIL ALL *bam FILES MERGED

	VERSION		0.01

	HISTORY

		0.01 BASIC VERSION MERGING BAM FILES

    INPUTS

        1. 'N' NUMBER OF COMMA SEPARATED LIST OF 'INPUTFILES' (E.G., '/Full/path/to/directory')

        2. SPECIES NAME TO DETERMINE CORRECT SAMTOOLS INDEXES

    OUTPUTS

        1. 'N' NUMBER OF SNP FILES IN INCREASING ORDER OF READ HIT QUANTITY

		2. ONE FOR EACH CHROMOSOME

    USAGE

    ./cumulativeBam.pl \
	<--inputdirs String>  <--outputdir String> \
	<--reference String> <--label Integer> \
    [--keep] [--queue String] [--maxjobs Integer] [--cpus Integer ][--help]

    --inputdirs	:   Comma-separated list of directories containing
	                         chr* subdirs
    --outputdir :   Create this directory and write output files to it
	--zipped    :   Input file is compressed (gzip, zip or bz2)
    --species   :   Name of the reference species (e.g., 'human', 'mouse')
    --label     :   Name to used to submit maxjobs to cluster
    --keep      :   Keep intermediate files
    --queue     :   Cluster queue options
    --maxjobs   :   Max. number of concurrent cluster maxjobs (DEFAULT: 30)
    --cpus      :   Max. number of cpus per job (DEFAULT: 4)
    --help      :   print help info

    EXAMPLES


perl /nethome/bioinfo/apps/agua/0.5/bin/apps/snp/cumulativeBam.pl \
--inputdirs /ngs/bioinfo/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/autochr22/eland/1,/ngs/bioinfo/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/autochr22/eland/2 \
--zipped gzip \
--filename hit.sam.gz \
--outputdir /ngs/bioinfo/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/autochr22/eland/cumulative-1 \
--zipped \
--species human \
--label cumulative-1 \
--walltime 48 \
--cluster LSF \
--queue small \
--cpus 1 \
--maxjobs 1000 \
--clean

\
--stdout /ngs/bioinfo/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/autochr22/eland/cumulative-1/cumulativeBam-stdout.txt


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
my $bindir;
my $zipped;
my $clean;
my $stdout;


# SPECIFIC
my $start;
my $label;
my $params;
my $species;
my $samtools_index;
my $keep;

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
    'bindir=s' 		=> \$bindir,
    'clean'        	=> \$clean,
    'stdout=s' 		=> \$stdout,

	#### SPECIFIC
    'start=s'       => \$start,
    'label=s'       => \$label,
    'params=s'      => \$params,
    'species=s'     => \$species,
    'samtoolsindex=s' => \$samtools_index,
    'keep'        	=> \$keep,

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
	print "cumulativeBam.pl    Printing STDOUT to stdoutfile:\n\n$stdout\n\n";

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
die "species not defined (Use --help for usage)\n" if not defined $species;
die "label not defined (Use --help for usage)\n" if not defined $label;
die "filename not defined (Use --help for usage)\n" if not defined $filename;
$zipped = "gzip" if $filename =~ /\.gz$/;
$zipped = "zip" if $filename =~ /\.zip$/;
$zipped = "bz2" if $filename =~ /\.bz2$/;

#### MAKE OUTPUT DIR IF NOT EXISTS
print "cumulativeBam.pl    outputdir is a file: $outputdir\n" and exit if -f $outputdir;
File::Path::mkpath($outputdir) if not -d $outputdir;
print "cumulativeBam.pl    Can't create output directory: $outputdir\n" and exit if not -d $outputdir;


#### CHECK INPUT DIRS
my @indirs = split ",", $inputdirs;
foreach my $indir ( @indirs )
{   
    print "cumulativeBam.pl    Can't find inputdir: $indir\n" and exit if not -d $indir;
}

#### DEBUG
print "cumulativeBam.pl    inputdirs: $inputdirs\n";
print "cumulativeBam.pl    outputdir: $outputdir\n";
print "cumulativeBam.pl    cluster: $cluster\n";
print "cumulativeBam.pl    params: $params\n" if defined $params;

#### IF NOT DEFINED samtools_index, SET IT BASED ON SPECIES
if ( not defined $samtools_index )
{
	my $parameter = "SAMTOOLS" . uc($species);
	$samtools_index = $conf->getKeyValue("data", $parameter);
}
print "cumulativeBam.pl    samtools index: $samtools_index\n";
die "Can't find samtools index directory: $samtools_index\n" if not -d $samtools_index;

#### RETRIEVE COMMAND
my $command = "$0 @arguments";

#### INSTANTIATE cumulativeBam
my $snp = SNP->new(
	{
		#### INPUTS (FROM USER)
		inputdirs  	=> \@indirs,
		outputdir   => $outputdir,
		filename   	=> $filename,
		binlevel   	=> $binlevel,
		bindir   	=> $bindir,
		zipped		=> $zipped,
		clean     	=> $clean,
		label       => $label,

		#### SNP-SPECIFIC INPUTS
		start       => $start,
		params      => $params,
		species     => $species,
		samtools	=> $samtools,
		samtoolsindex	=> $samtools_index,

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
print "cumulativeBam.pl    Doing cumulativeBam()\n";
$snp->cumulativeBam();

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "cumulativeBam.pl    Run time: $runtime\n";
print "cumulativeBam.pl    Completed $0\n";
print "cumulativeBam.pl    ";
print Timer::datetime(), "\n";
print "cumulativeBam.pl    ****************************************\n\n\n";
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


