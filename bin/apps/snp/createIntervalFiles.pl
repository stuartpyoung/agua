#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     createIntervalFiles

    PURPOSE

		REALIGN *bam FILE READS AROUND INDEL OR MISMATCH REGIONS AND

		OUTPUT TO A 'CLEANED' *bam FILE

	NOTES

		GATK REFERENCE:

		Unlike most mappers, the GATK IndelRealigner walker uses the full alignment context to determine whether an appropriate alternate reference (i.e. indel) exists and updates SAMRecords accordingly.


		THE PROCESS CONSISTS OF TWO STEPS:

			1. CREATE INDEL INTERVALS FILE

			Emit intervals for the Local Indel Realigner to target for cleaning.

			time /usr/local/java/bin/java \
			-jar /nethome/syoung/base/apps/gatk/1.0.4705/GenomeAnalysisTK.jar \
			-I aln.bam \
			-R $REF \
			-T RealignerTargetCreator \
			-U \
			-S SILENT \
			-o realignTargets.intervals

			2. REALIGN READS IN MERGED INTERVALS

			Perform local realignment of reads based on misalignments due to the presence of indels. 
			time /usr/local/java/bin/java \
			-jar /nethome/syoung/base/apps/gatk/1.0.4705/GenomeAnalysisTK.jar \
			-I aln.bam \
			-R $REF \
			-T IndelRealigner \
			-targetIntervals realignTargets.intervals \
			--out realigned.bam \
			-U \
			-S SILENT

	VERSION		0.01

	HISTORY

		0.01 CREATE TARGETS AND GENERATE CLEANED *bam FILE USING GATK

    INPUTS

        1. SORTED *bam FILE

        2. OUTPUT DIRECTORY LOCATION

    OUTPUTS

        1. 'TARGETS' FILE CONTAINING REGIONS WITH INDELS AND MISMATCHES

		2. CLEANED *bam FILE CONTAINING READS REALIGNED IN TARGET REGIONS

    USAGE

    ./createIntervalFiles.pl \
	<--inputdirs String>  <--outputdir String> \
	<--reference String> <--label Integer> \
    [--keep] [--queue String] [--maxjobs Integer] [--cpus Integer ][--help]

    --inputdirs	:   Comma-separated list of directories containing
	                         chr* subdirs
    --outputdir :   Create this directory and write output files to it
    --species   :   Name of the reference species (e.g., 'human', 'mouse')
    --label     :   Name to used to submit maxjobs to cluster
    --keep      :   Keep intermediate files
    --queue     :   Cluster queue options
    --maxjobs   :   Max. number of concurrent cluster maxjobs (DEFAULT: 30)
    --cpus      :   Max. number of cpus per job (DEFAULT: 4)
    --help      :   print help info

    EXAMPLES


perl /nethome/bioinfo/apps/agua/0.5/bin/apps/snp/createIntervalFiles.pl \
--binlevel 2 \
--filename hit.bam \
--inputdirs \
/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/maq/1,\
/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/maq/2,\
--outputdir /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/maq/realign \
--species human \
--label eland-realign \
--walltime 48 \
--cluster LSF \
--queue small \
--cpus 1 \
--maxjobs 1000 \
--clean \
--stdout /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/maq/realign/createIntervalFiles.binlevel2.out


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
my $gatk = $conf->getKeyValue("applications", 'GATK');
my $java = $conf->getKeyValue("applications", 'JAVA');
my $cluster = $conf->getKeyValue("agua", 'CLUSTERTYPE');
my $qstat = $conf->getKeyValue("cluster", 'QSTAT');
my $qsub = $conf->getKeyValue("cluster", 'QSUB'); #### /usr/local/bin/qsub

#### SET msub TO qsub FOR ARRAY JOB SUBMISSION
$qsub = "/usr/local/bin/qsub";

#### GET OPTIONS
# SPECIFIC
my $inputdirs;
my $outputdir;
my $window = 5;
my $mode = "pileup";
my $referencedir;
my $species;
my $filename;
my $binlevel;
my $bindir;
my $params;
my $label;

# GENERAL
my $clean;
my $stdout;
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
    'window=s'		=> \$window,
    'mode=s'   		=> \$mode,
    'referencedir=s'=> \$referencedir,
    'filename=s' 	=> \$filename,
    'binlevel=s' 	=> \$binlevel,
    'bindir=s' 		=> \$bindir,
    'clean'        	=> \$clean,
    'stdout=s' 		=> \$stdout,

	#### SPECIFIC
    'label=s'       => \$label,
    'params=s'      => \$params,
    'species=s'     => \$species,
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
	print "createIntervalFiles.pl    Printing STDOUT to stdoutfile:\n\n$stdout\n\n";

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
die "mode not supported: $mode (must be 'pileup' or 'gatk')" if $mode !~ /^(pileup|gatk)$/;

#### MAKE OUTPUT DIR IF NOT EXISTS
print "createIntervalFiles.pl    outputdir is a file: $outputdir\n" and exit if -f $outputdir;
File::Path::mkpath($outputdir) if not -d $outputdir;
print "createIntervalFiles.pl    Can't create output directory: $outputdir\n" and exit if not -d $outputdir;

#### IF NOT DEFINED referencedir, SET IT BASED ON SPECIES
if ( not defined $referencedir )
{
	my $basedir = $conf->getKeyValue("data", uc($species));
	my $build = $conf->getKeyValue("data", uc($species) . "LATESTBUILD");
	$referencedir = "$basedir/$build/fasta";
}
print "cumulativeSnp.pl    referencedir: $referencedir\n";
die "Can't find referencedir: $referencedir\n" if not -d $referencedir;


#### CHECK INPUT DIRS
my @indirs = split ",", $inputdirs;
foreach my $indir ( @indirs )
{   
    print "createIntervalFiles.pl    Can't find inputdir: $indir\n" and exit if not -d $indir;
}

#### DEBUG
print "createIntervalFiles.pl    inputdirs: $inputdirs\n";
print "createIntervalFiles.pl    outputdir: $outputdir\n";
print "createIntervalFiles.pl    cluster: $cluster\n";
print "createIntervalFiles.pl    params: $params\n" if defined $params;


#### RETRIEVE COMMAND
my $command = "$0 @arguments";

#### INSTANTIATE createIntervalFiles
my $snp = SNP->new(
	{
		#### SPECIFIC INPUTS
		inputdirs  	=> \@indirs,
		outputdir   => $outputdir,
		window   	=> $window,
		referencedir=> $referencedir,
		filename   	=> $filename,
		params   	=> $params,
		binlevel   	=> $binlevel,
		bindir   	=> $bindir,
		label       => $label,
		params      => $params,
		species     => $species,

		#### GENERAL INPUTS
		clean     	=> $clean,
		gatk		=> $gatk,
		java		=> $java,
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
print "createIntervalFiles.pl    Doing createIntervalFiles()\n";
$snp->createIntervalFiles() if $mode =~ /^pileup$/;
$snp->createGatkIntervals() if $mode =~ /^gatk$/;

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "createIntervalFiles.pl    Run time: $runtime\n";
print "createIntervalFiles.pl    Completed $0\n";
print "createIntervalFiles.pl    ";
print Timer::datetime(), "\n";
print "createIntervalFiles.pl    ****************************************\n\n\n";
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


