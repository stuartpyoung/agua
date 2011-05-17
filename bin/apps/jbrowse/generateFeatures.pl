#!/usr/bin/perl -w

#### DEBUG
$DEBUG = 1;

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     generateFeatures

	PURPOSE

		1. GENERATE JBROWSE FEATURES FOR ALL CHROMOSOMES IN PARALLEL

            USING flatfile-to-json.pl

		2. SUPPORTS GFF, GFF2, BED AND BAM INPUT FILES

        3. PREPARE INPUT FILES BY USING columnSplit.pl TO GENERATE

            ONE INPUT FILE IN EACH chr* DIRECTORY 

	INPUTS

		1. DIRECTORY CONTAINING INPUT GFF OR BAM FILES

			OR

			A GFF OR BAM FILE

		2. refSeqs.js FILE CONTAINING REFERENCE GENOME INFO

		3. A WORKING JBROWSE INSTALLATION

		4. OUTPUT DIRECTORY TO CREATE JBROWSE data DIRECTORY

	OUTPUTS

		1. CREATE DIRECTORIES FOR OUTPUT FILES
fs		
		2. COPY refSeqs.js TO EACH SUB DIRECTORY

	NOTES

		1. CREATE DIRECTORIES FOR OUTPUT FILES

		2. COPY refSeqs.js TO EACH SUB DIRECTORY

		3. RUN flatfile-to-json.pl AS ARRAY JOB

    USAGE

    ./generateFeatures.pl 
		<--inputdir String> \
		<--outputdir String> \
		<--filename String> \
		<--filetype String> \
		<--label String> \
		<--key String> \
		<--species String> \
		<--build String> \
		[--jbrowse String] \
		[--maxjobs=i] \
		[--cpus=i] \
		[--cluster String] \
		[--queue String] \
		[--sleep=i] \
		[--verbose] \
		[--tempdir String] \
		[--cleanup String] \
		[--help]

	REQUIRED JBROWSE-SPECIFIC ARGUMENTS

	inputdir 	:	Directory with chromosome subdirs containing input file
	outputdir 	:	JBrowse files will be printed to chromosome subdirs in this folder
    filename 	:	Name of input file (e.g., out.bam). If not defined, process all files. 
    filetype 	:	Type of input file (bam or gff)
    label 		:	Label for experiment (e.g., bowtie)
    key 		:	Key for features in JBrowse
    species 	:	Reference species (e.g., human, mouse, rat)
    build 		:	Reference species build (e.g., hg19, mm9, rn4)

	OPTIONAL (mostly CLUSTER) ARGUMENTS
    jbrowse 	:	Location of JBrowse installation
	maxjobs		:	Maximum number of jobs to be run concurrently
    cpus		:	Number of cpus per job
    cluster 	:	Type of cluster (e.g., LSF, PBS)
    queue 		:	Name of queue (e.g., gsmall)
    sleep		:	Length of sleep between queue polls
    verbose		:	Print verbose/debug information
    tempdir 	:	Use this temporary directory to write data on execution host
    cleanup 	:	Remove all script files and logs after completion of runs
    help    	:	Print this help information

	NOTES

		CALCULATE LOCATION OF refSeqs.js FILE USING SPECIES AND BUILD

		AND SAMTOOLS:

			1. CONF FILE CONTAINS SPECIES REFERENCE LOCATION, E.G.:

				HUMAN 	/nethome/bioinfo/data/sequence/chromosomes/human

			2. refSeqs.js FILE IS PLACED IN <BUILD>/jbrowse SUBFOLDER, E.G.:

				/nethome/bioinfo/data/sequence/chromosomes/human/hg19/jbrowse/refSeqs.js

	EXAMPLES

/nethome/bioinfo/apps/agua/0.5/bin/apps/generateFeatures.pl \
--build hg19 \
--cluster PBS \
--filename out.bam \
--filetype bam \
--inputdir /nethome/syoung/.agua/Project1/Workflow1/bowtie \
--key bowtie \
--label bowtie \
--outputdir /nethome/syoung/.agua/Project1/Workflow1/jbrowse \
--queue gsmall \
--species human \
--username syoung


=cut

use strict;

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
use Agua::View;
use Timer;
use Util;
use Conf::Agua;

##### STORE ARGUMENTS TO PRINT TO FILE LATER
my @arguments = @ARGV;
unshift @arguments, $0;

#### FLUSH BUFFER
$| =1;

#### GET CONF
my $conffile = "$Bin/../../../conf/default.conf";
my $conf = Conf::Agua->new(inputfile=>$conffile);
my $jbrowse = $conf->getKeyValue("data", 'JBROWSE');
my $qstat = $conf->getKeyValue("cluster", 'QSTAT');
my $qsub = $conf->getKeyValue("cluster", 'QSUB');
my $htmlroot = $conf->getKeyValue("agua", 'HTMLROOT');

#### DEFAULT MAX FILE LINES (MULTIPLE OF 4 BECAUSE EACH FASTQ RECORD IS 4 LINES)
my $maxlines = 4000000;

#### GET OPTIONS
my $username;
my $project;
my $workflow;
#my $outputdir;
my $inputdir;
my $filename;
my $filetype;
my $chunksize;
my $compress;
my $sortmem;
my $label;
my $refseqfile;
my $configfile;
my $key;
my $species;
my $build;

#### CLUSTER OPTIONS
my $tempdir;
#my $tempdir = "/tmp";
my $cluster;
my $queue = "gsmall";
my $maxjobs = 100;
my $cpus = 1;
my $sleep = 5;
my $parallel;
my $verbose;
my $dot = 1;
my $cleanup = undef;

my $help;
print "generateFeatures.pl    Use option --help for usage instructions.\n" and exit if not GetOptions (	
	#### JBROWSE
    'username=s'      => \$username,
    'project=s'       => \$project,
    'workflow=s'      => \$workflow,
    'inputdir=s'	=> \$inputdir,
    #'outputdir=s'	=> \$outputdir,
    'filename=s' 	=> \$filename,
    'filetype=s' 	=> \$filetype,
    'label=s' 		=> \$label,
    'refseqfile=s' 	=> \$refseqfile,
    'configfile=s' 	=> \$configfile,
    'chunksize=i' 	=> \$chunksize,
    'sortmem=i' 	=> \$sortmem,
    'compress' 	    => \$compress,
    'key=s' 		=> \$key,
    'species=s' 	=> \$species,
    'build=s' 		=> \$build,
    'jbrowse=s' 	=> \$jbrowse,
    'htmlroot=s' 	=> \$htmlroot,

	#### CLUSTER
    'maxjobs=i' 	=> \$maxjobs,
    'cpus=i'        => \$cpus,
    'cluster=s' 	=> \$cluster,
    'queue=s' 		=> \$queue,
    'sleep=i' 		=> \$sleep,
    'verbose' 		=> \$verbose,
    'tempdir=s' 	=> \$tempdir,
    'cleanup=s' 	=> \$cleanup,
    'help' 			=> \$help
);

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "username not defined (Use --help for usage)\n" if not defined $username;
die "project not defined (Use --help for usage)\n" if not defined $project;
die "workflow not defined (Use --help for usage)\n" if not defined $workflow;
#die "outputdir not defined (Use --help for usage)\n" if not defined $outputdir;
die "inputdir not defined (Use --help for usage)\n" if not defined $inputdir;
die "filetype not defined (Use --help for usage)\n" if not defined $filetype;
die "jbrowse not defined (Use --help for usage)\n" if not defined $jbrowse;
die "species not defined (Use --help for usage)\n" if not defined $species;
die "build not defined (Use --help for usage)\n" if not defined $build;

#### IF NOT DEFINED samtools_index, SET IT BASED ON SPECIES
my $speciesdir = $conf->getKeyValue("data", uc($species));
die "Can't find samtools index directory: $speciesdir\n" if not -d $speciesdir;
$refseqfile = "$speciesdir/$build/jbrowse/refSeqs.js" if not defined $refseqfile;
print "generateFeatures.pl    Can't find refseqfile: $refseqfile\n" and exit if not -f $refseqfile;
print "generateFeatures.pl    refseqfile is empty: $refseqfile\n" and exit if -z $refseqfile;

#### CREATE OUTPUT DIR IF NOT EXISTS
#File::Path::mkpath($outputdir) if not -d $outputdir;

#### CHECK FILETYPE IS SUPPORTED
print "filetype not supported (bam|gff): $filetype\n" and exit if $filetype !~ /^(bam|gff)$/;

#### DEBUG
print "generateFeatures.pl    inputdir: $inputdir\n";
print "generateFeatures.pl    filetype: $filetype\n";
print "generateFeatures.pl    refseqfile: $refseqfile\n";

#### INSTANTIATE JBrowse OBJECT
my $viewObject = Agua::View->new(
	{
        #### VIEW
		project	=> 	$project,
		workflow	=> 	$workflow,
		inputdir	=> 	$inputdir,
		username	=>	$username,
		htmlroot	=>	$htmlroot,    

		#### JBROWSE
		inputdir	=> 	$inputdir,
		#outputdir 	=> 	$outputdir,
		refseqfile 	=> 	$refseqfile,
		configfile 	=> 	$configfile,
		filetype	=>	$filetype,
		filename	=>	$filename,
		label		=>	$label,
		key			=>	$key,
		jbrowse		=>	$jbrowse,
		species		=>	$species,
		build		=>	$build,

		#### CLUSTER
		cluster 	=> $cluster,
		queue 		=> $queue,
		maxjobs 	=> $maxjobs,
		cpus        => $cpus,
		sleep 		=> $sleep,
		tempdir 	=> $tempdir,
		dot 		=> $dot,
		verbose 	=> $verbose,
		qstat 		=> $qstat,
		cleanup     => $cleanup,
		qsub        => $qsub,

		conf        => $conf,
		command 	=> \@arguments
	}
);

#### GENERATE FEATURES
$viewObject->generateFeatures();

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "generateFeatures.pl    Run time: $runtime\n";
print "generateFeatures.pl    Completed $0\n";
print "generateFeatures.pl    ";
print Timer::current_datetime(), "\n";
print "generateFeatures.pl    ****************************************\n\n\n";
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

