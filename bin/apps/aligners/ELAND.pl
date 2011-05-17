#!/usr/bin/perl -w

#### DEBUG
$DEBUG = 1;

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     ELAND

    VERSION         0.01

    PURPOSE

        WRAPPER SCRIPT FOR RUNNING ELAND ASSEMBLY

    INPUT

        1. ASSEMBLY DIRECTORY

        2. FASTA OR FASTQ INPUT FILES

    OUTPUT

        1. ELAND OUTPUT FILES IN ASSEMBLY DIRECTORIES

		2. CREATES ONE ASSEMBLY SUB-DIRECTORY PER CHROMOSOME IN

			THE USER-SPECIFIED OUTPUT DIRECTORY

	NOTES

		Running ELAND as a standalone program
		http://rulai.cshl.edu/xuan/Doc/Solexa/docs/Whole%20genome%20alignments%20using%20ELAND.html

		The eland_extended and eland_pair modes introduced in Pipeline version 0.3 share a common 'export' file format. The intention of this file is to combine all the important information for a lane into one file (or two files in the case of a inputtype read run). These scripts:

		Bring in the base quality value information
		Use the base quality values to pick the best alignment (taking into account read-pairing if relevant)
		Give a quality score to the alignment
		You do not get these files just by running the ELAND executable, the export files are generated by running several scripts to postprocess the raw ELAND output.
		There seems to be an increasingly common requirement to run this set of scripts in a 'standalone' fashion. This can be done via the script Pipeline/Eland/ELAND_standalone.pl. Minimal usage:

		./Pipeline/Eland/ELAND_standalone.pl -if read1.fastq -if read2.fastq -it fastq
		-eg /lustre/data01/Mondas_software/Genomes/E_coli_ELAND		



		ELAND_standalone.pl

			Usage: ./ELAND_standalone.pl options

			--input-file|if 			must specify at least one file, specify two for inputtype analysis
			--input-type|it 			type of input file: fastq, fasta, export or qseq
			--read-length|rl			will work it out from input file(s) if not specified
			--seed-length|sl			length of read used for ELAND alignment 
										- default is to min of read-length and 32
			--eland-genome|eg			full path of a squashed genome directory
			--output-prefix|op			will produce a set of output files whose names prefixed by this
										- defaults to 'reanalysis'  
			--pipeline-dir|pd			path of pipeline installation to call 
										- by default, looks in same dir as executable is found)
			--base-quality|bq			in fasta mode assume all bases have this quality
										- default is set to 30
			--pair-params|pp			inputtype read analysis parameters to pass to pickBestPair
										- multiple arguments must go in quotes
										- defaults to '--circular' - treats all chromosomes as circular


			NB: For a inputtype read analysis, both reads must share same input-type, read-length and seed-length

			NB: If you want to e.g. produce a set of export files for submission then it is probably best to re-run GERALD part of the pipeline to get a 'proper' analysis.

			NB: ELAND_standalone.pl uses Solexa's fastq encoding scheme (ASCII character = quality value + 64) rather than the more standard ASCII character = quality value + 32


			See: Data analysis - documentation : Whole genome alignments using ELAND
			Copyright (c) Illumina Inc. 2007. All rights reserved.
			Author: Anthony J. Cox
			http://bisr.osumc.edu/docs/Whole%20genome%20alignments%20using%20ELAND.html


			IMPORTANT CHANGES IN 1.4:

			- The format of the Firecrest output has changed. The intensity and noise
			  files are now generated cycle by cycle (same format as IPAR)
			- Bustard supports the binary data format generated by RTA (with the --CIF
			- The pipeline uses the data from the file RunInfo.xml (normally generated
			  by SCS) to identify the boundaries of the reads (including index reads).
			- Improved the estimation of the alignment scores of longer reads.
			- phageAlign produces export files


			NOTE: Running ELAND_standalone.pl is better than running the ELAND executable, which just gives the raw match output for all sequences. It does not take account of quality-filtering nor of base quality values. 

			IMPORTANT CHANGES IN 1.3:

			- A new build system is used. The installation still involves installing the
			  prerequisites and then typing "make" and "make install" in the top-level
			  pipeline folder. The executables can now be found in the bin/ directory
			  (e.g. bin/goat_pipeline.py).
			- The Bustard output formats have changed; a new file format called "qseq.txt"
			  is used to store read IDs, sequence and quality information as well as
			  filter information. The old file formats can be produced optionally with the
			  "--with-seq", "--with-prb", "--with-siq2", "--with-qval" options.
			- For base-call auto-calibration, the option "--with-qval" needs to be
			  specified at the goat_pipeline.py or bustard.py command line.
			- The Gerald analysis modes "expression" and "eland" are deprecated. They are
			  replaced by "eland_tag" and "eland_extended" respectively.
			- The "CONTAM_FILE" feature in PhageAlign mode is deprecated.

    USAGE

    ./ELAND.pl <--inputfiles String> <--outputdir String> <--referencedir String> <--lines Integer> [--convert] [--clean] [--help]

		--inputfiles				:   /full/path/to/readfile.fastq
									if inputtype-end: /../read_1.fastq,/../read_2.fastq
		--rundir				:	/full/path/to/pipeline/rundir
		--inputtype				:   type of input file: fastq, fasta, export or qseq
		--outputdir				:   /full/path/to/output_directory
		--referencedir			:   /full/path/to/squashed_genome_files
		--seedlength			:   Length of read used for ELAND alignment (Default: min. of seedlength and 32)
		--reads					:   Number of reads per sub-file
		--quality				:   Set quality value of all bases in fasta mode (Default: 30)
		--pairparams			:   inputtype read analysis parameters to pass to pickBestPair
										- multiple arguments must go in quotes
										- defaults to '--circular' - treats all chromosomes as circular
		--outputprefix			:	Prefix output files with this (Default: 'reanalysis')
		--clean					:   Clean run (remove old splitfile)
		--queue					:   Cluster queue options
		--jobs					:   Max. number of concurrent cluster jobs
		--help                 	:   print help info

    EXAMPLES

cd /p/NGS/syoung/base/pipeline/SRA/eland/1

/nethome/bioinfo/apps/agua/0.4/bin/apps/ELANDc.pl \
--inputfiles /p/NGS/syoung/base/pipeline/SRA/NA18507/samples/reads_ELAND_1.1.fastq,/p/NGS/syoung/base/pipeline/SRA/NA18507/samples/reads_ELAND_2.1.fastq \
--inputtype fastq \
--referencedir /q/nextGen/GA/refSeqs/chromosomes/hs_chr_sq \
--outputdir /p/NGS/syoung/base/pipeline/SRA/eland/1 \
--splitfile /p/NGS/syoung/base/pipeline/SRA/eland/1/reads.1.fastq.split \
--queue "-q psmall" \
--jobs 40

	\
	&> /p/NGS/syoung/base/pipeline/SRA/eland/1/outerr.txt


WHICH CALLS ELAND_standalone.pl

/nethome/syoung/base/apps/GAPipeline-1.4.0/bin/ELAND_standalone.pl \
--input-type fastq \
--eland-genome /q/nextGen/GA/refSeqs/chromosomes/hs_chr_sq \
--input-file /p/NGS/syoung/base/pipeline/SRA/eland/1/1/reads_ELAND_1.1.fastq \
--input-file /p/NGS/syoung/base/pipeline/SRA/eland/1/1/reads_ELAND_2.1.fastq \
... etc.

=cut

use strict;

#### EXTERNAL MODULES
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use Getopt::Long;
use FindBin qw($Bin);
use File::Path;

#### USE LIBRARY
use lib "$Bin/../../../lib";	
use lib "$Bin/../../../lib/external";	

#### INTERNAL MODULES
use ELAND;
use Sampler;
use Timer;
use Util;
use Conf::Agua;

##### STORE ARGUMENTS FOR PRINTING TO USAGE FILE LATER
my @arguments = @ARGV;
unshift @arguments, $0;

#### FLUSH BUFFER
$| =1;

#### SET ELAND LOCATION
my $conf = Conf::Agua->new(inputfile=>"$Bin/../../../conf/default.conf");
my $casava = $conf->getKeyValue("applications", 'CASAVA');
my $samtools = $conf->getKeyValue("applications", 'SAMTOOLS');
my $qstat = $conf->getKeyValue("cluster", 'QSTAT');
my $qsub = $conf->getKeyValue("cluster", 'QSUB'); #### /usr/local/bin/qsub

#### DEFAULT MAX FILE LINES (MULTIPLE OF 4 BECAUSE EACH FASTQ RECORD IS 4 LINES)
my $maxlines = 4000000;

#### GET OPTIONS
# GENERAL
my $stdout;
my $inputfiles;
my $matefiles;
my $rundir;
my $outputdir;
my $referencedir;
my $splitfile;
my $clean;
my $label;
my $readhits;
my $check;

# ELAND-SPECIFIC
my $seedlength;
my $inputtype;
my $reads;
my $quality;
my $pairparams;

# SAMTOOLS-SPECIFIC
my $samtools_index;
my $species;

# CLUSTER 
my $cluster;
my $maxjobs = 30;
my $walltime = 24; #### WALL TIME IN HOURS (INTEGER)
my $cpus = 1;
my $sleep = 5;
my $queue;
my $verbose;
my $tempdir;
my $parallel;
my $dot = 1;

my $help;
if ( not GetOptions (

	#### GENERAL
    'inputfiles=s' 	=> \$inputfiles,
    'matefiles=s' 	=> \$matefiles,
    'outputdir=s' 	=> \$outputdir,
    'referencedir=s' => \$referencedir,
    'splitfile=s' 	=> \$splitfile,
    'clean' 		=> \$clean,
    'label=s' 		=> \$label,
    'stdout=s' 		=> \$stdout,
    'readhits' 		=> \$readhits,
    'check' 		=> \$check,

	#### ELAND
    'reads=i' 		=> \$reads,
    'seedlength=i' 	=> \$seedlength,	
    'inputtype=s' 	=> \$inputtype,
    'quality=i' 	=> \$quality,
    'pairparams=s' 	=> \$pairparams,

	#### SAMTOOLS
    'samtoolsindex=s' => \$samtools_index,
    'species=s'     => \$species,

	#### CLUSTER
    'maxjobs=i' 	=> \$maxjobs,
    'cpus=i'        => \$cpus,
    'cluster=s' 	=> \$cluster,
    'queue=s' 		=> \$queue,
    'walltime=i'    => \$walltime,
    'sleep=i' 		=> \$sleep,
    'verbose' 		=> \$verbose,
    'tempdir=s' 	=> \$tempdir,
    'help' 			=> \$help
) )
{ print "Use option --help for usage instructions.\n";  exit;    };

#### SET NUMBER OF LINES IF reads OPTION WAS USED
$maxlines = $reads * 4 if defined $reads;

#### MAKE OUTPUT DIR IF NOT EXISTS
File::Path::mkpath($outputdir) if not -d $outputdir;
print "Could not create output directory: $outputdir\n" if not -d $outputdir;

#### PRINT TO STDOUT IF DEFINED stdout
if ( defined $stdout )
{
	my ($stdout_path) = $stdout =~ /^(.+)\/[^\/]+$/;
	File::Path::mkpath($stdout_path) if not -d $stdout_path;
	open(STDOUT, ">$stdout") or die "Can't open STDOUT file: $stdout\n" if defined $stdout;
}

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "Input file not defined (Use --help for usage)\n" if not defined $inputfiles;
die "Input type not defined (Use --help for usage)\n" if not defined $inputtype;
die "Output directory not defined (Use --help for usage)\n" if not defined $outputdir;
die "Reference directory not defined (Use --help for usage)\n" if not defined $referencedir;
die "Label not defined (Use --help for usage)\n" if not defined $label;
die "Input type not supported (must be fastq, fasta, export or qseq): $inputtype\n" if $inputtype !~ /^(fastq|fasta|export|qseq)$/i;
die "species not defined (Use --help for usage)\n" if not defined $species;

#### IF NOT DEFINED samtools_index, SET IT BASED ON SPECIES
if ( not defined $samtools_index )
{
	my $parameter = "SAMTOOLS" . uc($species);
	$samtools_index = $conf->getKeyValue("data", $parameter);
}
die "Can't find samtools index directory: $samtools_index\n" if not -d $samtools_index;

if ( defined $clean )
{
	my $message = "\nYou specified the option '--clean'.\nAre you sure you want to continue and overwrite any existing splitfiles?\n(Y to continue, N to exit)\n\n";
	print $message;

	exit if not Util::yes($message);
}

#### DEBUG
print "ELAND.pl    inputfiles: $inputfiles\n";
print "ELAND.pl    matefiles: $matefiles\n" if defined $matefiles;
print "outputdir: $outputdir\n";
print "referencedir: $referencedir\n";


my $runEland = ELAND->new(
	{
		inputfiles 	=> $inputfiles,
		matefiles 	=> $matefiles,
		referencedir => $referencedir,
		outputdir 	=> $outputdir,
		maxlines 	=> $maxlines,
		splitfile 	=> $splitfile,
		clean 		=> $clean,
		label 		=> $label,
		verbose 	=> $verbose,
		readhits 	=> $readhits,

		#### ELAND
		casava		=> $casava,
		seedlength 	=> $seedlength,
		inputtype 	=> $inputtype,
		quality	 	=> $quality,
		pairparams 	=> $pairparams,

		#### SAMTOOLS
		samtools	=> $samtools,
		samtoolsindex	=> $samtools_index,

		#### CLUSTER
		maxjobs 	=> $maxjobs,
		cpus        => $cpus,
		cluster 	=> $cluster,
		queue 		=> $queue,
		walltime 	=> $walltime,
		qstat 		=> $qstat,
		qsub 		=> $qsub,
		sleep 		=> $sleep,
		tempdir 	=> $tempdir,
		dot 		=> $dot,

		command 	=>	\@arguments
	}
);

if ( not defined $check )
{
	$runEland->run();		
}
else
{
	$runEland->check();	
}

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "\nRun time: $runtime\n";
print "Completed $0\n";
print Timer::datetime(), "\n";
print "****************************************\n\n\n";

#### PRINT TO STDOUT IF DEFINED stdout
close(STDOUT) or die "Can't close STDOUT file\n";
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


