#!/usr/bin/perl -w

#### DEBUG
$DEBUG = 1;

#### TIME
my $time = time();

=head2

	APPLICATION     trimReads

    PURPOSE

        1. RUN filetools.pl trimRead FOR ALL .fastq/.fastq.gz FILES IN A DIRECTORY

    INPUT

        1. DIRECTORY CONTAINING .fastq OR fastq.gz FILES

    OUTPUT

        1. FOR EACH INPUT FILE:

			1. A TRIMMED FASTQ FILE CONTAINING READS OF THE SPECIFIED LENGTH

	USAGE

    ./trimReads.pl  <--inputdir String> <--outputdir String> [--cluster Integer] [--maxjobs Integer] [--limit Integer] [--compress] [--fixed] [--id_length Integer] [--sequence_length Integer] [--compress String] [--dot Integer] [-h]

	--inputdir		:   /Full/path/to/inputdir 
	--outputdir		:   /Full/path/to/outputdir 
	--type			:   Format of input file: fastq|fasta
	--min			:   Skip read if average quality < min
	--max			:   Skip read if average quality >= max
    --length		:   Only count this length of bases (after skip if defined)
	--compress		:   Compress with gzip or zip
	--cluster		:   Run on cluster with 'clusterSubmit.pl'
	--maxjobs		:   Maximum number of concurrent maxjobs
	--dot			:   Print counter every 'dot' number of records
	--help			:   print help info

    EXAMPLES


/nethome/bioinfo/apps/agua/0.5/bin/apps/readprep/trimReads.pl \
--inputdir /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/SRX000600,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/SRX000601,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/SRX000602,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/SRX000603,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/SRX001539,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/SRX001540 \
--outputdir /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/SRX000600,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/SRX000601,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/SRX000602,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/SRX000603,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/SRX001539,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/SRX001540 \
--compress gzip \
--bins -50,-40,-30,-20,-10,0,1,2,3,4,5,6,7,8,9,10,15,20,30,40,50,60,70,80,90,100,200,300 \
--type sanger \
--cluster LSF \
--maxjobs 600 \
--queue small \
--clean \
&> /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/samples/100M/trimReads.out


=cut

use strict;

#### USE LIBRARY
use FindBin qw($Bin);
use lib "$Bin/../../../lib";
use lib "$Bin/../../../lib/external";

#### FLUSH BUFFER
$| = 1;

#### INTERNAL MODULES
use Monitor;
use ReadPrep;
use SolexaUtil;
use Timer;
use Util;
use Conf::Agua;

#### INITIALISE SolexaUtil OBJECT
my $solexa = SolexaUtil->new();

#### EXTERNAL MODULES
use Term::ANSIColor qw(:constants);
use Getopt::Long;
use IO::Pipe;
use File::Path;
use Data::Dumper;

#### SAVE ARGUMENTS
my @arguments = @ARGV;

#### GET OPTIONS
my $inputdir;
my $outputdir;
#my $outputfile;
my $dot = 1000000;
my $compress;
my $zipped;
my $paired;
my $bins;
my $format;
my $type;
my $clean;
my $min;
my $max;
my $skip;
my $length;
my $cluster;
my $maxjobs = 30;
my $sleep = 5;
my $queue;
my $cleanup;
my $help;
if ( not GetOptions (
    'inputdir=s' 	=> \$inputdir,
    'outputdir=s' 	=> \$outputdir,
    'zipped=s' 		=> \$zipped,
    'compress=s' 	=> \$compress,
    'paired' 		=> \$paired,
    'bins=s' 		=> \$bins,
    'format=s' 		=> \$format,
    'type=s' 		=> \$type,
    'clean' 		=> \$clean,
    'min=i' 		=> \$min,
    'max=i' 		=> \$max,
    'skip=i' 		=> \$skip,
    'length=i' 		=> \$length,
    'dot=i' 		=> \$dot,
    'cluster=s' 	=> \$cluster,
    'maxjobs=i' 	=> \$maxjobs,
    'queue=s' 		=> \$queue,
    'cleanup' 		=> \$cleanup,
    'sleep=i' 		=> \$sleep,
    'help' 			=> \$help
) )
{ print "Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "inputdir not defined (Use --help for usage)\n" if not defined $inputdir;
die "outputdir not defined (Use --help for usage)\n" if not defined $outputdir;
#die "outputfile not defined (Use --help for usage)\n" if not defined $outputfile;
die "zipped must be 'gzip' or 'zip'\n" if defined $zipped and $zipped !~ /^(gzip|zip)$/;
die "compress must be 'gzip' or 'zip'\n" if defined $compress and $compress !~ /^(gzip|zip)$/;
die "maxjobs not defined\n" if defined $cluster and not defined $maxjobs;
die "skip not numeric\n" if defined $skip and not $skip =~ /^\d+$/;
die "length not numeric\n" if defined $length and not $length =~ /^\d+$/;
die "can't specify min and max at the same time\n" if defined $min and defined $max;

#### SET EXECUTABLE
my $executable = "$Bin/trimRead.pl";
print "executable: $executable\n";

#### RUN fastq2fasta.pl ON EVERY FILE
my $commands = [];  #### FOR RUNNING LOCALLY
my $jobs = [];		#### FOR RUNNING ON CLUSTER
my $clusterObject = ReadPrep->new(
	{
		cluster 	=> $cluster,
		queue 		=> $queue,
		cleanup 	=> $cleanup,
		maxjobs     => $maxjobs,
		sleep       => $sleep,
		dot         => $dot
	}
);

#### COLLECT JOBS OR COMMANDS
my $counter = 0;
my @indirs = split ",", $inputdir;
my @outdirs = split ",", $outputdir;
for ( my $i = 0; $i < $#indirs + 1; $i++ )
{
	my $indir = $indirs[$i];
	my $outdir = $outdirs[$i];
	print "trimReads.pl    indir: $indir\n";
	print "trimReads.pl    outdir: $outdir\n";

	my $files = Util::files($indir);
	$files = Util::by_suffix($files, "\.fastq.gz") if defined $zipped and $zipped eq "gzip";
	$files = Util::by_suffix($files, "\.fastq.zip") if defined $zipped and $zipped eq "zip";
	$files = Util::by_suffix($files, "\.fastq") if not defined $zipped;
	die "No files in input directory: $inputdir\n" if not @$files or scalar(@$files) == 0;

	foreach my $inputfile ( @$files )
	{
		$counter++;

		#### SKIP IF IS MATE FILE AND paired IS DEFINED
		next if $inputfile =~ /_2\./ and defined $paired;

		print "trimReads.pl    inputfile: $inputfile\n";

		#### SKIP IF OUTPUT FILE EXISTS ALREADY AND clean NOT SPECIFIED
		my $outputfile = "$outdir/$inputfile";
		$outputfile =~ s/\.gz$//;
		$outputfile =~ s/\.zip$//;
		print "trimReads.pl    skipping because found outputfile: $outputfile\n"
			and next if -f $outputfile and not defined $clean;
		print "trimReads.pl    outputfile: $outputfile\n";

		#### EXAMPLE: 
		#/nethome/bioinfo/apps/agua/0.5/bin/apps/readprep/trimRead.pl \
		#--inputfile /scratch/syoung/base/pipeline/SRA/NA12878/SRR001115.fastq \
		#--outputfile /scratch/syoung/base/pipeline/SRA/NA12878/SRR001115-35bp.fastq \
		#--type fastq \
		#--length 35 

		#### SET COMMAND
		my $command = qq{$executable --inputfile $indir/$inputfile};
		$command .= qq{ --outputfile $outputfile};
		$command .= qq{ --format $format};
		$command .= qq{ --type $type};
		$command .= qq{ --length $length};
		$command .= qq{ --compress $compress} if defined $compress;
		$command .= qq{ --paired $paired} if defined $paired;
		$command .= qq{ --min $min} if defined $min;
		$command .= qq{ --max $max} if defined $max;

#exit;

		#### SET LABEL
		my $filename = $inputfile;
		$filename =~ s/\.gz$//;
		$filename =~ s/\.zip$//;
		$filename =~ s/\.fastq$//;	
		my $label = "trimReads-$counter-$filename";

		#### SET JOB
		my $job = $clusterObject->setJob([$command], $label, $outdir);
		push @$jobs, $job;
	}
}

#### RUN JOBS
$clusterObject->runJobs($jobs, "trimReads");

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "\nRun time: $runtime\n";
print "Completed $0\n";
print Util::datetime(), "\n";
print "****************************************\n\n\n";
exit;

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#									SUBROUTINES
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#=head2
#
#	SUBROUTINE		files
#	
#	PURPOSE
#	
#		GENERATE ARRAY OF HASHES CONTAINING
#		
#		INPUTFILE AND OUTPUT FILES
#
#=cut
#
#sub files
#{
#	my $inputdir		=	shift;
#	my $compress		=	shift;
#	
#	#### COLLECT ALL THE INPUT FILES
#	my $filearray = [];
#	my @indirs = split ",", $inputdir;
#	for ( my $i = 0; $i < $#indirs + 1; $i++ )
#	{
#		my $indir = $indirs[$i];
#		my $files = Util::files($indir);
#		$files = Util::by_suffix($files, "\.fastq.gz") if defined $compress and $compress eq "gzip";
#		$files = Util::by_suffix($files, "\.fastq.zip") if defined $compress and $compress eq "zip";
#		$files = Util::by_suffix($files, "\.fastq") if not defined $compress;
#		die "No files in input directory: $inputdir\n" if not @$files or scalar(@$files) == 0;
#	
#		foreach my $file ( @$files )
#		{
#			push @$filearray, "$indir/$file";
#		}
#	}
#
#	return $filearray;
#}
#


sub usage
{
	print `perldoc $0`;
    exit;
}
