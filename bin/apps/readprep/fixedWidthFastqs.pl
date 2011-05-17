#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();

=head2

	APPLICATION     fixedWidthFastqs

    PURPOSE

        1. RUN fastq2fasta.pl FOR ALL .fastq/.fastq.gz FILES IN A DIRECTORY

        2. PRINT FASTA FILES TO USER-DEFINED OUTPUT DIRECTORY

        3. (OPTIONAL) COMPRESS FASTA FILES

    INPUT

        1. DIRECTORY CONTAINING .fastq OR fastq.gz FILES

        2. OUTPUT DIRECTORY

    OUTPUT

        1. FASTA FILES IN OUTPUT DIRECTORY

            mkdir -p /data/NGS/syoung/base/pipeline/SRA/SRA000271/chromosomes/NA18507/eland

        2. OUTPUT FROM fastq2fasta.pl

            /nethome/syoung/base/bin/nextgen/fastq2fasta.pl \
            --inputfile /mihg/data/NGS/syoung/base/pipeline/SRA/NA18507/SRX000600/SRR002271_1.fastq.gz \
            --outputfile /mihg/data/NGS/syoung/base/pipeline/SRA/SRA000271/chromosomes/NA18507/SRR002271_1.fasta \
            --dot 100000

	USAGE

    ./fixedWidthFastqs.pl  <--inputdir String> <--outputdir String> [--cluster Integer] [--maxjobs Integer] [--limit Integer] [--compress] [--fixed] [--id_length Integer] [--length Integer] [--compress String] [--dot Integer] [-h]

        --inputdir		       :   /Full/path/to/inputdir 
        --outputdir		       :   /Full/path/to/outputdir 
        --zipped               :   Input file is compressed (gzip, zip or bz2)
        --compress             :   Compress with gzip or zip
        --fixed                :   Print FASTA record of fixed length
        --id_length            :   Fixed id length
        --length      :   Fixed sequence length
        --cluster              :   Run on cluster with 'clusterSubmit.pl'
        --maxjobs                 :   Maximum number of concurrent maxjobs
        --dot                  :   Print counter every 'dot' number of records
        --help                 :   print help info

    EXAMPLES

/nethome/bioinfo/apps/agua/0.5/bin/apps/readprep/fixedWidthFastqs.pl \
--inputdir /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/SRX000600,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/SRX000601,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/SRX000602,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/SRX000603,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/SRX001539,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/SRX001540 \
--outputdir /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/SRX000600/fasta,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/SRX000601/fasta,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/SRX000602/fasta,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/SRX000603/fasta,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/SRX001539/fasta,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/SRX001540/fasta \
--compress gzip \
--cluster LSF \
--maxjobs 100 \
--queue large \
--compress gzip \
--id 60 \
--length 26 \
--dot 10000000
--dot 10 

=cut

use strict;

#### USE LIBRARY
use FindBin qw($Bin);
use lib "$Bin/../../../lib";
use lib "$Bin/../../../lib/external";

#### FLUSH BUFFER
$| = 1;

#### INTERNAL MODULES
use Timer;
use Util;
use Conf::Agua;
use ReadPrep;
use Monitor;


#### EXTERNAL MODULES
use Term::ANSIColor qw(:constants);
use Getopt::Long;
use IO::Pipe;
use File::Path;
use Data::Dumper;

#### SET DEFAULT SEPARATOR
my $DEFAULT_SEPARATOR = "\\s+";

#### SAVE ARGUMENTS
my @arguments = @ARGV;

#### GET OPTIONS
my $inputdir;
my $readfile;
my $outputdir;
my $dot = 1000000;
my $zipped;
my $compress;
my $cluster;
my $maxjobs = 30;
my $sleep = 5;
my $queue;
my $cleanup;
my $quality;
my $id;
my $length;
my $help;
if ( not GetOptions (
    'inputdir=s' 	=> \$inputdir,
    'outputdir=s' 	=> \$outputdir,
    'zipped=s' 		=> \$zipped,
    'compress=s' 	=> \$compress,
    'dot=i' 		=> \$dot,
    'cluster=s' 	=> \$cluster,
    'maxjobs=i' 	=> \$maxjobs,
    'queue=s' 		=> \$queue,
    'cleanup' 		=> \$cleanup,
    'sleep=i' 		=> \$sleep,
    'quality' 		=> \$quality,
    'id_length=i' 	=> \$id,
    'length=i' 		=> \$length,
    'help' 			=> \$help
) )
{ print "Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "Input directory not defined (Use --help for usage)\n" if not defined $inputdir;
die "Output directory not defined (Use --help for usage)\n" if not defined $outputdir;
die "compress must be 'gzip', 'zip' or 'bz2'\n" if defined $compress and $compress !~ /^(gzip|zip|bz2)$/;
die "zipped must be 'gzip', 'zip' or 'bz2'\n" if defined $zipped and $zipped !~ /^(gzip|zip|bz2)$/;
die "Jobs not defined\n" if defined $cluster and not defined $maxjobs;
die "Id length not defined\n" if not defined $id;
die "Sequence length not defined\n" if not defined $length;

#### SET EXECUTABLE
my $executable = "$Bin/fixedWidthFastq.pl";
print "executable: $executable\n";

#### GET DIRECTORIES
my $files = files($inputdir, $outputdir, $zipped);
print "No. files: ", scalar(@$files), "\n";

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
print "clusterObject:\n";
print Dumper $clusterObject;

#### COLLECT JOBS OR COMMANDS
my $counter = 0;
foreach my $file ( @$files )
{
	$counter++;

    #### SET FILES
    my $inputfile = $file->{inputfile};
    my $outputfile = $file->{outputfile};

    my $outdir = $file->{outdir};

	print "fixedWidthFastqs.pl    quitting because outdir is a file: $outdir\n" and exit if -f $outdir;
	File::Path::mkpath($outdir) if not -d $outdir;
	print "fixedWidthFastqs.pl    Can't create outdir: $outdir\n" if not -d $outdir;

	#### SET COMMAND
    my $command = qq{$executable --inputfile $inputfile};
	$command .= qq{ --outputfile $outputfile};
	$command .= " --compress $compress" if defined $compress;
    $command .= " --id $id " if defined $id;
    $command .= " --length $length " if defined $length;
    $command .= " --dot $dot " if defined $dot;

	#### SET LABEL
	my ($filename) = $inputfile =~ /([^\/]+)$/;
	$filename =~ s/\.gz$//;
	$filename =~ s/\.zip$//;
	$filename =~ s/\.bz2$//;
	$filename =~ s/\.fastq$//;
	my $label = "fixedWidthFastqs-$counter-$filename";

	#### SET JOB
	my $job = $clusterObject->setJob([$command], $label, $outdir);
	push @$jobs, $job;

    print "Ran inputfile: $inputfile\n";
}

#### RUN JOBS 
$clusterObject->runJobs($jobs, "fixedWidthFastqs");

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


=head2

	SUBROUTINE		files

	PURPOSE

		GENERATE ARRAY OF HASHES CONTAINING

		INPUTFILE AND OUTPUT FILES

=cut

sub files
{
	my $inputdir		=	shift;
	my $outputdir		=	shift;
	my $zipped			=	shift;

	#### ARRAY OF HASHES CONTAINING inputfile AND outputfile
	my $filearray = [];

	my @indirs = split ",", $inputdir;
	my @outdirs = split ",", $outputdir;

	#### COLLECT ALL THE *.fastq.gz INPUT FILES AND THEIR
	#### COMPLEMENTARY OUTPUT FILES
	for ( my $i = 0; $i < $#indirs + 1; $i++ )
	{
		my $indir = $indirs[$i];
		my $outdir = $outdirs[$i];

		#### CREATE OUTPUT DIRECTORY
		File::Path::mkpath($outdir) or die "Can't create output directory: $outdir\n" if not -d $outdir;

		my $files = Util::files($indir);
		die "No files in input directory: $inputdir\n" if not @$files or scalar(@$files) == 0;
		$files = Util::by_suffix($files, "\.fastq.gz") if defined $zipped and $zipped eq "gzip";
		$files = Util::by_suffix($files, "\.fastq.zip") if defined $zipped and $zipped eq "zip";
		$files = Util::by_suffix($files, "\.fastq.bz2") if defined $zipped and $zipped eq "bz2";
		$files = Util::by_suffix($files, "\.fastq") if not defined $zipped;

		#### SANITY CHECK
		die "No sequence files in input directory: $indir\n" if not @$files or scalar(@$files) == 0;

		foreach my $file ( @$files )
		{
			#### SET FILES
			my $inputfile = "$indir/$file";
			my $outputfile = "$outdir/$file";
			$outputfile =~ s/\.fastq.gz$/.fastq/;
			$outputfile =~ s/\.fastq.bz2$/.fastq/;
			$outputfile =~ s/\.fastq.zip$/.fastq/;

			my $hash = {};
			$hash->{inputfile} = $inputfile;
			$hash->{outputfile} = $outputfile;
			$hash->{outdir} = $outdir;
			push @$filearray, $hash;
		}
	}

	return $filearray;
}


sub usage
{
	print `perldoc $0`;
    exit;
}

