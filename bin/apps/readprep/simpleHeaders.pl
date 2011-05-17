#!/usr/bin/perl -w

#### DEBUG
$DEBUG = 1;

#### TIME
my $time = time();

=head2

	APPLICATION     simpleHeaders

    PURPOSE

        1. RUN simpleHeader.pl FOR ALL .fastq/.fastq.gz FILES IN A DIRECTORY

		2. ADD MATE NUMBERS TO ARGUMENTS FOR simpleHeader.pl IF paired IS DEFINED

    INPUT

        1. DIRECTORY CONTAINING .fastq OR fastq.gz FILES

    OUTPUT

        1. FOR EACH INPUT FILE:

			1. TAB-SEPARATED *.avqual.tsv FILE

			2. TAB-SEPARATED FILE *.qualstats.tsv FILE

	USAGE

    ./simpleHeaders.pl  <--inputdir String> <--outputdir String> [--cluster Integer] [--maxjobs Integer] [--limit Integer] [--compress] [--fixed] [--id_length Integer] [--sequence_length Integer] [--compress String] [--dot Integer] [-h]

        --inputdir	:   /Full/path/to/inputdir 
        --outputdir	:   /Full/path/to/outputdir 
        --paired	:   Add mate number to arguments for simpleHeader.pl
        --barcode   :   Bar code value for sequence file
        --cleanup	:   Remove all run scripts after completion
        --dot       :   Print counter every 'dot' number of records
        --help      :   print help info

    EXAMPLES

/nethome/bioinfo/apps/agua/0.5/bin/apps/readprep/simpleHeaders.pl \
--inputdir /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/samples/100M \
--outputdir /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/samples/100M/simpleheader \
--paired \
--cluster LSF \
--maxjobs 600 \
--queue large \
&> /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/samples/100M/simpleheaders/simpleHeaders.out


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
my $barcode;
my $paired;
my $zipped;
my $clean;
my $dot = 1000000;
my $cluster;
my $maxjobs = 30;
my $sleep = 5;
my $queue;
my $cleanup;
my $help;
if ( not GetOptions (
    'inputdir=s' 	=> \$inputdir,
    'outputdir=s' 	=> \$outputdir,
    'barcode' 		=> \$barcode,
    'paired' 		=> \$paired,
    'zipped=s' 		=> \$zipped,
    'clean' 		=> \$clean,
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
die "maxjobs not defined\n" if defined $cluster and not defined $maxjobs;
die "zipped must be 'gzip' or 'zip'\n" if defined $zipped and $zipped !~ /^(gzip|zip)$/;

#### CREATE OUTPUT DIRECTORY
File::Path::mkpath($outputdir) if not -d $outputdir;
print "Can't create outputdir: $outputdir\n" if not -d $outputdir;

#### SET EXECUTABLE
my $executable = "$Bin/simpleHeader.pl";
print "executable: $executable\n";

#### GET DIRECTORIES
my $files = files($inputdir, $zipped);
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

#### COLLECT JOBS OR COMMANDS
my $counter = 0;
foreach my $inputfile ( @$files )
{
	$counter++;
	print "simpleHeaders.pl    inputfile: $inputfile\n";

	#### SET OUTPUT AND REJECT FILES
	my ($filename) = $inputfile =~ /^.+?\/([^\/]+)$/;
	$filename =~ s/\.gz$//;
	$filename =~ s/\.zip$//;
	$filename =~ s/\.fastq$//;

	my $outputfile = "$outputdir/$filename.sequence.txt";
	my $rejectfile = "$outputdir/$filename.reject.txt";

	#### SKIP IF OUTPUT FILE EXISTS ALREADY AND clean NOT SPECIFIED
	print "simpleHeaders.pl    skipping because found outputfile: $outputfile\n"
		and next if -f $outputfile and not defined $clean;
	print "simpleHeaders.pl    outputfile: $outputfile\n";

	my $matenumber;
	($matenumber) = $inputfile =~ /_((1|2))\./ if defined $paired;

	#### SET COMMAND
    my $command = qq{$executable --inputfile $inputfile};
	$command .= qq{ --outputfile $outputfile};
	$command .= qq{ --rejectfile $rejectfile};
	$command .= qq{ --barcode} if defined $barcode;
	$command .= qq{ --matenumber $matenumber} if defined $matenumber;

	#### SET LABEL
	my $label = "simpleHeaders-$counter-$filename";

	#### SET JOB
	my $job = $clusterObject->setJob([$command], $label, $outputdir);
	push @$jobs, $job;
}

#### RUN JOBS
$clusterObject->runJobs($jobs, "simpleHeaders");

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
	my $zipped		=	shift;

	#### COLLECT ALL THE INPUT FILES
	my $filearray = [];
	my @indirs = split ",", $inputdir;
	for ( my $i = 0; $i < $#indirs + 1; $i++ )
	{
		my $indir = $indirs[$i];
		my $files = Util::files($indir);
		$files = Util::by_suffix($files, "\.fastq.gz") if defined $zipped and $zipped eq "gzip";
		$files = Util::by_suffix($files, "\.fastq.zip") if defined $zipped and $zipped eq "zip";
		$files = Util::by_suffix($files, "\.fastq") if not defined $zipped;
		die "No files in input directory: $inputdir\n" if not @$files or scalar(@$files) == 0;

		foreach my $file ( @$files )
		{
			push @$filearray, "$indir/$file";
		}
	}

	return $filearray;
}



sub usage
{
	print `perldoc $0`;
    exit;
}
