#!/usr/bin/perl -w

#### DEBUG
$DEBUG = 1;

#### TIME
my $time = time();

=head2

	APPLICATION     printSubfiles

    PURPOSE

        1. PRINT FASTA, QUAL AND FASTQ FILES GIVEN .list FILES

			FOR EACH INPUT READ FILE

    INPUTS

        1. DIRECTORY CONTAINING .fasta OR fasta.gz FILES AND .list FILES

    OUTPUTS

        1. FASTA, QUAL OR FASTQ FILES IN OUTPUT DIRECTORY

	USAGE

    ./printSubfiles.pl  <--inputdirs String> <--outputdir String> [--cluster Integer] [--maxjobs Integer] [--limit Integer] [--compress] [--fixed] [--id_length Integer] [--sequence_length Integer] [--compress String] [--dot Integer] [-h]

        --inputdirs		       :   /Full/path/to/inputdirs 
        --outputdir		       :   /Full/path/to/outputdir 
        --compress             :   Compress with gzip or zip
        --fixed                :   Print FASTA record of fixed length
        --id_length            :   Fixed id length
        --sequence_length      :   Fixed sequence length
        --cluster              :   Run on cluster with 'clusterSubmit.pl'
        --maxjobs                 :   Maximum number of concurrent maxjobs
        --dot                  :   Print counter every 'dot' number of records
        --help                 :   print help info

    EXAMPLES



/nethome/bioinfo/apps/agua/0.5/bin/apps/readprep/printSubfiles.pl \
--inputdirs /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/SRX000600/fasta,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/SRX000601/fasta,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/SRX000602/fasta,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/SRX000603/fasta,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/SRX001539/fasta,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/SRX001540/fasta \
--compress gzip \
--mode samples \
--size 10000000 \
&> /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/samples/100M/sampleReads.out

=cut

use strict;

#### USE LIBRARY
use FindBin qw($Bin);
use lib "$Bin/../../../lib";
use lib "$Bin/../../../lib/external";

#### FLUSH BUFFER
$| = 1;

#### INTERNAL MODULES
use Sampler;
use Timer;
use Util;
use Conf::Agua;
use ReadPrep;

#### EXTERNAL MODULES
use Term::ANSIColor qw(:constants);
use Getopt::Long;
use IO::Pipe;
use Data::Dumper;

#### SET DEFAULT SEPARATOR
my $DEFAULT_SEPARATOR = "\\s+";

#### SAVE ARGUMENTS
my @arguments = @ARGV;

#### GET OPTIONS
my $inputdirs;
my $readfile;
#my $paired;
my $dot = 1000000;
my $compress;
my $tempdir;
my $mode;
my $cluster;
my $maxjobs = 100;
my $sleep = 5;
my $queue = '';
my $size;
my $help;
if ( not GetOptions (
    'inputdirs=s' => \$inputdirs,
    #'paired' => \$paired,
    'compress=s' => \$compress,
    'tempdir=s' => \$tempdir,
    'mode=s' => \$mode,
    'size=i' => \$size,
    'dot=i' => \$dot,
    'cluster=s' => \$cluster,
    'maxjobs=i' => \$maxjobs,
    'queue=s' => \$queue,
    'sleep=i' => \$sleep,
    'help' => \$help
) )
{ print "Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "Input directory not defined (Use --help for usage)\n" if not defined $inputdirs;
#die "Output directory not defined (Use --help for usage)\n" if not defined $outputdir;
die "Compress type must be 'gzip' or 'zip'\n" if defined $compress and $compress !~ /^(gzip|zip)$/;
die "Jobs not defined\n" if defined $cluster and not defined $maxjobs;

#### SET EXECUTABLE 
my $executable = "$Bin/printSubfile.pl";

#### CREATE OUTPUT DIRECTORY
#mkdir($outputdir) or die "Can't create output directory: $outputdir\n" if not -d $outputdir;

my $directories;
@$directories = split ",", $inputdirs;


#### INITIALISE Cluster OBJECT
my $clusterObject = ReadPrep->new(
	{
		cluster 	=> $cluster,
		queue 		=> $queue,
		maxjobs     => $maxjobs,
		sleep       => $sleep,
		dot         => $dot
	}
) if defined $cluster;


#### COLLECT JOBS
my $commands = [];  #### FOR RUNNING LOCALLY
my $jobs = [];		#### FOR RUNNING ON CLUSTER
my $counter = 0;
foreach my $inputdir ( @$directories )
{
	print "inputdir: $inputdir\n";

	my $args =
	{
		'inputdir'	=>	$inputdir,
		'mode'		=>	$mode,
		'compress'	=>	$compress,
		'size'		=>	$size,
		'queue'		=>	$queue,
		'dot'		=>	$dot
	};

	my $infofiles = Sampler::infofiles($args);

	#### SET OUTPUT DIR FOR CLUSTER JOBS
	my ($outdir) = $inputdir =~ /^(.+?)\/[^\/]+$/;
	print "printSubfiles.pl    outdir: $outdir\n";

	foreach my $infofile ( @$infofiles )
	{
		#### GET INFO FILE NAME AND THEN INPUT FILE NAME
		my @keys = keys %$infofile;
		my $infofilename = $keys[0];
		my $width = $infofile->{$infofilename}->{totallength};
		my $inputfile = Sampler::inputfile($infofilename, $args);

		#### SET COMMAND
		my $command = qq{$executable --inputfile $inputfile};
		$command .= qq{ --mode $mode};
		$command .= qq{ --width $width};
		$command .= qq{ --compress $compress } if defined $compress;
		$command .= qq{ --tempdir $tempdir } if defined $tempdir;
		$command .= qq{ --dot $dot } if defined $dot;
#exit;		

		#### SET LABEL
		my ($filename) = $inputfile =~ /([^\/]+)$/;
		$filename =~ s/\.gz$//;
		$filename =~ s/\.zip$//;
		$filename =~ s/\.fastq$//;
		my $label = "printSubfiles-$counter-$filename";

		#### SET JOB
		my $job = $clusterObject->setJob([$command], $label, $outdir);
		push @$jobs, $job;

		$counter++;
	}
}

#### RUN JOBS (OR COMMANDS IF LOCAL)
print "printSubfiles.pl    Running " , scalar(@$commands), " jobs \n";
$clusterObject->runJobs($jobs, "printListFiles");
print "printSubfiles.pl    Finished generating read list files\n";




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


sub usage
{
	print `perldoc $0`;
    exit;
}

