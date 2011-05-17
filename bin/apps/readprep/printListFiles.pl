#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();

=head2

	APPLICATION     printListFiles

    PURPOSE

        1. PRINT *.list FILES CONTAINING THE LOCATION OF A FIXED-WIDTH

			READ ON EACH LINE

    INPUT

		1. FASTA OR FASTA.GZ FILE

    OUTPUT

		OPTION: total

			1. TOTAL READS IN ALL READ FILES IN INPUT DIRECTORY

		OPTION: sample

			1. size NUMBER OF FASTA FILES, RECORDS RANDOMLY SELECTED FROM

				FILES IN THE INPUT DIRECTORY

    USAGE

    ./printListFiles.pl <--inputdir String> <--mode String> [-h]

    --inputdir             	:   /full/path/to/inputdir
	--outputdir				:	/full/path/to/outputdir
	--size					:	Desired no. reads per file
    --mode                 	:   Functional mode (total, samples)
    --compress				:   Input files are compressed (gzip|zip)
    --help					:   print help info

    EXAMPLES

/nethome/bioinfo/apps/agua/0.5/bin/apps/readprep/printListFiles.pl \
--inputdir /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/SRX000600/fasta,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/SRX000601/fasta,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/SRX000602/fasta,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/SRX000603/fasta,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/SRX001539/fasta,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/SRX001540/fasta \
--size 100000000 \
--compress gzip \
--cluster LSF \
--maxjobs 100 \
--queue large


=cut

use strict;

#### FLUSH BUFFER
$| = 1;

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
use Sampler;
use Timer;
use Util;
use Conf::Agua;
use ReadPrep;

#### GET OPTIONS
use Getopt::Long;
my $inputdir;
my $compress;
my $cluster;
my $maxjobs = 30;
my $sleep = 5;
my $queue;
my $size;
my $dot = 1000000;
my $help;
GetOptions (
    'inputdir=s' => \$inputdir,
    'compress=s' => \$compress,
    'size=i' => \$size,
    'dot=i' => \$dot,
    'cluster=s' => \$cluster,
    'maxjobs=i' => \$maxjobs,
    'queue=s' => \$queue,
    'sleep=i' => \$sleep,

    'help' => \$help             
) or die "No modes specified. Use --help for usage\n";
usage() if defined $help;

#### CHECK INPUTS
die "Input directory not defined (use --help for usage)\n" if not defined $inputdir;
die "Size not defined (use --help for usage)\n" if not defined $size;

#### SET ARGS
my $args =
{
	'inputdir'	=>	$inputdir,
	'compress'	=>	$compress,
	'size'		=>	$size,
	'dot'		=>	$dot
};


#### CREATE LIST FILES USING THE INFO FILES AND
#### RETURN THE TOTAL READ COUNTS FOR EACH DIRECTORY
my $samples = Sampler::listfiles($args);


#### SET EXECUTABLE 
my $executable = "$Bin/printListFile.pl";
print "printListFiles.pl    executable: $executable\n";


#### PRINT READ SUBFILES BASED ON READ LOCATIONS IN LIST FILES
print "printListFiles.pl    Generating read subfiles using locations in list files...\n";
my $directories;
@$directories = split ",", $inputdir;
print "printListFiles.pl    directories: \n";
print join "\n", @$directories;
print "\n";
print "printListFiles.pl    No. directories: ", scalar(@$directories), "\n";


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
print "cluster:\n";
print Dumper $clusterObject;

#### COLLECT JOBS
my $jobs = [];		#### FOR RUNNING ON CLUSTER
my $counter = 0;
foreach my $inputdir ( @$directories )
{
	print "printListFiles.pl    inputdir: $inputdir\n";
	my $infofiles 		=	$samples->{$inputdir}->{infofiles};

	#### FRACTION THIS DIRECTORY MAKES UP OF THE TOTAL READS
	print "printListFiles.pl    Doing read subfiles for inputdir: $inputdir\n";
	my $fraction 		=	$samples->{$inputdir}->{fraction};
	print "printListFiles.pl    This inputdir contains this fraction of the total reads: $fraction\n";
	print "printListFiles.pl    Size (desired no. reads per file): $size\n";

	my $total 		=	$samples->{$inputdir}->{total};
	print "printListFiles.pl    total reads in inputdir: $total\n";

	#### SET OUTPUT DIR FOR CLUSTER JOBS
	my ($outdir) = $inputdir =~ /^(.+?)\/[^\/]+$/;
	print "printListFiles.pl    outdir: $outdir\n";

	foreach my $info ( @$infofiles )
	{
		my @keys = keys %$info;
		my $infofile = $keys[0];

		#### SKIP MAKING MATE PAIR INFOFILE
		next if $infofile =~ /_2\./;

		#### INCREMENT JOB COUNTER
		$counter++;

		#### SET COMMAND
		my $command = qq{$executable --inputfile $infofile --dot 100000};
		$command .= " --compress $compress" if defined $compress;
		$command .= " --dot $dot" if defined $dot;
		$command .= " --size $size" if defined $size;
		$command .= " --fraction $fraction" if defined $fraction;
		$command .= " --total $total" if defined $total;

		#### SET JOB
		my $label = "printListFiles-" . $counter;
		my $job = $clusterObject->setJob([$command], $label, $outdir);
		push @$jobs, $job;
	}
}
print "printListFiles.pl    Collected ", scalar(@$jobs), " jobs\n";


#### RUN JOBS (OR COMMANDS IF LOCAL)
print "printListFiles.pl    Running jobs\n";
$clusterObject->runJobs($jobs, "printListFiles");
print "printListFiles.pl    Finished generating read list files\n";


#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "\nRun time: $runtime\n";
print "Completed $0\n";
print Timer::datetime(), "\n";
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




