#!/usr/bin/perl -w

#### DEBUG
$DEBUG = 1;

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     splitFiles

    VERSION         0.01

    PURPOSE

        SPLIT FASTQ INPUT FILES INTO SUBFILES, EACH IN ITS OWN SUBDIR

    INPUT

        1. OUTPUT DIRECTORY

        2. FASTA OR FASTQ INPUT FILES

    OUTPUT

        1. PRINT SUBFILES IN SUBDIRECTORIES OF OUTPUT DIRECTORY

	EXAMPLE

/nethome/syoung/0.5/bin/apps/aligners/splitFiles.pl \
--inputfiles /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/100M/simpleheader/yoruba1-1.reads_1.sequence.txt \
--matefiles /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/100M/simpleheader/yoruba1-1.reads_2.sequence.txt \
--outputdir /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland/1 \
--label eland-1 \
--maxjobs 2000 \
--walltime 24 \
--queue small \
--cluster LSF \
--cpus 1 \
--reads 250000


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
use Cluster;
use Sampler;
use Timer;
use Util;
use Conf::Agua;

##### STORE ARGUMENTS FOR PRINTING TO USAGE FILE LATER
my @arguments = @ARGV;
unshift @arguments, $0;

#### FLUSH BUFFER
$| =1;

#### DEFAULT MAX FILE LINES (MULTIPLE OF 4 BECAUSE EACH FASTQ RECORD IS 4 LINES)

#### GET OPTIONS
# GENERAL
my $stdout;
my $inputfiles;
my $matefiles;
my $outputdir;
my $splitfile;
my $clean;
my $maxlines;
my $reads;
my $label;
my $check;


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
    'maxlines=i' 	=> \$maxlines,
    'reads=i' 		=> \$reads,
    'splitfile=s' 	=> \$splitfile,
    'clean' 		=> \$clean,
    'label=s' 		=> \$label,
    'stdout=s' 		=> \$stdout,
    'check' 		=> \$check,

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

#### CHECK INPUTS
die "inputfiles not defined (Use --help for usage)\n" if not defined $inputfiles;
die "outputdir not defined (Use --help for usage)\n" if not defined $outputdir;
die "label not defined (Use --help for usage)\n" if not defined $label;
die "neither reads nor maxlines is defined (Use --help for usage)\n" if not defined $reads and not defined $maxlines;

#### reads OVERRIDES maxlines
$maxlines = $reads * 4 if defined $reads;

#### MAKE OUTPUT DIR IF NOT EXISTS
File::Path::mkpath($outputdir) if not -d $outputdir;
print "Could not create output directory: $outputdir\n" if not -d $outputdir;

#### SET DEFAULT SPLITFILE IF NOT DEFINED
$splitfile = "$outputdir/splitfile.txt" if not defined $splitfile;

#### PRINT TO STDOUT IF DEFINED stdout
if ( defined $stdout )
{
	my ($stdout_path) = $stdout =~ /^(.+)\/[^\/]+$/;
	File::Path::mkpath($stdout_path) if not -d $stdout_path;
	open(STDOUT, ">$stdout") or die "Can't open STDOUT file: $stdout\n" if defined $stdout;
}

#### PRINT HELP
if ( defined $help )	{	usage();	}


#### DEBUG
print "splitFiles.pl    inputfiles: $inputfiles\n";
print "splitFiles.pl    matefiles: $matefiles\n" if defined $matefiles;
print "outputdir: $outputdir\n";

#my $clusterObject = Cluster->new(
#	{
#		verbose 	=> $verbose,
#		
#		#### CLUSTER
#		maxjobs 	=> $maxjobs,
#		cpus        => $cpus,
#		cluster 	=> $cluster,
#		queue 		=> $queue,
#		walltime 	=> $walltime,
#		sleep 		=> $sleep,
#		tempdir 	=> $tempdir,
#		dot 		=> $dot,
#				
#		command 	=>	\@arguments
#	}
#);
##print Dumper $clusterObject;
#
#$clusterObject->doSplitfiles($splitfile, $label, $inputfiles, $matefiles, $outputdir, $maxlines, $clean);		

#### SET SUFFIX FOR Sampler::splitFiles FROM
#### SUFFIX OF LAST COMMA-SEPARATED INPUT FILE
my ($suffix) = $inputfiles =~ /(\.[^\.]{1,5})$/;

print "splitFiles.pl    splitfile: $splitfile\n";
print "splitFiles.pl    label: $label\n";
print "splitFiles.pl    outputdir: $outputdir\n";
print "splitFiles.pl    inputfiles: $inputfiles\n";
print "splitFiles.pl    matefiles: $matefiles\n";
print "splitFiles.pl    maxlines: $maxlines\n";
print "splitFiles.pl    suffix: $suffix\n";
print "splitFiles.pl    clean: $clean\n" if defined $clean;

Sampler::splitFiles(
	{
		'inputfiles' 	=> 	$inputfiles,
		'matefiles' 	=> 	$matefiles,
		'lines'			=> 	$maxlines,
		'splitfile'		=>	$splitfile,
		'outputdir'		=>	$outputdir,
		'label'			=>	$label,
		'suffix'		=> 	$suffix
	}
);	


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





