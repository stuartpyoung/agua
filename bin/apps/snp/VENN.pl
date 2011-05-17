#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     VENN

    PURPOSE

		WRAPPER SCRIPT FOR DOING SNP ANNOTATION AND VERIFICATION

    INPUT

			1. SNP CALLS IN pileup FORMAT IN PARALLEL FOR EACH REFERENCE

			2. TAKES AS INPUT THE OUTPUTS OF BOWTIE.pm, TOPHAT.pm, AND ELAND.pm

			3. LOCATES pileup FILES USING THE STANDARD DIRECTORY ARCHITECTURE

				LAID OUT IN PARENT Cluster.pm

    OUTPUT

        1. MODIFIED pileup FORMAT FILE WITH APPENDED COLUMNS CONTAINING

			SNP ANNOTATION INFORMATION

	NOTES

			INHERITS Cluster.pm IN ORDER TO:

				1. ACCESS THE DIRECTORY ARCHITECTURE INFO TO LOCATE INPUT FILES

				2. RUN SCRIPTS ON AN HPC CLUSTER

    USAGE

    ./VENN.pl <--inputfiles String> <--matefiles String> <--outputdir String>
        <--referencedir String> [--splitfile String] [--reads Integer] [--convert]
        [--clean] [--queue String] [--maxjobs Integer] [--cpus Integer] [--help]

    --outputdir     :   Directory with one subdirectory per reference chromosome
						containing an out.sam or out.bam alignment output file
    --referencedir  :   Location of directory containing chr*.fa reference files
    --inputdirs		:   Tab-separated locations of directories containing input files
    --inputfile     :   Name of input file, e.g., "hit.snp"
    --outputfile    :   Name of output file, e.g., "hit.sav"
    --species       :   Name of the reference species (e.g., 'mouse')
    --queue         :   Cluster queue name
    --cluster       :   Cluster type (LSF|PBS)
    --help          :   print help info


	EXAMPLES




=cut

use strict;

#### EXTERNAL MODULES
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use Getopt::Long;
use FindBin qw($Bin);

#### USE LIBRARY
use lib "$Bin/../../../lib";	
use lib "$Bin/../../../lib/external";	

#### INTERNAL MODULES
use SNP;
use Timer;
use Util;
use Conf::Agua;

##### STORE ARGUMENTS TO PRINT TO FILE LATER
my @arguments = @ARGV;
unshift @arguments, $0;

#### FLUSH BUFFER
$| =1;

#### SET VENN LOCATION
my $conf = Conf::Agua->new(inputfile=>"$Bin/../../../conf/default.conf");
my $qstat = $conf->getKeyValue("cluster", 'QSTAT');
my $qsub = $conf->getKeyValue("cluster", 'QSUB'); #### /usr/local/bin/qsub
my $samtools = $conf->getKeyValue("applications", 'SAMTOOLS');

#### GET OPTIONS
#### GENERAL
my $outputdir;
my $bamfile;
my $querydir;
my $targetdir;
my $querylabel;
my $targetlabel;
my $suffix;
my $queryindex;
my $targetindex;
my $binlevel;
my $filename;
my $stdout;

#### SPECIFIC
my $verbose;

#### CLUSTER OPTIONS
my $tempdir;
my $cluster;
my $walltime;
my $queue = "gsmall";
my $maxjobs = 30;
my $cpus = 1;
my $sleep = 5;
my $parallel;
my $dot = 1;

my $help;
print "VENN.pl    Use option --help for usage instructions.\n" and exit if not GetOptions (

	#### GENERAL VENN VARIABLES
    'outputdir=s'	=> \$outputdir,
    'bamfile=s'	=> \$bamfile,
    'querydir=s'	=> \$querydir,
    'targetdir=s'	=> \$targetdir,
    'querylabel=s'	=> \$querylabel,
    'targetlabel=s'	=> \$targetlabel,
    'suffix=s'		=> \$suffix,
    'queryindex=s'	=> \$queryindex,
    'targetindex=s'	=> \$targetindex,
    'binlevel=s' 	=> \$binlevel,
    'filename=s' 	=> \$filename,
    'stdout=s' 		=> \$stdout,

	#### CLUSTER
    'maxjobs=i' 	=> \$maxjobs,
    'cpus=i'        => \$cpus,
    'walltime=i' 	=> \$walltime,
    'cluster=s' 	=> \$cluster,
    'queue=s' 		=> \$queue,
    'sleep=i' 		=> \$sleep,
    'verbose' 		=> \$verbose,
    'tempdir=s' 	=> \$tempdir,
    'help' 			=> \$help
);

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "outputdir not defined (Use --help for usage)\n" if not defined $outputdir;
die "bamfile not defined (Use --help for usage)\n" if not defined $bamfile;
die "querydir not defined (Use --help for usage)\n" if not defined $querydir;
die "targetdir not defined (Use --help for usage)\n" if not defined $targetdir;
die "querylabel not defined (Use --help for usage)\n" if not defined $querylabel;
die "targetlabel not defined (Use --help for usage)\n" if not defined $targetlabel;
die "suffix not defined (Use --help for usage)\n" if not defined $suffix;
die "queryindex not defined (Use --help for usage)\n" if not defined $queryindex;
die "targetindex not defined (Use --help for usage)\n" if not defined $targetindex;
print "Can't find outputdir: $outputdir\n" if not -d $outputdir;

#### PRINT TO STDOUT IF DEFINED stdout
print "Printing STDOUT to file:\n\n$stdout\n\n" if defined $stdout;
open(STDOUT, ">$stdout") or die "Can't redirect STDOUT to file: $stdout\n" if defined $stdout;
open(STDERR, ">>$stdout") or die "Can't redirect STDERR to file: $stdout\n" if defined $stdout;

#### INSTANTIATE OBJECT
my $snp = SNP->new(
	{
		outputdir 	=> $outputdir,
		bamfile 	=> $bamfile,
		querydir 	=> $querydir,
		targetdir 	=> $targetdir,
		querylabel 	=> $querylabel,
		targetlabel => $targetlabel,
		suffix 		=> $suffix,
		queryindex 	=> $queryindex,
		targetindex => $targetindex,
		binlevel   	=> $binlevel,
		filename   	=> $filename,
		querylabel  => $querylabel,
		targetlabel => $targetlabel,
		verbose 	=> $verbose,

		#### SAMTOOLS
		samtools 	=> $samtools,

		#### CLUSTER
		cluster 	=> $cluster,
		queue 		=> $queue,
		walltime	=> $walltime,
		maxjobs 	=> $maxjobs,
		cpus        => $cpus,
		qstat 		=> $qstat,
		qsub 		=> $qsub,
		sleep 		=> $sleep,
		tempdir 	=> $tempdir,
		dot 		=> $dot,

		command 	=>	\@arguments
	}
);

#### DO IT
$snp->snpToVenn();

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "VENN.pl    Run time: $runtime\n";
print "VENN.pl    Completed $0\n";
print "VENN.pl    ";
print Timer::current_datetime(), "\n";
print "VENN.pl    ****************************************\n\n\n";
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

