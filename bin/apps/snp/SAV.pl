#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     SAV

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

    ./SAV.pl <--inputfiles String> <--matefiles String> <--outputdir String>
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


/nethome/bioinfo/apps/agua/0.4/bin/apps/SAV.pl \
--dbsnp snp130 \
--species human \
--queue large \
--outputdir /scratch/syoung/base/pipeline/SRA/NA18507/maq/maq2 \
--referencedir /nethome/bioinfo/data/sequence/chromosomes/human/hg19/fasta \
--cluster LSF


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
use SAV;
use Timer;
use Util;
use Conf::Agua;

##### STORE ARGUMENTS TO PRINT TO FILE LATER
my @arguments = @ARGV;
unshift @arguments, $0;

#### FLUSH BUFFER
$| =1;

#### SET SAV LOCATION
my $conf = Conf::Agua->new(inputfile=>"$Bin/../../../conf/default.conf");
my $qstat = $conf->getKeyValue("cluster", 'QSTAT');
my $qsub = $conf->getKeyValue("cluster", 'QSUB'); #### /usr/local/bin/qsub
my $samtools = $conf->getKeyValue("applications", 'SAMTOOLS');

#### DEFAULT MAX FILE LINES (MULTIPLE OF 4 BECAUSE EACH FASTQ RECORD IS 4 LINES)
my $maxlines = 4000000;

#### GET OPTIONS
#### GENERAL
my $outputdir;
my $inputdirs;
my $binlevel;
my $filename;
my $referencedir;
my $dbdir;
my $tempfile;
my $stdout;

#### SAV-SPECIFIC
my $dbsnp;
my $verbose;
my $species;
my $chunksize = 50;

#### CLUSTER OPTIONS
my $tempdir;
my $cluster;
my $walltime;
my $queue;
my $maxjobs = 100;
my $cpus = 1;
my $sleep = 5;
my $parallel;
my $dot = 1;

my $help;
print "SAV.pl    Use option --help for usage instructions.\n" and exit if not GetOptions (

	#### GENERAL SAV VARIABLES
    'outputdir=s'	=> \$outputdir,
    'referencedir=s'=> \$referencedir,
    'inputdirs=s'	=> \$inputdirs,
    'binlevel=s' 	=> \$binlevel,
    'filename=s' 	=> \$filename,
    'dbdir=s'		=> \$dbdir,
    'tempfile=s'	=> \$tempfile,
    'dbsnp=s' 		=> \$dbsnp,
    'stdout=s' 		=> \$stdout,
    'species=s'     => \$species,
    'chunksize=i'     => \$chunksize,

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
die "inputdirs not defined (Use --help for usage)\n" if not defined $inputdirs;
die "referencedir not defined (Use --help for usage)\n" if not defined $referencedir;
die "dbsnp not defined (Use --help for usage)\n" if not defined $dbsnp;
die "tempfile not defined (Use --help for usage)\n" if not defined $tempfile;
die "species not defined (Use --help for usage)\n" if not defined $species;
print "Can't find outputdir: $outputdir\n" if not -d $outputdir;
print "Can't find dbdir: $dbdir\n" if defined $dbdir and not -f $dbdir;

#### SET DEFAULT dbdir
$dbdir = "$Bin/../t/04-Filter-SNP/dbfile" if not defined $dbdir;

#### PRINT TO STDOUT IF DEFINED stdout
print "Printing STDOUT to file:\n\n$stdout\n\n" if defined $stdout;
open(STDOUT, ">$stdout") or die "Can't redirect STDOUT to file: $stdout\n" if defined $stdout;
open(STDERR, ">>$stdout") or die "Can't redirect STDERR to file: $stdout\n" if defined $stdout;

#### DEBUG
print "SAV.pl    outputdir: $outputdir\n";
print "SAV.pl    referencedir: $referencedir\n";

#### SET INPUT DIRS
my $indirs;
@$indirs = split ",", $inputdirs;

#### INSTANTIATE SAV OBJECT
my $sav = SAV->new(
	{
		referencedir=> $referencedir,
		outputdir 	=> $outputdir,
		inputdirs 	=> $indirs,
		binlevel   	=> $binlevel,
		filename   	=> $filename,
		dbdir 		=> $dbdir,
		tempfile 	=> $tempfile,
		dbsnp 		=> $dbsnp,
		verbose 	=> $verbose,
		chunksize	=> $chunksize,

		#### SAMTOOLS	
		species		=> $species,
		samtools	=> $samtools,

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

#### RUN SAV AND DO SNP CALLS
$sav->run();

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "SAV.pl    Run time: $runtime\n";
print "SAV.pl    Completed $0\n";
print "SAV.pl    ";
print Timer::current_datetime(), "\n";
print "SAV.pl    ****************************************\n\n\n";
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

