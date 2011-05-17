#!/usr/bin/perl -w

#### DEBUG
$DEBUG = 1;

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     gzipSamfile

    PURPOSE

        WRAPPER SCRIPT FOR RUNNING gzipSamfile:

            1. RUN ONE LANE AT A TIME (MATE PAIRED)

            2. OUTPUT TO A USER-SPECIFIED DIRECTORY

	VERSION		0.02

	HISTORY

		LATER: ADD SPLIT FILES, CONVERSION OF SAM TO BAM
			MERGING OF BAM, AND RECONVERSION TO SAM

		0.02 ADDED CHUNK-BY-CHROMOSOME IF referencedir SPECIFIED
			AND RUN CUFFLINKS COMMAND

		0.01 BASIC VERSION

    INPUT

        1. LIST OF 'INPUTFILES' (s_*_1_sequence.txt)

        2. LIST OF 'MATEFILES' (s_*_2_sequence.txt)

        3. OUTPUT DIRECTORY CONTAINING THE SEQUENCES

    OUTPUT

        1. 

    NOTES

        1. IN THE CASE OF MULTIPLE INPUT FILES, THE READS IN ALL

            FILES MUST BE THE SAME LENGTH BECAUSE A SINGLE s_*_pair.xml

            FILE WILL BE USED FOR THE READS MERGED INTO ONE FILE

    USAGE

    ./gzipSamfile.pl <--inputfiles String> [--matesfiles String] <--distance String> <--params String> <--outputdir String> <--reference String> <--label Integer>
    [--keep] [--queue String] [--maxjobs Integer] [--cpus Integer ][--help]

    --inputfiles          :   Comma-separated '*sequence.txt' file names (e.g., 's_1_1_sequence.txt,s_2_1_sequence.txt')
    --matefiles           :   Comma-separated '*sequence.txt' mate file names
    --params              :   Additional parameters to be passed to bowtie
    --outputdir           :   Create this directory and write output files to it
    --reference           :   Location of indexed reference file
                              (NB: Use file stub, e.g., /some/path/hg19)
    --species             :   Name of the reference species (e.g., 'human', 'mouse')
    --label               :   Name to used to submit maxjobs to cluster
    --keep                :   Keep intermediate files
    --queue               :   Cluster queue options
    --maxjobs                :   Max. number of concurrent cluster maxjobs (DEFAULT: 30)
    --cpus                :   Max. number of cpus per job (DEFAULT: 4)
    --help                :   print help info

    EXAMPLES

perl /nethome/bioinfo/apps/agua/0.5/bin/apps/converters/gzipSamfile.pl \
--reference chr22 \
--outputdir /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/maq/1 \
--samfile hit.sam \
--label maq-gzipSamfile \
--keep \
--queue small \
--cluster LSF \
--cpus 1 \
--maxjobs 1000 \
--clean


=cut

use strict;

#### EXTERNAL MODULES
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use Getopt::Long;
use FindBin qw($Bin);
use File::Path;
use File::Copy;

#### USE LIBRARY
use lib "$Bin/../../../lib";	
use lib "$Bin/../../../lib/external";	

#### INTERNAL MODULES
use Converter;
use Timer;
use Util;
use Conf::Agua;

##### STORE ARGUMENTS FOR PRINTING TO USAGE FILE LATER
my @arguments = @ARGV;
unshift @arguments, $0;

#### FLUSH BUFFER
$| =1;

#### GET CONF 
my $conf = Conf::Agua->new(inputfile=>"$Bin/../../../conf/default.conf");
my $samtools = $conf->getKeyValue("applications", 'SAMTOOLS');
my $qstat = $conf->getKeyValue("cluster", 'QSTAT');
my $qsub = $conf->getKeyValue("cluster", 'QSUB'); #### /usr/local/bin/qsub

#### SET msub TO qsub FOR ARRAY JOB SUBMISSION
$qsub = "/usr/local/bin/qsub";
#exit;

#### GET OPTIONS
# GENERAL
my $outputdir;
my $reference;
my $samfile;
my $species;
my $samtools_index;
my $clean;
my $stdout;
my $label = "gzipSamfile";
my $keep;

# CLUSTER
my $cluster;
my $maxjobs = 30;
my $sleep = 5;
my $verbose;
my $tempdir;
my $cpus = 1;
my $queue;
my $walltime = 24; #### WALL TIME IN HOURS (INTEGER)
my $dot = 1;
my $cleanup;

my $help;
if ( not GetOptions (

	#### GENERAL
    'outputdir=s'   => \$outputdir,
	'samfile=s'		=> \$samfile,
	'reference=s'	=> \$reference,
    'stdout=s' 		=> \$stdout,
    'label=s'       => \$label,

	#### CLUSTER
    'keep'        	=> \$keep,
    'clean'        	=> \$clean,
    'maxjobs=i'     => \$maxjobs,
    'cpus=i'        => \$cpus,
    'cluster=s' 	=> \$cluster,
    'queue=s'       => \$queue,
    'walltime=i'    => \$walltime,
    'sleep=s' 		=> \$sleep,
    'cleanup'       => \$cleanup,
    'verbose' 		=> \$verbose,
    'tempdir=s' 	=> \$tempdir,
    'help'          => \$help
) )
{ print "Use option --help for usage instructions.\n";  exit;    };

#### PRINT TO STDOUT IF DEFINED stdout
if ( defined $stdout )
{
	print "gzipSamfile.pl    Printing STDOUT to stdoutfile:\n\n$stdout\n\n";

	my ($stdout_path) = $stdout =~ /^(.+)\/[^\/]+$/;
	File::Path::mkpath($stdout_path) if not -d $stdout_path;
	open(STDOUT, ">$stdout") or die "Can't open STDOUT file: $stdout\n" if defined $stdout;
}

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "outputdir not defined (Use --help for usage)\n" if not defined $outputdir;
die "reference not defined (Use --help for usage)\n" if not defined $reference;


#### DEBUG
print "gzipSamfile.pl    outputdir: $outputdir\n";
print "gzipSamfile.pl    reference: $reference\n";
print "gzipSamfile.pl    cluster: $cluster\n" if defined $cluster;

#### INSTANTIATE gzipSamfile
my $converter = Converter->new(
	{
		#### INPUTS (FROM USER)
		clean     	=> $clean,
		label       => $label,
		samtools	=> $samtools,
		samtoolsindex	=> $samtools_index,

		#### CLUSTER (PRESETS AND USER)
		cluster 	=> $cluster,
		queue 		=> $queue,
		walltime 	=> $walltime,
		qstat 		=> $qstat,
		maxjobs     => $maxjobs,
		cpus        => $cpus,
		sleep       => $sleep,
		cleanup     => $cleanup,
		verbose 	=> $verbose,
		tempdir 	=> $tempdir,
		dot         => $dot,
		qsub        => $qsub,

		command 	=>	\@arguments
	}
);

my $references;
@$references = split ",", $reference;

#### RUN RNA PAIRED TRANSCRIPTOME ANALYSIS
print "gzipSamfile.pl    Doing run()\n";
$converter->gzipFile($outputdir, $references, $samfile);


#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "gzipSamfile.pl    Run time: $runtime\n";
print "gzipSamfile.pl    Completed $0\n";
print "gzipSamfile.pl    ";
print Timer::datetime(), "\n";
print "gzipSamfile.pl    ****************************************\n\n\n";
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


