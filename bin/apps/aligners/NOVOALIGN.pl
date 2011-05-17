#!/usr/bin/perl -w

#### DEBUG
$DEBUG = 1;

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     NOVOALIGN

    PURPOSE

        WRAPPER SCRIPT FOR RUNNING NOVOALIGN:

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

    ./NOVOALIGN.pl <--inputfiles String> [--matesfiles String] <--distance String> <--params String> <--outputdir String> <--reference String> <--label Integer>
    [--keep] [--queue String] [--maxjobs Integer] [--cpus Integer ][--help]

    --inputfiles          :   Comma-separated '*sequence.txt' file names (e.g., 's_1_1_sequence.txt,s_2_1_sequence.txt')
    --matefiles           :   Comma-separated '*sequence.txt' mate file names
    --params              :   Additional parameters to be passed to novoalign
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

/nethome/bioinfo/apps/agua/0.5/bin/apps/aligners/NOVOALIGN.pl \
--inputfiles /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/2000bp/100M/simpleheader/yoruba1-1.reads_1.sequence.txt \
--matefiles /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/2000bp/100M/simpleheader/yoruba1-1.reads_2.sequence.txt \
--reference /nethome/bioinfo/data/sequence/chromosomes/human/hg19/novoalign/chr22 \
--outputdir /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/2000bp/autochr22/novoalign \
--distance 2000 \
--deviation 300 \
--species human \
--label novoalign-autochr22 \
--maxjobs 1000 \
--queue large \
--cluster LSF \
--walltime 24 \
--cpus 1 \
--reads 1000000 \
--stdout /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/2000bp/autochr22/novoalign/autochr22-1.out 


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
use Aligner::NOVOALIGN;
use FileTools;
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
my $novoalign = $conf->getKeyValue("applications", 'NOVOALIGN');
my $samtools = $conf->getKeyValue("applications", 'SAMTOOLS');
my $cluster = $conf->getKeyValue("agua", 'CLUSTERTYPE');
my $qstat = $conf->getKeyValue("cluster", 'QSTAT');
my $qsub = $conf->getKeyValue("cluster", 'QSUB'); #### /usr/local/bin/qsub

#### SET msub TO qsub FOR ARRAY JOB SUBMISSION
$qsub = "/usr/local/bin/qsub";

#### GET OPTIONS
my $maxlines = 4000000; #### MULTIPLE OF 4 BECAUSE EACH FASTQ RECORD IS 4 LINES
my $threshold;
my $inputfiles;
my $matefiles;
my $outputdir;
my $reference;
my $reads;
my $splitfile;
my $clean;
my $chunks;
my $stdout;

# NOVOALIGN-SPECIFIC
my $distance;
my $deviation;
my $params;
my $gtf;
my $species;
my $samtools_index;
my $label;
my $keep;

# CLUSTER
my $maxjobs = 30;
my $sleep = 5;
my $verbose;
my $tempdir;
my $cpus = 1;
my $queue = "-q gsmall";
my $walltime = 24; #### WALL TIME IN HOURS (INTEGER)
my $dot = 1;
my $cleanup;

my $help;
if ( not GetOptions (

	#### GENERAL
    'threshold=i'  	=> \$threshold,
    'inputfiles=s'  => \$inputfiles,
    'matefiles=s'   => \$matefiles,
    'outputdir=s'   => \$outputdir,
    'reference=s' 	=> \$reference,
    'reads=i' 		=> \$reads,
    'splitfile=s' 	=> \$splitfile,
    'clean'        	=> \$clean,
    'chunks=s'      => \$chunks,
    'stdout=s' 		=> \$stdout,

	#### NOVOALIGN-SPECIFIC
    'distance=i' 	=> \$distance,
    'deviation=i' 	=> \$deviation,
    'params=s'      => \$params,
    'gtf=s'        	=> \$gtf,
    'species=s'     => \$species,
    'samtoolsindex=s' => \$samtools_index,
    'label=s'       => \$label,
    'keep'        	=> \$keep,

	#### CLUSTER
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

#### MAKE OUTPUT DIR IF NOT EXISTS
print "NOVOALIGN.pl    Output directory is a file: $outputdir\n" and exit if -f $outputdir;
File::Path::mkpath($outputdir) if not -d $outputdir;
print "NOVOALIGN.pl    Can't create output directory: $outputdir\n" if not -d $outputdir;

#### PRINT TO STDOUT IF DEFINED stdout
if ( defined $stdout )
{
	print "NOVOALIGN.pl    Printing STDOUT to stdoutfile:\n\n$stdout\n\n";

	my ($stdout_path) = $stdout =~ /^(.+)\/[^\/]+$/;
	File::Path::mkpath($stdout_path) if not -d $stdout_path;
	open(STDOUT, ">$stdout") or die "Can't open STDOUT file: $stdout\n" if defined $stdout;
}

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
print "inputfiles not defined (option --help for usage)\n" and exit  if not defined $inputfiles;
print "outputdir not defined (option --help for usage)\n" and exit  if not defined $outputdir;
print "reference not defined (option --help for usage)\n" and exit  if not defined $reference;
print "species not defined (option --help for usage)\n" and exit  if not defined $species;
print "label not defined (option --help for usage)\n" and exit  if not defined $label;
print "deviation not defined for paired end insert (option --help for usage)\n" and exit if defined $matefiles and not defined $deviation;
print "distance not defined for paired end insert (option --help for usage)\n" and exit if defined $matefiles and not defined $distance;

#### DEBUG
print "NOVOALIGN.pl    inputfiles: $inputfiles\n";
print "NOVOALIGN.pl    outputdir: $outputdir\n";
print "NOVOALIGN.pl    reference: $reference\n";
print "NOVOALIGN.pl    cluster: $cluster\n";
print "NOVOALIGN.pl    params: $params\n" if defined $params;

#### SET NUMBER OF LINES IF reads OPTION WAS USED
$maxlines = $reads * 4 if defined $reads;

#### IF NOT DEFINED samtools_index, SET IT BASED ON SPECIES
if ( not defined $samtools_index )
{
	my $parameter = "SAMTOOLS" . uc($species);
	$samtools_index = $conf->getKeyValue("data", $parameter);
}
die "Can't find samtools index directory: $samtools_index\n" if not -d $samtools_index;

#### CONVERT INPUTFILES AND MATEFILES INTO ARRAYS
my (@ins, @mates);
@ins = split ",", $inputfiles;
@mates = split ",", $matefiles if defined $matefiles;
@mates = () if not defined $matefiles;
print "NOVOALIGN.pl    no. files: ", $#ins + 1, "\n";
print "NOVOALIGN.pl    ins: @ins\n";
print "NOVOALIGN.pl    mates: @mates\n";


##### CHECK LENGTH OF ARRAYS MATCHES AND ORDER OF FILES 
#if ( defined $matefiles and $matefiles )
#{
#    my $filetool = FileTools->new();
#    my $mates_paired = $filetool->matesPaired(\@ins, \@mates);
#}


#### INSTANTIATE NOVOALIGN
my $runNovoalign = Aligner::NOVOALIGN->new(
	{
		#### INPUTS (FROM USER)
		threshold  	=> $threshold,
		inputfiles  => $inputfiles,
		matefiles   => $matefiles,
		distance   	=> $distance,
		deviation   => $deviation,
		params      => $params,
		outputdir   => $outputdir,
		referencedir   => $reference,
		splitfile   => $splitfile,
		maxlines    => $maxlines,
		chunks    	=> $chunks,
		clean     	=> $clean,
		label       => $label,
		keep		=> $keep,

		#### EXECUTABLES (PRESETS FROM CONF)
		novoalign   => $novoalign,
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


#### RUN
print "NOVOALIGN.pl    Doing run()\n";
$runNovoalign->run();


#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "NOVOALIGN.pl    Run time: $runtime\n";
print "NOVOALIGN.pl    Completed $0\n";
print "NOVOALIGN.pl    ";
print Timer::datetime(), "\n";
print "NOVOALIGN.pl    ****************************************\n\n\n";
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


