#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     TOPHAT

    PURPOSE

        WRAPPER SCRIPT FOR RUNNING TOPHAT:

            1. RUN ONE LANE AT A TIME (MATE PAIRED)

            2. OUTPUT TO A USER-SPECIFIED DIRECTORY

	VERSION		0.03

	HISTORY

		LATER: ADD SPLIT FILES, CONVERSION OF SAM TO BAM
			MERGING OF BAM, AND RECONVERSION TO SAM

		0.03  ADDED PYRAMID MERGE AND IMPROVED USAGE STATS

		0.02 ADDED CHUNK-BY-CHROMOSOME IF referencedirdir SPECIFIED
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

    ./TOPHAT.pl <--inputfiles String> [--matesfiles String] <--distance String> <--params String> <--outputdir String> <--referencedir String> <--label Integer>
    [--keep] [--queue String] [--maxjobs Integer] [--cpus Integer ][--help]

    --inputfiles          :   Comma-separated '*sequence.txt' file names (e.g., 's_1_1_sequence.txt,s_2_1_sequence.txt')
    --matefiles           :   Comma-separated '*sequence.txt' mate file names
    --params              :   Additional parameters to be passed to tophat
    --outputdir           :   Create this directory and write output files to it
    --referencedir           :   Location of indexed referencedir file
                              (NB: Use file stub, e.g., /some/path/hg19)
    --species             :   Name of the referencedir species (e.g., 'human', 'mouse')
    --label               :   Name to used to submit maxjobs to cluster
    --keep                :   Keep intermediate files
    --queue               :   Cluster queue options
    --maxjobs                :   Max. number of concurrent cluster maxjobs (DEFAULT: 30)
    --cpus                :   Max. number of cpus per job (DEFAULT: 4)
    --help                :   print help info

    EXAMPLES

screen -S analysis1
mkdir -p /scratch/syoung/base/pipeline/bixby/run1/ln/tophat-multi/analysis1/

cd /scratch/syoung/base/pipeline/bixby/run1/ln/tophat-multi/analysis1/

perl /nethome/bioinfo/apps/agua/0.4/bin/apps/TOPHAT.pl \
--inputfiles /nethome/bioinfo/data/sequence/labs/bixby/run1/s_1_sequence.txt,/nethome/bioinfo/data/sequence/labs/bixby/run1/s_2_sequence.txt \
--params "--solexa1.3-quals" \
--distance 200 \
--queue "-q ngs" \
--gtf /nethome/bioinfo/data/sequence/chromosomes/mouse/mm9/cufflinks/mm9-knownGene.gtf \
--referencedir /nethome/bioinfo/data/sequence/chromosomes/mouse/mm9/bowtie \
--outputdir /scratch/syoung/base/pipeline/bixby/run1/ln/tophat-multi/analysis1 \
--label bixln \
--keep \
--cluster LSF \
--cpus 4 \
--maxjobs 1000 \
--species mouse \
--clean




TEST WITH SHORT FILE:

cd /nethome/syoung/base/pipeline/bixby/run1/ln/tophat-multi/

perl /nethome/bioinfo/apps/agua/0.4/bin/apps/TOPHAT.pl \
--inputfiles /nethome/bioinfo/data/sequence/labs/bixby/run1/s_1-100000_sequence.txt \
--params "--mate-inner-dist 200 --solexa1.3-quals" \
--distance 200 \
--queue "-q ngs" \
--gtf /nethome/bioinfo/data/sequence/chromosomes/mouse/mm9/cufflinks/mm9-knownGene.gtf \
--referencedir /nethome/bioinfo/data/sequence/chromosomes/mouse/mm9/bowtie \
--outputdir /nethome/syoung/base/pipeline/bixby/run1/cspg/tophat-multi \
--label bixln \
--keep \
--cluster LSF \
--cpus 4

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
use Expression::TOPHAT;
use Util::FileTools;
use Timer;
use Util;
use Conf::Agua;

##### STORE ARGUMENTS TO PRINT TO FILE LATER
my @arguments = @ARGV;
print "TOPHAT.pl    arguments: @arguments\n";

#### FLUSH BUFFER
$| =1;

#### GET CONF 
my $conf = Conf::Agua->new(inputfile=>"$Bin/../../../conf/default.conf");
my $tophat = $conf->getKeyValue("applications", 'TOPHAT');
my $bowtie = $conf->getKeyValue("applications", 'BOWTIE');
my $cufflinks = $conf->getKeyValue("applications", 'CUFFLINKS');
my $samtools = $conf->getKeyValue("applications", 'SAMTOOLS');
my $cluster = $conf->getKeyValue("agua", 'CLUSTERTYPE');
my $qstat = $conf->getKeyValue("cluster", 'QSTAT');
my $qsub = $conf->getKeyValue("cluster", 'QSUB'); #### /usr/local/bin/qsub
print "TOPHAT.pl    tophat: $tophat\n";
print "TOPHAT.pl    bowtie: $bowtie\n";

#### GET OPTIONS
# GENERAL
my $maxlines = 4000000; #### MULTIPLE OF 4 BECAUSE EACH FASTQ RECORD IS 4 LINES
my $inputfiles;
my $matefiles;
my $outputdir;
my $referencedir;
my $reads;
my $splitfile;
my $clean;
my $chunks;

# TOPHAT-SPECIFIC
my $distance;
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
my $tempdir;  #= "/tmp";
my $cpus = 1;
my $queue;
my $dot = 1;
my $cleanup;
my $walltime = 24; #### WALL TIME IN HOURS (INTEGER)

my $help;
if ( not GetOptions (

	#### GENERAL
    'inputfiles=s'  => \$inputfiles,
    'matefiles=s'   => \$matefiles,
    'outputdir=s'   => \$outputdir,
    'referencedir=s'=> \$referencedir,
    'reads=i' 		=> \$reads,
    'splitfile=s' 	=> \$splitfile,
    'clean'        	=> \$clean,
    'chunks=s'      => \$chunks,

	#### TOPHAT-SPECIFIC
    'distance=i' 	=> \$distance,
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

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
#die "Distance not defined (option --help for usage)\n" if not defined $distance;
die "Inputfiles not defined (option --help for usage)\n" if not defined $inputfiles;
die "Output directory not defined (Use --help for usage)\n" if not defined $outputdir;
die "referencedir not defined (Use --help for usage)\n" if not defined $referencedir;
die "Species not defined (Use --help for usage)\n" if not defined $species;
die "Label not defined (Use --help for usage)\n" if not defined $label;

#### DEBUG
print "TOPHAT.pl    inputfiles: $inputfiles\n";
print "TOPHAT.pl    outputdir: $outputdir\n";
print "TOPHAT.pl    referencedir: $referencedir\n";
print "TOPHAT.pl    cluster: $cluster\n";

print "TOPHAT.pl    SGE_CELL:\n";
print `echo \$SGE_CELL`;
print "TOPHAT.pl    SGE_QMASTER_PORT:\n";
print `echo \$SGE_QMASTER_PORT`;
print "TOPHAT.pl    SGE_EXECD_PORT:\n";
print `echo \$SGE_EXECD_PORT`;
print "TOPHAT.pl    ENV{'SGE_CELL'}: $ENV{'SGE_CELL'}\n";
print "TOPHAT.pl    ENV{'SGE_QMASTER_PORT'}: $ENV{'SGE_QMASTER_PORT'}\n";
print "TOPHAT.pl    ENV{'SGE_EXECD_PORT'}: $ENV{'SGE_EXECD_PORT'}\n";

print "TOPHAT.pl    env | grep SGE:\n";
print `env | grep SGE`;

#### SET NUMBER OF LINES IF reads OPTION WAS USED
$maxlines = $reads * 4 if defined $reads;

#### IF NOT DEFINED samtools_index, SET IT BASED ON SPECIES
if ( not defined $samtools_index )
{
	my $parameter = "SAMTOOLS" . uc($species);
	$samtools_index = $conf->getKeyValue("data", $parameter);
}
print "TOPHAT.pl    samtools index: $samtools_index\n";
die "Can't find samtools index directory: $samtools_index\n" if not -d $samtools_index;

#### CONVERT INPUTFILES AND MATEFILES INTO ARRAYS
my (@ins, @mates);
@ins = split ",", $inputfiles;
@mates = split ",", $matefiles if defined $matefiles;
@mates = () if not defined $matefiles;
print "TOPHAT.pl    no. files: ", $#ins + 1, "\n";
print "TOPHAT.pl    ins: @ins\n";
print "TOPHAT.pl    mates: @mates\n";

#### GET ARGS IF queue DEFINED
my $args = '';
if ( defined $queue )
{
    if ( $queue =~ /^(.*)-q\s*(\S+)(.*)$/ )
    {
        $args = $1 . $3;
        $queue = $2;
    }
    print "TOPHAT.pl    queue: $queue\n";
    print "TOPHAT.pl    args: $args\n";    
}

#### CHECK ALL INPUT FILES ARE *sequence.txt FILES
foreach my $infile ( @ins )
{   
    print "TOPHAT.pl    Input file is not a *_sequence.txt file: $infile\n" and exit if $infile !~ /sequence\.txt$/;
}

#### CHECK MATES IF AVAILABLE
if ( @mates and $#mates > 0 )
{
    foreach my $matefile ( @mates )
    {   
        print "TOPHAT.pl    Mate file is not a *_sequence.txt file: $matefile\n" and exit if $matefile !~ /sequence\.txt$/;
    }
}


##### CHECK LENGTH OF ARRAYS MATCHES AND ORDER OF FILES 
#if ( defined $matefiles and $matefiles )
#{
#    my $filetool = FileTools->new();
#    my $mates_paired = $filetool->matesPaired(\@ins, \@mates);
#}

#### MAKE OUTPUT DIR IF NOT EXISTS
print "TOPHAT.pl    Output directory is a file: $outputdir\n" and exit if -f $outputdir;
File::Path::mkpath($outputdir) if not -d $outputdir;
print "TOPHAT.pl    Can't create output directory: $outputdir\n" if not -d $outputdir;

#### INSTANTIATE TOPHAT
my $runTophat = Expression::TOPHAT->new(
	{
		#### INPUTS (FROM USER)
		conf		=>	$conf,
		inputfiles  => 	$inputfiles,
		matefiles   => 	$matefiles,
		outputdir   => 	$outputdir,
		referencedir=> 	$referencedir,
		splitfile   => 	$splitfile,
		maxlines    => 	$maxlines,
		chunks    	=> 	$chunks,
		clean     	=> 	$clean,

		distance    => 	$distance,
		params      => 	$params,
		label       => 	$label,
		gtf       	=> 	$gtf,

		#### EXECUTABLES (FROM CONF)
		tophat      => 	$tophat,
		bowtie      => 	$bowtie,
		cufflinks   => 	$cufflinks,
		samtools	=> 	$samtools,
		samtoolsindex=> $samtools_index,

		#### CLUSTER (PRESETS AND USER)
		cluster 	=> 	$cluster,
		queue 		=> 	$queue,
		walltime 	=> 	$walltime,
		qstat 		=> 	$qstat,
		maxjobs     => 	$maxjobs,
		cpus        => 	$cpus,
		sleep       => 	$sleep,
		cleanup     => 	$cleanup,
		verbose 	=> 	$verbose,
		tempdir 	=> 	$tempdir,
		dot         => 	$dot,
		qsub        => 	$qsub
	}
);

#### RUN RNA PAIRED TRANSCRIPTOME ANALYSIS
print "TOPHAT.pl    Doing run()\n";
$runTophat->run();

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "TOPHAT.pl    Run time: $runtime\n";
print "TOPHAT.pl    Completed $0\n";
print "TOPHAT.pl    ";
print Timer::datetime(), "\n";
print "TOPHAT.pl    ****************************************\n\n\n";
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


