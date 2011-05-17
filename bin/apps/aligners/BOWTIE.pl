#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     BOWTIE

    PURPOSE

        WRAPPER SCRIPT FOR RUNNING BOWTIE:

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

    ./BOWTIE.pl <--inputfiles String> [--matesfiles String] <--distance String> <--params String> <--outputdir String> <--reference String> <--label Integer> \
    [--queue String] [--maxjobs Integer] [--cpus Integer ] [--cluster String] [--walltime Integer] [--sleep Integer] [--tempdir String] [--cleanup] [--keep] [--verbose] [--help]

    --inputfiles    :   Comma-separated filenames (e.g., 's_1_1_sequence.txt,s_2_1_sequence.txt')
    --matefiles     :   Comma-separated '*sequence.txt' mate file names
    --params        :   Additional parameters to be passed to bowtie
    --outputdir     :   Create this directory and write output files to it
    --reference     :   Location of indexed reference file stub (e.g., /some/path/hg19)
    --species       :   Name of the reference species (e.g., 'human', 'mouse')
    --label         :   Submit jobs to cluster using this label
    --keep          :   Keep intermediate files
    --maxjobs       :   Max. number of concurrent cluster maxjobs (DEFAULT: 30)
    --cpus          :   Max. number of cpus per job (DEFAULT: 4)
    --cluster       :   Type of cluster job scheduler (e.g., 'PBS', 'SGE', 'LSF')
    --queue         :   Cluster queue options
	--walltime		:	Terminate job if run time exceeds this number of hours
    --cleanup       :   Delete *stdout.txt and *sh files once runs finished
    --sleep         :   Pause for this number of seconds waiting for *stdout.txt files
	--verbose       :	Print debug information
	--tempdir       :	Use this temporary output directory on execution host
    --help          :   print help info

    EXAMPLES

screen -S analysis1
mkdir -p /scratch/syoung/base/pipeline/bixby/run1/ln/bowtie-multi/analysis1/

cd /scratch/syoung/base/pipeline/bixby/run1/ln/bowtie-multi/analysis1/

perl /nethome/bioinfo/apps/agua/0.4/bin/apps/BOWTIE.pl \
--inputfiles /nethome/bioinfo/data/sequence/labs/bixby/run1/s_1_sequence.txt,/nethome/bioinfo/data/sequence/labs/bixby/run1/s_2_sequence.txt \
--params "--solexa1.3-quals" \
--distance 200 \
--queue "-q ngs" \
--gtf /nethome/bioinfo/data/sequence/chromosomes/mouse/mm9/cufflinks/mm9-knownGene.gtf \
--reference /nethome/bioinfo/data/sequence/chromosomes/mouse/mm9/bowtie \
--outputdir /scratch/syoung/base/pipeline/bixby/run1/ln/bowtie-multi/analysis1 \
--label bixln \
--keep \
--cluster LSF \
--cpus 4 \
--maxjobs 1000 \
--species mouse \
--clean




TEST WITH SHORT FILE:

cd /nethome/syoung/base/pipeline/bixby/run1/ln/bowtie-multi/

perl /nethome/bioinfo/apps/agua/0.5/bin/apps/BOWTIE.pl \
--inputfiles /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/samples/100M/1/min3length26-1.reads_1.sequence.txt \
--matefiles /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/samples/100M/1/min3length26-1.reads_2.sequence.txt \
--distance 200 \
--queue large \
--reference /nethome/bioinfo/data/sequence/chromosomes/human/hg19/bowtie/chr22 \
--outputdir /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/samples/100M/1/bowtie \
--keep \
--cluster LSF \
--cpus 4 \
--maxjobs 1000 \
--species human \
--label sample1-bowtie \
--reads 1000000 \
&> /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/samples/100M/1/bowtie/chr22/sample1-bowtie.out 



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
use BOWTIE;
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
my $bowtie = $conf->getKeyValue("applications", 'BOWTIE');
my $samtools = $conf->getKeyValue("applications", 'SAMTOOLS');
my $cluster = $conf->getKeyValue("agua", 'CLUSTERTYPE');
my $qstat = $conf->getKeyValue("cluster", 'QSTAT');
my $qsub = $conf->getKeyValue("cluster", 'QSUB'); #### /usr/local/bin/qsub

#### SET msub TO qsub FOR ARRAY JOB SUBMISSION
#$qsub = "/usr/local/bin/qsub";

#### GET OPTIONS
# GENERAL
my $maxlines = 4000000; #### MULTIPLE OF 4 BECAUSE EACH FASTQ RECORD IS 4 LINES
my $inputfiles;
my $matefiles;
my $outputdir;
my $reference;
my $reads;
my $splitfile;
my $clean;
my $chunks;
my $stdout;
my $check;

# BOWTIE-SPECIFIC
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
my $tempdir;
my $cpus = 1;
my $queue = "-q gsmall";
my $walltime = 24; #### WALL TIME IN HOURS (INTEGER)
my $dot = 1;
my $cleanup;

my $help;
if ( not GetOptions (

	#### GENERAL
    'inputfiles=s'  => \$inputfiles,
    'matefiles=s'   => \$matefiles,
    'outputdir=s'   => \$outputdir,
    'reference=s' 	=> \$reference,
    'reads=i' 		=> \$reads,
    'splitfile=s' 	=> \$splitfile,
    'clean'        	=> \$clean,
    'chunks=s'      => \$chunks,
    'stdout=s' 		=> \$stdout,
    'check' 		=> \$check,

	#### BOWTIE-SPECIFIC
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

#### MAKE OUTPUT DIR IF NOT EXISTS
print "BOWTIE.pl    Output directory is a file: $outputdir\n" and exit if -f $outputdir;
File::Path::mkpath($outputdir) if not -d $outputdir;
print "BOWTIE.pl    Can't create output directory: $outputdir\n" if not -d $outputdir;

#### PRINT TO STDOUT IF DEFINED stdout
if ( defined $stdout )
{
	print "BOWTIE.pl    Printing STDOUT to stdoutfile:\n\n$stdout\n\n";

	my ($stdout_path) = $stdout =~ /^(.+)\/[^\/]+$/;
	File::Path::mkpath($stdout_path) if not -d $stdout_path;
	open(STDOUT, ">$stdout") or die "Can't open STDOUT file: $stdout\n" if defined $stdout;
}

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "Distance not defined (option --help for usage)\n" if not defined $distance;
die "Inputfiles not defined (option --help for usage)\n" if not defined $inputfiles;
die "Output directory not defined (Use --help for usage)\n" if not defined $outputdir;
die "Reference not defined (Use --help for usage)\n" if not defined $reference;
die "Species not defined (Use --help for usage)\n" if not defined $species;
die "Label not defined (Use --help for usage)\n" if not defined $label;



#### DEBUG
print "BOWTIE.pl    inputfiles: $inputfiles\n";
print "BOWTIE.pl    outputdir: $outputdir\n";
print "BOWTIE.pl    reference: $reference\n";
print "BOWTIE.pl    cluster: $cluster\n";
print "BOWTIE.pl    params: $params\n" if defined $params;

#### SET NUMBER OF LINES IF reads OPTION WAS USED
$maxlines = $reads * 4 if defined $reads;

#### IF NOT DEFINED samtools_index, SET IT BASED ON SPECIES
if ( not defined $samtools_index )
{
	my $parameter = "SAMTOOLS" . uc($species);
	$samtools_index = $conf->getKeyValue("data", $parameter);
}
print "BOWTIE.pl    samtools index: $samtools_index\n";
die "Can't find samtools index directory: $samtools_index\n" if not -d $samtools_index;

#### CONVERT INPUTFILES AND MATEFILES INTO ARRAYS
my (@ins, @mates);
@ins = split ",", $inputfiles;
@mates = split ",", $matefiles if defined $matefiles;
@mates = () if not defined $matefiles;
print "BOWTIE.pl    no. files: ", $#ins + 1, "\n";
print "BOWTIE.pl    ins: @ins\n";
print "BOWTIE.pl    mates: @mates\n";

#### GET ARGS IF queue DEFINED
my $args = '';
if ( defined $queue )
{
    if ( $queue =~ /^(.*)-q\s*(\S+)(.*)$/ )
    {
        $args = $1 . $3;
        $queue = $2;
    }
    print "BOWTIE.pl    queue: $queue\n";
    print "BOWTIE.pl    args: $args\n";    
}

#### CHECK ALL INPUT FILES ARE *sequence.txt FILES
foreach my $infile ( @ins )
{   
    print "BOWTIE.pl    Input file is not a *_sequence.txt file: $infile\n" and exit if $infile !~ /sequence\.txt$/;
}

#### CHECK MATES IF AVAILABLE
if ( @mates and $#mates > 0 )
{
    foreach my $matefile ( @mates )
    {   
        print "BOWTIE.pl    Mate file is not a *_sequence.txt file: $matefile\n" and exit if $matefile !~ /sequence\.txt$/;
    }
}


##### CHECK LENGTH OF ARRAYS MATCHES AND ORDER OF FILES 
#if ( defined $matefiles and $matefiles )
#{
#    my $filetool = FileTools->new();
#    my $mates_paired = $filetool->matesPaired(\@ins, \@mates);
#}

#### INSTANTIATE BOWTIE
my $runBowtie = BOWTIE->new(
	{
		#### INPUTS (FROM USER)
		inputfiles  => $inputfiles,
		matefiles   => $matefiles,
		outputdir   => $outputdir,
		referencedir   => $reference,
		splitfile   => $splitfile,
		maxlines    => $maxlines,
		chunks    	=> $chunks,
		clean     	=> $clean,

		distance    => $distance,
		params      => $params,
		label       => $label,
		gtf       	=> $gtf,

		#### EXECUTABLES (FROM CONF)
		bowtie      => $bowtie,
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

#### RUN RNA PAIRED TRANSCRIPTOME ANALYSIS
if ( not defined $check )
{
	print "BOWTIE.pl    Doing run()\n";
	$runBowtie->run();
}
else
{
	print "BOWTIE.pl    Doing check()\n";
	$runBowtie->check();	
}

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "BOWTIE.pl    Run time: $runtime\n";
print "BOWTIE.pl    Completed $0\n";
print "BOWTIE.pl    ";
print Timer::datetime(), "\n";
print "BOWTIE.pl    ****************************************\n\n\n";
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


