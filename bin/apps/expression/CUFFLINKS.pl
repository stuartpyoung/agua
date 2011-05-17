#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     CUFFLINKS

    PURPOSE

        WRAPPER SCRIPT FOR RUNNING CUFFLINKS:

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

    ./CUFFLINKS.pl <--inputfiles String> [--matesfiles String] <--distance String> <--params String> <--outputdir String> <--reference String> <--label Integer>
    [--keep] [--queue String] [--maxjobs Integer] [--cpus Integer ][--help]

    --outputdir  :   Create this directory and write output files to it
    --reference  :   Location of directory containing reference subdirectories
    --distance   :   Average inner distance between mate pairs
    --gtf		 :   Location GTF reference file to ouput expression by gene/region
    --params     :   Additional parameters to be passed to CUFFLINKS
    --label      :   Name to used to submit jobs to cluster
    --queue      :   Cluster queue options
    --maxjobs    :   Max. number of concurrent cluster maxjobs (DEFAULT: 30)
	--chunks     :	 Any combination of hyphen- or comma-separated chunk numbers
    --cpus       :   Max. number of cpus per job (DEFAULT: 4)
    --help       :   print help info

    EXAMPLES

=over

 perl /nethome/bioinfo/apps/agua/0.4/bin/apps/CUFFLINKS.pl \
 --distance 200 \
 --queue "-q large" \
 --gtf /nethome/bioinfo/data/sequence/chromosomes/mouse/mm9/cufflinks/mm9-knownGene.gtf \
 --reference /nethome/bioinfo/data/sequence/chromosomes/mouse/mm9/bowtie \
 --species mouse \
 --inputdir /nethome/syoung/base/pipeline/bixby/run1/tophat/analysis1 \
 --outputdir /nethome/syoung/base/pipeline/bixby/run1/tophat/analysis5 \
 --cluster PBS \
 --cpus 4 \
 --maxjobs 1000 \
 --label run1+2 \
 --splitfile /nethome/syoung/base/pipeline/bixby/run1/tophat/analysis1/splitfiles \
 --chunks 1-10

=back

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
use CUFFLINKS;
use FileTools;
use Timer;
use Util;
use Conf::Agua;

##### STORE ARGUMENTS TO PRINT TO FILE LATER
my @arguments = @ARGV;

#### FLUSH BUFFER
$| = 1;

#### GET CONF 
my $conf = Conf::Agua->new(inputfile=>"$Bin/../../../conf/default.conf");
my $cufflinks = $conf->getKeyValue("applications", 'CUFFLINKS');
my $bowtie = $conf->getKeyValue("applications", 'BOWTIE');
my $samtools = $conf->getKeyValue("applications", 'SAMTOOLS');
my $cluster = $conf->getKeyValue("agua", 'CLUSTERTYPE');
my $qstat = $conf->getKeyValue("cluster", 'QSTAT');
my $qsub = $conf->getKeyValue("cluster", 'QSUB'); #### /usr/local/bin/qsub

#### SET msub TO qsub FOR ARRAY JOB SUBMISSION
$qsub = "/usr/local/bin/qsub";


#### GET OPTIONS
# GENERAL
my $maxlines = 4000000; #### MULTIPLE OF 4 BECAUSE EACH FASTQ RECORD IS 4 LINES
my $matefiles;
my $inputdir;
my $outputdir;
my $reference;
my $reads;
my $splitfile;
my $clean;
my $chunks;

# CUFFLINKS-SPECIFIC
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
my $queue = "-q gsmall";
my $dot = 1;
my $cleanup;

my $help;
if ( not GetOptions (

	#### GENERAL
    'outputdir=s'   => \$outputdir,
    'inputdir=s'   	=> \$inputdir,
    'reference=s' 	=> \$reference,
    'reads=i' 		=> \$reads,
    'splitfile=s' 	=> \$splitfile,
    'chunks=s'      => \$chunks,

	#### CUFFLINKS
    'distance=i' 	=> \$distance,
    'params=s'      => \$params,
    'gtf=s'        	=> \$gtf,
    'species=s'     => \$species,
    'samtoolsindex=s' => \$samtools_index,
    'label=s'       => \$label,
    'keep'        	=> \$keep,

	#### CLUSTER
    'maxjobs=i'     => \$maxjobs,
    'cluster=s' 	=> \$cluster,
    'queue=s'       => \$queue,
    'sleep=s' 		=> \$sleep,
    'verbose' 		=> \$verbose,
    'cpus=i'        => \$cpus,
    'help'          => \$help
) )
{ print "Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "Distance not defined (option --help for usage)\n" if not defined $distance;
die "Input directory not defined (Use --help for usage)\n" if not defined $inputdir;
die "Output directory not defined (Use --help for usage)\n" if not defined $outputdir;
die "Reference not defined (Use --help for usage)\n" if not defined $reference;
die "Species not defined (Use --help for usage)\n" if not defined $species;
die "Label not defined (Use --help for usage)\n" if not defined $label;

#### DEBUG

#### SET NUMBER OF LINES IF reads OPTION WAS USED
$maxlines = $reads * 4 if defined $reads;

#### IF NOT DEFINED samtools_index, SET IT BASED ON SPECIES
if ( not defined $samtools_index )
{
	my $parameter = "SAMTOOLS" . uc($species);
	$samtools_index = $conf->getKeyValue("data", $parameter);
}
die "Can't find samtools index directory: $samtools_index\n" if not -d $samtools_index;

#### GET ARGS IF queue DEFINED
my $args = '';
if ( defined $queue )
{
    if ( $queue =~ /^(.*)-q\s*(\S+)(.*)$/ )
    {
        $args = $1 . $3;
        $queue = $2;
    }
}

#### MAKE OUTPUT DIR IF NOT EXISTS
print "CUFFLINKS.pl    Output directory is a file: $outputdir\n" and exit if -f $outputdir;
File::Path::mkpath($outputdir) if not -d $outputdir;
print "CUFFLINKS.pl    Can't create output directory: $outputdir\n" if not -d $outputdir;

#### INSTANTIATE CUFFLINKS
my $runCufflinks = CUFFLINKS->new(
	{
		#### INPUTS (FROM USER)
		inputdir   	=> $inputdir,
		outputdir   => $outputdir,
		reference   => $reference,
		splitfile   => $splitfile,
		maxlines    => $maxlines,
		chunks    	=> $chunks,

		distance    => $distance,
		params      => $params,
		label       => $label,
		gtf       	=> $gtf,

		#### EXECUTABLES (FROM CONF)
		cufflinks	=> $cufflinks,
		bowtie      => $bowtie,
		cufflinks   => $cufflinks,
		samtools	=> $samtools,
		samtoolsindex	=> $samtools_index,

		#### CLUSTER (PRESETS AND USER)
		cluster 	=> $cluster,
		queue 		=> $queue,
		qstat 		=> $qstat,
		maxjobs     => $maxjobs,
		cpus        => $cpus,
		sleep       => $sleep,
		verbose 	=> $verbose,
		dot         => $dot,
		qsub        => $qsub
	}
);

#### RUN RNA PAIRED TRANSCRIPTOME ANALYSIS
$runCufflinks->run();

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "CUFFLINKS.pl    Run time: $runtime\n";
print "CUFFLINKS.pl    Completed $0\n";
print "CUFFLINKS.pl    ";
print Timer::datetime(), "\n";
print "CUFFLINKS.pl    ****************************************\n\n\n";
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


