#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     MAQ

    PURPOSE

        WRAPPER SCRIPT FOR RUNNING MAQ ASSEMBLY AND SNP PREDICTION

    INPUT

        1. ASSEMBLY DIRECTORY

        2. FASTA OR FASTQ INPUT FILES

    OUTPUT

        1. MAQ OUTPUT FILES IN ASSEMBLY DIRECTORY

    USAGE

    ./MAQ.pl <--inputfiles String> <--matefiles String> <--outputdir String>
        <--referencedir String> [--splitfile String] [--reads Integer] [--convert]
        [--clean] [--queue String] [--maxjobs Integer] [--cpus Integer] [--help]

    --inputfiles      :   Single FASTQ sequence file
    --matefiles       :   Single FASTQ mate pair file
    --outputdir       :   Create this directory and write output files to it
    --referencedir    :   Location of squashed genome reference files
    --species         :   Name of the reference species (e.g., 'mouse')
    --splitfile       :    Location of file containing list of split input files
    --reads           :   Number of reads per sub-file
    --convert         :   Convert from Solexa to Sanger FASTQ ('pre1.3' or 'post-1.3')
							(NB: MUST USE for all Solexa data)
    --clean           :   Clean run (remove old splitfile)
    --queue           :   Cluster queue options
    --maxjobs            :   Max. number of concurrent cluster maxjobs
    --cpus            :   Max. number of cpus per job
    --help            :   print help info


	EXAMPLES


cd /scratch/syoung/base/pipeline/SRA/NA18507/maq/maq1

/nethome/bioinfo/apps/agua/0.4/bin/apps/MAQ.pl \
--reads 500000 \
--cluster LSF \
--queue large \
--outputdir /scratch/syoung/base/pipeline/SRA/NA18507/maq/maq1 \
--inputfile /scratch/syoung/base/pipeline/SRA/NA18507/samples/reads_1.1.fastq \
--matefile /scratch/syoung/base/pipeline/SRA/NA18507/samples/reads_2.1.fastq \
--referencedir /nethome/bioinfo/data/sequence/chromosomes/human/hg19/maq \
--maxjobs 5000 \
--parallel \
--label maq1 \
--species human

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
use MAQ;
use Timer;
use Util;
use Conf::Agua;

##### STORE ARGUMENTS TO PRINT TO FILE LATER
my @arguments = @ARGV;
unshift @arguments, $0;

#### FLUSH BUFFER
$| =1;

#### SET MAQ LOCATION
my $conf = Conf::Agua->new(inputfile=>"$Bin/../../../conf/default.conf");
my $maq = $conf->getKeyValue("agua", 'MAQ');
my $qstat = $conf->getKeyValue("cluster", 'QSTAT');
my $qsub = "/usr/local/bin/qsub";	#### USE QSUB FOR ARRAY JOBS (MSUB CAN'T)
#my $qsub = $conf->getKeyValue("cluster", 'QSUB'); #### /usr/local/bin/qsub
my $samtools = $conf->getKeyValue("applications", 'SAMTOOLS');

#### DEFAULT MAX FILE LINES (MULTIPLE OF 4 BECAUSE EACH FASTQ RECORD IS 4 LINES)
my $maxlines = 4000000;

#### GET OPTIONS
#### GENERAL
my $stdout;
my $inputfiles;
my $matefiles;
my $readfile;
my $outputdir;
my $reads;
my $referencedir;
my $splitfile;

#### MAQ-SPECIFIC
my $clean;			#### ONLY GENERATE SPLIT FILES IF THEY DON'T EXIST
my $label;
my $convert;		#### CONVERT FROM SOLEXA TO SANGER QUALITIES
my $verbose;
my $solexa = "post-1.3";
my $params;

# SAMTOOLS-SPECIFIC
my $samtools_index;
my $species;

#### CLUSTER OPTIONS
my $tempdir;
my $cluster = "PBS";
my $queue = "gsmall";
my $maxjobs = 30;
my $cpus = 1;
my $sleep = 5;
my $parallel;
my $dot = 1;
my $walltime = 24; #### WALL TIME IN HOURS (INTEGER)

my $help;
print "MAQ.pl    Use option --help for usage instructions.\n" and exit if not GetOptions (

	#### GENERAL
    'inputfiles=s' 	=> \$inputfiles,
    'matefiles=s' 	=> \$matefiles,	#### PAIRED END MATE
    'outputdir=s'	=> \$outputdir,
    'referencedir=s'=> \$referencedir,
    'splitfile=s' 	=> \$splitfile,
    'clean' 		=> \$clean,
    'label=s' 		=> \$label,
    'stdout=s' 		=> \$stdout,

	#### MAQ
    'reads=i' 		=> \$reads,
    'convert=s' 	=> \$convert,	#### CONVERT pre-/post-GA Pipeline v1.3 TO SANGER FASTQ
    'params=s'      => \$params,

	#### SAMTOOLS
    'samtoolsindex=s' => \$samtools_index,
    'species=s'     => \$species,

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
);

#### SET NUMBER OF LINES IF reads OPTION WAS USED
$maxlines = $reads * 4 if defined $reads;

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "inputfiles not defined (Use --help for usage)\n" if not defined $inputfiles;
die "outputdir not defined (Use --help for usage)\n" if not defined $outputdir;
die "referencedir not defined (Use --help for usage)\n" if not defined $referencedir;
die "label not defined (Use --help for usage)\n" if not defined $label;
die "species not defined (Use --help for usage)\n" if not defined $species;

#### IF NOT DEFINED samtools_index, SET IT BASED ON SPECIES
if ( not defined $samtools_index )
{
	my $parameter = "SAMTOOLS" . uc($species);
	$samtools_index = $conf->getKeyValue("data", $parameter);
}
die "Can't find samtools index directory: $samtools_index\n" if not -d $samtools_index;

#### MAKE OUTPUT DIR IF NOT EXISTS
File::Path::mkpath($outputdir) if not -d $outputdir;
print "Could not create output directory: $outputdir\n" if not -d $outputdir;

#### PRINT TO STDOUT IF DEFINED stdout
if ( defined $stdout )
{
	my ($stdout_path) = $stdout =~ /^(.+)\/[^\/]+$/;
	File::Path::mkpath($stdout_path) if not -d $stdout_path;
	open(STDOUT, ">$stdout") or die "Can't open STDOUT file: $stdout\n" if defined $stdout;
}


#### PRINT TO STDOUT IF DEFINED stdout
print "Printing STDOUT to file:\n\n$stdout\n\n" if defined $stdout;
open(STDOUT, ">$stdout") or die "Can't open STDOUT file: $stdout\n" if defined $stdout;
open(STDERR, ">>$stdout") or die "Can't open STDOUT file: $stdout\n" if defined $stdout;

#### DEBUG
print "MAQ.pl    inputfiles: $inputfiles\n";
print "MAQ.pl    matefiles: $matefiles\n" if defined $matefiles;
print "MAQ.pl    outputdir: $outputdir\n";
print "MAQ.pl    referencedir: $referencedir\n";
print "MAQ.pl    convert: $convert\n" if defined $convert;

#### CHECK INPUTS
print "--convert type not supported: $convert (must be 'post-1.3' or 'pre1.3')\n" and exit if ( defined $convert and $convert !~ /^(post|pre)-1.3$/ );


#### INSTANTIATE MAQ OBJECT
my $runMaq = MAQ->new(
	{

		inputfiles 	=> $inputfiles,
		matefiles 	=> $matefiles,
		referencedir=> $referencedir,
		outputdir 	=> $outputdir,
		maxlines 	=> $maxlines,
		splitfile 	=> $splitfile,
		clean 		=> $clean,
		label 		=> $label,
		verbose 	=> $verbose,

		#### MAQ
		maq 		=> $maq,
		convert 	=> $convert,
		params      => $params,

		#### SAMTOOLS
		samtools	=> $samtools,
		samtoolsindex	=> $samtools_index,

		#### CLUSTER
		cluster 	=> $cluster,
		queue 		=> $queue,
		walltime 	=> $walltime,
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

#### RUN MAQ AND DO SNP CALLS
$runMaq->run();

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "MAQ.pl    Run time: $runtime\n";
print "MAQ.pl    Completed $0\n";
print "MAQ.pl    ";
print Timer::current_datetime(), "\n";
print "MAQ.pl    ****************************************\n\n\n";
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

