#!/usr/bin/perl -w

#### DEBUG
$DEBUG = 1;

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     samToSnp

    PURPOSE

        CONVERT SAM ALIGNMENTS INTO SNP PREDICTIONS

	VERSION		0.01

	HISTORY

		0.01 BASIC VERSION, ACCEPTS *.sam.gz GZIP FILES

    INPUT

        1. A GZIPPED *.sam.gz FILE

    OUTPUT

        1. A TSV *.snp FILE IN THE SPECIFIED OUTPUT DIRECTORY

    USAGE

    ./samToSnp.pl <--inputfiles String> [--matesfiles String] <--distance String> <--species String> <--outputdir String> <--reference String> <--label Integer>
    [--keep] [--queue String] [--maxjobs Integer] [--cpus Integer ][--help]

    --samfile		:   Name of gzipped SAM file (e.g., hit.sam.gz)
    --outputdir     :   Create this directory and write output files to it
    --species       :   Name of the reference species (e.g., 'human', 'mouse')
    --reference     :   Comma-separated names of references (e.g., 'chr22')
    --label         :   Name to used to submit maxjobs to cluster
    --keep          :   Keep intermediate files
    --cluster		:   Cluster type (e.g., 'LSF', 'PGS' or 'SGE')
    --queue         :   Cluster queue options
    --maxjobs       :   Max. number of concurrent cluster maxjobs (DEFAULT: 30)
    --walltime      :   Max. number of hours to spend running each job
    --cpus          :   Max. number of cpus per job (DEFAULT: 4)
    --help          :   print help info


TO DO INFO FOR:

	'species=s'		=> \$species,
    'samtoolsindex=s' => \$samtools_index,
    'stdout=s' 		=> \$stdout,
    'label=s'       => \$label,

	#### CLUSTER
    'keep'        	=> \$keep,
    'clean'        	=> \$clean,
    'sleep=s' 		=> \$sleep,
    'cleanup'       => \$cleanup,
    'verbose' 		=> \$verbose,
    'tempdir=s' 	=> \$tempdir,
    'help'          => \$help




    EXAMPLES

perl /nethome/bioinfo/apps/agua/0.5/bin/apps/snp/samToSnp.pl \
--reference chr22 \
--outputdir /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/maq/%REPLICATE% \
--species human \
--samfile hit.sam \
--label msq-samToSnp \
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
my $cluster;
my $outputdir;
my $reference;
my $samfile;
my $species;
my $samtools_index;
my $clean;
my $stdout;
my $label = "samToSnp";
my $keep;

# CLUSTER
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
    'reference=s' 	=> \$reference,
	'samfile=s'		=> \$samfile,
	'species=s'		=> \$species,
    'samtoolsindex=s' => \$samtools_index,
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

#### MAKE OUTPUT DIR IF NOT EXISTS
print "samToSnp.pl    Output directory is a file: $outputdir\n" and exit if -f $outputdir;
File::Path::mkpath($outputdir) if not -d $outputdir;
print "samToSnp.pl    Can't create output directory: $outputdir\n" if not -d $outputdir;
die "samToSnp.pl    reference not defined (Use --help for usage)\n" if not defined $reference;

die "samToSnp.pl    neither species nor samtoolsindex are defined (Use --help for usage)\n" if not defined $species and not defined $samtools_index;

#### IF NOT DEFINED samtools_index, SET IT BASED ON SPECIES AND REFERENCE
if ( not defined $samtools_index )
{
	my $parameter = "SAMTOOLS" . uc($species);
	$samtools_index = $conf->getKeyValue("data", $parameter);
}
die "samToSnp.pl    Can't find samtools index directory: $samtools_index\n" if not -d $samtools_index;


#### PRINT TO STDOUT IF DEFINED stdout
if ( defined $stdout )
{
	print "samToSnp.pl    Printing STDOUT to stdoutfile:\n\n$stdout\n\n";

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
print "samToSnp.pl    outputdir: $outputdir\n";
print "samToSnp.pl    reference: $reference\n";
print "samToSnp.pl    cluster: $cluster\n";

#### INSTANTIATE samToSnp
my $converter = Converter->new(
	{
		#### INPUTS (FROM USER)
		outputdir   => $outputdir,
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
print "samToSnp.pl    Doing run()\n";
$converter->samToSnp($outputdir, $references, $samfile);


#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "samToSnp.pl    Run time: $runtime\n";
print "samToSnp.pl    Completed $0\n";
print "samToSnp.pl    ";
print Timer::datetime(), "\n";
print "samToSnp.pl    ****************************************\n\n\n";
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


