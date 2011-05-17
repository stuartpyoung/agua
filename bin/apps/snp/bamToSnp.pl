#!/usr/bin/perl -w

#### DEBUG
$DEBUG = 1;

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     bamToSnp

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

    ./bamToSnp.pl <--inputfiles String> [--matesfiles String] <--distance String> <--species String> <--outputdir String> <--reference String> <--label Integer>
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

perl /nethome/bioinfo/apps/agua/0.5/bin/apps/snp/bamToSnp.pl \
--reference chr22 \
--outputdir /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/maq/%REPLICATE% \
--species human \
--samfile hit.sam \
--label msq-bamToSnp \
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
use SNP;
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
my $chunksize;
my $inputdirs;
my $outputdir;
my $maxdepth;
my $filename;
my $binlevel;
my $bindir;
my $species;
my $samtools_index;
my $label = "bamToSnp";

# SPECIFIC
my $clean;
my $stdout;
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
    'chunksize=i'  	=> \$chunksize,
    'inputdirs=s'  	=> \$inputdirs,
    'outputdir=s'   => \$outputdir,
    'maxdepth=s'   => \$maxdepth,
    'filename=s' 	=> \$filename,
    'binlevel=s' 	=> \$binlevel,
    'bindir=s' 		=> \$bindir,
    'species=s'     => \$species,
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

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "queue not defined (option --help for usage)\n" if not defined $queue;
die "inputdirs not defined (option --help for usage)\n" if not defined $inputdirs;
die "outputdir not defined (Use --help for usage)\n" if not defined $outputdir;
die "species not defined (Use --help for usage)\n" if not defined $species;
die "label not defined (Use --help for usage)\n" if not defined $label;
die "filename not defined (Use --help for usage)\n" if not defined $filename;

#### DEBUG
print "bamToSnp.pl    outputdir: $outputdir\n";
print "bamToSnp.pl    cluster: $cluster\n";

#### MAKE OUTPUT DIR IF NOT EXISTS
print "bamToSnp.pl    Output directory is a file: $outputdir\n" and exit if -f $outputdir;
File::Path::mkpath($outputdir) if not -d $outputdir;
print "bamToSnp.pl    Can't create output directory: $outputdir\n" if not -d $outputdir;

die "bamToSnp.pl    neither species nor samtoolsindex are defined (Use --help for usage)\n" if not defined $species and not defined $samtools_index;

#### IF NOT DEFINED samtools_index, SET IT BASED ON SPECIES AND REFERENCE
if ( not defined $samtools_index )
{
	my $parameter = "SAMTOOLS" . uc($species);
	$samtools_index = $conf->getKeyValue("data", $parameter);
}
die "bamToSnp.pl    Can't find samtools index directory: $samtools_index\n" if not -d $samtools_index;

##### PRINT TO STDOUT IF DEFINED stdout
#if ( defined $stdout )
#{
#
#	my ($stdout_path) = $stdout =~ /^(.+)\/[^\/]+$/;
#	File::Path::mkpath($stdout_path) if not -d $stdout_path;
#	open(STDOUT, ">$stdout") or die "Can't open STDOUT file: $stdout\n" if defined $stdout;
#}

#### CHECK INPUT DIRS
my @indirs = split ",", $inputdirs;
foreach my $indir ( @indirs )
{   
    print "cumulativeSnp.pl    Can't find inputdir: $indir\n" and exit if not -d $indir;
}

#### INSTANTIATE bamToSnp
my $snp = SNP->new(
	{
		#### INPUTS (FROM USER)
		chunksize  	=> $chunksize,
		inputdirs  	=> \@indirs,
		outputdir   => $outputdir,
		maxdepth   	=> $maxdepth,
		filename 	=> $filename,
		binlevel 	=> $binlevel,
		bindir 		=> $bindir,

		#### SPECIFIC
		species     => $species,
		samtools	=> $samtools,
		samtoolsindex	=> $samtools_index,

		#### CLUSTER (PRESETS AND USER)
		clean     	=> $clean,
		label       => $label,
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
print "bamToSnp.pl    Doing run()\n";

#$snp->bamToSnp($outputdir, $references, $samfile);
#### STARTING WITH THE FIRST DIRECTORY AND GOING IN ORDER
####		
if ( not defined $binlevel )
{
	$snp->bamSnps();
}
else
{
	$snp->binBamToSnp();
}

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "bamToSnp.pl    Run time: $runtime\n";
print "bamToSnp.pl    Completed $0\n";
print "bamToSnp.pl    ";
print Timer::datetime(), "\n";
print "bamToSnp.pl    ****************************************\n\n\n";
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


