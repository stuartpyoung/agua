#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     samToBam

    PURPOSE

        WRAPPER SCRIPT FOR RUNNING samToBam ASSEMBLY AND SNP PREDICTION

    INPUT

        1. ASSEMBLY DIRECTORY

        2. FASTA OR FASTQ INPUT FILES

    OUTPUT

        1. samToBam OUTPUT FILES IN ASSEMBLY DIRECTORY

    USAGE

    ./samToBam.pl <--inputfiles String> <--outputfiles String> <--outputdir String>
        <--reference String> [--splitfile String] [--reads Integer] [--convert]
        [--clean] [--queue String] [--maxjobs Integer] [--cpus Integer] [--help]

    --inputfiles      :   Comma-separated list of input file paths (/full/path/to/file)
    --outputfiles     :   Comma-separated list of output files (/full/path/to/file)
    --species         :   Use the samtools index for this species (e.g., 'mouse', 'human') 
    --reference       :   SAM file alignments are based on this reference (e.g., 'chr22')
    --samtoolsindex   :   Location of the directory containing the samtools index files
    --cluster         :   Type of job scheduler (e.g., 'LSF', 'PBS, 'SGE')
    --clean           :   Clean up shell scripts, error files, etc. after jobs are done
    --queue           :   Cluster queue options
    --maxjobs         :   Max. number of concurrent cluster maxjobs
    --cpus            :   Max. number of cpus per job
    --help            :   print help info


	EXAMPLES


/nethome/bioinfo/apps/agua/0.5/bin/apps/converters/samToBam.pl \
--inputfile /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/autochr22/eland/1/chr22/hit.sam \
--outputfile /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/autochr22/eland/1/chr22/hit.bam \
--reference chr22 \
--species human \
--cluster LSF \
--queue large


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
use Converter;
use Timer;
use Util;
use Conf::Agua;

##### STORE ARGUMENTS TO PRINT TO FILE LATER
my @arguments = @ARGV;
unshift @arguments, $0;

#### FLUSH BUFFER
$| =1;

#### SET samToBam LOCATION
my $conf = Conf::Agua->new(inputfile=>"$Bin/../../../conf/default.conf");
my $qstat = $conf->getKeyValue("cluster", 'QSTAT');
my $qsub = "/usr/local/bin/qsub";	#### USE QSUB FOR ARRAY JOBS (MSUB CAN'T)
my $samtools = $conf->getKeyValue("applications", 'SAMTOOLS');

#####################
#### GET OPTIONS ####
#####################

#### GENERAL
my $stdout;
my $inputfiles;
my $outputfiles;

#### samToBam-SPECIFIC
my $clean;			#### ONLY GENERATE SPLIT FILES IF THEY DON'T EXIST
my $label;
my $verbose;

# SAMTOOLS-SPECIFIC
my $samtools_index;
my $species;
my $reference;

#### CLUSTER OPTIONS
my $tempdir;
my $cluster = "PBS";
my $queue = "gsmall";
my $maxjobs = 1000;
my $cpus = 1;
my $sleep = 5;
my $parallel;
my $dot = 1;
my $walltime = 24; #### WALL TIME IN HOURS (INTEGER)

my $help;
print "samToBam.pl    Use option --help for usage instructions.\n" and exit if not GetOptions (

	#### GENERAL
    'inputfiles=s' 	=> \$inputfiles,
    'outputfiles=s' => \$outputfiles,	#### PAIRED END MATE
    'clean' 		=> \$clean,
    'label=s' 		=> \$label,
    'stdout=s' 		=> \$stdout,

	#### SAMTOOLS
    'samtoolsindex=s' => \$samtools_index,
    'species=s'     => \$species,
    'reference=s'   => \$reference,

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

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "samToBam.pl    inputfiles not defined (Use --help for usage)\n" if not defined $inputfiles;
die "samToBam.pl    reference not defined (Use --help for usage)\n" if not defined $reference;
die "samToBam.pl    neither species nor samtoolsindex are defined (Use --help for usage)\n" if not defined $species and not defined $samtools_index;

#### IF NOT DEFINED samtools_index, SET IT BASED ON SPECIES AND REFERENCE
if ( not defined $samtools_index )
{
	my $parameter = "SAMTOOLS" . uc($species);
	$samtools_index = $conf->getKeyValue("data", $parameter);
}
die "samToBam.pl    Can't find samtools index directory: $samtools_index\n" if not -d $samtools_index;

#### CHECK INFILES AND OUTFILES
my ($samfiles, $bamfiles);
@$samfiles = split ",", $inputfiles;
@$bamfiles = split ",", $outputfiles if defined $outputfiles;
print "samToBam.pl    number of inputfiles (", scalar(@$samfiles), ") not the same as number of outputfiles (", scalar(@$bamfiles), ") (Use --help for usage)\n" and exit if scalar(@$samfiles) != scalar(@$bamfiles);

#### SET DEFAULT LABEL
$label = "samToBam" if not defined $label;

#### PRINT TO STDOUT IF DEFINED stdout
print "samToBam.pl    Printing STDOUT to file:\n\n$stdout\n\n" if defined $stdout;
open(STDOUT, ">$stdout") or die "samToBam.pl    Can't open STDOUT file: $stdout\n" if defined $stdout;
open(STDERR, ">>$stdout") or die "samToBam.pl    Can't open STDOUT file: $stdout\n" if defined $stdout;

#### DEBUG

#### INSTANTIATE samToBam OBJECT
my $converter = Converter->new(
	{
		inputfiles 	=> $inputfiles,
		outputfiles => $outputfiles,
		clean 		=> $clean,
		label 		=> $label,
		verbose 	=> $verbose,

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

#### COLLECT JOBS
my $jobs = [];
for ( my $i = 0; $i < @$samfiles; $i++ )
{
	my ($samfile, $bamfile);
	$samfile = $$samfiles[$i];
	$bamfile = $$bamfiles[$i] if defined $outputfiles;
	$bamfile = $samfile if not defined $outputfiles;
	$bamfile =~ s/\.sam$// if not defined $outputfiles;
	$bamfile .= ".bam" if not defined $outputfiles;

	#### SET OUTDIR
	my ($outdir) = $bamfile =~ /^(.+?)\/[^\/]+$/;
	File::Path::mkpath($outdir) if not -d $outdir;
	print "samToBam.pl    Could not create output directory: $outdir\n" and exit if not -d $outdir;

	#### SET COMMAND
	my $command = "$samtools/samtools view -bt $samtools_index/$reference.fai -o $bamfile $samfile";

	### SET LABEL
	my $joblabel = "samToBam";
	$joblabel = "$label-$reference-$i";

	#### SET JOB
	my $job = $converter->setJob( [ $command ], $label, $outdir);
	push @$jobs, $job;
}

#### RUN JOBS
$converter->runJobs( $jobs, $label);


#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "samToBam.pl    Run time: $runtime\n";
print "samToBam.pl    Completed $0\n";
print "samToBam.pl    ";
print Timer::current_datetime(), "\n";
print "samToBam.pl    ****************************************\n\n\n";
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

