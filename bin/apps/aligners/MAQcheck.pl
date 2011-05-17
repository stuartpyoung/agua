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
my $qsub = $conf->getKeyValue("cluster", 'QSUB'); #### /usr/local/bin/qsub

#### GET OPTIONS
my $min;
my $max;
my $paired;
my $replicates;
my $outputdir;
my $referencedir;
my $stdout;
my $label;
my $params;

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
    'min=s' 		=> \$min,
    'max=s' 		=> \$max,
    'paired'   		=> \$paired,
    'replicates=s' 	=> \$replicates,
    'outputdir=s'	=> \$outputdir,
    'referencedir=s'=> \$referencedir,
    'label=s' 		=> \$label,
    'stdout=s' 		=> \$stdout,
    'params=s'      => \$params,

	#### CLUSTER
    'maxjobs=i' 	=> \$maxjobs,
    'cpus=i'        => \$cpus,
    'cluster=s' 	=> \$cluster,
    'queue=s' 		=> \$queue,
    'walltime=i'    => \$walltime,
    'sleep=i' 		=> \$sleep,
    'tempdir=s' 	=> \$tempdir,
    'help' 			=> \$help
);

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "min not defined (Use --help for usage)\n" if not defined $min;
die "max not defined (Use --help for usage)\n" if not defined $max;
die "outputdir not defined (Use --help for usage)\n" if not defined $outputdir;
die "referencedir not defined (Use --help for usage)\n" if not defined $referencedir;
die "label not defined (Use --help for usage)\n" if not defined $label;
die "replicates not defined (Use --help for usage)\n" if not defined $replicates;

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

#### INSTANTIATE MAQ OBJECT
my $runMaq = MAQ->new(
	{
		#### INPUTS (FROM USER)
		min 		=> $min,
		max 		=> $max,
		paired   	=> $paired,
		replicates 	=> $replicates,
		referencedir=> $referencedir,
		outputdir 	=> $outputdir,
		label 		=> $label,
		params      => $params,

		#### EXECUTABLES (FROM CONF)
		maq 		=> $maq,

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

my ($completed, $sublabels, $missingfiles, $dubiousfiles);
($completed, $label, $sublabels, $missingfiles, $dubiousfiles) = $runMaq->check();	

#### SEND JOB COMPLETION SIGNAL
print "\n------------------------------------------------------------\n";
print "---[completed $label: $completed $sublabels]---";
if ( scalar(@$missingfiles) > 0 )
{
	print "\n";
	print scalar(@$missingfiles);
	print " missing file:\n" if scalar(@$missingfiles) == 1;
	print " missing files:\n" if scalar(@$missingfiles) != 1;
	print @$missingfiles if scalar(@$missingfiles) > 0;	
}
if ( scalar(@$dubiousfiles) > 0 )
{
	print "\n";
	print scalar(@$dubiousfiles);
	print " dubious file:\n" if scalar(@$dubiousfiles) == 1;
	print " dubious files:\n" if scalar(@$dubiousfiles) != 1;
	print "lowerbound\tupperbound\tsizemultiple\taveragesize\tfilesize\tlocation\n";
	print @$dubiousfiles if scalar(@$dubiousfiles) > 0;	
}
print "\n------------------------------------------------------------\n";

#### PRINT RUN TIME

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

