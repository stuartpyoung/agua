#!/usr/bin/perl -w

#### DEBUG
$DEBUG = 1;

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     sortBam

    PURPOSE

        RUN SAMTOOLS sort ON BAM FILE

    INPUT

        1. BAM FILE TO BE SORTED

    OUTPUT

        1. SORTED BAM FILE

    USAGE

    ./sortBam.pl <--inputfiles String> <--outputfiles String> [--clean]
		[--cluster String] [--queue String] [--maxjobs Integer] [--cpus Integer] [--help]

    --inputfiles      :   Comma-separated list of input file paths (/full/path/to/file)
    --outputfiles     :   Comma-separated list of output files (/full/path/to/file)
    --cluster         :   Type of job scheduler (e.g., 'LSF', 'PBS, 'SGE')
    --clean           :   Clean up shell scripts, error files, etc. after jobs are done
    --queue           :   Cluster queue options
    --maxjobs         :   Max. number of concurrent cluster maxjobs
    --cpus            :   Max. number of cpus per job
    --help            :   print help info


	EXAMPLES


/nethome/bioinfo/apps/agua/0.5/bin/apps/converters/sortBam.pl \
--inputfile /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/autochr22/eland/1/chr22/hit.Bam \
--outputfile /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/autochr22/eland/1/chr22/hit.sorted.bam \
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

#### SET sortBam LOCATION
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

#### sortBam-SPECIFIC
my $clean;			#### ONLY GENERATE SPLIT FILES IF THEY DON'T EXIST
my $label;
my $verbose;

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
print "sortBam.pl    Use option --help for usage instructions.\n" and exit if not GetOptions (

	#### GENERAL
    'inputfiles=s' 	=> \$inputfiles,
    'outputfiles=s' => \$outputfiles,	#### PAIRED END MATE
    'clean' 		=> \$clean,
    'label=s' 		=> \$label,
    'stdout=s' 		=> \$stdout,

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
die "sortBam.pl    inputfiles not defined (Use --help for usage)\n" if not defined $inputfiles;

#### CHECK INFILES AND OUTFILES
my ($infiles, $outfiles);
@$infiles = split ",", $inputfiles;
@$outfiles = split ",", $outputfiles if defined $outputfiles;
print "sortBam.pl    number of inputfiles (", scalar(@$infiles), ") not the same as number of outputfiles (", scalar(@$outfiles), ") (Use --help for usage)\n" and exit if scalar(@$infiles) != scalar(@$outfiles);

#### SET DEFAULT LABEL
$label = "sortBam" if not defined $label;

#### PRINT TO STDOUT IF DEFINED stdout
print "sortBam.pl    Printing STDOUT to file:\n\n$stdout\n\n" if defined $stdout;
open(STDOUT, ">$stdout") or die "sortBam.pl    Can't open STDOUT file: $stdout\n" if defined $stdout;
open(STDERR, ">>$stdout") or die "sortBam.pl    Can't open STDOUT file: $stdout\n" if defined $stdout;

#### DEBUG

#### INSTANTIATE sortBam OBJECT
my $converter = Converter->new(
	{
		inputfiles 	=> $inputfiles,
		outputfiles => $outputfiles,
		clean 		=> $clean,
		label 		=> $label,
		verbose 	=> $verbose,

		#### SAMTOOLS
		samtools	=> $samtools,

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
for ( my $i = 0; $i < @$infiles; $i++ )
{
	my ($infile, $outfile);
	$infile = $$infiles[$i];
	$outfile = $$outfiles[$i] if defined $outputfiles;
	$outfile = $infile if not defined $outputfiles;

	#### SET SORTED FILE STUB
	my $sortedfile_stub = $outfile;
	$sortedfile_stub =~ s/\.bam$//;

	#### SET SORTED FILE STUB AS '.temp' IF INPUT FILE IS SAME AS OUTPUT FILE
	my $tempfile = "$sortedfile_stub.temp";
	$sortedfile_stub = $tempfile if $infile eq $outfile;

	#### SET OUTDIR
	my ($outdir) = $outfile =~ /^(.+?)\/[^\/]+$/;
	File::Path::mkpath($outdir) if not -d $outdir;
	print "sortBam.pl    Could not create output directory: $outdir\n" and exit if not -d $outdir;

	#### SET COMMAND
	my $commands = [];
	my $command = "$samtools/samtools sort $infile $sortedfile_stub";
	push @$commands, "echo '$command'";
	push @$commands, $command;


	#### MOVE .temp FILE TO OUTPUT FILE IF INPUT FILE = OUTPUT FILE
	my $move = "mv -f $tempfile.bam $outfile" ;
	if ( $infile eq $outfile )
	{
		push @$commands, "date" ;
		push @$commands, "echo '$move'";
		push @$commands, $move;
		push @$commands, "echo 'completed move'";
		push @$commands, "date";
	}


	### SET LABEL
	my $joblabel = "sortBam";
	$joblabel = "$label-$i";

	#### SET JOB
	my $job = $converter->setJob( $commands, $label, $outdir);
	push @$jobs, $job;
}

#### RUN JOBS
$converter->runJobs( $jobs, $label);


#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "sortBam.pl    Run time: $runtime\n";
print "sortBam.pl    Completed $0\n";
print "sortBam.pl    ";
print Timer::current_datetime(), "\n";
print "sortBam.pl    ****************************************\n\n\n";
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

