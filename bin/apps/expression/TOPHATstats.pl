#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     TOPHATstats

    PURPOSE

        EXTRACT STATISTICS FROM A TOPHAT RUN

    INPUT

        1. INPUT DIRECTORIES OF TOPHAT ANALYSES 

        2. SPLITFILE SHOWING LOCATIONS OF INPUT FILES

		3. DIRECTORY CONTAINING REFERENCE SEQUENCE FILES

    OUTPUT

        1. A 'TOPHAT.stats' OUTPUT FILE IN EACH INPUT DIRECTORY

    USAGE

    ./TOPHATstats.pl <--inputdirs String> <--matefiles String> <--outputdir String>
        <--referencedir String> [--help]

    --inputdirs		:   Single FASTQ sequence file
    --outputdir		:   Create this directory and write output files to it
    --referencedir	:   Location of squashed genome reference files
    --chromosomes	:   Location of chromosome sizes file
    --random		:   Use this flag to include chr*_random.fa reference hits
    --cluster		:   Use cluster type ('PBS' or 'LSF') (default: run locally)
    --clean			:   Recount hits in all Tophat output files
    --help     		:   print help info


	NOTES

		TOPHAT OUTPUT accepted_hits.sam FILE CONTAINS ONLY READS WITH HITS (NO UNMATCHED READS):

		head /scratch/syoung/base/pipeline/bixby/run1/tophat/analysis13/chr1/accepted_hits.sam

			HWI-EAS185:3:13:1397:1556#0     16      chr1    3005531 0       52M     *       0       0       TGTGCGAGTTCCTCTGGGATTGGGTTACCTCACTCAGGATGATGCCCTCCAG aaaa_a\a`_`aWV`aaba]aabbbaW_^aa]aaaaababbaabZ`a`_aabNM:i:2
			HWI-EAS185:3:29:1842:1863#0     16      chr1    3005656 0       52M     *       0       0       CCACAATTTATGTATCCATTCCTCTGTTGAGGGGCATCTGGGTTCATTCCAG VYUVUTV^RSUVZY]YXX^M[^TY]____Ya`aa_`a`ababaaaa_]]aabNM:i:2


		CHECK FOR LINES WITH NO MATCH:

		grep -n "0   \*      \*" /scratch/syoung/base/pipeline/bixby/run1/tophat/analysis13/chr1/accepted_hits.sam
			<NO RESULT>

	EXAMPLES

chmod 755 /nethome/bioinfo/apps/agua/0.4/bin/apps/TOPHATstats.pl

perl /nethome/bioinfo/apps/agua/0.4/bin/apps/TOPHATstats.pl \
--inputdirs /scratch/syoung/base/pipeline/bixby/run1/tophat/analysis13,/scratch/syoung/base/pipeline/bixby/run1/tophat/analysis14,/scratch/syoung/base/pipeline/bixby/run1/tophat/analysis15,/scratch/syoung/base/pipeline/bixby/run1/tophat/analysis16 \
--referencedir /nethome/bioinfo/data/sequence/chromosomes/mouse/mm9/bowtie \
--chromosomes /nethome/bioinfo/data/sequence/chromosomes/mouse/mm9/bowtie/chromosome-sizes.txt \


perl /nethome/bioinfo/apps/agua/0.4/bin/apps/TOPHATstats.pl \
--inputdirs /scratch/syoung/base/pipeline/bixby/run1/tophat/analysis13 \
--referencedir /nethome/bioinfo/data/sequence/chromosomes/mouse/mm9/bowtie \
--chromosomes /nethome/bioinfo/data/sequence/chromosomes/mouse/mm9/bowtie/chromosome-sizes.txt \



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
use lib "$Bin/../../../lib/external/lib64/perl5/site_perl/5.8.8/x86_64-linux-thread-multi";

#### INTERNAL MODULES
use DBase::SQLite;
use TOPHAT;
use Timer;
use Util;
use Conf::Agua;

##### STORE ARGUMENTS TO PRINT TO FILE LATER
my @arguments = @ARGV;

#### FLUSH BUFFER
$| =1;

#### SET TOPHATstats LOCATION
my $conf = Conf::Agua->new(inputfile=>"$Bin/../../../conf/default.conf");

#### GET OPTIONS
my $inputdirs;
my $referencedir;
my $splitfiles;
my $chromosomes;
my $random;
my $sqlite;
my $clean;

# CLUSTER
my $cluster;
my $maxjobs = 500;
my $cpus = 1;
my $sleep = 5;
my $verbose;
my $tempdir;  #= "/tmp";
my $queue;
my $dot = 1;
my $cleanup;

my $help;
print "TOPHATstats.pl    Use option --help for usage instructions.\n" and exit if not GetOptions (
    'inputdirs=s' 	=> \$inputdirs,
    'referencedir=s'=> \$referencedir,
    'splitfiles=s'	=> \$splitfiles,
    'chromosomes=s'	=> \$chromosomes,
    'random'		=> \$random,
    'sqlite=s'		=> \$sqlite,
    'clean'			=> \$clean,

	#### CLUSTER
    'cluster=s' 	=> \$cluster,
    'maxjobs=i'     => \$maxjobs,
    'cpus=i'        => \$cpus,
    'queue=s'       => \$queue,
    'sleep=s' 		=> \$sleep,
    'verbose' 		=> \$verbose,
    'tempdir=s' 	=> \$tempdir,
    'cleanup'       => \$cleanup,

    'help' 			=> \$help
);

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "inputdirs not defined (Use --help for usage)\n" if not defined $inputdirs;
die "referencedir file not defined (Use --help for usage)\n" if not defined $referencedir;

#### GET INPUT DIRS
my $indirs;
@$indirs = split ",", $inputdirs;
foreach my $indir ( @$indirs )
{
	print "Can't find input directory: $indir\n" and exit if not -d $indir;
}
print "TOPHATstats.pl    inputdirs: $inputdirs\n";
print "TOPHATstats.pl    referencedir: $referencedir\n";

#### SET SPLITFILES IF DEFINED
my $splits = [];
@$splits = split ",", $splitfiles if defined $splitfiles;

#### GET REFERENCE FILES
chdir($referencedir) or die "Can't change to reference directory: $referencedir\n";
my $references;
@$references = <*.fa*>;
print "TOPHATstats.pl    Quitting because no files found in directory: $referencedir\n" and exit if scalar(@$references) == 0;

#### REMOVE 'RANDOM' FROM REFERENCES
if ( not defined $random )
{
	print "TOPHATstats.pl    Removing chr*_random.fa references\n";
	for ( my $i = 0; $i < @$references; $i++ )
	{
		if ( $$references[$i] =~ /random/ )
		{
			splice @$references, $i, 1;
			$i--;
		}
	}
}

#### SORT BY NUMBER
@$references = Util::sort_naturally($references);
foreach my $reference ( @$references )	{	$reference =~ s/\.fa$//;	}

#### DEBUG
@$references = reverse @$references;
print "TOPHATstats.pl    references: @$references\n";

#### SET CLUSTER OBJECT
print "TOPHATstats.pl    SETTING cluster object\n";
my $tophat = TOPHAT->new(
	{
		#### CLUSTER
		cluster 	=> 	$cluster,
		queue		=>	$queue,
		maxjobs     => $maxjobs,
		cpus        => $cpus,
		sleep       => $sleep,
		cleanup     => $cleanup,
		tempdir 	=> $tempdir
	}
);

for ( my $i = 0; $i < @$indirs; $i++ ){
	$tophat->stats($$indirs[$i], $$splits[$i], $references, $clean, $sqlite, $cluster);	
}

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "TOPHATstats.pl    Run time: $runtime\n";
print "TOPHATstats.pl    Completed $0\n";
print "TOPHATstats.pl    ";
print Timer::current_datetime(), "\n";
print "TOPHATstats.pl    ****************************************\n\n\n";
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


