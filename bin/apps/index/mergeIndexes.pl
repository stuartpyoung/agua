#!/usr/bin/perl -w

#### DEBUG


#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     mergeIndexes

    VERSION         0.01

    PURPOSE

        FAST RETRIEVAL OF READ HIT POSITIONS GIVEN THE READ NAME

    INPUT

        1. COMMA-SEPARATED LIST OF SAM-FORMAT HIT FILES (HITS ONLY)

    OUTPUT

        1. DB INDEX FILE CONTAINING B-TREE INDEX OF INPUT READS

    USAGE

    ./mergeIndexes.pl <--inputfile String> [--help]

		--inputfiles	:   Comma-separated list of SAM-format hit files
		--outputfile	    :   Print B-tree index to this file
		--help          :   print help info

    EXAMPLES

/nethome/bioinfo/apps/agua/0.4/bin/apps/mergeIndexes.pl \
--inputfiles /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/autochr22/bowtie/1/chr22/hit-200.sam \
--outputfile /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/autochr22/bowtie/1/chr22/index.db

=cut


use strict;

#### EXTERNAL MODULES
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use Getopt::Long;
use File::Path;
use FindBin qw($Bin);

#### USE LIBRARY
use lib "$Bin/../../../lib";	

#### INTERNAL MODULES
use IndexRead;
use Timer;
use Util;
use Conf::Agua;

#### SET MAQ LOCATION
my $conf = Conf::Agua->new(inputfile=>"$Bin/../../../conf/default.conf");
my $samtools = $conf->getKeyValue("applications", 'SAMTOOLS');
$samtools = "/nethome/apps/ccsngs/apps/samtools/0.1.8";

#### GET OPTIONS
my $inputfiles;
my $outputfile;
my $help;
if ( not GetOptions (

	#### GENERAL
    'inputfiles=s' 	=> \$inputfiles,
    'outputfile=s' 	=> \$outputfile,
    'help' 			=> \$help
) )
{
    print "mergeIndexes.pl    Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
print "mergeIndexes.pl    inputfiles not defined (use --help option)\n" and exit if not defined $inputfiles;
print "mergeIndexes.pl    outputfile not defined (use --help option)\n" and exit if not defined $outputfile;

my $infiles;
@$infiles = split ",", $inputfiles;
foreach my $infile ( @$infiles )
{
    print "mergeIndexes.pl    Can't find inputfile: $infile\n" and exit if not -f $infile;
}

my $indexer = IndexRead->new(
    {
        inputfiles  =>  $infiles,
        outputfile   =>  $outputfile
    }
);

$indexer->mergeIndexes();

print "mergeIndexes.pl    outputfile printed:\n\n$outputfile\n\n";


#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "mergeIndexes.pl    Run time: $runtime\n";
print "mergeIndexes.pl    Completed $0\n";
print "mergeIndexes.pl    ";
print Timer::datetime(), "\n";
print "mergeIndexes.pl    ****************************************\n\n\n";
exit;
