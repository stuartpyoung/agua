#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     indexReads

    VERSION         0.01

    PURPOSE

        Index SAM files by read name

    INPUT

        1. COMMA-SEPARATED LIST OF SAM-FORMAT HIT FILES (HITS ONLY)

    OUTPUT

        1. DB INDEX FILE CONTAINING B-TREE INDEX OF INPUT READS

    USAGE

    ./indexReads.pl <--inputfile String> [--help]

		--samfiles	:   Comma-separated list of SAM-format hit files
		--dbfile	:   Print B-tree index to this file
		--help      :   print help info

    EXAMPLES

/nethome/bioinfo/apps/agua/0.4/bin/apps/indexReads.pl \
--samfiles /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/autochr22/bowtie/1/chr22/hit-200.sam \
--dbfile /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/autochr22/bowtie/1/chr22/index.db

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
my $samfiles;
my $dbfile;
my $help;
if ( not GetOptions (

	#### GENERAL
    'samfiles=s' 	=> \$samfiles,
    'dbfile=s' 	=> \$dbfile,
    'help' 			=> \$help
) )
{
    print "indexReads.pl    Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
print "indexReads.pl    samfiles not defined (use --help option)\n" and exit if not defined $samfiles;
print "indexReads.pl    dbfile not defined (use --help option)\n" and exit if not defined $dbfile;

my $infiles;
@$infiles = split ",", $samfiles;
foreach my $infile ( @$infiles )
{
    print "indexReads.pl    Can't find inputfile: $infile\n" and exit if not -f $infile;
}

my $indexer = IndexRead->new(
    {
		samtools	=> $samtools,
        inputfiles  =>  $infiles,
        dbfile   	=>  $dbfile
    }
);

my $count = $indexer->buildIndex();
print "indexReads.pl    Reads indexed: $count\n";

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "indexReads.pl    Run time: $runtime\n";
print "indexReads.pl    Completed $0\n";
print "indexReads.pl    ";
print Timer::datetime(), "\n";
print "indexReads.pl    ****************************************\n\n\n";
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


