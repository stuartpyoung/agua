#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     indexSnp

    VERSION         0.01

    PURPOSE

        FAST RETRIEVAL OF READ HIT POSITIONS GIVEN THE READ NAME

    INPUT

        1. COMMA-SEPARATED LIST OF SAM-FORMAT HIT FILES (HITS ONLY)

    OUTPUT

        1. DB INDEX FILE CONTAINING B-TREE INDEX OF INPUT READS

    USAGE

    ./indexSnp.pl <--inputfile String> [--help]

		--inputfiles	:   Comma-separated list of SAM-format hit files
		--dbfile	    :   Print B-tree index to this file
		--help          :   print help info

    EXAMPLES

/nethome/bioinfo/apps/agua/0.4/bin/apps/indexSnp.pl \
--inputfiles /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/autochr22/bowtie/1/chr22/hit-200.sam \
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
use DBIndex::Snp;
use Timer;
use Util;
use Conf::Agua;

#### SET MAQ LOCATION
my $conf = Conf::Agua->new(inputfile=>"$Bin/../../../conf/default.conf");
my $samtools = $conf->getKeyValue("applications", 'SAMTOOLS');
$samtools = "/nethome/apps/ccsngs/apps/samtools/0.1.8";

#### GET OPTIONS
my $inputfiles;
my $dbfile;
my $help;
if ( not GetOptions (

	#### GENERAL
    'inputfiles=s' 	=> \$inputfiles,
    'dbfile=s' 	=> \$dbfile,
    'help' 			=> \$help
) )
{
    print "indexSnp.pl    Use option --help for usage instructions.\n" and exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
print "indexSnp.pl    inputfiles not defined (use --help option)\n" and exit if not defined $inputfiles;
print "indexSnp.pl    dbfile not defined (use --help option)\n" and exit if not defined $dbfile;

my $infiles;
@$infiles = split ",", $inputfiles;
foreach my $infile ( @$infiles )
{
    print "indexSnp.pl    Can't find inputfile: $infile\n" and exit if not -f $infile;
}

my $indexer = DBIndex::Snp->new(
    {
		samtools	=> $samtools,
        inputfiles  =>  $infiles,
        dbfile   	=>  $dbfile
    }
);


my $count = $indexer->buildIndex();
print "indexSnp.pl    Reads indexed: $count\n";

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "indexSnp.pl    Run time: $runtime\n";
print "indexSnp.pl    Completed $0\n";
print "indexSnp.pl    ";
print Timer::datetime(), "\n";
print "indexSnp.pl    ****************************************\n\n\n";
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


