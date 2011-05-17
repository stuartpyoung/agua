#!/usr/bin/perl -w

=head2

    APPLICATION     lookupRead

    VERSION         0.01

    PURPOSE

        FAST RETRIEVAL OF READ HIT POSITIONS GIVEN THE READ NAME

    INPUT

        1. COMMA-SEPARATED LIST OF SAM-FORMAT HIT FILES (HITS ONLY)

    OUTPUT

        1. DB INDEX FILE CONTAINING B-TREE INDEX OF INPUT READS

    USAGE

    ./lookupRead.pl <--inputfile String> [--help]

		--start	:   Comma-separated list of SAM-format hit files
		--stop    	:   Output stop: 'count', 'array' or 'hash'
		--dbfile	:   Use this B-tree index database file
		--help      :   print help info

    EXAMPLES

/nethome/bioinfo/apps/agua/0.5/bin/apps/aligners/lookupRead.pl \
--dbfile /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/autochr22/bowtie/1/chr22/index.db \
--start "SRR005718.15139750:3:267:564:782#0"

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
my $dbfile;
my $start;
my $stop;
my $help;
if ( not GetOptions (

	#### GENERAL
    'dbfile=s' 	    => \$dbfile,
    'start=i' 	=> \$start,
    'stop=i' 	    => \$stop,
    'help' 			=> \$help
) )
{
    print "lookupRead.pl    Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
print "lookupRead.pl    dbfile not defined (use --help option)\n" and exit if not defined $dbfile;
print "lookupRead.pl    Can't find dbfile: $dbfile\n" and exit if not -f $dbfile;


my $indexer = IndexRead->new(
    {
        dbfile   =>  $dbfile
    }
);

$indexer->getReads($start);


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


