#!/usr/bin/perl -w

#### TIME
my $time = time();

=head2

    APPLICATION     numberReads

    VERSION         0.01

    PURPOSE

        FAST RETRIEVAL OF READ HIT POSITIONS GIVEN THE READ NAME

    INPUT

        1. COMMA-SEPARATED LIST OF SAM-FORMAT HIT FILES (HITS ONLY)

    OUTPUT

        1. DB INDEX FILE CONTAINING B-TREE INDEX OF INPUT READS

    USAGE

    ./numberReads.pl <--inputfile String> [--help]

		--dbfile	:   Use this B-tree index database file
		--help      :   print help info

    EXAMPLES

/nethome/bioinfo/apps/agua/0.5/bin/apps/aligners/numberReads.pl \
--dbfile /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/autochr22/bowtie/1/chr22/hit-200.db

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
my $readname;
my $dbfile;
my $mode;
my $help;
if ( not GetOptions (

	#### GENERAL
    'dbfile=s' 	    => \$dbfile,
    'help' 			=> \$help
) )
{
    print "numberReads.pl    Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
print "numberReads.pl    dbfile not defined (use --help option)\n" and exit if not defined $dbfile;
print "numberReads.pl    Can't find dbfile: $dbfile\n" and exit if not -f $dbfile;

my $indexer = IndexRead->new(
    {
        dbfile   =>  $dbfile
    }
);

my $output = $indexer->numberReads($readname, $mode);
print "$output\n";

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "indexVenn.pl    Run time: $runtime\n";
print "indexVenn.pl    Completed $0\n";
print "indexVenn.pl    ";
print Timer::datetime(), "\n";
print "indexVenn.pl    ****************************************\n\n\n";
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
