#!/usr/bin/perl -w

#### DEBUG
$DEBUG = 1;

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     samToTsv

    PURPOSE

        EXTRACT STATISTICS FROM A TOPHAT RUN

    INPUT

        1. INPUT DIRECTORIES OF TOPHAT ANALYSES 

        2. SPLITFILE SHOWING LOCATIONS OF INPUT FILES

		3. DIRECTORY CONTAINING REFERENCE SEQUENCE FILES

    OUTPUT

        1. A 'TOPHAT.stats' OUTPUT FILE IN EACH INPUT DIRECTORY

    USAGE

    ./samToTsv.pl <--samfile String> <--matefiles String> <--outputdir String>
        <--tsvfile String> [--help]

    --samfile		:   Location of input .sam file
    --tsvfile	:   Location of output .tsv file
	--help     		:   print this help info

	EXAMPLES

chmod 755 /nethome/bioinfo/apps/agua/0.4/bin/apps/samToTsv.pl

perl /nethome/bioinfo/apps/agua/0.4/bin/apps/samToTsv.pl \
--samfile /scratch/syoung/base/pipeline/bixby/run1/tophat/analysis1/chrY_random/1/accepted_hits.sam


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
use Timer;
use TOPHAT;
use Util;
use Conf::Agua;

##### STORE ARGUMENTS TO PRINT TO FILE LATER
my @arguments = @ARGV;

#### FLUSH BUFFER
$| =1;

#### SET samToTsv LOCATION
my $conf = Conf::Agua->new(inputfile=>"$Bin/../../../conf/default.conf");

#### GET OPTIONS
my $samfile;
my $tsvfile;
my $tempdir;
my $help;
print "samToTsv.pl    Use option --help for usage instructions.\n" and exit if not GetOptions (
    'samfile=s' 	=> \$samfile,
    'tsvfile=s'		=> \$tsvfile,
    'tempdir=s'		=> \$tempdir,
    'help' 			=> \$help
);

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "samfile not defined (Use --help for usage)\n" if not defined $samfile;

if ( not defined $tsvfile )
{
	$tsvfile = $samfile;
	$tsvfile =~ s/\.sam$/.tsv/;
}

#### CONVERT SAM TO TSV
my $tophat = TOPHAT->new();
$tophat->samToTsv($samfile, $tsvfile);

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "samToTsv.pl    Run time: $runtime\n";
print "samToTsv.pl    Completed $0\n";
print "samToTsv.pl    ";
print Timer::current_datetime(), "\n";
print "samToTsv.pl    ****************************************\n\n\n";
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


