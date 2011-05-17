#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     removeTop

    PURPOSE

        PARSE OUT THE HITS ONLY FROM A SAM FILE

    USAGE

    ./removeTop.pl <--inputfile String> [--outputfile String] [--help]

    --inputfile		:   Full path to input file
    --lines 		:   No. of lines to remove
    --outputfile	:   (optional) Full path to output file
	--topfile		:	Print removed lines to this file (DEFAULT =  <inputfile>.top)
    --help      	:   print help info

    EXAMPLES

chmod 755 /nethome/bioinfo/apps/agua/0.4/bin/apps/removeTop.pl

/nethome/bioinfo/apps/agua/0.5/bin/apps/removeTop.pl \
--inputfile /nethome/bioinfo/apps/agua/0.5/bin/apps/t/Cluster/bowtie/chrY/1/out.sam
--lines 3


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
use Cluster;
use Timer;
use Util;
use Conf::Agua;

##### STORE ARGUMENTS TO PRINT TO FILE LATER
my @arguments = @ARGV;

#### FLUSH BUFFER
$| =1;


#### GET OPTIONS
my $inputfile;
my $lines;
my $outputfile;
my $topfile;
my $help;
if ( not GetOptions (
    'inputfile=s'   => \$inputfile,
    'lines=s'   	=> \$lines,
    'outputfile=s'	=> \$outputfile,
    'topfile=s'   	=> \$topfile,
    'help'          => \$help
) )
{ print "Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "inputfile not defined (Use --help for usage)\n" if not defined $inputfile;

my $runCluster = Cluster->new();
$runCluster->removeTop($inputfile, $lines, $outputfile, $topfile);

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "removeTop.pl    Run time: $runtime\n";
print "removeTop.pl    Completed $0\n";
print "removeTop.pl    ";
print Timer::datetime(), "\n";
print "removeTop.pl    ****************************************\n\n\n";
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


