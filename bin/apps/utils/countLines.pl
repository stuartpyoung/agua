#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();

=head2

    APPLICATION     countLines

    VERSION         0.01

    PURPOSE

		1. COUNT THE NUMBER OF LINES IN A FILE

        2. PRINT NUMBER TO OUTPUT FILE IF PROVIDED, TO STDOUT OTHERWISE

    INPUT

        1. NON-BINARY FILE

    OUTPUT

        1. OPTIONAL OUTPUT FILE CONTAINING LINE COUNT

    USAGE

    ./countLines.pl <--inputfile String> [--outputfile String] [--help]

		--inputfile		:   Input non-binary file
		--outputfile	:   Print output to this file (default: STDOUT)
		--help          :   print help info

    EXAMPLES

/nethome/bioinfo/apps/agua/0.5/bin/apps/aligners/countLines.pl \
--inputfile /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland2/1/chr22/hit.bam \
--outputfile /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland2/1/chr22/hit.bam.idxstats 


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
use Timer;
use Util;
use Conf::Agua;

#### SET MAQ LOCATION
my $conf = Conf::Agua->new(inputfile=>"$Bin/../../../conf/default.conf");
my $samtools = $conf->getKeyValue("applications", 'SAMTOOLS');

#### GET OPTIONS
my $inputfile;
my $outputfile;
my $help;
if ( not GetOptions (

	#### GENERAL
    'inputfile=s' 	=> \$inputfile,
    'outputfile=s' 	=> \$outputfile,
    'help' 			=> \$help
) )
{ print "countLines.pl    Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
print "countLines.pl    inputfile not defined (use --help option)\n" and exit if not defined $inputfile;
print "countLines.pl    Can't find inputfile: $inputfile\n" and exit if not -f $inputfile;

#### OPEN FILEPUTFILE
open(FILE, $inputfile) or die "Can't open inputfile: $inputfile\n";
$/ = "\n";
my $counter = 0;
while(<FILE>)	{	$counter++;	}
close(FILE);

#### PRINT TO STDOUT IF OUTPUT FILE NOT DEFINED 
if ( not defined $outputfile )
{
	print "$counter\n" and exit;
}
#### OTHERWISE, PRINT TO OUTPUT FILE
else
{
	open(OUT, ">$outputfile") or die "Can't open outputfile: $outputfile\n";
	print OUT "$counter\n"; 
	close(OUT) or die "Can't close outputfile: $outputfile\n";
}

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "\nRun time: $runtime\n";
print "Completed $0\n";
print Timer::datetime(), "\n";
print "****************************************\n\n\n";
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

