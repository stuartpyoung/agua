#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     createIntervals

    PURPOSE

		Generate an interval file conforming to the GATK RealignerTargetCreator

		output file, suitable for input into the GATK IndelRealigner

	VERSION		0.01

	HISTORY

		0.01 BASIC INTERVAL CALCULATION WITH WINDOW

    INPUTS

        1. SORTED (BY CHROMOSOMAL POSITION) *snp FILE

        2. OUTPUT FILE LOCATION

    OUTPUTS

        1. OUTPUT FILE WITH THE FOLLOWING FORMAT:

				chr22:16050994-16051000
				chr22:16051217-16051222
				chr22:16051270-16051277
				chr22:16051546-16051555
				chr22:16051609-16051617
				chr22:16051662-16051724
				...

    USAGE

    ./createIntervals.pl <--inputfile String>  <--outputfile String> [--help]

    --inputfile	:   Input *bam file location
    --outputfile:   Comma-separated output file
    --window	:   Include <window> bases to left and right of SNP/indel
					(Default: 5)
    --help      :   Print help info

    EXAMPLES

perl /nethome/bioinfo/apps/agua/0.5/bin/apps/snp/createIntervals.pl \
--inputfile /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/maq/2/chr22/hit.bam \
--outputfile /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/maq/2/chr22/hit.bam.range


=cut

use strict;

#### EXTERNAL MODULES
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use Getopt::Long;
use FindBin qw($Bin);
use File::Path;
use File::Copy;

#### USE LIBRARY
use lib "$Bin/../../../lib";	
use lib "$Bin/../../../lib/external";	

#### INTERNAL MODULES
use SNP;
use Timer;
use Util;
use Conf::Agua;

##### STORE ARGUMENTS TO PRINT TO FILE LATER
my @arguments = @ARGV;

#### FLUSH BUFFER
$| =1;

#### GET CONF 
my $conf = Conf::Agua->new(inputfile=>"$Bin/../../../conf/default.conf");
my $samtools = $conf->getKeyValue("applications", 'SAMTOOLS');

#### GET OPTIONS
# SPECIFIC
my $inputfile;
my $outputfile;
my $window = 5;
my $help;
if ( not GetOptions (

	#### GENERAL
    'inputfile=s'	=> \$inputfile,
    'outputfile=s'	=> \$outputfile,
    'window=s'		=> \$window,
    'samtools=s'	=> \$samtools,
    'help'			=> \$help

) )
{ print "Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "inputfile not defined (option --help for usage)\n" if not defined $inputfile;
die "outputfile not defined (Use --help for usage)\n" if not defined $outputfile;

#### DEBUG
print "createIntervals.pl    inputfile: $inputfile\n";
print "createIntervals.pl    outputfile: $outputfile\n";

#### RETRIEVE COMMAND
my $command = "$0 @arguments";

#### INSTANTIATE 
my $snp = SNP->new();

#### GET HIT RANGE
print "createIntervals.pl    Doing createIntervals()\n";
$snp->createIntervals($inputfile, $outputfile, $window);

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "createIntervals.pl    Run time: $runtime\n";
print "createIntervals.pl    Completed $0\n";
print "createIntervals.pl    ";
print Timer::datetime(), "\n";
print "createIntervals.pl    ****************************************\n\n\n";
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


