#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     

    PURPOSE

		DETERMINE THE RANGE (START AND STOP) OF READ ALIGNMENT POSITIONS

		ON THE REFERENCE SEQUENCE IN A *bam FILE

	VERSION		0.01

	HISTORY

		0.01 BASIC RANGE CALCULATION

    INPUTS

        1. SORTED *bam FILE

        2. OUTPUT FILE LOCATION

    OUTPUTS

        1. TAB-SEPARATED VALUES OUTPUT FILE WITH THE FOLLOWING FORMAT:

			#READS	#START	#STOP		<= HEADER LINE
			343040	1500031	51389018	<= DATA LINE

    USAGE

    ./hitRange.pl <--inputfile String>  <--outputfile String> [--help]

    --inputfile	:   Input *bam file location
    --outputfile:   Comma-separated output file
    --help      :   Print help info

    EXAMPLES

perl /nethome/bioinfo/apps/agua/0.5/bin/apps/binner/hitRange.pl \
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
use UCSCBin;
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
my $help;
if ( not GetOptions (

	#### GENERAL
    'inputfile=s'	=> \$inputfile,
    'outputfile=s'	=> \$outputfile,
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
print "hitRange.pl    inputfile: $inputfile\n";
print "hitRange.pl    outputfile: $outputfile\n";


#### RETRIEVE COMMAND
my $command = "$0 @arguments";

#### INSTANTIATE 
my $binner = UCSCBin->new({	samtools	=>	$samtools	});

#### GET HIT RANGE
print "hitRange.pl    Doing hitRange()\n";
my ($counter, $start, $stop) = $binner->hitRange($inputfile);

#### PRINT TO FILE
open(OUT, ">$outputfile") or die "Can't open outputfile: $outputfile\n";
print OUT "#reads\t#start\t#stop\n";
print OUT "$counter\t$start\t$stop\n";
close(OUT) or die "Can't close outputfile: $outputfile\n";

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "hitRange.pl    Run time: $runtime\n";
print "hitRange.pl    Completed $0\n";
print "hitRange.pl    ";
print Timer::datetime(), "\n";
print "hitRange.pl    ****************************************\n\n\n";
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


