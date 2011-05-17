#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();

=head2

	APPLICATION     printListFile

    PURPOSE

        1. SAMPLE READS FROM A LIBRARY 

    INPUT

		1. FASTA OR FASTA.GZ FILE

    OUTPUT

		OPTION: total

			1. TOTAL READS IN ALL READ FILES IN INPUT DIRECTORY

		OPTION: sample

			1. size NUMBER OF FASTA FILES, RECORDS RANDOMLY SELECTED FROM

				FILES IN THE INPUT DIRECTORY

    USAGE

    ./printListFile.pl <--inputfile String> <--mode String> [-h]

    --inputfile     :   /full/path/to/inputfile
	--outputdir		:	/full/path/to/outputdir
	--fraction		:	Files in inputfile directory / total files
	--size			:	Desired no. reads per file
    --mode          :   Functional mode (total, samples)
    --compress		:   Input files are compressed (gzip|zip)
    --help			:   print help info

    EXAMPLES

/nethome/syoung/base/bin/comparison/printListFile.pl \
--inputfile /mihg/data/NGS/syoung/base/pipeline/SRA/test/fasta,/mihg/data/NGS/syoung/base/pipeline/SRA/test2/fasta \
--size 10000


/nethome/syoung/base/bin/comparison/printListFile.pl \
--inputfile /mihg/data/NGS/syoung/base/pipeline/SRA/test/fasta \
--size 10000


/nethome/bioinfo/apps/agua/0.5/bin/apps/readprep/printListFile.pl \
--inputfile /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/SRX000600/fasta,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/SRX000601/fasta,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/SRX000602/fasta,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/SRX000603/fasta,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/SRX001539/fasta,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/SRX001540/fasta \
--fraction 
--size 100000000 \
--compress gzip


=cut

use strict;

#### FLUSH BUFFER
$| = 1;

#### EXTERNAL MODULES
use Data::Dumper;
use File::Path;
use FindBin qw($Bin);
use Getopt::Long;
use Term::ANSIColor qw(:constants);

#### USE LIBRARY
use lib "$Bin/../../../lib";	
use lib "$Bin/../../../lib/external";	

#### INTERNAL MODULES
use Sampler;
use Timer;
use Util;
use Conf::Agua;

#### GET OPTIONS
use Getopt::Long;
my $inputfile;
my $size;
my $compress;
my $paired;
my $fraction;
my $total;
my $dot = 1000000;
my $help;
GetOptions (
    'inputfile=s' => \$inputfile,
    'size=i' => \$size,
    'compress=s' => \$compress,
    'paired' => \$paired,
    'fraction=s' => \$fraction,
    'total=s' => \$total,
    'dot=i' => \$dot,
    'help' => \$help             
) or die "No modes specified. Use --help for usage\n";
usage() if defined $help;


#### CHECK INPUTS
die "Input directory not defined (use --help for usage)\n" if not defined $inputfile;
die "Size not defined (use --help for usage)\n" if not defined $size;

#### SET ARGS
my $args =
{
	'inputfile'	=>	$inputfile,
	'compress'	=>	$compress,
	'size'		=>	$size,
	'paired'	=>	$paired,
	'dot'		=>	$dot,
    'fraction' 	=>	$fraction,
    'total' 	=>	$total,
    'dot' 		=> 	$dot,
};

#### CREATE LIST FILES USING THE INFO FILES AND TOTAL READ COUNTS FOR EACH DIRECTORY
Sampler::printListFile($args);

#### REPORT OUTPUT FILES PRINTED
print "Completed samples\n";	


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
	print `perldoc $0`;

	exit;
}




