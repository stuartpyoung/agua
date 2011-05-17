#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();

=head2

	APPLICATION     fastqQualityHeader

    PURPOSE

		1. ADD QUALITY HEADERS (FROM SEQUENCE HEADERS) TO FASTQ-FORMAT READS

    INPUTS

        1. FASTQ FILE WITH EMPTY QUALITY HEADERS

    OUTPUTS

		1. FASTQ FILE WITH NON-EMPTY QUALITY HEADERS

    USAGE

    ./fastqQualityHeader.pl <--inputfile String> <--outputfile String> <--rejectfile String> <paired Boolean> [-h]

    --inputfile		:   /full/path/to/inputfile with FASTQ reads lacking quality headers
    --outputfile	:   /full/path/to/outputfile with quality headers added
	--dot			:   Print progress count per multiple of this integer
    --help			:   print help info

    EXAMPLES

/nethome/bioinfo/apps/agua/0.5/bin/apps/readprep/fastqQualityHeader.pl \
--inputfile /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/SRX001540/SRR005735_2.noQualityHeader.fastq.gz \
--outputfile /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/SRX001540/SRR005735_2.fastq.gz

=cut

use strict;


#### EXTERNAL MODULES
use FindBin qw($Bin);

#### USE LIBRARY
use lib "$Bin/../../../lib";	
use lib "$Bin/../../../lib/external";	

#### INTERNAL MODULES
use BinData;
use Timer;
use SolexaUtil;
use Util;
use Conf::Agua;

#### INITIALISE SolexaUtil OBJECT
my $solexa = SolexaUtil->new();

#### EXTERNAL MODULES
use Data::Dumper;
use Term::ANSIColor qw(:constants);
use Getopt::Long;

#### GET OPTIONS
use Getopt::Long;
my $inputfile;
my $outputfile;
my $paired;
my $dot = 1000000;
my $help;
GetOptions (
    'inputfile=s' 	=> \$inputfile,
    'outputfile=s' 	=> \$outputfile,
    'paired' 		=> \$paired,
    'dot=i' 		=> \$dot,
    'help' 			=> \$help             
) or die "No options specified. Use --help for usage\n";
usage() if defined $help;

#### CHECK INPUTS
die "inputfile is not defined (use --help for options)\n" if not defined $inputfile;
die "outputfile is not defined (use --help for options)\n" if not defined $outputfile;

#### OPEN INPUT FILE
if( $inputfile =~ /\.gz$/ or $inputfile =~ /\.zip$/ )
{
	my $pipe_command = "zcat $inputfile |";
	open(FILE, $pipe_command) or die "Can't open inputfile: $inputfile\n"
}
else
{
	open(FILE, $inputfile) or die "Can't open input file: $inputfile\n";	
}


#### OPEN OUTFILE
open(OUTFILE, ">$outputfile") or die "Can't open outputfile: $outputfile\n";

$/ = "\n";
my $counter = 0;
while ( <FILE> )
{    
	print "$counter\n" if $counter % $dot == 0;
	$counter++;

	#### GET FOUR LINES
	my $sequence_header = $_;
	my $sequence = <FILE>;
	<FILE>;
	my $quality = <FILE>;
	my $quality_header = $sequence_header;
	$quality_header =~ s/^@/+/;

	print "sequence_header not defined \n" and last if not defined $sequence_header;
	print "sequence not defined \n" and last if not defined $sequence;
	print "quality_header not defined \n" and last if not defined $quality_header;
	print "quality_header not defined \n" and last if not defined $quality_header;


	print OUTFILE "$sequence_header$sequence$quality_header$quality";
}
close(FILE) or die "Can't close inputfile: $inputfile\n" if $inputfile !~ /\.zip$/ and $inputfile !~ /\.gz$/;
close(OUTFILE) or die "Can't close outputfile: $outputfile\n";

#### REPORT COMPLETED
print "fastqQualityHeader.pl    outputfile printed:\n\n$outputfile\n\n";

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "\nRun time: $runtime\n";
print "Completed $0\n";
print Util::datetime(), "\n";
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
