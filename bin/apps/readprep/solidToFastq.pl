#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();

=head2

	APPLICATION     solidToFastq

	VERSION		0.02

	HISTORY
				0.02 FIXED COLORSPACE TO BASESPACE MAPPING

				0.01 BASIC MAPPING WITH QUALITY CONVERTION

    PURPOSE

		1. CONVERT csfasta AND qual TO fastq

		2. SKIP ANY CSFASTA LINES THAT CONTAIN '.' NUMBERS (I.E., 'N' BASES)

    INPUTS

        1. INPUT FILE

    OUTPUTS

		1. OUTPUT FILE OF CORRECTLY FORMATTED READS

		2. REJECTS FILE OF READS THAT WERE NOT ABLE TO BE FORMATTED CORRECTLY

    USAGE

    ./solidToFastq.pl <--inputfile String> <format String> [-h]

    --inputfile		:   /full/path/to/inputfile
    --outputfile	:   /full/path/to/outputfile
    --label     	:   Label, e.g. name of experiment or sample
    --barcode		:   Barcode number (optional)
    --matenumber	:   Mate file number (optional, 1 or 2)
    --dot			:   Print progress count per multiple of this integer
    --help			:   print help info

    EXAMPLES

SINGLE READ
-----------


NO BARCODING

	/nethome/bioinfo/apps/agua/0.5/bin/apps/readprep/solidToFastq.pl \
	--inputfile /FULL/PATH/TO/myfile.csfasta \
	--outputfile /FULL/PATH/TO/myfile.fastq \
	--label myexperiment



PAIRED ENDS
-----------

NO BARCODING, MATE NUMBER IS 1

	/nethome/bioinfo/apps/agua/0.5/bin/apps/readprep/solidToFastq.pl \
	--inputfile /FULL/PATH/TO/myfile.csfasta \
	--outputfile /FULL/PATH/TO/myfile_1.fastq \
	--label myexperiment \
	--matenumber 1

NO BARCODING, MATE NUMBER IS 2 (COMPLEMENTARY TO ABOVE EXAMPLE)

	/nethome/bioinfo/apps/agua/0.5/bin/apps/readprep/solidToFastq.pl \
	--inputfile /FULL/PATH/TO/myfile.csfasta \
	--outputfile /FULL/PATH/TO/myfile_2.fastq \
	--label myexperiment \
	--matenumber 2

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
my $label;
my $barcode = 0;
my $matenumber;
my $dot = 1000000;
my $help;
GetOptions (
    'inputfile=s' 	=> \$inputfile,
    'label=s' 		=> \$label,
    'outputfile=s' 	=> \$outputfile,
    'barcode=i' 	=> \$barcode,
    'matenumber=i' 	=> \$matenumber,
    'dot=i' 		=> \$dot,
    'help' 			=> \$help             
) or die "No options specified. Use --help for usage\n";
usage() if defined $help;


my $mappings = mappings();

#### CHECK INPUTS
die "inputfile is not defined (use --help for options)\n" if not defined $inputfile;
die "label is not defined (use --help for options)\n" if not defined $label;
die "outputfile is not defined (use --help for options)\n" if not defined $outputfile;

my $qualfile = $inputfile;
$qualfile =~ s/\.gz$//;
$qualfile =~ s/\.csfasta$//;
$qualfile .= "_QV.qual";
$qualfile .= ".gz" if $inputfile =~ /\.gz$/;
die "Can't find inputfile: $inputfile\n" if not -f $inputfile;
die "Can't find qualfile: $qualfile\n" if not -f $qualfile;
die "Can't find inputfile: $inputfile\n" if -z $inputfile;
die "Can't find qualfile: $qualfile\n" if -z $qualfile;

#### OPEN QUAL FILE
if( $qualfile =~ /\.gz$/ or $qualfile =~ /\.zip$/ )
{
	my $pipe_command = "zcat $qualfile |";
	open(QUAL, $pipe_command) or die "Can't open qualfile: $qualfile\n"
}
else
{
	open(QUAL, $qualfile) or die "Can't open input file: $qualfile\n";	
}

#### OPEN CSFASTA FILE
if( $inputfile =~ /\.gz$/ or $inputfile =~ /\.zip$/ )
{
	my $pipe_command = "zcat $inputfile |";
	open(FILE, $pipe_command) or die "Can't open inputfile: $inputfile\n"
}
else
{
	open(FILE, $inputfile) or die "Can't open input file: $inputfile\n";	
}

#### GET ASCII
my $ascii = $solexa->_ascii();
%$ascii = reverse(%$ascii);

#### OPEN OUTFILE
open(OUTFILE, ">$outputfile") or die "Can't open outputfile: $outputfile\n";
$/ = "\n";
my $counter = 0;
while ( <FILE> )
{
	print "$counter\n" if $counter % $dot == 0;
	$counter++;

	my $line = $_;
	while ( $line =~ /^#/ )
	{
		$line = <FILE>;
		<QUAL>;
	}

	my $sequence_header = $line;
	$sequence_header =~ s/^>//;
	$sequence_header =~ s/\_F3\s*$//;
	$sequence_header =~ s/\_/:/g;
	$sequence_header = "$label:0:$sequence_header#$barcode";
	$sequence_header .= "/$matenumber" if defined $matenumber;
	my $sequence = <FILE>;

	#### CONVERT FROM COLORSPACE TO BASESPACE
	$sequence = solidToFastq($mappings, $sequence);

	my $quality_header = <QUAL>;
	$quality_header =~ s/^>//;
	my $quality = <QUAL>;

	#### SKIP IF 'N' IN SEQUENCE
	next if not defined $sequence;

	my @ascii = split " ", $quality;
	foreach my $ascii ( @ascii )	{	$ascii = 0 if $ascii < 0; $ascii += 33; }
	my $symbolic_quality = pack("C*", @ascii);

	print OUTFILE "\@$sequence_header\n$sequence\n+$sequence_header\n$symbolic_quality\n";
}

close(FILE) or die "Can't close inputfile: $inputfile\n" if $inputfile !~ /\.zip$/ and $inputfile !~ /\.gz$/;
close(QUAL) or die "Can't close qualfile: $qualfile\n" if $qualfile !~ /\.zip$/ and $qualfile !~ /\.gz$/;
close(OUTFILE) or die "Can't close outputfile: $outputfile\n";

#### REPORT COMPLETED
print "solidToFastq.pl    outputfile printed:\n\n$outputfile\n\n";

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


=head2

	SUBROUTINE		solidToFastq

	PURPOSE

		CONVERT A COLORSPACE SEQUENCE INTO BASESPACE

	INPUTS

		1. MAPPING OF NUMBER AND FIRST NUCLEOTIDE TO SECOND NUCLEOTIDE

		2. COLORSPACE SEQUENCE STARTING WITH 'T' OR OTHER BASE

	OUTPUTS

		1. BASE SEQUENCE STARTING WITH THE NEXT BASE AFTER THE 'T'

=cut

sub solidToFastq
{
	my $mappings 	=	shift;
	my $csfasta		=	shift;

	my $bases = "";
	my @array = split "", $csfasta;
	my $first = shift @array; #### I.E., 'T' AT THE BEGINNING OF THE CSFASTA SEQUENCE
	my ($number, $base);
	while ( $#array > 0 )
	{
		$number = shift @array;
		if ( $number eq "." )
		{
			return;
		}
		else
		{
			$base = $mappings->{$number}->{$first};
		}
		$first = $base;
		$bases .= $base;
	}

	return $bases;
}


	#my $sequence = "T120002200022023232222200200113101013";
	#my $expected = "GAAAAGAAAAGAAGCTAGAGAGGGAAACATGGTTGC";


=head2

	SUBROUTINE		mappings

	PURPOSE

		LOAD THE COLORSPACE TO BASESPACE MAPPINGS HASH

=cut

sub mappings
{
	my $mappings = {};
	my @seconds = split " ", <DATA>;
	while ( <DATA> )
	{
		my ($first, @array) = split " ", $_;
		for ( my $i = 0; $i < $#array + 1; $i++ )
		{
			$mappings->{$array[$i]}->{$first} = $seconds[$i];
		}
	}

	return $mappings;
}


sub usage
{
	print `perldoc $0`;

	exit;
}


__DATA__
    A   C   G   T
A   0   1   2   3
C   1   0   3   2
G   2   3   0   1
T   3   2   1   0


__END__

QUALITY TESTING 

# Title: WT_GUNEYPD_20100219_solid0398_WT_57_PD_04PM
>1_16_203_F3
12 27 5 -1 -1 6 30 8 -1 8 -1 23 -1 6 -1 -1 19 -1 -1 5 -1 -1 6 -1 4 16 25 5 -1 8 9 -1 5 -1 4 21 26 -1 -1 -1 -1 -1 -1 -1 9 11 8 -1 -1 4 
>1_16_249_F3
7 28 8 -1 -1 20 4 8 -1 7 -1 6 -1 7 -1 -1 5 -1 -1 7 -1 -1 8 -1 15 15 11 6 -1 6 6 -1 10 -1 6 9 5 -1 -1 -1 -1 -1 -1 -1 7 6 5 -1 -1 7 
>1_16_328_F3
25 29 8 -1 -1 5 11 8 -1 29 -1 7 -1 24 -1 -1 12 -1 -1 24 -1 -1 4 -1 26 26 6 13 -1 21 17 -1 4 -1 4 5 7 -1 -1 -1 -1 -1 -1 -1 11 9 8 -1 -1 6

@solidToFastq:0:1:16:203#0
ACNNACCNTNGNGNNCNNCNNGNCTGANGGNANGTANNNNNNNGTGNNC
+solidToFastq:0:1:16:203#0
<&!!'?)!)!8!'!!4!!&!!'!%1:&!)*!&!%6;!!!!!!!*,)!!%
@solidToFastq:0:1:16:249#0
CCNNGCCNTNGNANNGNNTNNANAGAGNCGNANGAGNNNNNNNGGGNNA
+solidToFastq:0:1:16:249#0
=)!!5%)!(!'!(!!&!!(!!)!00,'!''!+!'*&!!!!!!!('&!!(
@solidToFastq:0:1:16:328#0
TTNNCCCNANCNCNNCNNTNNGNATGCNCGNANAAGNNNNNNNAACNNC
+solidToFastq:0:1:16:328#0
>)!!&,)!>!(!9!!-!!9!!%!;;'.!62!%!%&(!!!!!!!,*)!!'



SEQUENCE TESTING

my $sequence = "T120002200022023232222200200113101013";
my $expected = "GAAAAGAAAAGAAGCTAGAGAGGGAAACATGGTTGC";

