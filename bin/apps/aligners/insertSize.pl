#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     insertSize

    PURPOSE

        CALCULATE THE AVERAGE INSERT SIZE BASED ON ELAND ALIGNMENTS OF PAIRED END READS

    USAGE

    ./insertSize.pl <--inputfile String> [--outputfile String] [--help]

    --inputfile		:   Location of ELAND first mate alignment file ('myTest_1_export.txt')
    --matefile		:   Location of ELAND second mate alignment file ('myTest_2_export.txt')
    --outputfile	:   Location to print insert size report file ('myTest.inserts.txt')
	--lines			:	Number of lines to use to estimate insert size (DEFAULT = 10000)
    --help      	:   print help info

	INPUT

		1. EXTRACT LINES FROM INPUT FASTQ PAIRED END FILES AND RUN ELAND_standalone.pl

		2. RUN insertSize.pl ON OUTPUT *export.txt FILES TO DETERMINE INSERT SIZES

	OUTPUT

		1. PRINT TO OUTPUT FILE:

			-	MIN, MAX AND AVERAGE INSERT SIZE

			-	DISTRIBUTION OF INSERT SIZES

    EXAMPLES

/nethome/bioinfo/apps/agua/0.5/bin/apps/insertSize.pl \
--inputfile /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/rerun/eland/SRX000600/chrY/read_1/reanalysis_export.txt \
--matefile /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/rerun/eland/SRX000600/chrY/read_2/reanalysis_export.txt \
--outputfile /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/rerun/eland/SRX000600/chrY/paired/insertsize.txt \
--lines 10000


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
use BinData;
use Timer;
use Util;
use Conf::Agua;

##### STORE ARGUMENTS TO PRINT TO FILE LATER
my @arguments = @ARGV;

#### FLUSH BUFFER
$| = 1;

#### GET OPTIONS
my $inputfile;
my $matefile;
my $outputfile;
my $lines = 10000;
my $help;
if ( not GetOptions (
    'inputfile=s'	=> \$inputfile,
    'matefile=s'	=> \$matefile,
    'outputfile=s'	=> \$outputfile,
    'lines=i'		=> \$lines,
    'help'			=> \$help
) )
{ print "Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "inputfile not defined (Use --help for usage)\n" if not defined $inputfile;
die "outputfile not defined (Use --help for usage)\n" if not defined $outputfile;

#### SET MATE FILE
$matefile = $inputfile if not defined $matefile;
$matefile =~ s/_1_export/_2_export/;
die "Can't set matefile correctly: $matefile" if $matefile eq $inputfile;

#### CHECK FILES PRESENT AND NON-EMPTY
die "Can't find inputfile: $inputfile\n" if not -f $inputfile;
die "Can't find matefile: $matefile\n" if not -f $matefile;
die "inputfile is empty: $inputfile\n" if -z $inputfile;
die "matefile is empty: $matefile\n" if -z $matefile;

#### OPEN FILES
open(FILE, $inputfile) or die "Can't open input file: $inputfile\n";
open(MATE, $matefile) or die "Can't open mate file: $matefile\n";

#### OPEN OUTFILE IF DEFINED
open(OUTFILE, ">$outputfile") or die "Can't open output file: $outputfile\n" if defined $outputfile;

#### LOAD MATE HITS INTO MEMORY
$/ = "\n";
my $counter = 0;
my $line = <MATE>;
my $mates = {};
while ( defined $line and $counter < ($lines * 10) )
{
	$counter++;
	#last if $counter >= 100;

	#### STORE ALIGNMENT POSITION IF MATCHED REFERENCE
	$line =~ /^(\S+\s*\t){10}\t(\S+)\t.+\t\d+/;
	my $match = $1;
	if ( defined $match and $match !~ /NM\s*/ )
	{
		my @array = split "\t", $line;
		my $id = $array[0];
		my $chromosome = $array[10];
		my $machine = $array[2];
		my $lane = $array[3];
		my $x = $array[4];
		my $y = $array[5];

		my $position = $array[12];

		$mates->{"$id$machine$lane$x$y"} = {
			chromosome => $chromosome,
			position => $position
		};
	}
	$line = <MATE>;	
}
close(MATE);

my @keys = keys %$mates;
print "insertSize.pl    Number of mates: ", $#keys + 1, "\n";

#### SEARCH FOR MATE FOR EACH READ
my $readline = <FILE>;
my $reads = 0;

print "insertSize.pl    Doing inputfile: $inputfile\n";
my $distances = [];
while ( defined $readline and $reads < $lines )
{
	$reads++;
	#last if $reads >= 10;

	#### STORE ALIGNMENT POSITION IF MATCHED REFERENCE
	$readline =~ /^(\S+\s*\t){10}\t(\S+)\t.+\t\d+/;
	my $match = $1;
	if ( defined $match and $match !~ /NM\s*/ )
	{
		my @array = split "\t", $readline;
		my $id = $array[0];
		my $chromosome = $array[10];
		my $position = $array[12];
		my $machine = $array[2];
		my $lane = $array[3];
		my $x = $array[4];
		my $y = $array[5];

		my $mate = $mates->{"$id$machine$lane$x$y"};
		if ( defined $mate and $chromosome eq $mate->{chromosome} )
		{
			my $distance = abs($position - $mate->{position});
			push @$distances, $distance;
		}
	}	
	$readline = <FILE>;	

	print "$reads\n" if 1000 % $reads == 0;
	#last if scalar(@$distances) >= 100;
}
close(FILE);
print "insertSize.pl    No. distances: ", scalar(@$distances), "\n";

#### GET MEDIAN AND MEAN
@$distances = sort { $a <=> $b } @$distances;
my $median = $$distances[scalar(@$distances)/2 - 1];
my $sum = 0;
foreach my $distance ( @$distances ) {	$sum += $distance; }
my $mean = $sum/scalar(@$distances);
$mean =~ s/\..+$//;

#### GET UNITS OF MEDIAN
my $units = length($median);

#### CALCULATE SAMPLE STANDARD DEVIATION
#### (USES THE MEDIAN INSTEAD OF THE MEAN)
#### I.E., s = sqrt( 1/(N-1) * sum 1..N [ (distance - median)**2] )
my $deviation = 0;
foreach my $distance ( @$distances )
{
	$deviation += ($distance - $median)**2;
}
$deviation = $deviation / (scalar(@$distances) - 1);
$deviation = $deviation**0.5;



#### ROUND TO TWO UNITS BELOW
my $divisor = 1 . "0" x ($units - 2);
my $center = (int($median/$divisor)) * $divisor;

#### SET BIN ARRAY TO STRADDLE MEDIAN
my $binarray = [];
my $offset = -10 * $divisor;
for ( my $i = 0; $i < 21; $i++ )
{
	$$binarray[$i] = $center + $offset;
	$offset += $divisor;
}

#### BIN AVERAGE QUALITY	
my $binner = BinData->new(
	{
		'BINS'	=>  $binarray
	}
);
$binner->add($distances);

#### PRINT BINNED DATA TO OUTPUT FILE
my $data_bins = $binner->get_data_bins();
my $bins_array = $binner->get_bins();

print OUTFILE "INSERT SIZE STATISTICS\n";
print OUTFILE "----------------------\n";
print OUTFILE "inputfile  : $inputfile\n";
print OUTFILE "matefile   : $matefile\n";
print OUTFILE "\n";
print OUTFILE "median: $median\n";
print OUTFILE "sample standard deviation: $deviation\n";
print OUTFILE "reads: $reads\n";
print OUTFILE "mean: $mean\n";
print OUTFILE "distance\tcount\n";
for ( my $i = 0; $i < @$data_bins; $i++ )
{
	my $binstop = $$bins_array[$i];
	print OUTFILE "$binstop\t$$data_bins[$i]\n";
}
close(OUTFILE);

print "\n";

#### PRINT OUTPUT FILE
print "outputfile printed:\n\n$outputfile\n\n";

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
	print GREEN;
	print `perldoc $0`;
	print RESET;
	exit;
}


