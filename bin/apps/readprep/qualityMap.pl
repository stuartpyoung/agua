#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();

=head2

	APPLICATION     qualityMap

    PURPOSE

        1. CALCULATE THE AVERAGE QUALITY FOR EACH DISTINCT SQUARE PORTION

			OF THE FLOWCELL

    INPUT

		1. *avqual.tsv.gz FILE PRODUCED BY averageQuality.pl

		2. NUMBER OF SQUARE REGIONS TO GENERATE (10, 20, ..., 100)

    OUTPUT

        1. TAB-SEPARATED *.qualmap.tsv FILE:

			total reads		XXXX
			total avg qual	XXXX
				X values (right bound)
			Y (lower bound)
			1		<number_reads>/<average_quality>
			2		<number_reads>/<average_quality>
			...

			WHERE 'READ*' IS THE FASTQ IDENTIFIER FOR THE READYY

		2. TAB-SEPARATED FILE *.qualstats.tsv FILE:

			TOTAL	<total_average_quality>
			10	<count_reads_with_average_quality=0-10>
			20	<count_reads_with_average_quality=11-20>
			...
			100	<count_reads_with_average_quality=90+>

    NOTES

		ILLUMINA READ HEADER FORMAT

		Illumina sequence identifiers
		http://en.wikipedia.org/wiki/FASTQ_format		

		Sequences from the Illumina software use a systematic identifier:

		@HWUSI-EAS100R:6:73:941:1973#0/1

		where

		HWUSI-EAS100R	the unique instrument name
		6				flowcell lane
		73				tile number within the flowcell lane
		941				'x'-coordinate of the cluster within the tile
		1973			'y'-coordinate of the cluster within the tile
		#0				index number for a multiplexed sample (0 for no indexing)
		/1				the member of a pair, /1 or /2 (paired-end or mate-pair reads only)
		Versions of the Illumina pipeline since 1.4 appear to use #NNNNNN instead of #0 for the multiplex ID, where NNNNNN is the sequence of the multiplex tag from the beginning of the read.



    USAGE

    ./qualityMap.pl <--inputfile String> <type String> [-h]

    --inputfile		:   /full/path/to/inputfile
    --compress		:   Optional compression: gzip|zip
    --type			:   Format of FASTQ file: solexa|sanger
    --help			:   print help info

    EXAMPLES

/nethome/bioinfo/apps/agua/0.5/bin/apps/readprep/qualityMap.pl \
--inputfile /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/SRX000600/SRR004186_1.avqual.gz \
--compress gzip \
--bins -50,-40,-30,-20,-10,0,1,2,3,4,5,6,7,8,9,10,15,20,30,40,50,60,70,80,90,100,200,300 \
--type sanger


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

#### FIXED WIDTH SETTINGS
#my $DEFAULT_ID_LENGTH = 20;
#my $DEFAULT_SEQUENCE_LENGTH = 40;

#### GET OPTIONS
use Getopt::Long;
my $inputfile;
my $type;
my $compress;
my $bins;
my $dot = 1000000;
my $help;
GetOptions (
    'inputfile=s' 	=> \$inputfile,
    'type=s' 		=> \$type,
    'compress=s' 	=> \$compress,
    'bins=s' 		=> \$bins,
    'dot=i' 		=> \$dot,
    'help' 			=> \$help             
) or die "No options specified. Use --help for usage\n";
usage() if defined $help;

#### CHECK INPUTS
die "inputfile not defined (use --help for usage)\n" if not defined $inputfile;
die "type not defined\n" if not defined $type;

#### 1. GENERATE AVERAGE QUALITIES FILE
#### OPEN OUTPUT FILE
my $outputfile = $inputfile;
$outputfile =~ s/\.avqual.gz$//;
$outputfile .= ".qualmap";
open(OUTFILE, ">$outputfile") or die "Can't open output file: $outputfile\n";

#### OPEN FILE AND SET RECORD SEPARATOR
if( $inputfile =~ /\.gz$/ or $inputfile =~ /\.zip$/ )
{
    my $pipe_command = "zcat $inputfile |";
    open(FILE, $pipe_command);
}
else
{
    open(FILE, $inputfile) or die "Can't open input file: $inputfile\n";
}


#### REPORT OUTPUT FILES PRINTED
print "qualityMap.pl    Output file printed:\n\n$outputfile\n\n";


#### 2. BIN AVERAGE QUALITIES
print "qualityMap.pl    Doing bins\n";
my $binarray = [ 0,10,20,30,40,50,60,70,80,90,100,500 ];
@$binarray = split ",", $bins if defined $bins and $bins;

#### CREATE BinData OBJECT
#my $min = -500;
#my $max = 500;
my $binner = BinData->new(
    {
        #'MIN'	=>  $min,
        #'MAX'	=>  $max,
        'BINS'	=>  $binarray
    }
);

#my $bins_array = $binner->get_bins();
#exit;

#### OPEN FILE AND SET RECORD SEPARATOR
if( $outputfile =~ /\.gz$/ or $outputfile =~ /\.zip$/ )
{
    my $pipe_command = "zcat $outputfile |";
    open(FILE, $pipe_command);
}
else
{
    open(FILE, $outputfile) or die "Can't open output file: $outputfile\n";
}

$counter = 0;
my $total = 0;
my $avqual;
while ( defined ($avqual = <FILE>) and defined $avqual )
{
	print "$counter\n" if $counter % $dot == 0;

	#$avqual =~ s/\n$//;

	my ($qual) = $avqual =~ /(\S+)\s*$/;

	next if not defined $qual;

	#### BIN AVERAGE QUALITY	
	$binner->add($qual);

	#### ADD TO TOTAL AVERAGE QUALITY
	$total += $qual;

	$counter++;
	#last if $counter >= 1000;
}

$total = $total / $counter;
print "qualityMap.pl    total: $total\n";

#### 3. PRINT AVERAGE QUALITY STATISTICS (TOTAL AVERAGE QUALITY, BINS)
print "qualityMap.pl    Printing statsfile\n";
my $statsfile = $outputfile;
$statsfile =~ s/\.avqual.tsv//;

#### OPEN OUTPUT FILE
open(STATSFILE, ">$statsfile") or die "Can't open output file: $statsfile\n";
print STATSFILE "total reads\t$counter\n";
print STATSFILE "total avg quality\t$total\n";

my $data_bins = $binner->get_data_bins();
my $bins_array = $binner->get_bins();
print STATSFILE "data bins\t$@$data_bins\n";
print "data bins: @$data_bins\n";

for ( my $i = 0; $i < @$data_bins; $i++ )
{
	my $binstop = $$bins_array[$i];
	print STATSFILE "$binstop\t$$data_bins[$i]\n";
	print "$binstop\t$$data_bins[$i]\n";
}

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

	SUBROUTINE		compress

	PURPOSE



=cut

sub compress
{
	my $filename	=	shift;
	my $compress	=	shift;


	my $compressfile = $filename;
	if ( $compress =~ /^gzip$/ )
	{
		$compressfile = $filename . ".gz";
		`rm -fr $compressfile`;
		my $command = "$compress $filename";
		print `$command`;
	}
	elsif ( $compress =~ /^zip$/ )
	{
		$compressfile = $filename . ".zip";
		`rm -fr $compressfile`;
		my $command = "$compress $filename";
		print `$command`;
	}	

	return $compressfile;
}



sub average
{
	my $array	=	shift;

	my $average = 0;
	foreach my $entry ( @$array )
	{
		$average += $entry if defined $entry;
	}

	return $average/scalar(@$array);
}

sub usage
{
	print `perldoc $0`;

	exit;
}
