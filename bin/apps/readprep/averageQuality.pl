#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();

=head2

	APPLICATION     averageQuality

    PURPOSE

        1. CALCULATE THE AVERAGE QUALITY VALUES FOR EACH READ IN A FASTQ FILE

			AND COLLECT THE AVERAGE QUALITY FOR AT EACH BASE ACROSS ALL

			READS IN THE FILE

		2. BIN THE VALUES TO GENERATE THE DISTRIBUTION OF AVERAGE QUALITY VALUES

			AND CALCULATE THE TOTAL AVERAGE QUALITY VALUE

    INPUT

		1. SOLEXA-FORMAT OR SANGER-FORMAT FASTQ FILE

    OUTPUT

        1. TAB-SEPARATED *.avqual FILE:

			READ1	<average_quality>
			READ2	<average_quality>
			READ3	<average_quality>
			...

			WHERE 'READ*' IS THE FASTQ IDENTIFIER FOR THE READYY

		2. TAB-SEPARATED FILE *.qualstats.tsv FILE:

			TOTAL	<total_average_quality>
			10	<count_reads_with_average_quality=0-10>
			20	<count_reads_with_average_quality=11-20>
			...
			100	<count_reads_with_average_quality=90+>

    NOTES

    USAGE

    ./averageQuality.pl <--inputfile String> <type String> [-h]

    --inputfile		:   /full/path/to/inputfile
    --compress		:   Optional compression: gzip|zip
    --type			:   Format of FASTQ file: solexa|sanger
    --min			:   Skip read if average quality < min
    --max			:   Skip read if average quality >= max
    --skip			:   Skip first N number of bases
    --length		:   Only count this length of bases (after skip if defined)
    --help			:   print help info

    EXAMPLES

/nethome/bioinfo/apps/agua/0.5/bin/apps/readprep/averageQuality.pl \
--inputfile /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/SRX000600/SRR004186_1.fastq.gz \
--inputfile /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/SRX000600/SRR004186_1.avqual \
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

#### INITIALISE SolexaUtil
my $solexa = SolexaUtil->new();

#### EXTERNAL MODULES
use Data::Dumper;
use Term::ANSIColor qw(:constants);
use Getopt::Long;

#### DEFAULTS
my $binarray = [ 0,1,2,3,4,5,6,7,8,9,10,15,20,30,40,50,60,70,80,90,100,200,300 ];

#### GET OPTIONS
use Getopt::Long;
my $inputfile;
my $outputfile;
my $type;
my $min;
my $max;
my $skip;
my $length;
my $clean;
my $compress;
my $bins;
my $dot = 1000000;
my $help;
GetOptions (
    'inputfile=s' 	=> \$inputfile,
    'outputfile=s' 	=> \$outputfile,
    'type=s' 		=> \$type,
    'clean' 		=> \$clean,
    'min=i' 		=> \$min,
    'max=i' 		=> \$max,
    'skip=i' 		=> \$skip,
    'length=i' 		=> \$length,
    'compress=s' 	=> \$compress,
    'bins=s' 		=> \$bins,
    'dot=i' 		=> \$dot,
    'help' 			=> \$help             
) or die "No options specified. Use --help for usage\n";
usage() if defined $help;

#### CHECK INPUTS
die "inputfile not defined (use --help for usage)\n" if not defined $inputfile;
die "type not defined\n" if not defined $type;
die "type must be 'sanger' or 'solexa'\n" if defined $type and $type !~ /^(sanger|solexa)$/;
die "compress must be 'gzip' or 'zip'\n" if defined $compress and $compress !~ /^(gzip|zip)$/;

#### SET OUTPUT FILE
if ( not defined $outputfile )
{
	$outputfile = $solexa->avqualfile($inputfile, $min, $max, $skip, $length);
}

#### SET BINARRAY IF NOT DEFAULT
@$binarray = split ",", $bins if defined $bins and $bins;

#### 1. CALCULATE AVERAGE QUALITIES IF FILE DOES NOT EXIST
#### OR clean IS SPECIFIED
my $outfile = $outputfile;
$outfile .= ".gz" if defined $compress and $compress eq "gzip";
$outfile .= ".zip" if defined $compress and $compress eq "zip";

average_qualities($inputfile, $outputfile, $min, $max, $skip, $length, $dot)
	if not -f $outfile or defined $clean;

#### COMPRESS IF DEFINED
compress($outputfile, $compress) if defined $compress;
$outputfile .= ".gz" if defined $compress and $compress eq "gzip";
$outputfile .= ".zip" if defined $compress and $compress eq "zip";

#### 2. BIN AVERAGE QUALITIES AND PRINT TO STATS FILE
print "averageQuality.pl    Doing bins: @$binarray\n";
bin_qualities($outputfile, $binarray, $dot);


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

	SUBROUTINE		bin_qualities

	PURPOSE

		CALCULATE QUALITY DISTRIBUTION STATS FROM *.avqual FILE

=cut

sub bin_qualities
{
	my $outputfile	=	shift;
	my $binarray	=	shift;
	my $dot			=	shift;


	#### CREATE BinData OBJECT
	my $binner = BinData->new(
		{
			'BINS'	=>  $binarray
		}
	);

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

	my $counter = 0;
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
	}

	#### CALCULATE TOTAL AVERAGE QUALITY
	my $total_average = $total / $counter;

	#### PRINT AVERAGE QUALITY STATISTICS (TOTAL AVERAGE QUALITY, BINS)
	print "averageQuality.pl    Printing statsfile\n";
	my $statsfile = $outputfile;
	$statsfile =~ s/\.gz//;
	$statsfile =~ s/\.avqual//;
	$statsfile .= ".qualstats";

	#### OPEN OUTPUT FILE
	open(STATSFILE, ">$statsfile") or die "Can't open output file: $statsfile\n";
	print STATSFILE "total reads\t$counter\n";
	print STATSFILE "total avg quality\t$total_average\n";

	my $data_bins = $binner->get_data_bins();
	my $bins_array = $binner->get_bins();
	print STATSFILE "bins\t@$bins_array\n";

	for ( my $i = 0; $i < @$data_bins; $i++ )
	{
		my $binstop = $$bins_array[$i];
		print STATSFILE "$binstop\t$$data_bins[$i]\n";
		print "$binstop\t$$data_bins[$i]\n";
	}
	close(FILE);
	close(STATSFILE);

	#### REPORT OUTPUT FILES PRINTED
	print "averageQuality.pl::average_qualities    statsfile printed:\n\n$statsfile\n\n";
}


=head2

	SUBROUTINE		average_qualities

	PURPOSE

		1. CALCULATE AVERAGE QUALITIES FOR ALL READS AND PRINT

			TO TSV OUTPUT FILE

		2. COLLECT THE AVERAGE QUALITY FOR AT EACH BASE ACROSS

			ALL READS IN THE FILE

=cut

sub average_qualities
{
	my $inputfile	=	shift;
	my $outputfile	=	shift;
	my $min			=	shift;
	my $skip		=	shift;
	my $length		=	shift;
	my $dot			=	shift;


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

	#### COLLECT BASE QUALITIES
	my @bases = ();

	#### READ THROUGH ONE SEQUENCE RECORD AT A TIME
	my $record;
	$/ = "\n";
	my $counter = 0;
	while ( defined ($record = <FILE>) and defined $record )
	{    
		#### GET FOUR LINES
		$record .= <FILE>;
		$record .= <FILE>;
		$record .= <FILE>;    

		print "averageQuality.pl::average_qualities         $counter\n"if defined $dot and $counter % $dot == 0;


		my ($fasta, $quality) = $solexa->fastq2fasta($record);
		if ( not defined $fasta )
		{
			print "averageQuality.pl::average_qualities         FASTA not defined for record: $record\n";
			$counter++;
			next;
		}

		my ($header) = $record =~ /^\@(.+)\n/;
		$header =~ s/\s+length=.+$//;

		my ($symbolic) = $quality =~ /(\S+)\s*$/;	
		my $numeric = $solexa->symbolic2numeric($symbolic);

		my @numbers = split " ", $numeric;

		if ( defined $skip )
		{
			for ( my $i = 0; $i < $skip; $i++ )
			{
				my $discarded = pop(@numbers);
			}
		}

		if ( defined $length )
		{
			#### SKIP IF READ IS SHORTER THAN LENGTH
			next if ($#numbers + 1) < $length;

			while ( $#numbers >= $length )
			{
				my $discarded = shift(@numbers);
			}
		}

		#### POPULATE ARRAY OF AVERAGE BASE QUALITIES IF ITS THE FIRST TIME
		if ( $#bases == -1 )
		{
			for ( my $i = 0; $i < $#numbers + 1; $i++ )
			{
				$bases[$i] = $numbers[$i];
			}
		}

		#### OR ADD QUALITIES IF BASE QUALITIES ARRAY ALREADY POPULATED
		else
		{
			for ( my $i = 0; $i < $#numbers + 1; $i++ )
			{
				$bases[$i] += $numbers[$i] if defined $numbers[$i];
			}
		}

		my $average = average(\@numbers);	

		#### DO MIN OR MAX	
		next if defined $min and $min and $average < $min;
		next if defined $max and $max and $average >= $max;

		print OUTFILE "$header\t$average\n";

		#### INCREMENT COUNTER
		$counter++;
	}
	close(FILE);
	close(OUTFILE);

	#### PRINT AVERAGE BASE QUALITIES
	my $basequalfile = $outputfile;
	$basequalfile =~ s/\.avqual/.basequal/;
	open(BASEQUAL, ">$basequalfile") or die "Can't open basequalfile: $basequalfile\n";
	print "total reads\t$counter\n";
	print BASEQUAL "total reads\t$counter\n";
	for ( my $i = 0; $i < $#bases + 1; $i++ )
	{
		$bases[$i] = $bases[$i] / $counter;
		print "$i\t$bases[$i]\n";
		print BASEQUAL "$i\t$bases[$i]\n";
	}
	close(BASEQUAL) or die "Can't close basequalfile: $basequalfile\n";	

	#### REPORT OUTPUT FILES PRINTED
	print "averageQuality.pl::average_qualities    outputfile printed:\n\n$outputfile\n\n";
	print "averageQuality.pl::average_qualities    basequalfile printed:\n\n$basequalfile\n\n";
}

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
