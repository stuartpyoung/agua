#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();

=head2

	APPLICATION     trimRead

    PURPOSE

		1. TRIM FASTA OR FASTQ READS AT 5' END AND/OR TRUNCATE

			TO SPECIFIED LENGTH

		2. IGNORE READS SHORTER THAN THE SPECIFIED LENGTH

		3. IGNORE READS WITH AVERAGE QUALITY ABOVE OR BELOW

			THE SPECIFIED VALUE

    INPUTS

        1. INPUT FILE

		2. OUTPUT FILE

		3. FILE TYPE (fasta|fastq)

		4. START TRIM LENGTH (E.G., 1 TO REMOVE FIRST READ)

		5. TRUNCATED LENGTH (AFTER START TRIM)

    OUTPUTS

        1. FILE CONTAINING SEQUENCES

	NOTES

		THE SRA TAKES FASTQ DATA IN ILLUMINA FASTQ FORMAT (I.E., OFFSET IS 64):

		http://www.ebi.ac.uk/ena/about/page.php?page=sra_data_format#1.5

		Illumina Genome Analyzer: FASTQ format

		We support Illumina Genome Analyzer submissions in a well defined FASTQ format described below.
		We accept FASTQ submissions with application reads in separate files, i.e. paired reads can be submitted in two files: one file for the forward and one file for the reverse reads. When this is done we recommend that no technical reads, including barcode reads, are submitted to SRA.
		Each spot in the FASTQ file consists of three components: spot name, base calls and quality scores. The spot names are composed of seven fields:
		@<instrument name>:<lane>:<tile>:<x-coordinate>:<y-coordinate>#<tag or index>/<pair index>
		<instrument name>: the unique instrument name
		<lane>: the flowcell lane
		<tile>: the tile within the flowcell lane
		<x-coordinate>: the x-coordinate of the cluster within the tile
		<y-coordinate>: the y-coordinate of the cluster within the tile
		<tag or index>: the tag or the tag index for multiplexed samples (optional or 0 for no indexing)
		<pair index>: the pair index: 1 or 2 (optional or 0 for single-ended spots)
		Example:
		@HWUSI-EAS100R:6:73:941:1973#0/1
		The quality encoding should use a Phred quality score ranging from 0 to 62 using ASCII characters from 64 to 126.
		A single spot in the FASTQ file should look like this:
		@HWUSI-EAS100R:6:73:941:1973#0/1
		GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
		+
		!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
		The FASTQ files should be individually compressed using gzip prior uploading them to EBI.
		More information about the FASTQ format is available from Wikipedia.


    USAGE

    ./trimRead.pl <--inputfile String> <format String> [-h]

    --inputfile		:   /full/path/to/inputfile
    --outputfile	:   /full/path/to/outputfile
    --paired		:   Do both paired ends at once (FASTQ only)
    --compress		:   Optional compression: gzip|zip
    --format		:   Format of file: fastq|fasta
    --type			:   Format of FASTQ file: sanger|solexa
    --min			:   Skip read if average quality < min
    --max			:   Skip read if average quality >= max
    --skip			:   Skip first N number of bases
    --length		:   Only count this length of bases (after skip if defined)
    --dot			:   Print progress count per multiple of this integer
    --help			:   print help info

    EXAMPLES


/nethome/bioinfo/apps/agua/0.5/bin/apps/readprep/trimRead.pl \
--inputfile /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/SRX000600/SRR004865_1.fastq.gz \
--outputfile /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/SRX000600/SRR004865_1.fastq \
--type solexa \
--format fastq \
--length 26 \
--min 3 \
--paired


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
my $format;
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
    'paired' 		=> \$paired,
    'format=s' 		=> \$format,
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
$skip = 0 if not defined $skip;
die "length not defined\n" if not defined $length;
die "format not defined\n" if not defined $format;
die "format must be 'fastq' or 'fasta'\n" if defined $format and $format !~ /^(fastq|fasta)$/;
die "type not defined for 'fastq' format file'\n" if defined $format and not defined $type;
die "type must be 'sanger' or 'solexa'\n" if defined $type and $type !~ /^(sanger|solexa)$/;
die "format must be 'fastq' or 'fasta'\n" if defined $format and $format !~ /^(fastq|fasta)$/;
die "compress must be 'gzip' or 'zip'\n" if defined $compress and $compress !~ /^(gzip|zip)$/;
die "can't use min if format is 'fasta'\n" if $format eq "fasta" and defined $min;
die "can't use max if format is 'fasta'\n" if $format eq "fasta" and defined $max;

#### PRINT INPUTS

#### OPEN FILE AND SET RECORD SEPARATOR
$/ = "\n";
if( $inputfile =~ /\.gz$/ or $inputfile =~ /\.zip$/ )
{
	my $pipe_command = "zcat $inputfile |";
	open(FILE, $pipe_command) or die "Can't open inputfile: $inputfile\n"
}
else
{
	open(FILE, $inputfile) or die "Can't open input file: $inputfile\n";	
}

#### OPEN MATE INPUT FILE IF paired IS DEFINED
my $matefile = $inputfile;
$matefile =~ s/_1\./_2\./;
if ( defined $paired )
{
	#### OPEN FILE AND SET RECORD SEPARATOR
	if( $inputfile =~ /\.gz$/ or $inputfile =~ /\.zip$/ )
	{
		my $pipe_command = "zcat $inputfile |";
		open(MATEFILE, $pipe_command) or die "Can't open matefile: $matefile\n";
	}
	else
	{
		open(MATEFILE, "$matefile") or die "Can't open matefile: $matefile";
	}
}

#### OPEN MATE OUTPUT FILE IF paired IS DEFINED
my $mateoutfile = $outputfile;
$mateoutfile =~ s/_1\./_2\./;
open(MATEOUTFILE, ">$mateoutfile") or die "Can't open mateoutfile: $mateoutfile" if defined $paired;

#### OPEN OUTFILE
open(OUTFILE, ">$outputfile") or die "Can't open outputfile: $outputfile\n";

#### TRIM READS IN FILE (AND MATE IF paired IS DEFINED)
my $counter = 0;
if ( $format eq "fastq" )
{
	while ( <FILE> )
	{
		print "$counter\n" if $counter % $dot == 0;

		my $sequence_header = $_;
		my $sequence = <FILE>;
		my $quality_header = <FILE>;
		my $quality = <FILE>;

		$sequence =~ s/\s+//g;
		$quality =~ s/\s+//g;

		#### REMOVE 'length=...' 
		$sequence_header =~ s/\s+length=.+$//;
		$quality_header =~ s/\s+length=.+$//;


		#### DO MATEFILE IF paired IS DEFINED		
		my $mate_sequence_header;
		my $mate_sequence;
		my $mate_quality_header;
		my $mate_quality;
		if ( defined $paired )
		{
			$mate_sequence_header = <MATEFILE>;
			$mate_sequence = <MATEFILE>;
			$mate_quality_header = <MATEFILE>;
			$mate_quality = <MATEFILE>;

			$mate_sequence =~ s/\s+//g;
			$mate_quality =~ s/\s+//g;

			#### REMOVE 'length=...' 
			$mate_sequence_header =~ s/\s+length=.+$//;
			$mate_quality_header =~ s/\s+length=.+$//;

		}

		#### CHECK LENGTH OF SEQUENCE AND QUALITY
		next if not defined $sequence or not defined $quality;
		next if length($quality) < ($length + $skip);
		next if length($sequence) < ($length + $skip);

		if ( $skip )
		{
			$sequence = substr($sequence, $skip);
			$quality = substr($quality, $skip);
			$mate_sequence = substr($mate_sequence, $skip) if defined $paired;
			$mate_quality = substr($mate_quality, $skip) if defined $paired;
		}
		$sequence = substr($sequence, 0, $length);
		$quality = substr($quality, 0, $length);

		$mate_sequence = substr($mate_sequence, 0, $length) if defined $paired;
		$mate_quality = substr($mate_quality, 0, $length) if defined $paired;

		if ( defined $min or defined $max )
		{
			my $numeric = $solexa->symbolic2numeric($quality, $type);

			my @numbers = split " ", $numeric;

			my $average = average(\@numbers);	

			#### DO MIN OR MAX	
			next if defined $min and $min and $average < $min;
			next if defined $max and $max and $average >= $max;	
		}

		print OUTFILE "$sequence_header$sequence\n$quality_header$quality\n";
		print MATEOUTFILE "$mate_sequence_header$mate_sequence\n$mate_quality_header$mate_quality\n" if defined $paired;

		$counter++;
	}
}
else
{
	while ( <FILE> )
	{
		my $sequence_header = $_;
		my $sequence = <FILE>;


		#### CHOP FIRST CHARACTER OF HEADER
		$sequence_header =~ s/^.//;

		#### CHECK LENGTH OF SEQUENCE AND QUALITY
		next if not defined $sequence;
		next if length($sequence) != $length;

		print OUTFILE "$sequence_header$sequence";
	}
}
close(FILE) or die "Can't close inputfile: $inputfile\n" if $inputfile !~ /\.zip$/ and $inputfile !~ /\.gz$/;
close(MATEFILE) or die "Can't close matefile: $matefile\n" if defined $paired and $matefile !~ /\.zip$/ and $matefile !~ /\.gz$/;
close(OUTFILE) or die "Can't close outputfile: $outputfile\n";
close(MATEOUTFILE) or die "Can't close mateoutfile: $mateoutfile\n" if defined $paired;

#### COMPRESS IF DEFINED
compress($outputfile, $compress) if defined $compress;
$outputfile .= ".gz" if defined $compress and $compress eq "gzip";
$outputfile .= ".zip" if defined $compress and $compress eq "zip";

#### COMPRESS IF DEFINED
compress($mateoutfile, $compress) if defined $compress and defined $paired;
$mateoutfile .= ".gz" if defined $compress and $compress eq "gzip" and defined $paired;
$mateoutfile .= ".zip" if defined $compress and $compress eq "zip" and defined $paired;

#### REPORT COMPLETED
print "trimRead.pl    outputfile printed:\n\n$outputfile\n\n";
print "trimRead.pl    mateoutfile printed:\n\n$mateoutfile\n\n";

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


=head2

	SUBROUTINE		compress

	PURPOSE

		RETURN AVERAGE OF INPUT ARRAY

=cut


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
