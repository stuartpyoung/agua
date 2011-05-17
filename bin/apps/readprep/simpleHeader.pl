#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();

=head2

	APPLICATION     simpleHeader

    PURPOSE

		1. CONVERT INCORRECT HEADER FORMATS TO THE SIMPLE HEADER FORMAT

			ACCEPTED BY ELAND AND OTHER ALIGNERS MATCHING THIS REGEX:

			[A-Z0-9\.]:[0-9]:[0-9]:[0-9]:[0-9](#[0-9]\/[12])?

			WHERE THE LAST PORTION DIFFERS DEPENDING ON WHETHER THE FILE

			CONTAINS A SINGLE OR PAIRED END READS

    INPUTS

        1. INPUT FILE

    OUTPUTS

		1. OUTPUT FILE OF CORRECTLY FORMATTED READS

		2. REJECTS FILE OF READS THAT WERE NOT ABLE TO BE FORMATTED CORRECTLY

    USAGE

    ./simpleHeader.pl <--inputfile String> <format String> [-h]

    --inputfile		:   /full/path/to/inputfile
    --outputfile	:   /full/path/to/outputfile
    --rejectfile	:   /full/path/to/rejectfile
    --barcode		:   Barcode number (optional)
    --matenumber	:   Mate file number (1 or 2)
    --dot			:   Print progress count per multiple of this integer
    --help			:   print help info

    EXAMPLES

/nethome/bioinfo/apps/agua/0.5/bin/apps/readprep/simpleHeader.pl \
--inputfile /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/samples/100M/min3length26-4.reads_1.fastq \
--outputfile /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/samples/100M/4/min3length26-4.reads_1.sequence.txt \
--rejectfile /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/samples/100M/4/min3length26-4.reads_1.rejects.txt \
--matenumber 1

	Run time: 00:09:44
	Completed /nethome/bioinfo/apps/agua/0.5/bin/apps/readprep/simpleHeader.pl
	8:27AM, 4 September 2010
	****************************************

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
my $rejectfile;
my $barcode;
my $matenumber;
my $dot = 1000000;
my $help;
GetOptions (
    'inputfile=s' 	=> \$inputfile,
    'outputfile=s' 	=> \$outputfile,
    'rejectfile=s' 	=> \$rejectfile,
    'barcode=i' 	=> \$barcode,
    'matenumber=i' 	=> \$matenumber,
    'dot=i' 		=> \$dot,
    'help' 			=> \$help             
) or die "No options specified. Use --help for usage\n";
usage() if defined $help;

#### CHECK INPUTS
die "inputfile is not defined (use --help for options)\n" if not defined $inputfile;
die "outputfile is not defined (use --help for options)\n" if not defined $outputfile;
die "rejectfile is not defined (use --help for options)\n" if not defined $rejectfile;
die "matenumber must be 1 or 2 (use --help for options)\n" if defined $matenumber and $matenumber !~ /^(1|2)$/;

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

#### OPEN REJECT FILE
open(REJECTFILE, ">$rejectfile") or die "Can't open rejectfile: $rejectfile";

#### OPEN OUTFILE
open(OUTFILE, ">$outputfile") or die "Can't open outputfile: $outputfile\n";
$/ = "\n";
while ( <FILE> )
{
    if ( $_ =~ /^(@[A-Z0-9\_\-\.]+)\s+\S*([0-9]+:[0-9]+:[0-9]+:[0-9]+)\s*/i 
		or $_ =~ /^(@[A-Z0-9\_\-\.]+)\s+\S*([0-9]+:[0-9]+:[0-9]+:[0-9]+#[0-9]+\/(1|2))\s*/i )
	{
		my $sequence_header = "$1:$2";
		my $sequence = <FILE>;
		my $quality_header = <FILE>;
		my $quality = <FILE>;
		if ( $quality_header =~  /^(\+[A-Z0-9\_\-\.]+)\s+\S*([0-9]+:[0-9]+:[0-9]+:[0-9]+)\s*/i 
		or $quality_header =~ /^(\+[A-Z0-9\_\-\.]+)\s+\S*([0-9]+:[0-9]+:[0-9]+:[0-9]+#[0-9]+\/(1|2))\s*/i )
		{
			my $new_quality_header = "$1:$2";
			if ( defined $matenumber )
			{
				#### SET BARCODE = 0 IF NOT DEFINED
				$barcode = 0 if not defined $barcode;

				#### ADD BARCODE AND MATE NUMBER
				$sequence_header .= "#$barcode/$matenumber";
				$new_quality_header .= "#$barcode/$matenumber";
			}
			elsif ( defined $barcode )
			{
				#### ADD BARCODE
				$sequence_header .= "#$barcode";
				$new_quality_header .= "#$barcode";
			}
			print OUTFILE "$sequence_header\n$sequence$new_quality_header\n$quality";
		}
		else
		{
			my $rest = <FILE>;
			$rest .= <FILE>;
			$rest .= <FILE>;
			print REJECTFILE $_;
			print REJECTFILE $rest;
		}
    }


	#### @HWI-EAS185_6_20GVYAAXX_Jia_cDNA2_mtDNA_Total_Small_Medium_read2_JH:6:1:589:870/1

    elsif ( $_ =~ /^(@[A-Z0-9\_\-\.]+):([0-9]+:[0-9]+:[0-9]+:[0-9]+)/i 
		or $_ =~ /^(@[A-Z0-9\_\-\.]+):([0-9]+:[0-9]+:[0-9]+:[0-9]+#[0-9]+\/(1|2))\s*/i )
	{
		my $sequence_header = "$1:$2";

		my $sequence = <FILE>;
		my $quality_header = <FILE>;
		my $quality = <FILE>;
		if ( $quality_header =~  /^(\+[A-Z0-9\_\-\.]+):([0-9]+:[0-9]+:[0-9]+:[0-9]+)\s*/i 
		or $quality_header =~ /^(\+[A-Z0-9\_\-\.]+):([0-9]+:[0-9]+:[0-9]+:[0-9]+#[0-9]+\/(1|2))\s*/i )
		{
			my $new_quality_header = "$1:$2";
			if ( defined $matenumber )
			{
				#### SET BARCODE = 0 IF NOT DEFINED
				$barcode = 0 if not defined $barcode;

				#### ADD BARCODE AND MATE NUMBER
				$sequence_header .= "#$barcode/$matenumber";
				$new_quality_header .= "#$barcode/$matenumber";
			}
			elsif ( defined $barcode )
			{
				#### ADD BARCODE
				$sequence_header .= "#$barcode";
				$new_quality_header .= "#$barcode";
			}
			print OUTFILE "$sequence_header\n$sequence$new_quality_header\n$quality";
		}
		else
		{
			my $rest = <FILE>;
			$rest .= <FILE>;
			$rest .= <FILE>;
			print REJECTFILE $_;
			print REJECTFILE $rest;
		}
	}
    else
    {
		print "simpleHeader.pl    Printing to reject file because no match: $_\n";

		my $rest = <FILE>;
		$rest .= <FILE>;
		$rest .= <FILE>;
        print REJECTFILE $_;
        print REJECTFILE $rest;
    }
}

close(FILE) or die "Can't close inputfile: $inputfile\n" if $inputfile !~ /\.zip$/ and $inputfile !~ /\.gz$/;
close(OUTFILE) or die "Can't close outputfile: $outputfile\n";
close(REJECTFILE) or die "Can't close rejectfile: $rejectfile\n";

##### COMPRESS IF DEFINED
#compress($outputfile, $compress) if defined $compress;
#$outputfile .= ".gz" if defined $compress and $compress eq "gzip";
#$outputfile .= ".zip" if defined $compress and $compress eq "zip";
#
##### COMPRESS IF DEFINED
#compress($rejectfile, $compress) if defined $compress and defined $paired;
#$rejectfile .= ".gz" if defined $compress and $compress eq "gzip" and defined $paired;
#$rejectfile .= ".zip" if defined $compress and $compress eq "zip" and defined $paired;

#### REPORT COMPLETED
print "simpleHeader.pl    outputfile printed:\n\n$outputfile\n\n";
print "simpleHeader.pl    rejectfile printed:\n\n$rejectfile\n\n";

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



sub usage
{
	print `perldoc $0`;

	exit;
}
