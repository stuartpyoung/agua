#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();

=head2


	************* NOT COMPLETED. INSTEAD USE checkMates.pl THEN simpleHeaders.pl **********


	APPLICATION     simpleHeader

    PURPOSE

		1. CONVERT INCORRECT HEADER FORMATS TO THE SIMPLE HEADER FORMAT

			ACCEPTED BY ELAND AND OTHER ALIGNERS MATCHING THIS REGEX:

			[A-Z0-9\.]:[0-9]:[0-9]:[0-9]:[0-9](#[0-9]\/[12])?

			WHERE THE LAST PORTION DIFFERS DEPENDING ON WHETHER THE $input_fh

			CONTAINS A SINGLE OR PAIRED END READS

    INPUTS

        1. INPUT $input_fh

    OUTPUTS

		1. OUTPUT $input_fh OF CORRECTLY FORMATTED READS

		2. REJECTS $input_fh OF READS THAT WERE NOT ABLE TO BE FORMATTED CORRECTLY

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
--inputfile /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/SRX000600/SRR004865_1.fastq.gz \
--outputfile /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/SRX000600/SRR004865_1.fastq \
--rejectfile /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/SRX000600/SRR004865_1.fastq \


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
my $barcode = 0;
my $paired;
my $dot = 1000000;
my $help;
GetOptions (
    'inputfile=s' 	=> \$inputfile,
    'outputfile=s' 	=> \$outputfile,
    'rejectfile=s' 	=> \$rejectfile,
    'barcode=i' 	=> \$barcode,
    'paired' 		=> \$paired,
    'dot=i' 		=> \$dot,
    'help' 			=> \$help             
) or die "No options specified. Use --help for usage\n";
usage() if defined $help;

#### CHECK INPUTS
die "inputfile is not defined (use --help for options)\n" if not defined $inputfile;
die "outputfile is not defined (use --help for options)\n" if not defined $outputfile;
die "rejectfile is not defined (use --help for options)\n" if not defined $rejectfile;

#### OPEN $input_fh AND SET RECORD SEPARATOR
$/ = "\n";
my ($input_fh, $output_fh, $mate_fh, $mateout_fh);
if( $inputfile =~ /\.gz$/ or $inputfile =~ /\.zip$/ )
{
	my $pipe_command = "zcat $inputfile |";
	open($input_fh, $pipe_command) or die "Can't open inputfile: $inputfile\n"
}
else
{
	open($input_fh, $inputfile) or die "Can't open input file: $inputfile\n";	
}

#### OPEN MATE INPUT $input_fh IF paired IS DEFINED
my $matefile = $inputfile;
$matefile =~ s/_1\./_2\./;
if ( defined $paired )
{
	#### OPEN $input_fh AND SET RECORD SEPARATOR
	if( $inputfile =~ /\.gz$/ or $inputfile =~ /\.zip$/ )
	{
		my $pipe_command = "zcat $inputfile |";
		open($mate_fh, $pipe_command) or die "Can't open matefile: $matefile\n";
	}
	else
	{
		open($mate_fh, "$matefile") or die "Can't open matefile: $matefile";
	}
}

#### OPEN MATE OUTPUT $input_fh IF paired IS DEFINED
my $mateoutfile = $outputfile;
$mateoutfile =~ s/_1\./_2\./;
open($mateout_fh, ">$mateoutfile") or die "Can't open mateoutfile: $mateoutfile" if defined $paired;

#### OPEN REJECT $input_fh
open($reject_fh, ">$rejectfile") or die "Can't open rejectfile: $rejectfile";

#### OPEN $output_fh
open($output_fh, ">$outputfile") or die "Can't open outputfile: $outputfile\n";
$/ = "\n";

#### PROCESS READS
my ($sequence_header, $sequence, $quality_header, $quality, $trash);
my ($mate_sequence_header, $mate_sequence, $mate_quality_header, $mate_quality);
while ( <$input_fh> )
{
    if ( $_ =~ /^(@[A-Z0-9\.]+)\s+\S*([0-9]+:[0-9]+:[0-9]+:[0-9]+)\s*/ 
		or $_ =~ /^(@[A-Z0-9\.]+)\s+\S*([0-9]+:[0-9]+:[0-9]+:[0-9]+#[0-9]+)\s*/
		or $_ =~ /^(@[A-Z0-9\.]+)\s+\S*([0-9]+:[0-9]+:[0-9]+:[0-9]+#[0-9]+\/(1|2))\s*/ )
	{
		$sequence_header = "$1:$2";
		$sequence = <$input_fh>;
		$trash = <$input_fh>;
		$quality = <$input_fh>;
		$quality_header = $sequence_header;
		$quality_header =~ s/^@/+/;

		my $new_quality_header = "$1:$2";
		if ( defined $paired )
		{
			#### SET BARCODE = 0 IF NOT DEFINED
			$barcode = 0 if not defined $barcode;

			#### ADD BARCODE AND MATE NUMBER
			my $matenumber = 1;
			$sequence_header .= "#$barcode/$matenumber";
			$new_quality_header .= "#$barcode/$matenumber";

			#### DO MATE
			$mate_sequence_header = <$mate_fh>;
			$mate_sequence_header =~ /^(@[A-Z0-9\.]+)\s+\S*([0-9]+:[0-9]+:[0-9]+:[0-9]+)\s*/
			or $mate_sequence_header =~ /^(@[A-Z0-9\.]+)\s+\S*([0-9]+:[0-9]+:[0-9]+:[0-9]+#[0-9]+)\s*/
			or $mate_sequence_header =~ /^(@[A-Z0-9\.]+)\s+\S*([0-9]+:[0-9]+:[0-9]+:[0-9]+#[0-9]+\/(1|2))\s*/;

			$mate_sequence = <$mate_fh>;
			$mate_quality_header = <$mate_fh>;
			$mate_quality = <$mate_fh>;

			#### REMOVE 'length=...' 
			$mate_sequence_header =~ s/\s+length=.+$//;
			$mate_quality_header =~ s/\s+length=.+$//;


			print $mateout_fh "$mate_sequence_header$mate_sequence\n$mate_quality_header$mate_quality\n";

		}


			elsif ( defined $barcode )
			{
				#### ADD BARCODE
				$sequence_header .= "#$barcode";
				$new_quality_header .= "#$barcode";
			}
			print $output_fh "$sequence_header\n$sequence$new_quality_header\n$quality";
		}
		else
		{
			my $rest = <$input_fh>;
			$rest .= <$input_fh>;
			$rest .= <$input_fh>;
			print $reject_fh $_;
			print $reject_fh $rest;

			#### DO MATE $input_fh
			if ( defined $paired )
			{
				$mate_sequence_header = <$mate_fh>;
				$mate_sequence = <$mate_fh>;
				$mate_quality_header = <$mate_fh>;
				$mate_quality = <$mate_fh>;
				print $reject_fh "$mate_sequence_header$mate_sequence\n$mate_quality_header$mate_quality\n";
			}
		}
    }
    else
    {
		my $rest = <$input_fh>;
		$rest .= <$input_fh>;
		$rest .= <$input_fh>;
        print $reject_fh $_;
        print $reject_fh $rest;

		#### DO MATE $input_fh
		if ( defined $paired )
		{
			$mate_sequence_header = <$mate_fh>;
			$mate_sequence = <$mate_fh>;
			$mate_quality_header = <$mate_fh>;
			$mate_quality = <$mate_fh>;
			print $reject_fh "$mate_sequence_header$mate_sequence\n$mate_quality_header$mate_quality\n";
		}
    }	
}

close($input_fh) or die "Can't close inputfile: $inputfile\n" if $inputfile !~ /\.zip$/ and $inputfile !~ /\.gz$/;
close($output_fh) or die "Can't close outputfile: $outputfile\n";
close($reject_fh) or die "Can't close rejectfile: $rejectfile\n";

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


__END__

#### PROCESS READS
my ($sequence_header, $sequence, $quality_header, $quality);
my ($mate_sequence_header, $mate_sequence, $mate_quality_header, $mate_quality);
while ( <$input_fh> )
{
    if ( $_ =~ /^(@[A-Z0-9\.]+)\s+\S*([0-9]+:[0-9]+:[0-9]+:[0-9]+)\s*/ 
		or $_ =~ /^(@[A-Z0-9\.]+)\s+\S*([0-9]+:[0-9]+:[0-9]+:[0-9]+#[0-9]+\/(1|2))\s*/ )
	{
		$sequence_header = "$1:$2";
		$sequence = <$input_fh>;
		$quality_header = <$input_fh>;
		$quality = <$input_fh>;
		if ( $quality_header =~  /^(\+[A-Z0-9\.]+)\s+\S*([0-9]+:[0-9]+:[0-9]+:[0-9]+)\s*/ 
		or $quality_header =~ /^(\+[A-Z0-9\.]+)\s+\S*([0-9]+:[0-9]+:[0-9]+:[0-9]+#[0-9]+\/(1|2))\s*/ )
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
			print $output_fh "$sequence_header\n$sequence$new_quality_header\n$quality";


			#### DO $mate_fh IF paired IS DEFINED		
			if ( defined $paired )
			{
				$mate_sequence_header = <$mate_fh>;
				$mate_sequence = <$mate_fh>;
				$mate_quality_header = <$mate_fh>;
				$mate_quality = <$mate_fh>;

				$mate_sequence =~ s/\s+//g;
				$mate_quality =~ s/\s+//g;

				#### REMOVE 'length=...' 
				$mate_sequence_header =~ s/\s+length=.+$//;
				$mate_quality_header =~ s/\s+length=.+$//;


				print $mateout_fh "$mate_sequence_header$mate_sequence\n$mate_quality_header$mate_quality\n";
			}
		}
		else
		{
			my $rest = <$input_fh>;
			$rest .= <$input_fh>;
			$rest .= <$input_fh>;
			print $reject_fh $_;
			print $reject_fh $rest;

			#### DO MATE $input_fh
			if ( defined $paired )
			{
				$mate_sequence_header = <$mate_fh>;
				$mate_sequence = <$mate_fh>;
				$mate_quality_header = <$mate_fh>;
				$mate_quality = <$mate_fh>;
				print $reject_fh "$mate_sequence_header$mate_sequence\n$mate_quality_header$mate_quality\n";
			}
		}
    }
    else
    {
		my $rest = <$input_fh>;
		$rest .= <$input_fh>;
		$rest .= <$input_fh>;
        print $reject_fh $_;
        print $reject_fh $rest;

		#### DO MATE $input_fh
		if ( defined $paired )
		{
			$mate_sequence_header = <$mate_fh>;
			$mate_sequence = <$mate_fh>;
			$mate_quality_header = <$mate_fh>;
			$mate_quality = <$mate_fh>;
			print $reject_fh "$mate_sequence_header$mate_sequence\n$mate_quality_header$mate_quality\n";
		}
    }	
}

close($input_fh) or die "Can't close inputfile: $inputfile\n" if $inputfile !~ /\.zip$/ and $inputfile !~ /\.gz$/;
close($output_fh) or die "Can't close outputfile: $outputfile\n";
close($reject_fh) or die "Can't close rejectfile: $rejectfile\n";

