#!/usr/bin/perl -w

#### DEBUG
$DEBUG = 1;

#### TIME
my $time = time();

=head2

	APPLICATION     checkHeader

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

    ./checkHeader.pl <--inputfile String> <format String> [-h]

    --inputfile		:   /full/path/to/inputfile
    --outputfile	:   /full/path/to/outputfile
    --rejectfile	:   /full/path/to/rejectfile
    --dot			:   Print progress count per multiple of this integer
    --help			:   print help info

    EXAMPLES

/nethome/bioinfo/apps/agua/0.5/bin/apps/readprep/checkHeader.pl \
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
my $dot = 1000000;
my $help;
GetOptions (
    'inputfile=s' 	=> \$inputfile,
    'outputfile=s' 	=> \$outputfile,
    'rejectfile=s' 	=> \$rejectfile,
    'dot=i' 		=> \$dot,
    'help' 			=> \$help             
) or die "No options specified. Use --help for usage\n";
usage() if defined $help;

#### CHECK INPUTS
die "inputfile is not defined (use --help for options)\n" if not defined $inputfile;
die "outputfile is not defined (use --help for options)\n" if not defined $outputfile;
die "rejectfile is not defined (use --help for options)\n" if not defined $rejectfile;

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
	print "Processing read: $_\n";
    if ( $_ =~ /^@[A-Z0-9\.]+:[0-9]+:[0-9]+:[0-9]+:[0-9]+\s*$/ 
		or $_ =~ /^@[A-Z0-9\.]+:[0-9]+:[0-9]+:[0-9]+:[0-9]+#[0-9]+\/(1|2)\s*$/ )
	{
		print "GOOD READ: $_\n";
		my $rest = <FILE>;
		$rest .= <FILE>;
		$rest .= <FILE>;
		print "rest: $rest\n";

        print OUTFILE $_;
        print OUTFILE $rest;

    }
    else
    {
		print "BAD READ: $_\n";
		my $rest = <FILE>;
		$rest .= <FILE>;
		$rest .= <FILE>;
		print "rest: $rest\n";

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
print "checkHeader.pl    outputfile printed:\n\n$outputfile\n\n";
print "checkHeader.pl    rejectfile printed:\n\n$rejectfile\n\n";

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
