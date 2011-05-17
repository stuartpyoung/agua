#!/usr/bin/perl -w

$DEBUG = 1;

#### TIME
my $time = time();

=head2

    APPLICATION     nospaces

    PURPOSE

        EXTRACT THE WHOLE SEQUENCE FROM A SINGLE-FASTA FILE TO AN OUTPUT

		FILE WITH NO LINE RETURNS OR SPACES IN THE SEQUENCE

		(IN PREPARATION FOR EXTRACTING BASES FROM THE SEQUENCE USING

		BASE POSITION INFORMATION TO DETERMINE THE STARTING BYTE IN

		POSITION IN THE FILE). IF A DIRECTORY IS SUPPLIED AS THE INPUT

		FILE, DO THE ABOVE FOR ALL FILES IN THE DIRECTORY.

    INPUT

        1. SINGLE-FASTA SEQUENCE FILE

        2. DESTINATION OUTPUT FILE LOCATION

    OUTPUT

        1. SINGLE-FASTA SEQUENCE FILE WITH NO SPACE OR LINE RETURNS IN THE

		SEQUENCE

	USAGE

    ./nospaces.pl  <--inputfile String> [--outputfile String] [-h]

        --inputfile            :   Location of input single-FASTA file
        --outputfile           :   Destination output file (DEFAULT: inputfile.nospace.suffix)
        --help                 :   print help info

    EXAMPLES

/nethome/bioinfo/apps/agua/0.4/bin/apps/utils/nospaces.pl \
--inputfile /nethome/bioinfo/data/sequence/chromosomes/human-fa/hg19/chrY.fa \
--outputfile /nethome/bioinfo/data/sequence/chromosomes/human-fa/hg19/chrY.nospace.fa

=cut

use strict;

#### EXTERNAL MODULES
use File::Path;
use FindBin qw($Bin);
use Getopt::Long;
use Term::ANSIColor qw(:constants);

use lib "$Bin/../../../lib";

#### INTERNAL MODULES
use Timer;

#### GET OPTIONS
my $inputfile;	
my $outputfile;
my $help;
GetOptions (
    'inputfile=s' => \$inputfile,
    'outputfile=s' => \$outputfile,
    'help' => \$help
) or die "Can't find options. Try $0 --help for usage information\n";

#### CHECK INPUTS
usage() if defined $help;
usage() if not defined $inputfile;
print "Can't find inputfile: $inputfile\n" and exit if not -f $inputfile
	and not -d $inputfile;

#### DO SINGLE FILE
if ( -f $inputfile )
{
	print "\nDOING single file\n";
	nospaces($inputfile, $outputfile);
}

#### OR DO ALL FILES IN A DIRECTORY
elsif ( -d $inputfile )
{
	print "\nDOING directory\n";

	#### SET OUTPUT DIR IF NOT DEFINED AND CREATE IF MISSING
	$outputfile = $inputfile if not defined $outputfile;
	File::Path::mkpath($outputfile) if not -d $outputfile;

	chdir($inputfile) or die "Can't change to reference directory: $inputfile\n";
	my @inputfiles = <*>;
	foreach my $infile ( @inputfiles )
	{
		my $outfile = set_outputfile($infile);
		nospaces("$inputfile/$infile", "$outputfile/$infile");
	}
}

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

=head2

	SUBROUTINE		nospaces

	PURPOSE

		EXTRACT THE WHOLE SEQUENCE FROM A SINGLE-FASTA FILE TO AN OUTPUT

		FILE WHILST REMOVING LINE RETURNS AND SPACES IN THE OUTPUT SEQUENCE

=cut

sub nospaces
{
	my $inputfile		=	shift;
	my $outputfile		=	shift;



	$outputfile = set_outputfile($inputfile) if not defined $outputfile;	

	#### OPEN FILES
	my $filehandle;
	open(INFILE, $inputfile) or die "Can't open inputfile: $inputfile\n";
	open(OUTFILE, ">$outputfile") or die "Can't open outputfile: $outputfile\n";

	#### CHECK INPUT FASTA HEADER
	$/ = "\n";
	my $header = <INFILE>;
	print "First input line is not a FASTA header in inputfile: $inputfile\n" and return if not $header =~ /^>/;
	print OUTFILE $header;

	#### PRINT LINES WITHOUT SPACE TO OUTPUT FILE
	while ( <INFILE> )
	{
		next if $_ =~ /^\s*$/;
		if  ( $_ =~ /^>/ )
		{
			print OUTFILE "\n$_";
		}
		else
		{
			$_ =~ s/\s+//g;
			print OUTFILE $_;
		}
	}
	close(INFILE);
	close(OUTFILE);
	print "Outputfile printed:\n\n$outputfile\n\n";
}

=head2

	SUBROUTINE		set_outputfile

	PURPOSE

		SET OUTPUT FILE BASED ON INPUT FILE

=cut

sub set_outputfile
{
	my $inputfile		=	shift;

	my $outputfile;

	$inputfile =~ /^(.+?)\.([^\/^\.]+)$/;
	my $filestub = $1;
	my $suffix = $2;

	if ( defined $filestub )
	{
		$outputfile = $filestub . ".nospace." . $suffix;
	}
	else
	{
		$outputfile = $inputfile . ".nospace";
	}

	return $outputfile;
}


sub usage
{
    print `perldoc $0`;
	exit;
}

