#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();

=head2

	APPLICATION     polyAFilter

    PURPOSE

		1. REMOVE READS THAT ARE ALL POLY-A

    INPUTS

        1. INPUT FILE

		2. MATE FILE

    OUTPUTS

		1. FILE OF READS NONE OF WHICH CONTAIN ALL POLY-A

		2. REJECTS FILE OF READS THAT CONTAIN ALL POLY-A

		3. CORRESPONDING MATE AND MATE REJECTS FILES IF PAIRED

    USAGE

    ./polyAFilter.pl <--inputfile String> <--outputfile String> <--rejectfile String> <paired Boolean> [-h]

    --inputfile		:   /full/path/to/inputfile
    --outputfile	:   /full/path/to/outputfile containing no poly-A reads
	--rejectfile    :   /full/path/to/rejectfile containing rejected reads
	--paired		:   Do poly-A filter on mate file reads at same time
	--dot			:   Print progress count per multiple of this integer
    --help			:   print help info

    EXAMPLES

cd /scratch/syoung/base/pipeline/catfilter

/nethome/bioinfo/apps/agua/0.5/bin/apps/readprep/polyAFilter.pl \
--inputfile /scratch/syoung/base/pipeline/catfilter/srx001539-8.100reads.reads_1.fastq \
--outputfile /scratch/syoung/base/pipeline/catfilter//srx001539-8.100reads.nopolya.reads_1.fastq \
--rejectfile /scratch/syoung/base/pipeline/catfilter//srx001539-8.100reads.polya.reads_1.fastq

    Run time: 00:00:03
    Completed /nethome/bioinfo/apps/agua/0.5/bin/apps/readprep/polyAFilter.pl
    12:38AM, 21 September 2010
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
my $paired;
my $dot = 1000000;
my $help;
GetOptions (
    'inputfile=s' 	=> \$inputfile,
    'outputfile=s' 	=> \$outputfile,
    'rejectfile=s' 	=> \$rejectfile,
    'paired' 		=> \$paired,
    'dot=i' 		=> \$dot,
    'help' 			=> \$help             
) or die "No options specified. Use --help for usage\n";
usage() if defined $help;

#### CHECK INPUTS
die "inputfile is not defined (use --help for options)\n" if not defined $inputfile;
die "outputfile is not defined (use --help for options)\n" if not defined $outputfile;
die "rejectfile is not defined (use --help for options)\n" if not defined $rejectfile;

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


#### OPEN OUTFILE
open(OUTFILE, ">$outputfile") or die "Can't open outputfile: $outputfile\n";

#### OPEN REJECT FILE
open(REJECT, ">$rejectfile") or die "Can't open rejectfile: $rejectfile\n";

#### OPEN MATE OUTPUT FILE IF paired IS DEFINED
my $mateoutfile = $outputfile;
$mateoutfile =~ s/_1\./_2\./;
open(MATEOUTFILE, ">$mateoutfile") or die "Can't open mateoutfile: $mateoutfile" if defined $paired;

#### OPEN MATE OUTPUT FILE IF paired IS DEFINED
my $mate_rejectfile = $rejectfile;
$mate_rejectfile =~ s/_1\./_2\./;
open(MATEREJECT, ">$mate_rejectfile") or die "Can't open mate_rejectfile: $mate_rejectfile" if defined $paired;


$/ = "\n";
my $counter = 0;
while ( <FILE> )
{    
	print "$counter\n" if $counter % $dot == 0;
	$counter++;

	#### GET FOUR LINES
	my $sequence_header = $_;
	my $sequence = <FILE>;
	my $quality_header = <FILE>;
	my $quality = <FILE>;

	print "sequence_header not defined \n" and last if not defined $sequence_header;
	print "sequence not defined \n" and last if not defined $sequence;
	print "quality_header not defined \n" and last if not defined $quality_header;
	print "quality_header not defined \n" and last if not defined $quality_header;



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

	}

	#### REJECT IF READ IS ALL 'A's
	if ( $sequence =~ /^A/ and $sequence =~ /^[A]+$/ 
		or (defined $paired and $sequence =~ /^A/ and $sequence =~ /^[A]+$/ ) ) 
	#if ( $sequence =~ /^[A]+$/ 
	#	or (defined $paired and $mate_sequence =~ /^[A]+$/) ) 
	{
		print REJECT "$sequence_header$sequence$quality_header$quality";
		print MATEREJECT "$mate_sequence_header$mate_sequence$mate_quality_header$mate_quality" if defined $paired;
	}
	else
	{
		print OUTFILE "$sequence_header$sequence$quality_header$quality";
		print MATEOUTFILE "$mate_sequence_header$mate_sequence$mate_quality_header$mate_quality" if defined $paired;
	}

	#last if $counter >= 100000;
}
close(FILE) or die "Can't close inputfile: $inputfile\n" if $inputfile !~ /\.zip$/ and $inputfile !~ /\.gz$/;
close(MATEFILE) or die "Can't close matefile: $matefile\n" if defined $paired and $matefile !~ /\.zip$/ and $matefile !~ /\.gz$/;
close(OUTFILE) or die "Can't close outputfile: $outputfile\n";
close(MATEOUTFILE) or die "Can't close mateoutfile: $mateoutfile\n" if defined $paired;
close(REJECT) or die "Can't close rejectfile: $rejectfile\n" if defined $paired;
close(MATEREJECT) or die "Can't close rejectfile: $rejectfile\n" if defined $paired;

#### REPORT COMPLETED
print "polyAFilter.pl    outputfile printed:\n\n$outputfile\n\n";
print "polyAFilter.pl    rejectfile printed:\n\n$rejectfile\n\n";

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
