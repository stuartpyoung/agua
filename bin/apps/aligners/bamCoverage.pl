#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();

=head2

    APPLICATION     bamCoverage

    VERSION         0.01

    PURPOSE

        GENERATE COVERAGE FOR BAM FILE

    INPUT

        1. *.bam BAM-FORMAT FILE

    OUTPUT

        1. *.bam.coverage COVERAGE FILE OF WINDOW POSITIONS AND READ COVERAGES

		0	12
		100	23
		200	45
		300	32
		...

    USAGE

    ./bamCoverage.pl <--inputfile String> [--help]

		--inputfile		:   Input BAM-format file
		--outputfile	:   Print output to this file (default: STDOUT)
		--begin			:   Begin calculating average from this position
		--end			:   Stop calculating average at this position
		--window		:   Window size (default: 1000bp)
		--help          :   print help info

    EXAMPLES

/nethome/bioinfo/apps/agua/0.4/bin/apps/bamCoverage.pl \
--inputfile /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/autochr22/bowtie/1/chr22/hit.sorted.bam \
--window 100000000


=cut

use strict;

#### EXTERNAL MODULES
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use Getopt::Long;
use File::Path;
use FindBin qw($Bin);

#### USE LIBRARY
use lib "$Bin/../../../lib";	

#### INTERNAL MODULES
use Timer;
use Util;
use Conf::Agua;

#### SET MAQ LOCATION
my $conf = Conf::Agua->new(inputfile=>"$Bin/../../../conf/default.conf");
my $samtools = $conf->getKeyValue("applications", 'SAMTOOLS');

#### GET OPTIONS
my $inputfile;
my $outputfile;
my $begin = 1;
my $end;
my $window;
my $help;
if ( not GetOptions (

	#### GENERAL
    'inputfile=s' 	=> \$inputfile,
    'outputfile=s' 	=> \$outputfile,
    'begin=i' 		=> \$begin,
    'end=i' 		=> \$end,
    'window=i' 		=> \$window,
    'help' 			=> \$help
) )
{ print "bamCoverage.pl    Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
print "bamCoverage.pl    inputfile not defined (use --help option)\n" and exit if not defined $inputfile;

#### CREATE *.bai INDEX FILE IF NOT PRESENT
my $indexfile = "$inputfile.bai";
if ( not -f $indexfile )
{
	my $command = "$samtools/samtools index $inputfile $indexfile";
	print `$command`;
}

#### OPEN OUTPUTFILE
my $outfh;
open($outfh, ">$outputfile") or die "Can't open outputfile: $outputfile" if defined $outputfile;

#### GET THE END OF THE REFERENCE IF NOT DEFINED
if ( not defined $end )
{
	my $command = "$samtools/samtools idxstats $inputfile";
	my $idxstats = `$command`;
	#### FORMAT:
	#### chr22   51304566        914819  0
	#### *       0       0       0
	($end) = $idxstats =~ /^\S+\s+(\S+)/ 
}

my $pileup;
open($pileup, "$samtools/samtools pileup -s $inputfile | cut -f 2,4 |");

#### SET START
my $sum = 0;
my ($position, $depth) = <$pileup> =~ /^(\S+)\s+(\S+)/;
$sum += $depth;

if ( defined $window )
{
	#### SCROLL UNTIL begin IS REACHED
	if ( $position < $begin )
	{
		$sum = 0;
		while ( $position < $begin )
		{
			($position, $depth) = <$pileup> =~ /^(\S+)\s+(\S+)/;
		}
		$sum += $depth;
	}

	#### ADD ONE TO START BECAUSE SAM FILE IS 1-INDEXED
	my $start = (int($position/$window) * $window) + 1;

	#### PRINT EMPTY VALUES UP TO FIRST HIT
	fill_in($begin, $start - 1, $window) if $position >= $window;

	my $stop = $start + $window - 1;

	$stop = $end if $stop > $end;
	my $counter = 0;
	while ( <$pileup> )
	{
		if ( $position >= $stop )
		{
			my $average = sprintf "%.6f", $sum/($stop - $start);
			print "$start\t$stop\t$average\n" if not defined $outputfile;
			print $outfh "$start\t$stop\t$average\n" if defined $outputfile;

			#### RESET VARIABLES
			$sum = 0;
			$counter = 0;
			$start = (int($position/$window) * $window) + 1;
			$stop = $start + $window - 1;
			$stop = $end if $stop > $end;
			last if $start > $stop;
		}
		else
		{
			($position, $depth) = $_ =~ /^(\S+)\s+(\S+)/;
			$sum += $depth;
		}

		$counter++;
	}

	#### DEAL WITH THE LAST WINDOW
	my $average = sprintf "%.6f", $sum/($stop - $start);
	print "$start\t$stop\t$average\n" if $start < $stop and not defined $outputfile;
	print $outfh "$start\t$stop\t$average\n" if $start < $stop and defined $outputfile;

	#### PRINT EMPTY VALUES UP TO FIRST HIT
	fill_in($stop, $end, $window) if $stop < $end;
}
else
{
	#### SCROLL UNTIL begin IS REACHED
	if ( $position < $begin )
	{
		$sum = 0;
		while ( $position < $begin )
		{
			($position, $depth) = <$pileup> =~ /^(\S+)\s+(\S+)/;
		}
		$sum += $depth;
	}

	while ( <$pileup> )
	{
		($position, $depth) = $_ =~ /^(\S+)\s+(\S+)/;
		last if $position > $end;
		$sum += $depth;
	}

	my $range = $end - $begin;
	my $average = sprintf "%.6f", $sum/$range;

	print $outfh "average\t$average\n" if defined $outputfile;
	print "average\t$average\n" if not defined $outputfile;
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

	SUBROUTINE		fill_in

	PURPOSE

		FILL IN EMPTY VALUES FOR STANDARDIZED OUTPUT:

			1. PRINT LEAD IN OF EMPTY BINS BEFORE FIRST HIT

			2. PRINT LEAD OUT OF EMPTY BINS AFTER LAST HIT

	INPUTS

		1. BEGIN POSITION (NB: BAM/SAM COORDINATES ARE 1-BASED)

		2. END POSITION

		3. WINDOW OVER WHICH TO SPACE OUTPUTS


	OUTPUTS

		1. PRINT POSITION AND 0 FOR EACH EMPTY WINDOW IN THE

			RANGE BEGIN --> END

=cut

sub fill_in
{
	my $begin	=	shift;
	my $end		=	shift;
	my $window	=	shift;
	my $outfh	=	shift;

	my $position = $begin;
	while ( $position < $end )
	{
		my $next_position = $position + $window - 1;
		print "$position\t$next_position\t0\n" if not defined $outfh;
		print $outfh "$position\t$next_position\t0\n" if defined $outfh;
		$position += $window;
	}
}



sub usage
{
	print GREEN;
	print `perldoc $0`;
	print RESET;
	exit;
}

