#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();

=head2

    APPLICATION     snpFinder

    VERSION         0.01

    PURPOSE

        VERIFY THE POSITIONS OF dbSNP ENTRIES IN THE HUMAN REFERENCE

    INPUT

        1. dbSNP BUILD SEQUENCE WITH 10BP MARGINS

		2. HUMAN REFERENCE SEQUENCE

    OUTPUT

        1. TABLE WITH THE FOLLOWING FIELDS

			Name	dbSNP_start	matched_start

    USAGE

    ./snpFinder.pl <--inputfile String> <--outputfile String> <--referencefile String> [--help]

		--inputfile				:   /full/path/to/readfile.fastq
		--outputfile				:   /full/path/to/output_directory
		--referencefile			:   /full/path/to/squashed_genome_files
		--help                 	:   print help info

    EXAMPLES


chmod 755 /nethome/bioinfo/apps/agua/0.4/bin/apps/snpFinder.pl

/nethome/bioinfo/apps/agua/0.4/bin/apps/snpFinder.pl \
--inputfile /nethome/syoung/base/pipeline/benchmark/data/duan/run12/s_1_1_sequence.fastq \
--referencefile /nethome/bioinfo/data/sequence/chromosomes/human-sq \
--outputfile /nethome/syoung/base/pipeline/benchmark/eland/

=cut

use strict;

#### EXTERNAL MODULES
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use Getopt::Long;
use FindBin qw($Bin);

#### USE LIBRARY
use lib "$Bin/../../../lib";	

#### INTERNAL MODULES
use Timer;
use Util;
use Conf::Agua;

#### GET OPTIONS
my $inputfile;
my $referencefile;
my $outputfile;
my $help;
if ( not GetOptions (
    'inputfile=s' 	=> \$inputfile,
    'referencefile=s' 	=> \$referencefile,
    'outputfile=s' 	=> \$outputfile,
    'help' 			=> \$help
) )
{ print "snpFinder.pl    Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "Input file not defined (Use --help for usage)\n" if not defined $inputfile;
die "Output directory not defined (Use --help for usage)\n" if not defined $outputfile;
die "Reference file not defined (Use --help for usage)\n" if not defined $referencefile;
print "snpFinder.pl    inputfile: $inputfile\n";
print "snpFinder.pl    referencefile: $referencefile\n" if defined $referencefile;
print "snpFinder.pl    outputfile: $outputfile\n";

#### OPEN OUTPUT FILE
open(OUTFILE, ">$outputfile") or die "Can't open output file: $outputfile\n";
print OUTFILE "name\tstart\tstop\tstrand\tposition\tdifference\treference_matched\tsequence\treference_base\tdbsnp_base\n";


#### STORE REFERENCE IN MEMORY
open(REF, $referencefile) or die "Can't open reference file: $referencefile\n";
$/ = "END OF FILE";
my $reference = <REF>;
close(REF);
$reference =~ s/^>[^\n]+\n//;


#### GO THROUGH ALL SNPS
$/ = "\n>";
open(INFILE, $inputfile) or die "Can't open inputfile: $inputfile\n";
while ( <INFILE> )
{
	$_ =~ s/\s*>$//;

	#### >hg19_snp130_rs9785927 range=chrY:150850-150870 5'pad=10 3'pad=10 strand=+ repeatMasking=none
	#### CTCTGATGGGTGGGCAGGTGA

	my ($name, $chromosome, $start, $stop, $strand, $sequence) = $_ =~ /_(rs\d+) range=([^:]+):(\d+)-(\d+).+?strand=(.).+?\n\s*(\S+)$/;
	next if not defined $sequence;
	next if length($sequence) != 21;

	my ($dbsnp_base) = $sequence =~ /^.{10,10}(.)/;
	$sequence =~ s/^(.{10,10})./$1\./;

	my $matched = 0;
	while ( $reference =~ /($sequence)/g ) {
		my $reference_matched = $1;
		my ($reference_base) = $reference =~ /^.{10,10}(.)/;

		$matched = 1;
		my $position = pos $reference;
		my $difference = $position - $start;
		print OUTFILE "$name\t$start\t$stop\t$strand\t$position\t$difference\t$reference_matched\t$sequence\t$reference_base\t$dbsnp_base\n";
    }

	if ( not $matched )
	{
		print OUTFILE "$name\t$start\t$stop\t$strand\n";
	}

}
close(OUTFILE);

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

	SUBROUTINE		set_parameter

	PURPOSE

		SET THE VALUE OF A PARAMETER IN arguments

=cut

sub set_parameter
{	
	my $arguments		=	shift;
	my $parameter			=	shift;
	my $value			=	shift;


	for ( my $i = 0; $i < @$arguments; $i++ )
	{
		if ( "--$parameter" eq $$arguments[$i] )
		{
			$$arguments[$i + 1] = $value;
			return $arguments;
		}	
	}

	return $arguments;
}



=head2

	SUBROUTINE		fill_in

	PURPOSE

		SUBSTITUTE counter FOR ONE OR MORE '%COUNTER%' STRINGS IN ALL ARGUMENTS

=cut

sub fill_in
{	
	my $arguments		=	shift;
	my $pattern			=	shift;
	my $value			=	shift;

	print "\n";

	foreach my $argument ( @$arguments )
	{
		$argument =~ s/$pattern/$value/ig;
	}

	return $arguments;
}


=head2

	SUBROUTINE		get_argument

	PURPOSE

		EXTRACT AN ARGUMENT FROM THE ARGUMENT ARRAY AND RETURN IT 

=cut

sub get_argument
{	
	my $arguments		=	shift;
	my $name			=	shift;

print "\n";


	my $argument;
	for ( my $i = 0; $i < @$arguments; $i++ )
	{
		if ( "--$name" eq $$arguments[$i] )
		{
			$argument = $$arguments[$i + 1];
			splice @$arguments, $i, 2;
			return $argument;
		}
	}

	return;
}



sub usage
{
	print GREEN;
	print `perldoc $0`;
	print RESET;
	exit;
}


