#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();

=head2

    APPLICATION     findMatch

    VERSION         0.01

    PURPOSE

        FIND THE POSITION OF AN EXACT MATCH OF THE SUPPLIED SEQUENCE

		IN A TARGET FASTA FILE CONTAINING ONE FASTA RECORD

    INPUT

        1. QUERY SEQUENCE

		2. TARGET TARGETERENCE FASTA FILE

    OUTPUT

        1. TABLE WITH THE FOLLOWING FIELDS

			Name	dbSNP_start	matched_start

    USAGE

    ./findMatch.pl <--query String> <--targetfile String> [--help]

		--query				:   /full/path/to/readfile.fastq
		--outputfile				:   /full/path/to/output_directory
		--targetfile			:   /full/path/to/squashed_genome_files
		--help                 	:   print help info

    EXAMPLES


chmod 755 /nethome/bioinfo/apps/agua/0.4/bin/apps/utils/findMatch.pl

/nethome/bioinfo/apps/agua/0.4/bin/apps/utils/findMatch.pl \
--query TTTGGACGCGAGGAGTAAATCAAATGGTACCAACTTCTACT \
--targetfile /nethome/bioinfo/data/sequence/chromosomes/human-fa/hg19/chrY.fa \
--outputfile /nethome/syoung/base/pipeline/snpfilter/matches/rs2521574-chrY.match

NO MATCHES AS EXPECTED. NOW TRIED FLANKS OF SNP:

/nethome/bioinfo/apps/agua/0.4/bin/apps/utils/findMatch.pl \
--query TTTGGACGCGAGGAGTAA \
--targetfile /nethome/bioinfo/data/sequence/chromosomes/human-fa/hg19/chrY.fa \
--outputfile /nethome/syoung/base/pipeline/snpfilter/matches/rs2521574-chrY.match



MULTI-FASTA

--targetfile /nethome/bioinfo/data/sequence/genes/human/refseq/40/refseqgene.genomic.fna \


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
my $query;
my $targetfile;
my $outputfile;
my $help;
if ( not GetOptions (
    'query=s' 		=> \$query,
    'targetfile=s' 	=> \$targetfile,
    'outputfile=s' 	=> \$outputfile,
    'help' 			=> \$help
) )
{ print "findMatch.pl    Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "query not defined (Use --help for usage)\n" if not defined $query;
die "outputfile not defined (Use --help for usage)\n" if not defined $outputfile;
die "targetfile not defined (Use --help for usage)\n" if not defined $targetfile;
print "findMatch.pl    query: $query\n";
print "findMatch.pl    targetfile: $targetfile\n";
print "findMatch.pl    outputfile: $outputfile\n";

#### OPEN OUTPUT FILE
open(OUTFILE, ">$outputfile") or die "Can't open output file: $outputfile\n";

#### STORE TARGETERENCE IN MEMORY
open(TARGET, $targetfile) or die "Can't open target file: $targetfile\n";
$/ = "END OF FILE";
my $target = <TARGET>;
close(TARGET);

#### SKIP PAST FASTA HEADER
$target =~ s/^>[^\n]+\n//;

#### REMOVE WHITESPACE
$target =~ s/\s+//g;

#### SET LENGTH TO PRINT IN CAPITAL LETTERS (UPPERCASE)
my $length = length($query);

#### PRINT ALL PERFECT MATCH POSITIONS
while ( $target =~ /($query)/g ) {

	my $end_position = pos $target;
	print $end_position - $length, "\n";
	my $prestring = lc(substr($target, $end_position - $length - 20, 20));
	my $substring = uc(substr($target, $end_position - $length, $length ));
	my $poststring = lc(substr($target, $end_position, 20));
	print "$prestring$substring$poststring\n";

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


