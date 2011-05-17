#!/usr/bin/perl -w

$DEBUG = 1;

#### TIME
my $time = time();

=head2

    APPLICATION     complement

    PURPOSE

        CONVERT A DNA SEQUENCE TO ITS COMPLEMENT

    INPUT

        1. A DNA SEQUENCE AS STDIN OR A FASTA FILE

    OUTPUT

        1. THE COMPLEMENT SEQUENCE AS STDOUT OR PRINTED TO 

			A FASTA FILE (IF --output OPTION IS SPECIFIED)

	USAGE

    ./complement.pl  <--input String> [--output String] [-h]

        --input            :   Location of input single-FASTA file
        --output           :   Destination output file (DEFAULT: input.nospace.suffix)
        --help             :   print help info

    EXAMPLES

chmod 755 /nethome/bioinfo/apps/agua/0.4/bin/apps/utils/complement.pl

/nethome/bioinfo/apps/agua/0.4/bin/apps/utils/complement.pl \
--input CTATCCGCAGGTCCAGGTACC



=cut

use strict;

#### EXTERNAL MODULES
use FindBin qw($Bin);
use Term::ANSIColor qw(:constants);
use Getopt::Long;

use lib "$Bin/../../../lib";

#### INTERNAL MODULES
use Timer;

#### GET OPTIONS
my $input;	
my $output;
my $help;
GetOptions (
    'input=s' => \$input,
    'output=s' => \$output,
    'help' => \$help
) or die "Can't find options. Try $0 --help for usage information\n";

#### CHECK INPUTS
usage() if defined $help;
usage() if not defined $input;

#### DO FILE IF FOUND
if ( -f $input )
{
	open(OUTFILE, ">$output") or die "Can't open output file: $output\n" if defined $output;

	open(FILE, $input) or die "Can't open input file: $input\n";
	$/ = "\\n>";
	foreach my $record ( <FILE> )
	{
		my ($header, $sequence) = $record =~ /^([^\n]+)\n(.+)$/;
		$sequence =~ s/>$//;
		print "Input file is not in FASTA format: $header\n$sequence\n" and exit if $header !~ /^>/ or not defined $header;
		print "No sequence in input file: $input\n" and exit if not defined $sequence;
		my $complement = complement($sequence);

		print OUTFILE "$header\n$complement" if defined $output;
		print "$header\n$complement\n" if not defined $output;
	}
	close(FILE);
	close(OUTFILE);	
}

#### OTHERWISE, DO INPUT FROM STDIN
else
{
	my $complement = complement($input);

	open(OUTFILE, ">$output") or die "Can't open output file: $output\n" if defined $output;
	print OUTFILE $complement if defined $output;
	close(OUTFILE) if defined $output;

	print "$complement\n", if not defined $output;
}

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#									SUBROUTINES
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

=head2

	SUBROUTINE		complement

	PURPOSE

		COMPLEMENT A DNA SEQUENCE

=cut

sub complement
{
	my $sequence	=	shift;

    my $complement = reverse($sequence);
    $complement =~ tr/ACGTacgt/TGCAtgca/;

	return $complement;
}


sub usage
{
    print `perldoc $0`;
	exit;
}

