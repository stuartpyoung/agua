#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();

=head2

    APPLICATION     sortfile

    VERSION         0.01

    PURPOSE

        SORT A FILE BY COLUMN USING THE LINUX sort COMMAND

    INPUT

        1. INPUT FILE

		2. (OPTIONAL) OUTPUT FILE

    OUTPUT

        1. SORTED INPUT FILE OR OUTPUT FILE (IF SPECIFIED)

    USAGE

    ./sortfile.pl <--inputfile String>	[--column Integer]
				[--outputfile String]
					[--outputfile String]	[--help]

		--inputfile			:   Location of input file
		--column			:   Column number (1,2,3, ...)
		--outputfile		:   (Optional) Location of output file
		--numeric			:   Perform numeric sort on column
		--reverse			:   Order in reverse (e.g., from highest to lowest number)
		--help              :   print help info

    EXAMPLES

chmod 755 /nethome/bioinfo/apps/agua/0.4/bin/apps/utils/sortfile.pl

/nethome/bioinfo/apps/agua/0.4/bin/apps/utils/sortfile.pl \
--inputfile /nethome/syoung/base/pipeline/snpfilter/sortfile/454HCDiffs-header-chr22.txt \
--outputfile /nethome/syoung/base/pipeline/snpfilter/sortfile/454HCDiffs-header-chr22-sorted.txt \
--column 2 \
--numeric 


/nethome/bioinfo/apps/agua/0.4/bin/apps/utils/sortfile.pl \
--inputfile /nethome/syoung/base/pipeline/snpfilter/sortfile/454HCDiffs-header-chr22.txt \
--outputfile /nethome/syoung/base/pipeline/snpfilter/sortfile/454HCDiffs-header-chr22-sorted-reverse.txt \
--column 2 \
--numeric \
--reverse


=cut

use strict;

#### EXTERNAL MODULES
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use Getopt::Long;
use FindBin qw($Bin);

#### USE LIBRARY
use lib "$Bin/../../../lib";	
use lib "$Bin/../../../lib/external";	

#### INTERNAL MODULES
use FileTools;
use Timer;

#### GET OPTIONS
my $column;
my $numeric;
my $reverse;
my $inputfile;
my $outputfile;
my $help;
if ( not GetOptions (
    'column=s' 		=> \$column,
    'numeric' 	=> \$numeric,
    'reverse' 	=> \$reverse,
    'inputfile=s' 	=> \$inputfile,
    'outputfile=s' 	=> \$outputfile,
    'help' 			=> \$help
) )
{ print "sortfile.pl    Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "column not defined (Use --help for usage)\n" if not defined $column;
die "outputfile not defined (Use --help for usage)\n" if not defined $outputfile;
die "inputfile not defined (Use --help for usage)\n" if not defined $inputfile;
print "sortfile.pl    column: $column\n";
print "sortfile.pl    inputfile: $inputfile\n" if defined $inputfile;
print "sortfile.pl    outputfile: $outputfile\n";

my $filetool = FileTools->new();
print "sortfile.pl    filetool: $filetool\n";

print "sortfile.pl    Doing filetool->sortFile(args)\n";
$filetool->sortFile(
	{
		inputfile	=>	$inputfile,
		outputfile	=>	$outputfile,
		column		=>	$column,
		numeric		=>	$numeric,
		reverse		=>	$reverse
	}	
);



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


