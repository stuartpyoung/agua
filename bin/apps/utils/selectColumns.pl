#!/usr/bin/perl -w

$DEBUG = 1;

=head2

    TEST        selectColumns

    PURPOSE

		REMOVE DUPLICATE LINES BASED ON USER-INPUT KEYS 


	EXAMPLES

./selectColumns.pl \
--columns 5,2,1,3,4 \
--inputfile /nethome/syoung/base/pipeline/dbsnp/snp130-chr2.tsv \
--outputfile /nethome/syoung/base/pipeline/dbsnp/snp130-chr2-uniq.tsv 


=cut

use strict;

#### USE FINDBIN
use FindBin qw($Bin);

#### USE LIBS
use lib "$Bin/../../../lib";
use lib "$Bin/../../../lib/external";

#### INTERNAL MODULES
use FileTools;

#### EXTERNAL MODULES
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use Getopt::Long;
use File::Path;
use File::Copy;

##### STORE ARGUMENTS TO PRINT TO FILE LATER
my @arguments = @ARGV;

#### GET OPTIONS
my $inputfile;
my $columns;
my $outputfile;
my $help;
if ( not GetOptions (
    'inputfile=s'  => \$inputfile,
    'columns=s'   => \$columns,
    'outputfile=s'   => \$outputfile,
    'help'          => \$help
) )
{ print "Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "inputfile not defined (option --help for usage)\n" if not defined $inputfile;
die "outputfile not defined (Use --help for usage)\n" if not defined $outputfile;
die "columns not defined (Use --help for usage)\n" if not defined $columns;

#### DEBUG
print "selectColumns.pl    inputfile: $inputfile\n";
print "selectColumns.pl    outputfile: $outputfile\n";
print "selectColumns.pl    columns: $columns\n";

#### CHECK COLS ARE NUMERIC
my @cols = split ",", $columns;
foreach my $column ( @cols )
{
	print "column '$column' is non-numeric. Exiting\n" and exit if $column !~ /^\d+$/;
}

#### INSTANTIATE FileTools OBJECT
my $filetool = FileTools->new();

#### GET INPUT AND OUTPUT FILES
my @infiles = split ",", $inputfile;
my @outfiles = split ",", $outputfile;

#### CHECK EQUAL NUMBER OF INPUT AND OUTPUT FILES
print "No. inputfiles (", $#infiles + 1 ,") does not equal number of outputfiles (", $#outfiles + 1 , ")\n" and exit if $#infiles != $#outfiles;


for ( my $i = 0; $i < $#infiles + 1; $i++ )
{
	$filetool->selectColumns(
		{
			inputfile	=>	$infiles[$i],
			outputfile	=>	$outfiles[$i],
			columns		=>	$columns
		}
	);

}


################################################################################
################################################################################
########################           SUBROUTINES          ########################
################################################################################
################################################################################




