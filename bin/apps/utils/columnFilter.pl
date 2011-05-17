#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();

=head2

	APPLICATION     columnFilter

    PURPOSE

        FILTER LINES IN A FILE BASED ON THE VALUE OF A USER-DEFINED COLUMN

    INPUT

        1. INPUT FILE 

        2. COLUMN NUMBER (1, 2, 3, ...)

        3. OPERATOR

		4. VALUE

    OUTPUT

        1. FILE CONTAINING LINES PASSING THE threshold THRESHOLD VALUE IN THE USER-DEFINED COLUMN

    USAGE

    ./columnFilter.pl  <--inputfile String> <--outputfile String> <--column Integer> <--threshold String> [--separator String] [-h]

        --inputfile		       :   /Full/path/to/inputfile 
        --outputfile		   :   /Full/path/to/outputfile 
        --operator             :   Operator for the column value (5 choices: ">=" or ">" or "<=" or "<" or "==" or "=")
        --column 	           :   Use data in this column
        --threshold            :   Threshold value for operator
        --separator            :   (optional) Column separator (default = white space)
        --help                 :   print help info

    EXAMPLES

/nethome/syoung/base/bin/utils/columnFilter.pl --inputfile /mihg/data/NGS/syoung/base/pipeline/SRA/SRA000271/NA18507.snp --outputfile /mihg/data/NGS/syoung/base/pipeline/SRA/SRA000271/NA18507-20pctl-filtered.snp --column 5 --operator ">" --threshold 81

/nethome/syoung/base/bin/utils/columnFilter.pl --inputfile /mihg/data/NGS/syoung/base/pipeline/SRA/SRA000271/NA18507-short.snp --outputfile /mihg/data/NGS/syoung/base/pipeline/SRA/SRA000271/NA18507-short-20pctl-filtered.snp --column 5 --operator ">" --threshold 26

=cut

use strict;

#### USE LIBRARY
use FindBin qw($Bin);
use lib "$Bin/../../../lib";

#### FLUSH BUFFER
$| = 1;

#### INTERNAL MODULES
use Timer;
use Util;
use Conf::Agua;

#### EXTERNAL MODULES
use Term::ANSIColor qw(:constants);
use Getopt::Long;

#### SET DEFAULT SEPARATOR
my $DEFAULT_SEPARATOR = "\\s+";

#### SAVE ARGUMENTS
my @arguments = @ARGV;

#### GET OPTIONS
my $inputfile;
my $outputfile;
my $column;
my $threshold;
my $operator;
my $separator;
my $help;
if ( not GetOptions (
    'inputfile=s' => \$inputfile,
    'outputfile=s' => \$outputfile,
    'column=i' => \$column,
    'threshold=s' => \$threshold,
    'operator=s' => \$operator,
    'separator=s' => \$separator,
    'help' => \$help
) )
{ usage(); exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### SET VALUE TO ZERO IF inputfile NOT DEFINED
#### OTHERWISE, SET VALUE TO THE inputfile (ORDINAL) NUMBER
#### WHERE 1 IS THE LAST inputfile
die "Input file not defined (Use --help for usage)\n" if not defined $inputfile;
die "Column number not defined (Use --help for usage)\n" if not defined $column;
die "Threshold not defined (Use --help for usage)\n" if not defined $threshold;
die "Output file not defined (Use --help for usage)\n" if not defined $outputfile;
die "Operator value not valid (Use --help for usage)\n" if not $operator =~ /^(<=|<|>=|>|==|=|eq|regex|match)$/;

#### SET DEFAULT SEPARATOR
print "Using default separator: whitespace\n" if not defined $separator;
$separator = $DEFAULT_SEPARATOR if not defined $separator;

open(FILE, $inputfile) or die "Can't open input file: $inputfile\n";
open(OUTFILE, ">$outputfile") or die "Can't open output file: $outputfile\n";

while ( <FILE> )
{
    next if $_ =~ /^\s*$/;

    my @elements = split $separator, $_;
    die "No elements in line: $_\n" if $#elements < 0;

    #### ARRAY IS ZER0-INDEXED, SO USE COLUMN MINUS 1
    my $value = $elements[$column - 1];
    die "Column $column value not defined\n"  if not defined $value;

    #### USE OPERATOR AND THRESHOLD VALUE TO FILTER LINE
    my $passed = 0;
    if ( $operator eq ">=")
    {
        $passed = 1 if $value >= $threshold;
    }
    elsif ( $operator eq ">" )
    {
        $passed = 1 if $value > $threshold;
    }
    elsif ( $operator eq "<=" )
    {
        $passed = 1 if $value > $threshold;
    }
    elsif ( $operator eq "<" )
    {
        $passed = 1 if $value > $threshold;
    }
    elsif ( $operator eq "==" or $operator eq "=" )
    {
        $passed = 1 if $value == $threshold;
    }
    elsif ( $operator eq "eq" )
    {
        $passed = 1 if $value eq $threshold;
    }
    elsif ( $operator eq "regex" )
    {
		use re 'eval';# EVALUATE $pattern AS REGULAR EXPRESSION
        $passed = 1 if $value =~ /$threshold/;
		no re 'eval';# STOP EVALUATING AS REGULAR EXPRESSION
    }
    else
    {
        print "Operator not supported: $operator\n" and exit;
    }

    if ( $passed )
    {
        print OUTFILE $_;
    }
    else
    {
    }
}
close(FILE);
close(OUTFILE);
print "Output file printed: \n\n$outputfile\n\n";

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

