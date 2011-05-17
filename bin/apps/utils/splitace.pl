#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();



my $delay = 30;
print "Sleeping $delay seconds...\n";
sleep($delay);


=head2

	APPLICATION     splitace



        **** DUMMY EXECUTABLE TO TEST WORKFLOW ****



    PURPOSE

        1. SPLIT A LARGE .ace FILE INTO SMALLER .ace FILES FOR EACH CONTIG

        2. WRITE THE INDIVIDUAL CONTIG .ace FILE IN NUMBERED SUBDIRECTORIES

            BASED ON THE NUMBER OF THE CONTIG (E.G., THE .ace FILES FOR CONTIGS

            1 TO 100 GO IN SUBDIRECTORY '100', THE .ace FILES FOR CONTIGS

            701-800 GO IN SUBDIRECTORY '800', ETC.)

    INPUT

        1. AN .ace FILE

    OUTPUT

        1. MULTIPLE .ace FILES - ONE FOR EACH CONTIG IN THE INPUT .ace FILE

    USAGE

    ./splitace.pl <-i inputfile> [-h]

    -i inputfile            :   /full/path/to/inputfile
    -i outputdir     :   /full/path/to/outputdir
    -h help                 :   print help info

    EXAMPLES

./splitace.pl -i /home/syoung/base/pipeline/run2-lane6-mtdna-velvet/data/s_6_1_sequence.ace -o /home/syoung/base/pipeline/run2-lane6-mtdna-velvet/assembly/acefiles

=cut

use strict;

#### USE LIBRARY
use FindBin qw($Bin);
use lib "$Bin/../../../lib";

#### INTERNAL MODULES
use Timer;
use Util;
use Conf::Agua;

#### EXTERNAL MODULES
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use Getopt::Long;

my @arguments = @ARGV;


#### GET OPTIONS
my $inputfile;
my $outputdir;
my $help;
if ( not GetOptions (
    'inputfile=s' => \$inputfile,
    'outputdir=s' => \$outputdir,
    'help' => \$help
) )
{ print "Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

print "Inputfile: $inputfile\n";
print "Outputdir: $outputdir\n";

#### CHECK INPUTFILE
if ( not defined $inputfile )
{
    print "Input file not defined (option --inputfile)\n";
    usage();
}
if ( not -f $inputfile )
{
    die "Could not find input file: $inputfile\n";
    usage();
}


#### CHECK OUTPUT DIRECTORY
if ( not defined $outputdir )    {    print "Output directory not defined (option --outputdir)\n"; usage();	}
if ( not -d $outputdir )
{
    mkdir($outputdir) or die "Can't make directory: $outputdir\n";
    if ( not -d $outputdir )
    {
        die "Could not create output directory: $outputdir\n";
    }
}

#### OPEN INPUT FILE
open(FILE, $inputfile) or die "Can't open input file: $inputfile\n";

#### SET RECORD DIVIDER
$/ = "\nCO ";

#### CHECK FIRST LINE TO MAKE SURE ITS AN .ace FILE 
my $line = <FILE>;
if ( $line !~ /^AS\s+\d+\s+\d+\s*$/ms )
{
    die "Input file does not have .ace-format ('AS <no. contigs> <no. bases>') first line: $line";
}

while ( <FILE> )
{

    my ($contig_number, $number_reads, $contig_length) = $_ =~ /^(\S+)\s+(\d+)\s+(\d+)/i;

	#### REMOVE ANY INFO AFTER THE FIRST "|" BAR IN THE CONTIG NUMBER
	$contig_number =~ s/\|.+$//;

	#### REMOVE ANY LEADING NON-NUMERIC SYMBOLS IN CONTIG NUMBER
	$contig_number =~ s/^(\D+)//;

    my $output_subdir = int($contig_number / 100) + 1;

    $output_subdir = $output_subdir . "00";
    if ( not $output_subdir )
    {
        $output_subdir = "000";
    }

    print "Output subdir: $output_subdir\n";
    my $output_dirpath = "$outputdir/$output_subdir";
    print "output_dirpath: $output_dirpath\n";

    if ( not -d $output_dirpath )
    {
        mkdir($output_dirpath) or die "Can't create directory: $output_dirpath\n";
        if ( not -d $output_dirpath )
        {
            die "Could not create output subdirectory: $output_dirpath";
        }
    }

	if ( $^O =~ /^MSWin32$/ )   {   $output_dirpath =~ s/\//\\/g;  }

    print "output_dirpath: $output_dirpath\n";


    #### SET OUTPUT FILE
    my $outputfile = "$output_dirpath/contig.$contig_number.ace";
	if ( $^O =~ /^MSWin32$/ )   {   $outputfile =~ s/\//\\/g;  }
    print "Outputfile: $outputfile\n";


    #### REMOVE ENDING 'CO' FROM OUTPUT
    $_ =~ s/CO\s*$//;

    #### OPEN OUTPUT FILE
    open(OUTFILE, ">$outputfile") or die "Can't open output file: $outputfile\n";
    print OUTFILE "AS $number_reads $contig_length\n\n";
    print OUTFILE "CO Contig";
    print OUTFILE $_;
    close(OUTFILE);
    #`cat $outputfile`;
}


#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "Run time: $runtime\n";
print "Completed $0\n";
print Util::datetime(), "\n";
print "****************************************\n\n\n";
exit;

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#									SUBROUTINES
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


sub usage
{
	print GREEN <<"EOF";

	APPLICATION     splitace



        **** DUMMY EXECUTABLE TO TEST WORKFLOW ****



    PURPOSE

        1. SPLIT A LARGE .ace FILE INTO SMALLER .ace FILES FOR EACH CONTIG

        2. WRITE THE INDIVIDUAL CONTIG .ace FILE IN NUMBERED SUBDIRECTORIES

            BASED ON THE NUMBER OF THE CONTIG (E.G., THE .ace FILES FOR CONTIGS

            1 TO 100 GO IN SUBDIRECTORY '100', THE .ace FILES FOR CONTIGS

            701-800 GO IN SUBDIRECTORY '800', ETC.)

    INPUT

        1. AN .ace FILE

    OUTPUT

        1. MULTIPLE .ace FILES - ONE FOR EACH CONTIG IN THE INPUT .ace FILE

    USAGE

    ./splitace.pl <-i inputfile> [-h]

    -i inputfile            :   /full/path/to/inputfile
    -i outputdir     :   /full/path/to/outputdir
    -h help                 :   print help info

    EXAMPLES

./splitace.pl -i /home/syoung/base/pipeline/run2-lane6-mtdna-velvet/data/s_6_1_sequence.ace -o /home/syoung/base/pipeline/run2-lane6-mtdna-velvet/assembly/acefiles

=cut

EOF

	print RESET;

	exit(1);
}

