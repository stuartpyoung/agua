#!/usr/bin/perl -w
use strict;

#### DEBUG

#### TIME
my $time = time();

=head2

    APPLICATION     reduceReads.pl

    PURPOSE

        REDUCE READS TO A USER-DEFINED LENGTH

    INPUT

        1. INPUT FASTA OR FASTQ FILE

    OUTPUT

        1. LIST OF DIFFERENT READ LENGTH AND NUMBER OF READS 

    USAGE

    ./reduceReads.pl <--inputfile String> <--outputfile String> <--type String> <--length Integer> [--help]

    --inputfile        	:   /Full/path/to/inputfile
    --inputfile        	:   /Full/path/to/outputfile
    --type				:   File type (fasta|fastq)
    --length			:   Desired length of output reads
    --help              :   Print help info

    EXAMPLES

perl reduceReads.pl --inputfile 7/run12+15-s_2_1.7.fastq --outputfile 7/run12+15-s_2_1.7-reduced.fastq --type fastq --length 37

=cut

#### FLUSH BUFFER
$| = 1;

#### EXTERNAL MODULES
use FindBin qw($Bin);
use Data::Dumper;
use Term::ANSIColor qw(:constants);
use Getopt::Long;

#### USE LIB
use lib "$Bin/../../lib";

#### INTERNAL MODULES  
use Timer;
use FileTools;

#### GET OPTIONS
my $inputfile;
my $outputfile;
my $type;
my $length;
my $dot = 10000;
my $help;
GetOptions (
    'inputfile=s' => \$inputfile,
    'outputfile=s' => \$outputfile,
    'type=s' => \$type,
    'length=s' => \$length,
    'dot=s' => \$dot,
    'help' => \$help) or die "No options specified. Try '--help'\n";

#### PRINT HELP IF REQUESTED
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "Input file not defined (use -h option)\n" if not defined $inputfile;
die "Type not defined (use -h option)" if not defined $type;

#### INSTANTIATE FileTools OBJECT
my $filetools = FileTools->new();
my $success = $filetools->reduceReads($inputfile, $outputfile, $type, $length, $dot);
if ( $success )
{
	print "Finished printing reduced reads (length $length)\n";
}

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "\nRun time: $runtime\n";
print "Completed $0\n";
print Timer::datetime(), "\n";
print "****************************************\n\n\n";
exit;

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### ####             SUBROUTINES                 #### #### #### ####  
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 


sub usage
{
	print GREEN;
    print `perldoc $0`;
	print RESET;

	exit;
}

