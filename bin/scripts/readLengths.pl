#!/usr/bin/perl -w
use strict;

#### DEBUG

#### TIME
my $time = time();

=head2

    APPLICATION     readLengths.pl

    PURPOSE

        USES FileTools::readLengths TO CHECK READ LENGTHS

			1. PRINT LIST OF DIFFERENT READ LENGTHS

				WITH NUMBER OF READS OF EACH LENGTH

    INPUT

        1. INPUT FASTA OR FASTQ FILE

    OUTPUT

        1. LIST OF DIFFERENT READ LENGTH AND NUMBER OF READS 

    USAGE

    ./readLengths.pl <--inputfile String> [--help]

    --inputfile         :   /Full/path/to/input/directory
	--type				:	fasta|fastq

    --help              :   Print help info

    EXAMPLES

perl readLengths.pl --inputfile /mihg/data/NGS/syoung/base/pipeline/run12+15/1.4.0/lane1/3/run12+15-s_1_1.3.bfq --type fastq


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
my $type;
my $help;
GetOptions (
    'inputfile=s' => \$inputfile,
    'type=s' => \$type,
    'help' => \$help) or die "No options specified. Try '--help'\n";

#### PRINT HELP IF REQUESTED
if ( defined $help )	{	usage();	}

#### INSTANTIATE FileTools OBJECT
my $filetools = FileTools->new();
my $readlengths = $filetools->readLengths($inputfile, $type);

print "Length\tSequences\n";
foreach my $key ( keys %$readlengths )
{
	print "$key\t$readlengths->{$key}\n";
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

