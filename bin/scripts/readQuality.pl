#!/usr/bin/perl -w
use strict;

#### DEBUG

#### TIME
my $time = time();

=head2

    APPLICATION     readQuality.pl

    PURPOSE

        ANALYSE READ QUALITY IN A FILE:

			1. AVERAGE READ QUALITY

    INPUT

        1. INPUT FASTA OR FASTQ FILE

    OUTPUT

        1. LIST OF DIFFERENT READ LENGTH AND NUMBER OF READS 

    USAGE

    ./readQuality.pl <--inputfile String> <--type String> [--help]

    --inputfile        	:   /Full/path/to/first_inputfile
    --type				:   File type (fasta|fastq)
    --help              :   Print help info

    EXAMPLES


perl readQuality.pl --inputfile 6/run12+15-s_2_1.6.fastq 




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

#### CHECK INPUTS
die "Input file not defined (use -h option)\n" if not defined $inputfile;
die "Type not defined (use -h option)" if not defined $type;

#### INSTANTIATE FileTools OBJECT
my $filetools = FileTools->new();
my $read_quality = $filetools->readQuality($inputfile, $type);

print "Read quality stats for $inputfile\n";
print "Number reads   : $read_quality->{number}\n";
print "Average quality: $read_quality->{average}\n";

#my $read_id = '';
#while ( $read_id !~ /^q$/i )
#{
#	$read_id = <STDIN>;
#	if ( $read_id =~ /^q$/i )
#	{
#		last;
#	}
#	else
#	{
#		my $average_quality = $read_quality->{reads}->{$read_id};
#	}
#}

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

