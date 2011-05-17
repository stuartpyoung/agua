#!/usr/bin/perl -w

#### DEBUG
$DEBUG = 1;

#### TIME
my $time = time();

=head2

	APPLICATION     replaceAtQual

	PURPOSE

		REPLACE THE '@' QUALITY VALUE WITH 'A' (ONE VALUE HIGHER)

		BECAUSE IT INTERFERES WITH MAQ CONVERSION FROM solexa TO fastq

    INPUT

        1. solexa-FORMAT FASTQ FILE

    OUTPUT

        1. sanger-FORMAT FASTQ FILE

	USAGE

    ./replaceAtQual.pl  <--inputfile String> <--outputfile String> [--dot Integer] [-h]

        --inputfile		       :   /Full/path/to/inputfile 
        --outputfile		   :   /Full/path/to/outputfile 
        --dot                  :   Print counter every 'dot' number of records
        --help                 :   print help info

    EXAMPLES


perl /nethome/bioinfo/apps/agua/0.4/bin/apps/utils/replaceAtQual.pl \
--inputfile /mihg/data/NGS/syoung/base/pipeline/SRA/NA18507/samples/reads_1.1.fastq \
--outputfile /mihg/data/NGS/syoung/base/pipeline/SRA/NA18507/samples/reads_sanger_1.1.fastq \
--dot 100000

=cut

use strict;

#### USE LIBRARY
use FindBin qw($Bin);
use lib "$Bin/../../../lib";

#### FLUSH BUFFER
$| = 1;

#### GET TEMP DIRECTORY
my $configfile = "$Bin/../../../../conf/default.conf";
my $conf = Conf::Agua->new(inputfile=>$configfile);
my $tempdir = $conf->getKeyValue("agua", 'EXECUTION_TEMPDIR');

#### INTERNAL MODULES
use Sampler;
use Timer;
use Util;
use Conf::Agua;

#### EXTERNAL MODULES
use Term::ANSIColor qw(:constants);
use Getopt::Long;
use IO::Pipe;

#### SET DEFAULT SEPARATOR
my $DEFAULT_SEPARATOR = "\\s+";

#### SAVE ARGUMENTS
my @arguments = @ARGV;

#### GET OPTIONS
my $inputfile;
my $outputfile;
my $dot = 10000;
my $help;
if ( not GetOptions (
    'inputfile=s' => \$inputfile,
    'outputfile=s' => \$outputfile,
    'dot=i' => \$dot,
    'help' => \$help
) )
{ print "Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "Input file not defined (Use --help for usage)\n" if not defined $inputfile;
#die "Output file not defined (Use --help for usage)\n" if not defined $outputfile;

#### GET THE INFO FILES AND TOTAL READ COUNTS FOR EACH DIRECTORY
my $args =
{
	'inputfile'		=>	$inputfile,
	'outputfile'	=>	$outputfile,
	'dot'			=>	$dot
};

#### SET FILES
Sampler::replaceAtQual($args);

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

