#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();

=head2

	APPLICATION     printSubfile

    PURPOSE

        1. PRINT FASTQ FILES GIVEN .list FILE AND SOURCE FILE

    INPUT

        1. DIRECTORY CONTAINING fasta.gz FILES AND .list FILES

    OUTPUT

        1. FASTA, QUAL OR FASTQ FILES IN OUTPUT DIRECTORY

	USAGE

    ./printSubfile.pl  <--inputfile String> <--outputdir String> <--mode String> <--width width> [--compress] [-h]

        --inputfile		       :   /Full/path/to/inputfile 
        --outputdir		       :   /Full/path/to/outputdir 
        --mode				   :   Type of output file to print (fasta|qual|fastq)
        --compress             :   Compress with gzip or zip
        --dot                  :   Print counter every 'dot' number of records
        --help                 :   print help info

    EXAMPLES


/home/syoung/base/bin/comparison/printSubfile.pl \
--inputfile /mihg/data/NGS/syoung/base/pipeline/SRA/test/fasta/test_1_1.fasta \
--mode fasta \
--width 57 \
--dot 100000


/nethome/syoung/base/bin/comparison/printSubfiles.pl \
--inputdir /mihg/data/NGS/syoung/base/pipeline/SRA/test/fasta,/mihg/data/NGS/syoung/base/pipeline/SRA/test2/fasta \
--outputdir /mihg/data/NGS/syoung/base/pipeline/SRA/test/fasta/samples \
--mode fastq \
--size 10000


/nethome/syoung/base/bin/comparison/printSubfiles.pl \
--inputdir /mihg/data/NGS/syoung/base/pipeline/SRA/test/fasta,/mihg/data/NGS/syoung/base/pipeline/SRA/test2/fasta \
--outputdir /mihg/data/NGS/syoung/base/pipeline/SRA/test/fasta/samples \
--mode fastq \
--size 10000

/nethome/syoung/base/bin/comparison/printSubfiles.pl \
--inputdir /mihg/data/NGS/syoung/base/pipeline/SRA/test/fasta,/mihg/data/NGS/syoung/base/pipeline/SRA/test2/fasta \
--mode fastq \
--size 10000



=cut

use strict;

#### USE LIBRARY
use FindBin qw($Bin);
use lib "$Bin/../../../lib";
use lib "$Bin/../../../lib/external";

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
#my $outputdir;
my $compress;
my $mode;
my $width;
my $tempdir;
#my $paired;
my $dot = 1000000;
my $help;
if ( not GetOptions (
    'inputfile=s' 	=> \$inputfile,
    'compress=s' 	=> \$compress,
    'tempdir=s' 	=> \$tempdir,
    'mode=s' 		=> \$mode,
    'width=i' 		=> \$width,
    'dot=i' 		=> \$dot,
    'help' 			=> \$help
) )
{ print "Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "inputfile not defined (Use --help for usage)\n" if not defined $inputfile;
die "compress must be 'gzip' or 'zip'\n" if defined $compress and $compress !~ /^(gzip|zip)$/;
die "width not defined\n" if not defined $width;

#### CREATE OUTPUT DIRECTORY
print "printSubfile.pl    inputfile: $inputfile\n";

#### GET THE INFO FILES AND TOTAL READ COUNTS FOR EACH DIRECTORY
my $args =
{
	'inputfile'	=>	$inputfile,
	'tempdir'	=>	$tempdir,
	'mode'		=>	$mode,
	'compress'	=>	$compress,
	'width'		=>	$width,
	'dot'		=>	$dot
};

#### SET FILES
Sampler::subfiles($args);

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

