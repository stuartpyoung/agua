#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     faToBfa

    PURPOSE

        CONVERT FROM .fa FASTA FORMAT TO .bfa BINARY FASTA FORMAT (MAQ)

    INPUT

        1. INPUT DIRECTORY CONTAINING .fa FILES

        2. OUTPUT DIRECTORY TO PRINT .bfa FILES

    OUTPUT

        1. MAQ .bfa FILES IN OUTPUT DIRECTORY

    USAGE

    ./faToBfa.pl <--inputdir String> <--outputdir String> [--help]

    --inputdir     :   Single FASTQ sequence file
    --outputdir    :   Create this directory and write output files to it
    --help         :   print help info

	EXAMPLES

/nethome/bioinfo/apps/agua/0.5/bin/apps/converters/faToBfa.pl \
--inputdir /nethome/bioinfo/data/sequence/chromosomes/mouse/mm9/fasta \
--outputdir /nethome/bioinfo/data/sequence/chromosomes/mouse/mm9/maq

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
use MAQ;
use Timer;
use Util;
use Conf::Agua;

##### STORE ARGUMENTS TO PRINT TO FILE LATER
my @arguments = @ARGV;
unshift @arguments, $0;

#### FLUSH BUFFER
$| =1;

#### SET MAQ LOCATION
my $conf = Conf::Agua->new(inputfile=>"$Bin/../../../../conf/default.conf");
my $maq = $conf->getKeyValue("agua", 'MAQ');

#### GET OPTIONS
my $inputdir;
my $outputdir;
my $help;
print "faToBfa.pl    Use option --help for usage instructions.\n" and exit if not GetOptions (	
    'inputdir=s' 	=> \$inputdir,
    'outputdir=s'	=> \$outputdir,
    'help' 			=> \$help
);

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "inputdir not defined (Use --help for usage)\n" if not defined $inputdir;
die "outputdir not defined (Use --help for usage)\n" if not defined $outputdir;

#### MAKE OUTPUT DIR IF NOT EXISTS
File::Path::mkpath($outputdir) if not -d $outputdir;
print "Could not create output directory: $outputdir\n" if not -d $outputdir;

#### INSTANTIATE MAQ OBJECT
my $runMaq = MAQ->new( { maq => $maq } );
$runMaq->faToBfa($inputdir, $outputdir);

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "faToBfa.pl    Run time: $runtime\n";
print "faToBfa.pl    Completed $0\n";
print "faToBfa.pl    ";
print Timer::current_datetime(), "\n";
print "faToBfa.pl    ****************************************\n\n\n";
exit;

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#									SUBROUTINES
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


sub usage
{
	print GREEN;
	print `perldoc $0`;
	print RESET;
	exit;
}

