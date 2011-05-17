#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();

=head2

	APPLICATION     moveFile

    PURPOSE

        MOVE FILES MATCHING A REGEX TO A SPECIFIED DIRECTORY

    INPUT

        1. DIRECTORY CONTAINING FILES

    OUTPUT

        1. FILES MATCHING REGEX MOVED TO OUTPUT DIRECTORY

    ./moveFile.pl  <--inputdir String> <--outputdir String> [--regex] [-h]

        --inputdir		:   /Full/path/to/inputdir 
        --outputdir		:   /Full/path/to/outputdir 
        --regex 		:   Files to be moved must match this Perl regular expression
		--test			:	Do a test run showing files to be moved but without moving
        --help			:   print help info

    EXAMPLES

moveFile.pl \
--inputdir /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/SRX000600,/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/SRX000601,/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/SRX000602,/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/SRX001539,/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/SRX001540 \
--outputdir /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/SRX000600/singles,/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/SRX000601/singles,/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/SRX000602/singles,/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/SRX001539/singles,/scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/SRX001540/singles \
--regex "\d+\.fastq.bz2"

=cut

use strict;

#### USE LIBRARY
use FindBin qw($Bin);
use lib "$Bin/../../../lib";
use lib "$Bin/../../../lib/external";

#### FLUSH BUFFER
$| = 1;

#### INTERNAL MODULES
use Timer;
use Util;
use Conf::Agua;
use ReadPrep;
use Monitor;


#### EXTERNAL MODULES
use Term::ANSIColor qw(:constants);
use Getopt::Long;
use File::Path;

#### GET OPTIONS
my $inputdir;
my $outputdir;
my $regex;
my $test;
my $help;
if ( not GetOptions (
    'inputdir=s' 	=> \$inputdir,
    'outputdir=s' 	=> \$outputdir,
    'regex=s' 		=> \$regex,
    'test' 			=> \$test,
    'help' 			=> \$help
) )
{ print "Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "Input directory not defined (Use --help for usage)\n" if not defined $inputdir;
die "Output directory not defined (Use --help for usage)\n" if not defined $outputdir;
die "Perl regular expression (option --regex) not defined\n" if not defined $regex;

#### GO THROUGH ALL INPUT DIRECTORIES
my @indirs = split ",", $inputdir;
my @outdirs = split ",", $outputdir;
for ( my $i = 0; $i < $#indirs + 1; $i++ )
{
	my $indir = $indirs[$i];
	my $outdir = $outdirs[$i];

	#### CREATE OUTPUT DIRECTORY
	File::Path::mkpath($outdir) or die "Can't create output directory: $outdir\n" if not -d $outdir;

	my $files = Util::files($indir);

	foreach my $file ( @$files )
	{
		next if $file !~ /$regex/;

		my $inputfile = "$indir/$file";
		my $outputfile = "$outdir/$file";
		my $command = "mv $inputfile $outputfile";
		print "moveFile.pl    command: $command\n";
		`$command` if not defined $test;
	}
}

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

