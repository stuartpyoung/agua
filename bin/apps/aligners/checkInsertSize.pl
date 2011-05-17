#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     checkInsertSize

    PURPOSE

        CALCULATE THE AVERAGE INSERT SIZE BASED ON ELAND ALIGNMENTS OF PAIRED END READS

    USAGE

    ./checkInsertSize.pl <--inputfile String> [--outputdir String] [--help]

    --inputfile		:   Location of ELAND first mate alignment file ('myTest_1_export.txt')
    --matefile		:   Location of ELAND second mate alignment file ('myTest_2_export.txt')
    --outputdir		:   Location to print insert size report file ('myTest.inserts.txt')
	--genome		:	Location of squashed ELAND reference genome
	--reads			:	Number of reads to use to estimate insert size (DEFAULT = 10000)
    --help      	:   print help info

	INPUT

		1. EXTRACT LINES FROM INPUT FASTQ PAIRED END FILES AND RUN ELAND_standalone.pl

		2. RUN insertSize.pl ON OUTPUT *export.txt FILES TO DETERMINE INSERT SIZES

	OUTPUT

		1. PRINT TO OUTPUT FILE:

			-	MIN, MAX AND AVERAGE INSERT SIZE

			-	DISTRIBUTION OF INSERT SIZES

    EXAMPLES

/nethome/bioinfo/apps/agua/0.5/bin/apps/checkInsertSize.pl \
--inputfile /scratch/syoung/base/pipeline/jvance/pdx/control/1/s_1_1.nfilter.sequence.txt \
--matefile /scratch/syoung/base/pipeline/jvance/pdx/control/1/s_1_2.nfilter.sequence.txt \
--genome /nethome/bioinfo/data/sequence/chromosomes/human/hg19/eland \
--outputdir /scratch/syoung/base/pipeline/jvance/pdx/control/1/eland \
--reads 1000000

=cut

use strict;

#### EXTERNAL MODULES
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use Getopt::Long;
use FindBin qw($Bin);
use File::Path;

#### USE LIBRARY
use lib "$Bin/../../../lib";
use lib "$Bin/../../../lib/external";

#### INTERNAL MODULES
use BinData;
use Timer;
use Util;
use Conf::Agua;

##### STORE ARGUMENTS TO PRINT TO FILE LATER
my @arguments = @ARGV;

#### GET DIRECTORY OF ELAND_standalone.pl
my $configfile = "$Bin/../../../conf/default.conf";
print "Can't find configfile: $configfile\n" if not -f $configfile;
print "configfile is empty: $configfile\n" if -z $configfile;
my $conf = Conf::Agua->new(inputfile=>$configfile);
my $casava = $conf->getKeyValue("applications", 'CASAVA');
print "casava is not defined in configfile: $configfile\n" if not defined $casava;

#### FLUSH BUFFER
$| = 1;

#### GET OPTIONS
my $inputfile;
my $matefile;
my $genome;
my $outputdir;
my $reads = 10000;
my $help;
if ( not GetOptions (
    'inputfile=s'	=> \$inputfile,
    'matefile=s'	=> \$matefile,
    'genome=s'		=> \$genome,
    'outputdir=s'	=> \$outputdir,
    'reads=i'		=> \$reads,
    'help'			=> \$help
) )
{ print "Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "inputfile not defined (Use --help for usage)\n" if not defined $inputfile;
die "outputdir not defined (Use --help for usage)\n" if not defined $outputdir;

#### SET MATE FILE
$matefile = $inputfile if not defined $matefile;
$matefile =~ s/_1\./_2\./;
die "Can't set matefile correctly: $matefile" if $matefile eq $inputfile;

#### CHECK FILES PRESENT AND NON-EMPTY
die "Can't find inputfile: $inputfile\n" if not -f $inputfile;
die "Can't find matefile: $matefile\n" if not -f $matefile;
die "inputfile is empty: $inputfile\n" if -z $inputfile;
die "matefile is empty: $matefile\n" if -z $matefile;

#### CHECK DIRECTORIES
die "Can't find genome directory: $genome\n" if not -d $genome;
File::Path::mkpath($outputdir) if not -d $outputdir;
print "Can't create outputdir: $outputdir\n" if not -d $outputdir;

#### CHECK insertSize.pl EXECUTABLE
my $executable = "$Bin/insertSize.pl";
die "Can't find executable: $executable\n" if not -f $executable;
die "executable is empty: $executable\n" if -z $executable;

#### SET FILES TO INPUT TO ELAND
my ($file) = $inputfile =~ /([^\/]+)$/;
my $infile = "$outputdir/$file";
my $mate = $infile;
$mate =~ s/_1\./_2\./;

#### CHOP FILES BY SPECIFIED NUMBER OF LINES
chopfile($inputfile, $infile, $reads);
chopfile($matefile, $mate, $reads);
print "infile: $infile\n";
print "mate: $mate\n";

#### SET ELAND COMMAND
chdir($outputdir) or die "Can't chdir to outputdir: $outputdir\n";

my $command = qq{time $casava/ELAND_standalone.pl \\
--input-type fastq \\
--eland-genome $genome \\
--input-file $infile \\
--input-file $mate};
print "command: $command\n";
print `$command`;


my $outputfile = "$outputdir/insertsize.txt";
my $run = qq{$executable \\
--inputfile $infile \\
--matefile $mate \\
--outputfile $outputfile \\
--lines $reads};
print "run: $run\n";
print `$run`;


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

=head2

	SUBROUTINE		chopfile

	PURPOSE

		PRINT A SPECIFIED NUMBER OF LINES FROM A FILE TO ANOTHER FILE

=cut

sub chopfile
{
	my $inputfile	=	shift;
	my $outputfile	=	shift;
	my $reads		=	shift;

	open(FILE, "$inputfile") or die "Can't open inputfile: $inputfile\n";
	open(OUT, ">$outputfile") or die "Can't open outputfile: $outputfile\n";

	my $lines_per_read = 4;
	my $counter = 0;
	$/ = "\n";
	my $line = <FILE>;
	while ( defined $line and $counter < ($reads * $lines_per_read))
	{
		print OUT $line;
		$line = <FILE>;
		$counter++;
	}
	close(FILE) or die "Can't close inputfile: $inputfile\n";
	close(OUT) or die "Can't close outputfile: $outputfile\n";
}



sub usage
{
	print GREEN;
	print `perldoc $0`;
	print RESET;
	exit;
}


