#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();

=head2

    APPLICATION     bamIdxstats

    VERSION         0.01

    PURPOSE

		1. INDEX *bam FILE IF NOT INDEXED

        2. PRINT *bam FILE idxstats TO OUTPUT FILE, INCLUDING NUMBER OF HITS/MISSES

    INPUT

        1. *.bam BAM-FORMAT FILE

    OUTPUT

        1. *.bam.idxstats FILE CONTAINING OUTPUT FROM SAMTOOLS idxstats

		0	12
		100	23
		200	45
		300	32
		...

    USAGE

    ./bamIdxstats.pl <--inputfile String> [--help]

		--inputfile		:   Input BAM-format file
		--outputfile	:   Print output to this file (default: STDOUT)
		--begin			:   Begin calculating average from this position
		--end			:   Stop calculating average at this position
		--window		:   Window size (default: 1000bp)
		--help          :   print help info

    EXAMPLES

/nethome/bioinfo/apps/agua/0.5/bin/apps/aligners/bamIdxstats.pl \
--inputfile /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland2/1/chr22/hit.bam \
--outputfile /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland2/1/chr22/hit.bam.idxstats 



=cut

use strict;

#### EXTERNAL MODULES
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use Getopt::Long;
use File::Path;
use FindBin qw($Bin);

#### USE LIBRARY
use lib "$Bin/../../../lib";	

#### INTERNAL MODULES
use Timer;
use Util;
use Conf::Agua;

#### SET MAQ LOCATION
my $conf = Conf::Agua->new(inputfile=>"$Bin/../../../conf/default.conf");
my $samtools = $conf->getKeyValue("applications", 'SAMTOOLS');

#### GET OPTIONS
my $inputfile;
my $outputfile;
my $help;
if ( not GetOptions (

	#### GENERAL
    'inputfile=s' 	=> \$inputfile,
    'outputfile=s' 	=> \$outputfile,
    'help' 			=> \$help
) )
{ print "bamIdxstats.pl    Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
print "bamIdxstats.pl    inputfile not defined (use --help option)\n" and exit if not defined $inputfile;

#### CREATE *.bai INDEX FILE IF NOT PRESENT
my $indexfile = "$inputfile.bai";
if ( not -f $indexfile )
{
	my $command = "$samtools/samtools index $inputfile $indexfile";
	print `$command`;
}

#### OPEN OUTPUTFILE
open(OUT, ">$outputfile") or die "Can't open outputfile: $outputfile" if defined $outputfile;

#### GET THE END OF THE REFERENCE IF NOT DEFINED
my $command = "$samtools/samtools idxstats $inputfile";
my $idxstats = `$command`;
#### FORMAT:
#### chr22   51304566        914819  0
#### *       0       0       0
my ($chromosome, $length, $hitcount, $misscount) = $idxstats =~ /^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;

#### PRINT TO OUTPUTFILE
print "$chromosome\t$length\t$hitcount\t$misscount\n"; 
print OUT "$chromosome\t$length\t$hitcount\t$misscount\n";
close(OUT) or die "Can't close outputfile: $outputfile\n";

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "\nRun time: $runtime\n";
print "Completed $0\n";
print Timer::datetime(), "\n";
print "****************************************\n\n\n";
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

