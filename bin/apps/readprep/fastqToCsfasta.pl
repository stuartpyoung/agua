#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();

=head2

	APPLICATION     fastqToCsfasta

	VERSION		0.02

	HISTORY

				0.01 BASIC CONVERSION

    PURPOSE

		1. CONVERT FROM FASTQ TO CSFASTA

    INPUTS

        1. INPUT FILE

    OUTPUTS

		1. CSFASTA FILE AND CORRESPONDING QUAL FILE

	NOTES

		MUCH CODE BORROWED FROM HERE:

		http://seqanswers.com/forums/showthread.php?t=2360

    USAGE

    ./fastqToCsfasta.pl <--inputfile String> <format String> [-h]

    --inputfile		:   /full/path/to/inputfile
    --outputfile	:   /full/path/to/outputfile
    --label     	:   Label, e.g. name of experiment or sample
    --help			:   print help info

    EXAMPLES


SINGLE READ
-----------


NO BARCODING

	/nethome/bioinfo/apps/agua/0.5/bin/apps/readprep/fastqToCsfasta.pl \
	--inputfile /FULL/PATH/TO/myfile.csfasta \
	--outputfile /FULL/PATH/TO/myfile.fastq \
	--label myexperiment



PAIRED ENDS
-----------

NO BARCODING, MATE NUMBER IS 1

	/nethome/bioinfo/apps/agua/0.5/bin/apps/readprep/fastqToCsfasta.pl \
	--inputfile /FULL/PATH/TO/myfile.csfasta \
	--outputfile /FULL/PATH/TO/myfile_1.fastq \
	--label myexperiment \
	--matenumber 1

NO BARCODING, MATE NUMBER IS 2 (COMPLEMENTARY TO ABOVE EXAMPLE)

	/nethome/bioinfo/apps/agua/0.5/bin/apps/readprep/fastqToCsfasta.pl \
	--inputfile /FULL/PATH/TO/myfile.csfasta \
	--outputfile /FULL/PATH/TO/myfile_2.fastq \
	--label myexperiment \
	--matenumber 2

=cut

use strict;


#### EXTERNAL MODULES
use FindBin qw($Bin);

#### USE LIBRARY
use lib "$Bin/../../../lib";	
use lib "$Bin/../../../lib/external";	

#### INTERNAL MODULES
use Timer;
use Util;
use Conf::Agua;

#### EXTERNAL MODULES
use Data::Dumper;
use Term::ANSIColor qw(:constants);
use Getopt::Long;

#### GET OPTIONS
use Getopt::Long;
my $inputfile;
my $outputfile;
my $label;
my $barcode = 0;
my $matenumber;
my $dot = 1000000;
my $help;
GetOptions (
    'inputfile=s' 	=> \$inputfile,
    'label=s' 		=> \$label,
    'outputfile=s' 	=> \$outputfile,
    'barcode=i' 	=> \$barcode,
    'matenumber=i' 	=> \$matenumber,
    'dot=i' 		=> \$dot,
    'help' 			=> \$help             
) or die "No options specified. Use --help for usage\n";
usage() if defined $help;

#### CHECK INPUTS
die "inputfile is not defined (use --help for options)\n" if not defined $inputfile;
die "outputfile is not defined (use --help for options)\n" if not defined $outputfile;
die "label is not defined (use --help for options)\n" if not defined $label;

my $qualfile = $outputfile;
$qualfile =~ s/\.csfasta$/_QV.qual/;

#### OPEN CSFASTA FILE
if( $inputfile =~ /\.gz$/ or $inputfile =~ /\.zip$/ )
{
	my $pipe_command = "zcat $inputfile |";
	open(FILE, $pipe_command) or die "Can't open inputfile: $inputfile\n"
}
else
{
	open(FILE, $inputfile) or die "Can't open input file: $inputfile\n";	
}

open(OUT, ">$outputfile") || die "Can't open outputfile: $outputfile\n";
open(QUAL, ">$qualfile") || die "Can't open qualfile: $qualfile\n";
my $state = 0;
my ($n, $r, $q) = ("", "", "");
while(defined(my $line = <FILE>)) {
    chomp($line);
    if(0 == $state) {
        &print_out(\*OUT, \*QUAL, $n, $r, $q);
        $n = $line;
        $n =~ s/^\@/>/;
    }
    elsif(1 == $state) {
        $r = $line;
    }
    elsif(3 == $state) {
        $q = $line;

        # convert back from SANGER phred
        my $tmp_q = "";
        for(my $i=0;$i<length($q);$i++) {
            my $Q = ord(substr($q, $i, 1)) - 33;
            die unless (0 < $Q);
            if(0 < $i) {
                $tmp_q .= " ";
            }
            $tmp_q .= "$Q";
        }
        $q = $tmp_q;
    }
    $state = ($state+1)%4;
}
&print_out(\*OUT, \*QUAL, $n, $r, $q);
close(OUT);
close(FILE);
close(QUAL);

sub print_out {
    my ($OUT, $QUAL, $n, $r, $q) = @_;

    if(0 < length($n)) {
        print $OUT "$n\n$r\n";
        print $QUAL "$n\n$q\n";
    }
}

#### REPORT COMPLETED
print "fastqToCsfasta.pl    outputfile printed:\n\n$outputfile\n\n";

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

