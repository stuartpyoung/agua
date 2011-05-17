#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     bowtieRefList

    PURPOSE

		FOR EACH *.fa FILE IN THE INPUT DIRECTORY:

		1. CREATE A chr*.fai FILE CONTAINING A SAMTOOLS INDEX OF THE

			chr*.fa FILE (REQUIRED FOR CONVERSION OF SAM TO BAM FILES)

		2. RUN EACH JOB IN SERIES LOCALLY OR IN PARALLEL ON THE CLUSTER	

    USAGE

    ./bowtieRefList.pl <--inputdir String> [--help]

    --inputdir           :   Location of directory containing *.fa files
    --help               :   print help info

    EXAMPLES

chmod 755 /nethome/bioinfo/apps/agua/0.4/bin/apps/bowtieRefList.pl

perl /nethome/bioinfo/apps/agua/0.4/bin/apps/bowtieRefList.pl \
--inputdir /nethome/bioinfo/data/sequence/chromosomes/rat/rn4/samtools \
--queue large

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
use Timer;
use TOPHAT;
use Util;
use Conf::Agua;

##### STORE ARGUMENTS TO PRINT TO FILE LATER
my @arguments = @ARGV;
print "bowtieRefList.pl    arguments: @arguments\n";

#### FLUSH BUFFER
$| =1;

#### GET CONF 
my $conf = Conf::Agua->new(inputfile=>"$Bin/../../../conf/default.conf");
my $samtools = $conf->getKeyValue("applications", 'SAMTOOLS');
print "bowtieRefList.pl    bowtie: $samtools\n";

#### GET OPTIONS
my $inputdir;
my $queue;
my $maxjobs = 30;

my $help;
if ( not GetOptions (
    'inputdir=s'   	=> \$inputdir,
    'queue=s'   	=> \$queue,
    'maxjobs=i'   	=> \$maxjobs,
    'help'          => \$help
) )
{ print "Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "inputdir not defined (Use --help for usage)\n" if not defined $inputdir;
die "queue not defined (Use --help for usage)\n" if not defined $queue;
die "maxjobs not defined (Use --help for usage)\n" if not defined $maxjobs;
print "bowtieRefList.pl    inputdir: $inputdir\n";

chdir($inputdir) or die "Can't change to inputdir directory: $inputdir\n";
my @files = <*fa>;

#### INSTANTIATE tophat OBJECT
my $tophat = TOPHAT->new(
	{
		cluster	=>	"LSF",
		queue	=> $queue,
		maxjobs	=> $maxjobs
	} );

#### TRUNCATE inputdir FILES TO CREATE CORRECT STUB IDENTIFIER
my $jobs = [];
foreach my $file ( @files )
{
	my $commands = [];

	#### CHANGE TO DIRECTORY
	push @$commands, "cd $inputdir";

	#### CREATE chr*.fa.fai INDEX FILE WITH SAMTOOLS faidx
	####
	push @$commands, "time $samtools/samtools faidx $file";

	#### REPLACE chr1 WITH chr1.fa IN INDEX FILES TO AVOID THIS ERROR WHEN GENERATING BAM FILE:
	####
	####    [sam_read1] reference 'chr1.fa' is recognized as '*'.
	####

	######## COPY TO 
	####	cp chr1.fa.fai chr1.fai
	my ($stub) = $file =~ /^(.+)\.fa$/;
	my $faifile = $stub . ".fai";
	push @$commands, "cp $file.fai $faifile";


	######## CONVERT WITH sed
	####	FROM:
	####	chr2    242951149       6       50      51
	####	TO:
	####	chr2.fa 242951149       6       50      51
	push @$commands, "sed -e 's/\(chr[A-Z0-9]*\)/\1.fa/' < $faifile > TMP; mv -f TMP  $faifile";


	my $label = "faidx-$stub";
	my $job = $tophat->setJob($commands, $stub, $inputdir);
	push @$jobs, $job;
	#exit;
}

print "Running jobs\n";
$tophat->runJobs($jobs, "FAIDX");
print "Completed jobs\n";


#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "bowtieRefList.pl    Run time: $runtime\n";
print "bowtieRefList.pl    Completed $0\n";
print "bowtieRefList.pl    ";
print Timer::datetime(), "\n";
print "bowtieRefList.pl    ****************************************\n\n\n";
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


__END__

cd /nethome/bioinfo/data/sequence/chromosomes/human-fa

/nethome/syoung/base/apps/samtools/0.1.6/samtools faidx chr1.fa
/nethome/syoung/base/apps/samtools/0.1.6/samtools faidx chr2.fa
/nethome/syoung/base/apps/samtools/0.1.6/samtools faidx chr3.fa
/nethome/syoung/base/apps/samtools/0.1.6/samtools faidx chr4.fa
/nethome/syoung/base/apps/samtools/0.1.6/samtools faidx chr5.fa
/nethome/syoung/base/apps/samtools/0.1.6/samtools faidx chr6.fa
/nethome/syoung/base/apps/samtools/0.1.6/samtools faidx chr7.fa
/nethome/syoung/base/apps/samtools/0.1.6/samtools faidx chr8.fa
/nethome/syoung/base/apps/samtools/0.1.6/samtools faidx chr9.fa
/nethome/syoung/base/apps/samtools/0.1.6/samtools faidx chr10.fa
/nethome/syoung/base/apps/samtools/0.1.6/samtools faidx chr11.fa
/nethome/syoung/base/apps/samtools/0.1.6/samtools faidx chr12.fa
/nethome/syoung/base/apps/samtools/0.1.6/samtools faidx chr13.fa
/nethome/syoung/base/apps/samtools/0.1.6/samtools faidx chr14.fa
/nethome/syoung/base/apps/samtools/0.1.6/samtools faidx chr15.fa
/nethome/syoung/base/apps/samtools/0.1.6/samtools faidx chr16.fa
/nethome/syoung/base/apps/samtools/0.1.6/samtools faidx chr17.fa
/nethome/syoung/base/apps/samtools/0.1.6/samtools faidx chr18.fa
/nethome/syoung/base/apps/samtools/0.1.6/samtools faidx chr19.fa
/nethome/syoung/base/apps/samtools/0.1.6/samtools faidx chr20.fa
/nethome/syoung/base/apps/samtools/0.1.6/samtools faidx chr21.fa
/nethome/syoung/base/apps/samtools/0.1.6/samtools faidx chr22.fa
/nethome/syoung/base/apps/samtools/0.1.6/samtools faidx chrX.fa
/nethome/syoung/base/apps/samtools/0.1.6/samtools faidx chrY.fa


WHICH CREATES FILES LIKE THIS:

cat /nethome/bioinfo/data/sequence/chromosomes/human-fa/chr1.fa.fai

    chr1.fa 247249719       6       50      51

cat /nethome/bioinfo/data/sequence/chromosomes/human-fa/chr2.fa.fai

    chr2    242951149       6       50      51


REPLACE chr1 WITH chr1.fa IN INDEX FILES TO AVOID THIS ERROR WHEN GENERATING BAM FILE:

    [sam_read1] reference 'chr1.fa' is recognized as '*'.


cd /nethome/bioinfo/data/sequence/chromosomes/human-fa
cp chr1.fa.fai chr1.fai
cp chr2.fa.fai chr2.fai
cp chr3.fa.fai chr3.fai
cp chr4.fa.fai chr4.fai
cp chr5.fa.fai chr5.fai
cp chr6.fa.fai chr6.fai
cp chr7.fa.fai chr7.fai
cp chr8.fa.fai chr8.fai
cp chr9.fa.fai chr9.fai
cp chr10.fa.fai chr10.fai
cp chr11.fa.fai chr11.fai
cp chr12.fa.fai chr12.fai
cp chr13.fa.fai chr13.fai
cp chr14.fa.fai chr14.fai
cp chr15.fa.fai chr15.fai
cp chr16.fa.fai chr16.fai
cp chr17.fa.fai chr17.fai
cp chr18.fa.fai chr18.fai
cp chr19.fa.fai chr19.fai
cp chr20.fa.fai chr20.fai
cp chr21.fa.fai chr21.fai
cp chr22.fa.fai chr22.fai
cp chrX.fa.fai chrX.fai
cp chrY.fa.fai chrY.fai


CONVERT WITH sed

cat chr2.fa.fai

    chr2    242951149       6       50      51

sed -e 's/\(chr[A-Z0-9]*\)/\1.fa/' chr2.fa.fai 

    chr2.fa 242951149       6       50      51


sed -e 's/\(chr[A-Z0-9]*\)/\1.fa/' < chr1.fa.fai > TMP; mv -f TMP  chr1.fa.fai
sed -e 's/\(chr[A-Z0-9]*\)/\1.fa/' < chr2.fa.fai > TMP; mv -f TMP  chr2.fa.fai
sed -e 's/\(chr[A-Z0-9]*\)/\1.fa/' < chr3.fa.fai > TMP; mv -f TMP  chr3.fa.fai
sed -e 's/\(chr[A-Z0-9]*\)/\1.fa/' < chr4.fa.fai > TMP; mv -f TMP  chr4.fa.fai
sed -e 's/\(chr[A-Z0-9]*\)/\1.fa/' < chr5.fa.fai > TMP; mv -f TMP  chr5.fa.fai
sed -e 's/\(chr[A-Z0-9]*\)/\1.fa/' < chr6.fa.fai > TMP; mv -f TMP  chr6.fa.fai
sed -e 's/\(chr[A-Z0-9]*\)/\1.fa/' < chr7.fa.fai > TMP; mv -f TMP  chr7.fa.fai
sed -e 's/\(chr[A-Z0-9]*\)/\1.fa/' < chr8.fa.fai > TMP; mv -f TMP  chr8.fa.fai
sed -e 's/\(chr[A-Z0-9]*\)/\1.fa/' < chr9.fa.fai > TMP; mv -f TMP  chr9.fa.fai
sed -e 's/\(chr[A-Z0-9]*\)/\1.fa/' < chr10.fa.fai > TMP; mv -f TMP  chr10.fa.fai
sed -e 's/\(chr[A-Z0-9]*\)/\1.fa/' < chr11.fa.fai > TMP; mv -f TMP  chr11.fa.fai
sed -e 's/\(chr[A-Z0-9]*\)/\1.fa/' < chr12.fa.fai > TMP; mv -f TMP  chr12.fa.fai
sed -e 's/\(chr[A-Z0-9]*\)/\1.fa/' < chr13.fa.fai > TMP; mv -f TMP  chr13.fa.fai
sed -e 's/\(chr[A-Z0-9]*\)/\1.fa/' < chr14.fa.fai > TMP; mv -f TMP  chr14.fa.fai
sed -e 's/\(chr[A-Z0-9]*\)/\1.fa/' < chr15.fa.fai > TMP; mv -f TMP  chr15.fa.fai
sed -e 's/\(chr[A-Z0-9]*\)/\1.fa/' < chr16.fa.fai > TMP; mv -f TMP  chr16.fa.fai
sed -e 's/\(chr[A-Z0-9]*\)/\1.fa/' < chr17.fa.fai > TMP; mv -f TMP  chr17.fa.fai
sed -e 's/\(chr[A-Z0-9]*\)/\1.fa/' < chr18.fa.fai > TMP; mv -f TMP  chr18.fa.fai
sed -e 's/\(chr[A-Z0-9]*\)/\1.fa/' < chr19.fa.fai > TMP; mv -f TMP  chr19.fa.fai
sed -e 's/\(chr[A-Z0-9]*\)/\1.fa/' < chr20.fa.fai > TMP; mv -f TMP  chr20.fa.fai
sed -e 's/\(chr[A-Z0-9]*\)/\1.fa/' < chr21.fa.fai > TMP; mv -f TMP  chr21.fa.fai
sed -e 's/\(chr[A-Z0-9]*\)/\1.fa/' < chr22.fa.fai > TMP; mv -f TMP  chr22.fa.fai
sed -e 's/\(chr[A-Z0-9]*\)/\1.fa/' < chrX.fa.fai > TMP; mv -f TMP  chrX.fa.fai
sed -e 's/\(chr[A-Z0-9]*\)/\1.fa/' < chrY.fa.fai > TMP; mv -f TMP  chrY.fa.fai



cat *fa.fai

    chr10.fa        135374737       7       50      51
    chr11.fa        134452384       7       50      51
    chr12.fa        132349534       7       50      51
    chr13.fa        114142980       7       50      51
    chr14.fa        106368585       7       50      51
    chr15.fa        100338915       7       50      51
    chr16.fa        88827254        7       50      51
    chr17.fa        78774742        7       50      51
    chr18.fa        76117153        7       50      51
    chr19.fa        63811651        7       50      51
    chr1.fa 247249719       6       50      51
    chr20.fa        62435964        7       50      51
    chr21.fa        46944323        7       50      51
    chr22.fa        49691432        7       50      51
    chr2.fa 242951149       6       50      51
    chr3.fa 199501827       6       50      51
    chr4.fa 191273063       6       50      51
    chr5.fa 180857866       6       50      51
    chr6.fa 170899992       6       50      51
    chr7.fa 158821424       6       50      51
    chr8.fa 146274826       6       50      51
    chr9.fa 140273252       6       50      51
    chrX.fa 154913754       6       50      51
    chrY.fa 57772954        6       50      51


