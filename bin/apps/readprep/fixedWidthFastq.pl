#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();

=head2

	APPLICATION     fixedWidthFastq

    PURPOSE

        1. CONVERT A SOLEXA-FORMAT FASTQ FILE OR ORDINARY FASTQ FILE INTO A FASTA FILE

    INPUT

		1. SOLEXA-FORMAT FASTQ (SYMBOLIC QUALITY) OR ORDINARY FASTQ FILE (NUMERIC QUALITY)

    OUTPUT

        1. FASTA FILE

    NOTES

    USAGE

    ./fixedWidthFastq.pl <--inputfile String> <--outputfile String> [--limit Integer] [--compress] [--fixed] [--id Integer] [--length Integer] [--compress String] [--dot Integer] [-h]

    --inputfile             :   /full/path/to/inputfile
    --outputfile            :   /full/path/to/outputfile
    --limit                 :   Limit output to this number of lines
    --compress              :   Compress with gzip or zip
    --fixed                 :   Print FASTA record of fixed length
    --id             :   Fixed id length
    --length       :   Fixed sequence length
    --dot                   :   Print counter every 'dot' number of records
    --help                  :   print help info

    EXAMPLES


/nethome/bioinfo/apps/agua/0.5/bin/apps/readprep/fixedWidthFastq.pl --inputfile \
/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/SRX000600/SRR004850_1.fastq.gz \
--outputfile /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/min3.length26/SRX000600/fasta/SRR004850_1.fasta \
--dot 100000 \
--compress gzip \
--id 50 \
--length 26



=cut

use strict;


#### EXTERNAL MODULES
use FindBin qw($Bin);

#### USE LIBRARY
use lib "$Bin/../../../lib";	
use lib "$Bin/../../../lib/external";	

#### INTERNAL MODULES
use Timer;
use SolexaUtil;
use Util;
use Conf::Agua;

#### INITIALISE SolexaUtil OBJECT
my $solexa = SolexaUtil->new();

#### EXTERNAL MODULES
use Data::Dumper;
use Term::ANSIColor qw(:constants);
use Getopt::Long;

#### FIXED WIDTH SETTINGS
#my $DEFAULT_ID_LENGTH = 20;
#my $DEFAULT_SEQUENCE_LENGTH = 40;

#### GET OPTIONS
use Getopt::Long;
my $inputfile;
my $outputfile;
my $limit;
my $compress;
my $dot = 1000000;
my $id;
my $length;
my $help;
GetOptions (
    'inputfile=s'	=> \$inputfile,
    'outputfile=s'	=> \$outputfile,
    'limit=i' 		=> \$limit,
    'compress=s' 	=> \$compress,
    'dot=i' 		=> \$dot,
    'id=i' 			=> \$id,
    'length=i' 		=> \$length,
    'help' 			=> \$help             
) or die "No options specified. Use --help for usage\n";

usage() if defined $help;

#### CHECK INPUTS
die "Input file not defined (use --help for usage)\n" if not defined $inputfile;
die "Output file not defined (use --help for usage)\n" if not defined $outputfile;
die "Compress type must be 'gzip' or 'zip'\n" if defined $compress and $compress !~ /^(gzip|zip)$/;
die "Id length not defined\n" if not defined $id;
die "Sequence length not defined\n" if not defined $length;

#### CHECK OUTPUT FILE
if ( not defined $outputfile )
{
    print "Output file not defined. (option -o)\n";
    usage();
}

#### OPEN OUTPUT FILE
open(OUTFILE, ">$outputfile") or die "Can't open output file: $outputfile\n";

#### OPEN FILE AND SET RECORD SEPARATOR
if( $inputfile =~ /(\.gz|\.zip)$/ )
{
    my $pipe_command = "zcat $inputfile |";
    open(FILE, $pipe_command);
}
elsif( $inputfile =~ /\.bz2$/ )
{
    my $pipe_command = "bzip2 -dc $inputfile |";
    open(FILE, $pipe_command);
}
else
{
    open(FILE, $inputfile) or die "Can't open input file: $inputfile\n";
}

#### READ THROUGH ONE SEQUENCE RECORD AT A TIME
my $record;
$/ = "\n";
my $counter = 0;
while ( <FILE> )
{    
	print "$counter\n" if $counter % $dot == 0;

	my $sequence_header = $_;
	my $sequence = <FILE>;

	#### SCRAP QUALITY HEADER IN CASE ITS EMPTY AND USE SEQUENCE HEADER INSTEAD
	<FILE>;
	my $quality = <FILE>;
	my $quality_header = $sequence_header;
	$quality_header =~ s/^@/+/;

	$sequence =~ s/\s+//g;
	$quality =~ s/\s+//g;

	#### REMOVE 'length=...' 
	$sequence_header =~ s/\s+length=.+$//;
	$quality_header =~ s/\s+length=.+$//;

	#### REMOVE WHITESPACE
	$sequence_header =~ s/\s+$//;
	$quality_header =~ s/\s+$//;

	#### WHITE SPACES LESS ONE SPACE BECAUSE OF THE ">"
    $sequence_header = substr($sequence_header, 0, $id);
    $sequence_header = $sequence_header . " " x ($id - length($sequence_header)); 

    $quality_header = substr($quality_header, 0, $id);
    $quality_header = $quality_header . " " x ($id - length($quality_header));


	#### CHECK LENGTH OF SEQUENCE AND QUALITY
	next if not defined $sequence or not defined $quality;

	$sequence = substr($sequence, 0, $length);
	$quality = substr($quality, 0, $length);

	print OUTFILE "$sequence_header\n$sequence\n$quality_header\n$quality\n";

    $counter++;
}
close(OUTFILE);

#### COMPRESS IF DEFINED
$outputfile = compress($outputfile, $compress) if defined $compress;

#### REPORT OUTPUT FILES PRINTED
print "Output files printed:\n\n$outputfile\n\n";

#### PRINT INFO FILE :
#### NUMBER_SEQUENCES   ID_LENGTH   SEQUENCE_LENGTH RECORD_LENGTH
#### ADD ONE TO FIXED WIDTH TO ACCOUNT FOR TWO LINE \n LINE ENDS
my $fixed_width = ($id + $length + 2) * 2;    
print "fixed width: $fixed_width\n";

my $infofile = $outputfile . ".info";
open(OUT, ">$infofile") or die "Can't open INFO file for output: $infofile\n";
print OUT "$counter\t$id\t$length\t$fixed_width";
close(OUT);
print "INFO file printed:\n\n$infofile\n\n";


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

	SUBROUTINE		compress

	PURPOSE



=cut

sub compress
{
	my $filename	=	shift;
	my $compress	=	shift;

	return $filename if not defined $compress;	
	print "fixedWidthFastq.pl    compress type not supported: $compress\n" if $compress !~ /^(gzip|zip|bz2)$/;

	my $compressfile = $filename;
	$compressfile = $filename . ".gz" if $compress =~ /^gzip$/;
	$compressfile = $filename . ".zip" if $compress =~ /^zip$/;
	$compressfile = $filename . ".bz2" if $compress =~ /^bz2$/;
	`rm -fr $compressfile`;
	my $command = "$compress $filename";
	print "$command\n";
	print `$command`;

	return $compressfile;
}


=head2

    SUBROUTINE      fixed_width

    PURPOSE

        MAKE RECORD OF FIXED LENGTH

=cut

sub fixed_width
{
    my $record    		=   shift;
	my $idlength		=	shift;
	my $sequencelength	=	shift;


    my ($id, $sequence) = $record =~ /^([^\n]+)\n(.+)$/;

    $id = substr($id, 0, $idlength);

	#### WHITE SPACES LESS ONE SPACE BECAUSE OF THE ">"
    $id = $id . " " x ($idlength - length($id)); 

    $sequence = substr($sequence, 0, $sequencelength);
    $sequence = $sequence . " " x ($sequencelength - length($sequence));

    return "$id\n$sequence";
}

sub usage
{
	print `perldoc $0`;

	exit;
}
