#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();

=head2

	APPLICATION     editTrackInfo

    PURPOSE

		EDIT trackInfo.js BY ADDING 'plugins/view/jbrowse/' TO URL ENTRIES:

		editTrackInfo.pl \
		--path "plugins/view/jbrowse/" \
		--inputfile /nethome/syoung/base/pipeline/jbrowse1/ucsc/chr1/data/trackInfo.js


		TO CHANGE THIS

			  "url" : "data/seq/{refseq}/",

		TO THIS

			  "url" : "plugsin/view/jbrowse/data/seq/{refseq}/",

		WITH THIS FILE FORMAT

		trackInfo = 
		[
		   {
			  "url" : "data/seq/{refseq}/",
			  "args" : {
				 "chunkSize" : 20000
			  },
			  "label" : "DNA",
			  "type" : "SequenceTrack",
			  "key" : "DNA"
		   },


		SO THAT LOADED URLS WILL CHANGE FROM

		http://localhost/agua/0.4/data/tracks/chr1/exon/trackData.json

		TO

		http://localhost/agua/0.4/plugins/view/jbrowse/data/tracks/chr1/exon/trackData.json

    INPUT

		1. ADDITIONAL FILE PATH

        2. INPUT trackInfo.js FILE

    OUTPUT

        1. EDITED trackInfo.js FILE

    USAGE

    ./editTrackInfo.pl  <--inputfile String> <--path String> [--help]

        --inputfile	        :   Name or full path to inputfile directory
        --path	       		:   Name or full path to path directory
        --help              :   print help info

    EXAMPLES

cd /nethome/syoung/base/pipeline/jbrowse1/ucsc/chr1

/data/agua/0.4/bin/apps/editTrackInfo.pl \
--inputfile /nethome/syoung/base/pipeline/jbrowse1/ucsc/chr1/data/trackInfo.js \
--path "plugins/view/jbrowse"

=cut


use strict;

#### USE LIBRARY
use FindBin qw($Bin);
use lib "$Bin/../../../lib";

#### FLUSH BUFFER
$| = 1;

#### INTERNAL MODULES
use Timer;
use Util;
use Conf::Agua;


#### SAVE ARGUMENTS
my @args = @ARGV;

#### EXTERNAL MODULES
use File::Copy;
#use JSON;
use Term::ANSIColor qw(:constants);
use Getopt::Long;
use Data::Dumper;

#### GET OPTIONS
my $inputfile;
my $path;
my $help;
if ( not GetOptions (
    'inputfile=s' => \$inputfile,
    'path=s' => \$path,
	'help' => \$help
) )
{	usage(); exit;	};

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "inputfile not defined (use --help option for usage)\n" if not defined $inputfile;
die "path not defined (use --help option for usage)\n" if not defined $path;

#### GET CHROMOSOME NAME AND LENGTH FROM CONFIG FILE
open(CONF, $inputfile) or die "Can't open config file: $inputfile\n";
$/ = undef;
my $contents = <CONF>;
close(CONF) or die "Can't close config file: $inputfile\n";

#### EDIT URL PATHS
$contents =~ s/"url"\s*:\s*"([^"]+)"/"url": "$path\/$1"/g;

#### PRINT TO TEMP FILE
my $tempfile = $inputfile . "-temp";
open(OUT, ">$tempfile") or die "Can't open temp file: $tempfile\n";
print OUT $contents;
close(OUT);

#### MOVE TO INPUT FILE
File::Copy::move($inputfile, "$inputfile.bak");
File::Copy::move($tempfile, $inputfile);

print "Edited input file printed: $inputfile\n";

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


sub editTrackInfo
{
	my $inputfile	=	shift;
	my $outputfile	=	shift;
	my $featurehash	=	shift;
	my $refseq		=	shift;

	#### GET REFERENCE SEQUENCE INFO FOR FIRST LINE OF OUTPUT FILE	
	my $name = $refseq->{name};
	my $start = $refseq->{start};
	my $end = $refseq->{end};

	#### ADD 1 TO START AND END FOR 1-INDEXED GFF FORMAT
	$start++;
	$end++;

	#### OPEN OUTPUT FILE AND PRINT RUN COMMAND, DATE AND REFERENCE SEQUENCE LINE
	open(OUTFILE, ">$outputfile") or die "Can't open output file: $outputfile\n";
	print OUTFILE "### $0 @args\n";
	print OUTFILE "####", Util::datetime(), "\n\n";
	print OUTFILE "$name\trefseq\trefseq\t$start\t$end\t.\t.\t.\tName=$name\n";

	#### OPEN INPUT FILE
	open(FILE, $inputfile) or die "Can't open input file: $inputfile\n";
	$/ = "\n";
	my $counter = 0;
	print "Doing input file: $inputfile\n";
	while ( <FILE> )
	{
		next if $_ =~ /^\s*$/;
		$counter++;
		if ( $counter % 10000 == 0 ) {	print "$counter\n";	}


		my ($start, $last) = $_ =~ /^(\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+)(.+)$/;

		my @elements = split ";", $last;
		$last = '';
		foreach my $element ( @elements )
		{
			next if $element =~ /^\s*$/;
			$element =~ s/^\s+//;
			$element =~ s/\s+$//;


			my ($key, $value) = $element =~ /^(\S+)\s+(.+)$/;

			if ( $key eq "transcript_id" )
			{
				$value =~ s/"//g;
				$last = "Name=$value";
			}
		}
		$last =~ s/;$//;

		my @fields = split " ", $start;

		#### SET FEATURE
		$fields[2] = $featurehash->{feature};

		#### CORRECT SCORE TO NO DECIMAL PLACES
		$fields[5] =~ s/\.[0]+//g;

		#### ADD LAST ENTRY TO FIELDS
		push(@fields, $last);

		my $line = join "\t", @fields;
		print OUTFILE "$line\n";



	}
	close(FILE);
	close(OUTFILE);
	print "Output file printed:\n\n$outputfile\n\n";


}




sub usage
{
	print `perldoc $0`;
    exit;
}

