#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     transcriptsVenn

    PURPOSE

        CREATE VENN DIAGRAM DATA FOR INTERSECTION OF GENE NAMES

		IN FIRST COLUMN OF TWO FILES

	VERSION		0.02

	HISTORY

    INPUTS

        1. QUERY INPUT FILE

        2. TARGET INPUT FILE

    OUTPUTS

        1. COUNTS OF QUERY-ONLY, JOIN AND TARGET-ONLY GENE NAMES 

    USAGE

    ./transcriptsVenn.pl <--queryfile String> <--targetfile String>
			[--queryonly String] [--targetonly String] [--both String]
			[--help]

    --queryfile  	:   Location of query input file
    --targetfile  	:   Location of target input file
    --targetonly  	:   Print to this file entries found in only in target file
    --queryonly  	:   Print to this file entries found in only in query file
    --both  		:   Print to this file entries found in both input files
    --xref  		:   Use this kgXref TSV table data file for knownGene annotations
    --help       	:   print help info

    EXAMPLES

perl /nethome/bioinfo/apps/agua/0.5/bin/apps/transcriptsVenn.pl \
--queryfile /scratch/syoung/base/pipeline/bixby/run1/tophat/analysis1/chrY \
--targetfile /scratch/syoung/base/pipeline/bixby/run1/tophat/analysis2/chrY \
--outputdir /scratch/syoung/base/pipeline/bixby/run1/tophat/venn \
--queryonly analysis1-only \
--targetonly analysis2-only \

=cut

use strict;

#### EXTERNAL MODULES
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use Getopt::Long;
use FindBin qw($Bin);
use File::Path;
use File::Copy;

#### USE LIBRARY
use lib "$Bin/../../../lib";	
use lib "$Bin/../../../lib/external";	

#### INTERNAL MODULES
use Venn;
use Timer;
use Util;
use Conf::Agua;

##### STORE ARGUMENTS TO PRINT TO FILE LATER
my @arguments = @ARGV;
print "transcriptsVenn.pl    arguments: @arguments\n";

#### GET OPTIONS
my $outputdir;
my $targetlabel;
my $querylabel;
my $targetfile;
my $queryfile;
my $queryonly;
my $targetonly;
my $both;
my $xref;
my $help;
if ( not GetOptions (
    'outputdir=s'  => \$outputdir,
    'targetlabel=s'  => \$targetlabel,
    'querylabel=s'   => \$querylabel,
    'targetfile=s'  => \$targetfile,
    'queryfile=s'   => \$queryfile,
    'queryonly=s'  	=> \$queryonly,
    'targetonly=s'  => \$targetonly,
    'both=s'  		=> \$both,
    'xref=s'  		=> \$xref
) )
{ print "Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "outputdir not defined (Use --help for usage)\n" if not defined $outputdir;
die "targetfile not defined (Use --help for usage)\n" if not defined $targetfile;
die "queryfile not defined (Use --help for usage)\n" if not defined $queryfile;
die "targetlabel not defined (Use --help for usage)\n" if not defined $targetlabel;
die "querylabel not defined (Use --help for usage)\n" if not defined $querylabel;

#### CREATE OUTPUT DIR IF NOT EXISTS
print "outputdir is a file: $outputdir\n" and exit if -f $outputdir;
File::Path::mkpath($outputdir) if not -d $outputdir;
print "Can't create outputdir: $outputdir\n" and exit if not -d $outputdir;

#### DEBUG
print "transcriptsVenn.pl    outputdir: $outputdir\n";
print "transcriptsVenn.pl    targetfile: $targetfile\n";
print "transcriptsVenn.pl    queryfile: $queryfile\n";

my $venn = Venn->new();

#### GENERATE TRANSCRIPT-ONLY INPUT FILES
my $query_transcript_file = $venn->filterTranscripts(	{ inputfile => $queryfile } );
my $target_transcript_file = $venn->filterTranscripts(	{ inputfile => $targetfile } );

#### SET OUTPUT FILES
($queryonly, $targetonly, $both) = $venn->setOutputFiles($outputdir, $querylabel, $targetlabel, $queryonly, $targetonly, $both);

#### DO QUERY AND BOTH
my $querycount = 0;
my $bothcount = 0;

my ($query_only_lines, $target_only_lines, $query_both_lines, $target_both_lines);

##### DO QUERY
($query_both_lines, $query_only_lines) = $venn->match($query_transcript_file, $target_transcript_file);


#### DO TARGET
($target_both_lines, $target_only_lines) = $venn->match($target_transcript_file, $query_transcript_file);

#


@$target_both_lines = sort @$target_both_lines;
@$query_both_lines = sort @$query_both_lines;

#

print "transcriptsVenn.pl    target both != query both\n" if @$target_both_lines != @$query_both_lines;


#### PRINT OUTPUT FILES
open(QUERYONLY, ">$queryonly") or die "Can't open queryonly file: $queryonly\n";
foreach my $line ( @$query_only_lines )	{	print QUERYONLY $line; }
close(QUERYONLY);
open(TARGETONLY, ">$targetonly") or die "Can't open targetonly file: $targetonly\n";
foreach my $line ( @$target_only_lines )	{	print TARGETONLY $line; }
close(TARGETONLY);
open(BOTH, ">$both") or die "Can't open both file: $both\n";
foreach my $line ( @$query_both_lines )	{	print BOTH $line; }
close(BOTH);



if ( defined $xref )
{
	open(XREF, $xref) or die "Can't open xref file: $xref\n";
	my $annolines;
	@$annolines = <XREF>;
	close(XREF);
	my $annohash = {};
	foreach my $anno ( @$annolines )
	{
		my ($id, $annotation) = $anno =~ /^(\S+)\s+(.+)$/;
		$annohash->{$id} = $annotation;
	}


#exit;
	#$query_only_lines = $venn->getAnnotation($query_only_lines, $annohash);
	#$query_both_lines = $venn->getAnnotation($query_both_lines, $annohash);
	#$target_only_lines = $venn->getAnnotation($target_only_lines, $annohash);


	#### PRINT OUTPUT FILES
	open(QUERYONLY, ">$queryonly-ANNOTATED") or die "Can't open queryonly file: $queryonly-ANNOTATED\n";
	foreach my $line ( @$query_only_lines )
	{
		my $summary = $venn->gtfSummary($line);
		my ($id) = $summary =~ /(\S+)$/;
		my $annotation = $annohash->{$id};
		$summary.= "\t$annotation" if defined $annotation;
#exit;
		print QUERYONLY $summary, "\n";
	}
	close(QUERYONLY);

	open(TARGETONLY, ">$targetonly-ANNOTATED") or die "Can't open targetonly file: $targetonly-ANNOTATED\n";
	foreach my $line ( @$target_only_lines )
	{
		my $summary = $venn->gtfSummary($line);
		my ($id) = $summary =~ /(\S+)$/;
		my $annotation = $annohash->{$id};
		$summary.= "\t$annotation" if defined $annotation;

		print TARGETONLY $summary, "\n";
	}
	close(TARGETONLY);

	open(BOTH, ">$both-ANNOTATED") or die "Can't open both file: $both-ANNOTATED\n";
	foreach my $line ( @$query_both_lines )
	{
		my $summary = $venn->gtfSummary($line);
		my ($id) = $summary =~ /(\S+)$/;
		my $annotation = $annohash->{$id};
		print BOTH $summary, "\n";
	}
	close(BOTH);
}

#### PRINT INFO
print "transcriptsVenn.pl    files printed:\n";
print "query_only: $queryonly\n";
print "target_only: $targetonly\n";
print "both: $both\n";
print "\n";


#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "transcriptsVenn.pl    Run time: $runtime\n";
print "transcriptsVenn.pl    Completed $0\n";
print "transcriptsVenn.pl    ";
print Timer::datetime(), "\n";
print "transcriptsVenn.pl    ****************************************\n\n\n";
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


