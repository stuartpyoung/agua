#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     simpleVenn

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

    ./simpleVenn.pl <--queryfile String> <--targetfile String>
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

perl /nethome/bioinfo/apps/agua/0.5/bin/apps/aligners/samVenn.pl \
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
use Venn::Index;
use Timer;
use Util;
use Conf::Agua;

##### STORE ARGUMENTS TO PRINT TO FILE LATER
my @arguments = @ARGV;
print "simpleVenn.pl    arguments: @arguments\n";

#### GET OPTIONS
my $outputdir;
my $targetlabel;
my $querylabel;
my $targetfile;
my $queryfile;
my $queryonly;
my $targetonly;
my $both;
my $help;
if ( not GetOptions (
    'outputdir=s'  => \$outputdir,
    'targetlabel=s'  => \$targetlabel,
    'querylabel=s'   => \$querylabel,
    'targetfile=s'  => \$targetfile,
    'queryfile=s'   => \$queryfile,
    'queryonly=s'  	=> \$queryonly,
    'targetonly=s'  => \$targetonly,
    'both=s'  		=> \$both
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

my $venn = Venn::Index->new(
	{
		queryfile => $queryfile,
		targetfile => $targetfile,
		outputdir => $outputdir,
		querylabel => $querylabel,
		targetlabel => $targetlabel
	}
);
$venn->simpleVenn();


#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "simpleVenn.pl    Run time: $runtime\n";
print "simpleVenn.pl    Completed $0\n";
print "simpleVenn.pl    ";
print Timer::datetime(), "\n";
print "simpleVenn.pl    ****************************************\n\n\n";
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


