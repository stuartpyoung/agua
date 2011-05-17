#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     snpVenn

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

    ./snpVenn.pl <--queryfile String> <--targetfile String>
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

perl /nethome/bioinfo/apps/agua/0.5/bin/apps/snpVenn.pl \
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
use Venn::Snp;
use Timer;
use Util;
use Conf::Agua;

##### STORE ARGUMENTS TO PRINT TO FILE LATER
my @arguments = @ARGV;

#### GET OPTIONS
my $outputdir;
my $prefix;
my $suffix;
my $targetfile;
my $queryfile;
my $targetlabel;
my $querylabel;
my $queryonly;
my $targetonly;
my $both;
my $stdout;
my $help;
if ( not GetOptions (
    'outputdir=s'	=> \$outputdir,
    'prefix=s'  	=> \$prefix,
    'suffix=s'  	=> \$suffix,
    'targetfile=s'  => \$targetfile,
    'queryfile=s'   => \$queryfile,
    'targetlabel=s' => \$targetlabel,
    'querylabel=s'  => \$querylabel,
    'queryonly=s'  	=> \$queryonly,
    'targetonly=s'  => \$targetonly,
    'stdout=s'  	=> \$stdout,
    'both=s'  		=> \$both
) )
{ print "Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "outputdir not defined (Use --help for usage)\n" if not defined $outputdir;
#die "prefix not defined (Use --help for usage)\n" if not defined $prefix;
die "targetfile not defined (Use --help for usage)\n" if not defined $targetfile;
die "queryfile not defined (Use --help for usage)\n" if not defined $queryfile;

#### PRINT TO STDOUT IF DEFINED stdout
if ( defined $stdout )
{
	print "Printing STDOUT to file:\n\n$stdout\n\n" if defined $stdout;
	my ($stdoutdir) = $stdout =~ /^(.+?)\/[^\/]+$/;
	print "stdoutdir: $stdoutdir\n";
	File::Path::mkpath($stdoutdir) if not -d $stdoutdir;
	print "Can't create stdoutdir: $stdoutdir\n" and exit if not -d $stdoutdir;
	open(STDOUT, ">$stdout") or die "Can't redirect STDOUT to file: $stdout\n" if defined $stdout;
	open(STDERR, ">>$stdout") or die "Can't redirect STDERR to file: $stdout\n" if defined $stdout;
}

#### CREATE OUTPUT DIR IF NOT EXISTS
print "outputdir is a file: $outputdir\n" and exit if -f $outputdir;
File::Path::mkpath($outputdir) if not -d $outputdir;
print "Can't create outputdir: $outputdir\n" and exit if not -d $outputdir;

#### DEBUG

my $venn = Venn::Snp->new(
	{
		queryfile 	=>	$queryfile,
		targetfile 	=>	$targetfile,
		querylabel 	=>	$querylabel,
		targetlabel =>	$targetlabel,
		outputdir 	=>	$outputdir,
		suffix		=>	$suffix,
		prefix		=>	$prefix
	}
);
$venn->compare();


#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "snpVenn.pl    Run time: $runtime\n";
print "snpVenn.pl    Completed $0\n";
print "snpVenn.pl    ";
print Timer::datetime(), "\n";
print "snpVenn.pl    ****************************************\n\n\n";
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


