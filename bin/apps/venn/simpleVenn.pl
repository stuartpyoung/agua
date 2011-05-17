#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     simpleVenn

    PURPOSE

        CREATE VENN FILES FOR TWO FILES CONTAINING ONE NUMERIC COLUMN EACH

	VERSION		0.02

	HISTORY

    INPUTS

        1. QUERY INPUT FILE

        2. TARGET INPUT FILE

    OUTPUTS

        1. PRINT 'A n B', 'A c B' and 'B c A' FILES

    USAGE

    ./simpleVenn.pl <--queryfile String> <--targetfile String>
			[--suffix String] [--targetonly String] [--both String]
			[--help]

    --queryfile  	:   Location of query input file
    --targetfile  	:   Location of target input file
    --outputdir 	:   Print output files to this directory
    --querylabel  	:   Label of query data
    --targetlabel  	:   Label of target data
    --suffix		:   Output files end in this suffix
    --help       	:   print help info

    EXAMPLES

perl /nethome/bioinfo/apps/agua/0.5/bin/apps/aligners/samVenn.pl \
--queryfile /scratch/syoung/base/pipeline/bixby/run1/tophat/analysis1/chrY \
--targetfile /scratch/syoung/base/pipeline/bixby/run1/tophat/analysis2/chrY \
--outputdir /scratch/syoung/base/pipeline/bixby/run1/tophat/venn \
--suffix analysis1-only \
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
use Venn::Simple;
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
my $suffix;
my $targetonly;
my $both;
my $help;
if ( not GetOptions (
    'outputdir=s'  	=> \$outputdir,
    'targetfile=s'  => \$targetfile,
    'queryfile=s'   => \$queryfile,
    'targetlabel=s' => \$targetlabel,
    'querylabel=s'  => \$querylabel,
    'suffix=s'  	=> \$suffix
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

my $venn = Venn::Simple->new(
	{
		queryfile 	=> $queryfile,
		targetfile 	=> $targetfile,
		outputdir 	=> $outputdir,
		querylabel	=> $querylabel,
		targetlabel => $targetlabel,
		suffix		=>	$suffix
	}
);
$venn->simple();


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


