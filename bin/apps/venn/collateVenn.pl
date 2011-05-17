#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     collateVenn

	PURPOSE

		COLLATE THE RESULTS OF A SERIES OF VENN COMPARISON

		IN TERMS OF LINE COUNTS PER OUTPUT FILE ('*-AND-*'

		AND '*-NOT-*' FILES) AND PRINT TO THE TABLE TO A

		FILE, E.G.:

		(COLUMN HEADINGS ADDED AFTERWARDS)

		Sample ID	Sample 	Different to	Shared with 	Sample 33  		
		(Billion 	SNPs 	Sample 33 (A)	Sample 33 (B)	ONLY (C)	A + B	B + C
		Reads)									
		0.1			24284	20720			3563	1		6904		24284	20467
		0.3			18250	14918			3331			17136		18250	20467
		0.4			21206	16801			4404			16063		21206	20467
		0.5			21206	16801			4404			16063		21206	20467
		0.6			22896	17261			5634			14833		22896	20467
		...
		3.2			20619	843				19775			692			20619	20467
		3.3			20468	0				20467			0			20468	20467

	VERSION		0.01

	HISTORY

		0.01 GENERIC COLLATOR BASED ON LINE COUNTS

    INPUTS

        1. OUTPUT FROM snpVenn.pl, samVenn.pl AND OTHER COMPARATORS:

			FILES CONTAINING ONE-RECORD PER LINE OF QUERY-ONLY,

			JOIN AND TARGET-ONLY RECORDS

    OUTPUTS

        1. A TABULAR FILE COMPARING AcB, AnB AND BcA

    USAGE

    ./collateVenn.pl <--queryfile String> <--targetfile String>
			[--queryfile String] [--targetfile String] 
			[--help]

    --queryfile  	:   Location of query input file
    --targetfile  	:   Location of target input file
    --targetonly  	:   Print to this file entries found in only in target file
    --queryonly  	:   Print to this file entries found in only in query file
    --both  		:   Print to this file entries found in both input files
    --xref  		:   Use this kgXref TSV table data file for knownGene annotations
    --help       	:   print help info

    EXAMPLES

/nethome/bioinfo/apps/agua/0.5/bin/apps/venn/collateVenn.pl \
--outputfile /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/autochr22/eland/cumulative/chr22/venn \
--queryfile hit-33-snp \
--targetfile hit-%REPLICATE%-snp \
--replicates 1-33 \

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

#### GET CONFIG INFO
my $conf = Conf::Agua->new(inputfile=>"$Bin/../../../conf/default.conf");
my $qstat = $conf->getKeyValue("cluster", 'QSTAT');
my $qsub = $conf->getKeyValue("cluster", 'QSUB'); #### /usr/local/bin/qsub
my $samtools = $conf->getKeyValue("applications", 'SAMTOOLS');

#### GET OPTIONS
my $binlevel;
my $bamfile;
my $queryfile;
my $targetfile;
my $queryindex;
my $targetindex;
my $querylabel;
my $targetlabel;
my $replicates;
my $outputfile;
my $suffix;
my $querydir;
my $targetdir;
my $help;
if ( not GetOptions (
    'binlevel=s'   	=> \$binlevel,
    'bamfile=s'   	=> \$bamfile,
    'outputfile=s'  => \$outputfile,
    'queryfile=s'   => \$queryfile,
    'targetfile=s'  => \$targetfile,
    'queryindex=s'   => \$queryindex,
    'targetindex=s'  => \$targetindex,
    'suffix=s'   	=> \$suffix,
    'querydir=s'   	=> \$querydir,
    'targetdir=s'   => \$targetdir,
    'replicates=s'  => \$replicates,
    'targetlabel=s'	=> \$targetlabel,
    'querylabel=s'	=> \$querylabel
) )
{ print "Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "suffix not defined (Use --help for usage)\n" if not defined $suffix;
die "outputfile not defined (Use --help for usage)\n" if not defined $outputfile;
die "replicates not defined (Use --help for usage)\n" if not defined $replicates;
die "targetlabel not defined (Use --help for usage)\n" if not defined $targetlabel;
die "querylabel not defined (Use --help for usage)\n" if not defined $querylabel;

#### DEBUG
print "collateVenn.pl    outputfile: $outputfile\n";
print "collateVenn.pl    replicates: $replicates\n";
print "collateVenn.pl    querylabel: $querylabel\n";
print "collateVenn.pl    targetlabel: $targetlabel\n";
print "collateVenn.pl    queryindex: $queryindex\n";
print "collateVenn.pl    targetindex: $targetindex\n";
print "collateVenn.pl    binlevel: $binlevel\n";
print "collateVenn.pl    bamfile: $bamfile\n";
print "collateVenn.pl    suffix: $suffix\n";

my $venn = Venn::Snp->new(	{
		binlevel 	=>	$binlevel,
		bamfile 	=>	$bamfile,
		querylabel 	=>	$querylabel,
		targetlabel =>	$targetlabel,
		queryindex 	=>	$queryindex,
		targetindex =>	$targetindex,
		outputfile 	=>	$outputfile,
		replicates 	=>	$replicates,
		binlevel	=>	$binlevel,
		bamfile		=>	$bamfile,
		suffix		=>	$suffix,
		querydir 	=>	$querydir,
		targetdir 	=>	$targetdir,

		#### SAMTOOLS
		samtools 	=> $samtools,
	}
);
$venn->collate();


#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "collateVenn.pl    Run time: $runtime\n";
print "collateVenn.pl    Completed $0 @arguments\n";
print "collateVenn.pl    ";
print Timer::datetime(), "\n";
print "collateVenn.pl    ****************************************\n\n\n";
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


