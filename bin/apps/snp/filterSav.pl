#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     filterSav

    PURPOSE

        PARSE A *.sav FILE AND SEPARATE OUT dbSNP AND NON-dbSNP ENTRIES

    INPUT

        1. LOCATION OF *.sav FILE

    USAGE

    ./Cluster.pl <--inputfile String> <--mode String> [--help]

    --inputfile     :   Full path to base alignment inputfile
	--help          :   Print help info

	EXAMPLES

/nethome/syoung/0.5/bin/apps/snp/filterSav.pl \
--inputfile /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/bowtie/cumulative/chr22/hit-1.sav

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
use Cluster;
use Timer;
use Util;
use Conf::Agua;

#### GET OPTIONS
my $inputfile;
my $stdout;
my $help;
print "Cluster.pl    Use option --help for usage instructions.\n" and exit if not GetOptions (

	#### GENERAL
    'inputfile=s' 	=> \$inputfile,
    'stdout=s' 		=> \$stdout,
    'help' 			=> \$help
);

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### PRINT TO STDOUT IF DEFINED stdout
print "Printing STDOUT to file:\n\n$stdout\n\n" if defined $stdout;
open(STDOUT, ">$stdout") or die "Can't redirect STDOUT to file: $stdout\n" if defined $stdout;
open(STDERR, ">>$stdout") or die "Can't redirect STDERR to file: $stdout\n" if defined $stdout;

#### CHECK INPUTS
die "inputfile not defined (Use --help for usage)\n" if not defined $inputfile;
die "Can't find inputfile: $inputfile\n" if not -f $inputfile;
open(FILE, "$inputfile") or die "Can't open inputfile: $inputfile\n";

#### OUTPUT FILES
my $dbsnp = $inputfile;
$dbsnp =~ s/\.sav$//;
$dbsnp .= "-dbsnp.sav";

my $nondbsnp = $inputfile;
$nondbsnp =~ s/\.sav$//;
$nondbsnp .= "-non-dbsnp.sav";

open(DBSNP, ">$dbsnp") or die "Can't open dbSNP file: $dbsnp\n";
open(NON, ">$nondbsnp") or die "Can't open non-dbSNP file: $nondbsnp\n";

while ( <FILE> )
{
	next if $_ =~ /^\s*$/ or $_ =~ /^#/;
	my @elements = split "\t", $_;
	my $rsnumber = $elements[10];
	if ( defined $rsnumber and $rsnumber )
	{
		print DBSNP $_, "\n" ;
	}
	else
	{
		print NON $_, "\n";
	}
}
close(DBSNP);
close(NON);

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "Cluster.pl    Run time: $runtime\n";
print "Cluster.pl    Completed $0\n";
print "Cluster.pl    ";
print Timer::current_datetime(), "\n";
print "Cluster.pl    ****************************************\n\n\n";
print "\n------------------------------------------------------------\n";
print "---[completed filtersav: completed]---";
print "\n------------------------------------------------------------\n";
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

