#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     cleanup

    PURPOSE

        CLEAN UP AFTER AN ALIGNMENT RUN, E.G.:

			1.	REMOVE INPUT FILE CHUNKS (SPLITFILES)

			2.  REMOVE ALIGNMENT SUBDIRS USED TO PERFORM

				INPUT FILE CHUNK ALIGNMENTS

			3. 	BUNDLE UP AND COMPRESS *.sh AND *stdout.txt FILES

			4. 	DELETE *.sh AND *stdout.txt FILES

    INPUT

        1. LOCATION OF OUTPUTDIR 

    USAGE

    ./Cluster.pl <--source String> <--mode String> [--help]

    --source     :   Full path to base alignment source
    --target	:   Comma-separated list of target (e.g., 'chr22,chrY')
    --mode       	:   Supported modes:
	                        subdirs - remove split file alignment subdirectories
	                        splitfiles - remove split file alignment subdirectories
	                        archive - tar and zip all *.sh and *.stdout.txt files
	                        delete - remove all *.sh and *.stdout.txt files
	--help          :   Print help info

	EXAMPLES

/nethome/bioinfo/apps/agua/0.5/bin/apps/aligners/cleanup.pl \
--mode subdirs \
--source /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/bowtie

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
my $source;
my $target;
my $mode;
my $help;
print "Cluster.pl    Use option --help for usage instructions.\n" and exit if not GetOptions (

	#### GENERAL
    'source=s' 	=> \$source,
    'target=s' => \$target,
    'mode=s' 	=> \$mode,
    'help' 		=> \$help
);

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "source not defined (Use --help for usage)\n" if not defined $source;
die "target not defined (Use --help for usage)\n" if not defined $target;
die "mode not defined (Use --help for usage)\n" if not defined $mode;

#### INSTANTIATE OBJECT AND RUN CLEANUP
my $cluster = Cluster->new();
$cluster->copySplitfiles($source, $target, $mode);

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "Cluster.pl    Run time: $runtime\n";
print "Cluster.pl    Completed $0\n";
print "Cluster.pl    ";
print Timer::current_datetime(), "\n";
print "Cluster.pl    ****************************************\n\n\n";
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

