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

    ./Cluster.pl <--directory String> <--mode String> [--help]

    --directory     :   Full path to base alignment directory
    --references	:   Comma-separated list of references (e.g., 'chr22,chrY')
    --mode       	:   Supported modes:
	                        subdirs - remove split file alignment subdirectories
	                        splitfiles - remove split file alignment subdirectories
	                        archive - tar and zip all *.sh and *.stdout.txt files
	                        delete - remove all *.sh and *.stdout.txt files
	--help          :   Print help info

	EXAMPLES

/nethome/bioinfo/apps/agua/0.5/bin/apps/aligners/cleanup.pl \
--mode subdirs \
--directory /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/bowtie

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
my $directory;
my $references;
my $mode;
my $help;
print "Cluster.pl    Use option --help for usage instructions.\n" and exit if not GetOptions (

	#### GENERAL
    'directory=s' 	=> \$directory,
    'references=s' 	=> \$references,
    'mode=s' 		=> \$mode,
    'help' 			=> \$help
);

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "directory not defined (Use --help for usage)\n" if not defined $directory;
die "mode not defined (Use --help for usage)\n" if not defined $mode;
my $modes;
@$modes = split ",", $mode;
foreach my $mode ( @$modes )
{
	die "mode not supported (subdirs|archive|delete) (Use --help for usage)\n" if $mode !~ /^(splitfiles|subdirs|archive|delete)$/;
	die "must specifiy references (e.g., 'chr22') for subdirs mode (Use --help for usage)\n" if not defined $references and $mode =~ /^subdirs$/;
}

#### IF DEFINED REFERENCES
my $refs;
@$refs = split ",", $references if defined $references;

#### INSTANTIATE OBJECT AND RUN CLEANUP
my $cleanup = Cluster->new();
foreach my $mode ( @$modes )
{
	$cleanup->deleteSplitfiles($directory, $refs) if $mode =~ /^splitfiles$/;
	$cleanup->deleteAlignmentSubdirs($directory, $refs) if $mode =~ /^subdirs$/;

	#### ADD ARCHIVE AND DELETE SUBS
	$cleanup->archiveMiscfiles($directory, $refs) if $mode =~ /^archive$/;
	$cleanup->deleteMiscfiles($directory, $refs) if $mode =~ /^delete$/;	
}

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

