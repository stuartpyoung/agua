#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     copyFiles

    PURPOSE

        COPY A LIST OF FILES FROM ONE DIRECTORY TO ANOTHER

    INPUT

		1. COMMA-SEPARATED LIST OF FILES TO BE COPIED

		2. LOCATION OF SOURCE DIRECTORY

        3. LOCATION OF TARGET DIRECTORY

    USAGE

    ./copyFiles.pl <--source String> <--target String> <--filename String> [--help]

    --source     :  Full path to base directory containing subdirs
    --target	 :	Full path to destination for subdirs
    --filename   :  Comma-separated list of files or directories
	--dir        :	Copy directories, not files (default: copy files only)
	--help       :  Print help info

	EXAMPLES

/nethome/bioinfo/apps/agua/0.5/bin/apps/aligners/copyFiles.pl \
--filename hit.bam \
--source /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland/1/chr22 \
--target /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/reeland/1/chr22

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
my $filename;
my $mode;
my $dir;
my $help;
print "Cluster.pl    Use option --help for usage instructions.\n" and exit if not GetOptions (

	#### GENERAL
    'source=s' 		=> \$source,
    'target=s' 		=> \$target,
    'filename=s' 	=> \$filename,
    'mode=s' 		=> \$mode,
    'dir' 			=> \$dir,
    'help' 			=> \$help
);

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "source not defined (Use --help for usage)\n" if not defined $source;
die "target not defined (Use --help for usage)\n" if not defined $target;
die "filename not defined (Use --help for usage)\n" if not defined $filename;
die "mode not supported (name|regex) (Use --help for usage)\n" if defined $mode and $mode !~ /^(name|regex)$/;

#### MAKE OUTPUT DIR IF NOT EXISTS
print "copyFiles.pl    target is a file: $target\n" and exit if -f $target;
File::Path::mkpath($target) if not -d $target;
print "copyFiles.pl    Can't create output directory: $target\n" and exit if not -d $target;

#### INSTANTIATE OBJECT AND RUN CLEANUP
my $movedir = Cluster->new();
$movedir->copyFiles($source, $target, $filename, $mode, $dir);

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

