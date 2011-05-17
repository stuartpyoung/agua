#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     samHits

    PURPOSE

        PARSE OUT THE HITS ONLY FROM A SAM FILE

    USAGE

    ./samHits.pl <--inputfile String> [--outputfile String] [--help]

    --inputfile		:   Location of input SAM format file
    --outputfile	:   Location of output file (default = STDOUT)
    --help      	:   print help info

    EXAMPLES

chmod 755 /nethome/bioinfo/apps/agua/0.4/bin/apps/samHits.pl

/nethome/bioinfo/apps/agua/0.5/bin/apps/samHits.pl \
--inputfile /nethome/bioinfo/apps/agua/0.5/bin/apps/t/Cluster/bowtie/chrY/out.sam


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
use Timer;
use Util;
use Conf::Agua;

##### STORE ARGUMENTS TO PRINT TO FILE LATER
my @arguments = @ARGV;

#### FLUSH BUFFER
$| =1;


#### GET OPTIONS
my $inputfile;
my $outputfile;
my $missfile;
my $help;
if ( not GetOptions (
    'inputfile=s'	=> \$inputfile,
    'outputfile=s'	=> \$outputfile,
    'missfile=s'	=> \$missfile,
    'help'			=> \$help
) )
{ print "Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "inputfile not defined (Use --help for usage)\n" if not defined $inputfile;
open(FILE, $inputfile) or die "Can't open input file: $inputfile\n";

my $writeover = 0;
if ( defined $outputfile and $inputfile eq $outputfile )
{
	$writeover = 1;
	$outputfile .= ".temp"
}


#### OPEN OUTFILE IF DEFINED
open(OUTFILE, ">$outputfile") or die "Can't open output file: $outputfile\n" if defined $outputfile;

#### OPEN MISSED FILE IF DEFINED
open(MISSED, ">$missfile") or die "Can't open missfile: $missfile\n" if defined $missfile;


#### IF HEADER PRESENT, SKIP LINES BEGINNING WITH '@'
$/ = "\n";
my $line = <FILE>;
while ( $line =~ /^\@/ )
{
	$line = <FILE>;
}


#### SEPARATE OUT HITS
#### PRINT MISSES TO MISSFILE IF DEFINED
if ( defined $missfile )
{
	while ( <FILE> )
	{
		#### NON-MATCHED ENTRY HAS '*' INSTEAD OF REFERENCE NAME IN 3RD FIELD
		print MISSED $_ if $_ =~ /^\S+\s*\t\S+\s*\t\*/;
		next if $_ =~ /^\S+\s*\t\S+\s*\t\*/;

		print OUTFILE "$_" if defined $outputfile;
		print $_ if not defined $outputfile;
	}	
}

#### OTHERWISE, IGNORE MISSES
else
{
	while ( <FILE> )
	{
		#### NON-MATCHED ENTRY HAS '*' INSTEAD OF REFERENCE NAME IN 3RD FIELD
		next if $_ =~ /^\S+\s*\t\S+\s*\t\*/;

		print OUTFILE "$_" if defined $outputfile;
		print $_ if not defined $outputfile;
	}
}
close(FILE);
close(OUTFILE) if defined $outputfile;


#### WRITE OVER ORIGINAL FILE IF writeover IS DEFINED
if ( $writeover )
{
	print `mv -f $outputfile $inputfile`;
}

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


