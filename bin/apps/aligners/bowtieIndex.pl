#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     bowtieIndex

    PURPOSE

        WRAPPER SCRIPT FOR RUNNING bowtie-build TO INDEX ALL 

		'.fa' inputdir FILES IN A SPECIFIED FOLDER

    USAGE

    ./bowtieIndex.plv <--inputdir String> [--help]

    --inputdir           :   Location of directory containing *.fa files
    --help               :   print help info

    EXAMPLES

chmod 755 /nethome/bioinfo/apps/agua/0.4/bin/apps/bowtieIndex.pl


perl /nethome/bioinfo/apps/agua/0.4/bin/apps/bowtieIndex.pl \
--inputdir /nethome/bioinfo/data/sequence/chromosomes/rat/rn4/bowtie


=cut

use strict;

#### EXTERNAL MODULES
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use Getopt::Long;
use FindBin qw($Bin);

#### USE LIBRARY
use lib "$Bin/../../../lib";

#### INTERNAL MODULES
use Timer;
use Util;
use Conf::Agua;

##### STORE ARGUMENTS TO PRINT TO FILE LATER
my @arguments = @ARGV;
print "bowtieIndex.pl    arguments: @arguments\n";

#### FLUSH BUFFER
$| =1;

#### GET CONF 
my $conf = Conf::Agua->new(inputfile=>"$Bin/../../../conf/default.conf");
my $bowtie = $conf->getKeyValue("applications", 'BOWTIE');
print "bowtieIndex.pl    bowtie: $bowtie\n";

#### GET OPTIONS
# GENERAL
my $inputdir;

my $help;
if ( not GetOptions (
    'inputdir=s'   => \$inputdir,
    'help'          => \$help
) )
{ print "Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "inputdir not defined (Use --help for usage)\n" if not defined $inputdir;
print "bowtieIndex.pl    inputdir: $inputdir\n";

chdir($inputdir) or die "Can't change to inputdir directory: $inputdir\n";
my @files = <*fa>;
#### TRUNCATE inputdir FILES TO CREATE CORRECT STUB IDENTIFIER
foreach my $file ( @files )
{
	my ($stub) = $file =~ /^(.+)\.fa$/;
	my $command = "time $bowtie/bowtie-build  $file $stub";
	print "command: $command\n";
	`$command`;
}

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "bowtieIndex.pl    Run time: $runtime\n";
print "bowtieIndex.pl    Completed $0\n";
print "bowtieIndex.pl    ";
print Timer::datetime(), "\n";
print "bowtieIndex.pl    ****************************************\n\n\n";
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


