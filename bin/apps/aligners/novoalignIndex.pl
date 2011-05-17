#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     novoalignIndex

    PURPOSE

        WRAPPER SCRIPT FOR RUNNING novoalign-build TO INDEX ALL 

		'.fa' inputdir FILES IN A SPECIFIED FOLDER

    USAGE

    ./novoalignIndex.pl <--inputdir String> [--help]

    --inputdir           :   Location of directory containing *.fa files
    --help               :   print help info

    EXAMPLES

chmod 755 /nethome/bioinfo/apps/agua/0.4/bin/apps/novoalignIndex.pl


/nethome/bioinfo/apps/agua/0.5/bin/apps/novoalignIndex.pl \
--inputdir /nethome/bioinfo/data/sequence/chromosomes/rat/rn4/novoalign


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
use NOVOALIGN;
use Timer;
use Util;
use Conf::Agua;

##### STORE ARGUMENTS TO PRINT TO FILE LATER
my @arguments = @ARGV;
print "novoalignIndex.pl    arguments: @arguments\n";

#### FLUSH BUFFER
$| =1;

#### GET CONF 
my $conf = Conf::Agua->new(inputfile=>"$Bin/../../../conf/default.conf");
my $novoalign = $conf->getKeyValue("applications", 'NOVOALIGN');
print "novoalignIndex.pl    novoalign: $novoalign\n";

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
print "novoalignIndex.pl    inputdir: $inputdir\n";

#### DO CONVERSION
my $novoalignObject = NOVOALIGN->new({	novoalign	=>	$novoalign });
$novoalignObject->convertReferences($inputdir);




#chdir($inputdir) or die "Can't change to inputdir directory: $inputdir\n";
#my @files = <*fa>;
##### TRUNCATE inputdir FILES TO CREATE CORRECT STUB IDENTIFIER
#foreach my $file ( @files )
#{
#	my ($stub) = $file =~ /^(.+)\.fa$/;
#	my $command = "time $novoalign/novoalign-build  $file $stub";
#	`$command`;
#}

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "novoalignIndex.pl    Run time: $runtime\n";
print "novoalignIndex.pl    Completed $0\n";
print "novoalignIndex.pl    ";
print Timer::datetime(), "\n";
print "novoalignIndex.pl    ****************************************\n\n\n";
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


