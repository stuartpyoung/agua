#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     maqIndex

    PURPOSE

        WRAPPER SCRIPT FOR RUNNING maq-build TO INDEX

		ALL '.fa' inputdir FILES IN A SPECIFIED FOLDER

    USAGE

    ./maqIndex.plv <--inputdir String> [--help]

    --inputdir           :   Location of directory containing *.fa files
    --help               :   print help info

    EXAMPLES

chmod 755 /nethome/bioinfo/apps/agua/0.4/bin/apps/maqIndex.pl


perl /nethome/bioinfo/apps/agua/0.4/bin/apps/maqIndex.pl \
--inputdir /nethome/bioinfo/data/sequence/chromosomes/rat/rn4/maq


=cut

use strict;

#### EXTERNAL MODULES
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use Getopt::Long;
use FindBin qw($Bin);
use File::Path;

#### USE LIBRARY
use lib "$Bin/../../../lib";

#### INTERNAL MODULES
use MAQ;
use Timer;
use Util;
use Conf::Agua;

##### STORE ARGUMENTS TO PRINT TO FILE LATER
my @arguments = @ARGV;
print "maqIndex.pl    arguments: @arguments\n";

#### FLUSH BUFFER
$| =1;

#### GET CONF 
my $conf = Conf::Agua->new(inputfile=>"$Bin/../../../conf/default.conf");
my $maq = $conf->getKeyValue("agua", 'MAQ');
print "maqIndex.pl    maq: $maq\n";

#### GET OPTIONS
my $inputdir;
my $outputdir;
my $help;
if ( not GetOptions (
    'inputdir=s'   	=> \$inputdir,
    'outputdir=s'   => \$outputdir,
    'help'          => \$help
) )
{ print "Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "inputdir not defined (Use --help for usage)\n" if not defined $inputdir;
die "outputdir not defined (Use --help for usage)\n" if not defined $outputdir;
print "outputdir is a file: $outputdir\n" if -f $outputdir;
File::Path::mkpath($outputdir) if not -d $outputdir;
print "Can't create outputdir: $outputdir\n" if not -d $outputdir;

#### INSTANIATE MAQ OBJECT
my $maqObject = MAQ->new(
	{
		maq => $maq 		
	}
);
$maqObject->faToBfa($inputdir, $outputdir);

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "maqIndex.pl    Run time: $runtime\n";
print "maqIndex.pl    Completed $0\n";
print "maqIndex.pl    ";
print Timer::datetime(), "\n";
print "maqIndex.pl    ****************************************\n\n\n";
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


