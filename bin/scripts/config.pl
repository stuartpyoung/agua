#!/usr/bin/perl -w
use strict;

$DEBUG = 1;

=head2

APPLICATION     config

PURPOSE

    1. CONFIGURE THE Agua DATABASE

    2. CONFIGURE DATA AND APPLICATION PATHS AND SETTINGS

        E.G., PATHS TO BASIC EXECUTABLES IN CONF FILE:

        [applications]
        STARCLUSTERDEFAULT      /data/apps/starcluster/110202bal/bin/starcluster
        BOWTIE                  /data/apps/bowtie/0.12.2
        CASAVA                  /data/apps/casava/1.6.0/bin
        CROSSMATCH              /data/apps/crossmatch/0.990329/cross_match
        CUFFLINKS               /data/apps/cufflinks/0.8.2
        ...

INPUT

    1. MODE OF ACTION, E.G., 'mysql', 'addKey', removeKey

OUTPUT

    MYSQL DATABASE CONFIGURATION AND EDITED CONFIG FILE            

USAGE

sudo ./config.pl <--mode String> [--key String] [--value String] [--help]

--mode      :   mysql | addKey | removeKey (initialise database or edit config file)
--database  :	Name of database
--key       :	Name of key
--value     :	Value for key
--help      :   Print help info

EXAMPLES

sudo config.pl --mode mysql --database 

=cut

#### FLUSH BUFFER
$| = 1;

#### USE LIB
use FindBin qw($Bin);
use lib "$Bin/../../lib";

#### EXTERNAL MODULES
use Getopt::Long;
use Data::Dumper;

#### INTERNAL MODULES
use Agua::Configure;
use Agua::DBaseFactory;
use Conf;

#### GET OPTIONS
my $mode;
my $database;
my $configfile;
my $logfile;
my $help;
GetOptions (
    'mode=s'        => \$mode,
    'database=s'    => \$database,
    'configfile=s'  => \$configfile,
    'logfile=s'     => \$logfile,
    'help'          => \$help
) or die "No options specified. Try '--help'\n";

#### ASK FOR USER INPUT FOR INSTALLATION DIRECTORIES
print "\n\nmode not defined (--mode). Using default mode: mysql\n\n"
    if not defined $mode;
$mode = "mysql" if not defined $mode;

my $configObject = Agua::Configure->new(
    {
        mode        =>  $mode,
        database    =>  $database,
        configfile  =>  $configfile,
        logfile     =>  $logfile
    }
);


#### RUN QUERY
no strict;
eval { $configObject->$mode() };
if ( $@ ){
	print "Error: $mode): $@\n";
}
print "\nCompleted $0\n";

sub usage {
    print `perldoc $0`;
}



