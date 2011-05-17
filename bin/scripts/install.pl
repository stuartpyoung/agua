#!/usr/bin/perl -w
use strict;

$DEBUG = 1;

=head2

APPLICATION     install

PURPOSE:

    1. INSTALL THE DEPENDENCIES FOR Agua

    2. CREATE THE REQUIRED DIRECTORY STRUCTURE

INPUT

    1. INSTALLATION DIRECTORY (DEFAULT: /agua)

    2. www DIRECTORY (DEFAULT: /var/www)

OUTPUT

    1. REQUIRED DIRECTORY STRUCTURE AND PERMISSIONS

        FOR PROPER RUNNING OF Agua

    2. RUNNING APACHE INSTALLATION CONFIGURED FOR Agua

    3. RUNNING MYSQL INSTALLATION AWAITING MYSQL DATABASE

        CONFIGURATION WITH CONFIG

USAGE

    sudo ./install.pl <--installdir String> <--wwwdir String> [--help]

    --installdir        :   /Full/path/to/input/directory
    --wwwdir			:	fasta|fastq
    --help              :   Print help info

EXAMPLES

sudo install.pl --installdir /path/to/installdir 

=cut

#### FLUSH BUFFER
$| = 1;

#### USE LIB
use FindBin qw($Bin);
use lib "$Bin/../../lib";

#### EXTERNAL MODULES
use Getopt::Long;

#### INTERNAL MODULES
use Agua::Installer;

#### GET OPTIONS
my $installdir = "/agua";
my $wwwdir = "/var/www";
my $logfile;
my $help;
GetOptions (
    'installdir=s' => \$installdir,
    'wwwdir=s' => \$wwwdir,
    'logfile=s' => \$logfile,
    'help' => \$help
) or die "No options specified. Try '--help'\n";

#### PRINT HELP IF REQUESTED
if ( defined $help )	{	usage();	}

print "install.pl    installdir: $installdir\n";
print "install.pl    wwwdir: $wwwdir\n";

my $installer = Agua::Installer->new(
    {
        installdir  =>  $installdir,
        wwwdir      =>  $wwwdir,
        logfile     =>  $logfile
    }
);

$installer->install();

print "Completed $0\n";


sub usage {
    print `perldoc $0`;
}

