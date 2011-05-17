#!/usr/bin/perl -w
use strict;

#### DEBUG


=head2

    NAME		loader

    PURPOSE

		LOAD A TABLE INTO A SQLITE DATABASE

    INPUT

		1. SQLITE DATABASE FILE

        2. SQL 'CREATE TABLE' FILE

        3. BSV (BAR-SEPARATED VALUE) FILE WITH DATA TO LOAD

        4. 'nodups' OPTION IF WANT NO DUPLICATES LOADED

    OUTPUT

		1. NEWLY CREATED TABLE IN DATABASE

    USAGE

		./loader.pl <--dbfile String> <--sqlfile String> <--bsvfile String> [--nodups] [-h] 

    --dbfile            :   /full/path/to/dbfile to be used or created if not present
    --sqlfile           :   /full/path/to/sqlfile containing 'CREATE DATABASE...' command
    --bsvfile           :   /full/path/to/bsvfile (bar-separated values) containing data to load
    --nodups            :   Load no duplicate lines
     --help             :   print this help message

	< option > denotes REQUIRED argument
	[ option ] denotes OPTIONAL argument

    EXAMPLE

perl loader.pl --dbfile C:/DATA/base/html/Bioptic0.2.5/bin/data.dbl --sqlfile C:/DATA/base/html/Bioptic0.2.5/bin/sql/groups.sql --bsvfile C:/DATA/base/html/Bioptic0.2.5/bin/sql/groups.bsv

./loader.pl --dbfile C:\DATA\base\html\Bioptic0.2.5\bin\data.dbl --sqlfile C:\DATA\base\html\Bioptic0.2.5\bin\sql\groups.sql --bsvfile C:\DATA\base\html\Bioptic0.2.5\bin\sql\groups.bsv


./loader.pl --dbfile /nethome/syoung/base/pipeline/nimblegen-run1/SID9637_exon_Map/ccds-readinfo/mapping/pairalign.dbl --sqlfile /nethome/syoung/base/bin/exome/sql/pairalign.sql --bsvfile /nethome/syoung/base/pipeline/nimblegen-run1/SID9637_exon_Map/ccds-readinfo/mapping/pairalign.bsv-short --nodups


lines /nethome/syoung/base/pipeline/nimblegen-run1/SID9637_exon_Map/ccds-readinfo/mapping/pairalign.bsv
399464




=cut

#### TIME
my $time = time();

#### USE LIBS
use FindBin qw($Bin);
use lib "$Bin/../../lib";

#### HACK FOR WINDOWS
use lib "E:/agua/lib/external";

#### INTERNAL MODULES
use DBaseFactory;
use Data::Dumper;
use Timer;
use Util;
use Conf::Agua;

#### EXTERNAL MODULES
use Getopt::Long;

#### GET OPTIONS
my $dbfile;	
my $sqlfile;	
my $bsvfile;
my $nodups;
my $help;
GetOptions ('dbfile=s' => \$dbfile,
            'sqlfile=s' => \$sqlfile,
            'bsvfile=s' => \$bsvfile,
            'nodups' => \$nodups,
            'help' => \$help) or die "No options specified. Try '--help'\n";

if ( defined $help )	{	usage();	}

#### FLUSH BUFFER
$| =1;

#### CHECK INPUTS
die "Database file not defined (option --dbfile)\n" if not defined $dbfile;
die "SQL file not defined (option --sqlfile)\n" if not defined $sqlfile; 
die "BSV file not defined (option --bsvfile)\n" if not defined $bsvfile;
die "Could not find SQL file (option --sqlfile)\n" if not -f $sqlfile;
die "Could not find BSV file (option --bsvfile)\n" if not -f $bsvfile;

print "DBFILE: $dbfile\n";

#### CREATE DB OBJECT USING DBASE FACTORY
my $dbobject = DBaseFactory->new( "SQLite", { 'DBFILE' => $dbfile } )
	or die "Can't open DB file '$dbfile': $!\n";
if ( not defined $dbobject )
{
    die "Database object not defined. Tried to access sqlite DB file: $dbfile\n";
}

#### SET TABLE
my ($table) = $sqlfile =~ /([^\/^\\]+)\.sql/;

#### DROP TABLE
print "Dropping table: $table\n";
$dbobject->drop_table($table);

#### CREATE PROJECTS TABLE
print "Creating table: $table\n";
$dbobject->create_table($sqlfile);

### LOAD DATA INTO TABLE
print "Loading data into table: $table\n";
if ( $nodups )
{
    print "Using nodups...\n";
    $dbobject->load_nodups($table, $bsvfile);    
}
else
{
    $dbobject->load($table, $bsvfile);
}

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "\nRun time: $runtime\n";
print "Completed $0\n";
print Util::datetime(), "\n";
print "****************************************\n\n\n";
exit;

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#									SUBROUTINES
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


sub usage
{
    print `perldoc $0`;

	exit;
}
