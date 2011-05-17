#!/usr/bin/perl -w
use strict;

#### DEBUG


=head2

    NAME		addindex

    PURPOSE

		CREATE AN INDEX IN AN SQLITE TABLE

    INPUT

		1. SQLITE DATABASE FILE

        2. TABLE, INDEX NAME AND FIELDS

    OUTPUT

		1. NEW INDEX IN DATABASE

    USAGE

		./addindex.pl <--dbfile String> <--table String> <--index String> <--fields String>  [-h] 

    --dbfile            :   /Full/path/to/dbfile to be used or created if not present
    --table             :   Name of table
    --index             :   Name of index
    --fields            :   Fields to be included in index
     --help             :   print this help message

	< option > denotes REQUIRED argument
	[ option ] denotes OPTIONAL argument

    EXAMPLE

cd /home/syoung/base/html/Bioptic0.2.5/bin

./addindex.pl \
--dbfile /nethome/syoung/base/pipeline/nimblegen-run1/SID9638_exon_Map/ccds-readinfo/mapping/pairalign.dbl \
--table pairalign \
--index gene \
--fields gene


=cut

#### TIME
my $time = time();

#### USE LIBS
use lib "../lib";

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
my $table;	
my $index;
my $fields;
my $help;
GetOptions ('dbfile=s' => \$dbfile,
            'table=s' => \$table,
            'index=s' => \$index,
            'fields=s' => \$fields,
            'help' => \$help) or die "No options specified. Try '--help'\n";

if ( defined $help )	{	usage();	}

#### FLUSH BUFFER
$| =1;

#### CHECK INPUTS
if ( not defined $dbfile)   {   print "Database file not defined (option --dbfile)\n";    usage();    }
if ( not defined $table)   {   print "Table not defined (option --table)\n";    usage();    }
if ( not defined $index)   {   print "Index not defined (option --index)\n";    usage();    }
if ( not defined $fields)   {   print "Fields not defined (option --fields)\n";    usage();    }

#### CHECK INPUT FILES
if ( not -f $dbfile )    {    print "Could not find input file (option --dbfile)\n"; usage();	}

#### CREATE DB OBJECT USING DBASE FACTORY
my $dbobject = DBaseFactory->new( "SQLite", { 'DBFILE' => $dbfile } ) or die "Can't open DB file '$dbfile': $!\n";

## GET ALL RECORDS IN TABLE
my $query = qq{CREATE INDEX $index ON $table ($fields)};
my $result = $dbobject->do($query);
print "Result: $result\n";

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
