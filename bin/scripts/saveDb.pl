#!/usr/bin/perl -w
use strict;

#### DEBUG
$DEBUG = 1;


=head2

    NAME		saveDb

    PURPOSE

		DUMP ALL THE TABLES IN A DATABASE TO .TSV FILES

    INPUT

		1. DATABASE NAME

		2. LOCATION TO PRINT .TSV FILES

    OUTPUT

		1. ONE .TSV FILE FOR EACH TABLE IN DATABASE

    USAGE

		./saveDb.pl <--db String> <--outputdir String> [-h] 

    --db 			:   Name of database
    --outputdir 	:   Location of output directory
	--help      	:   print this help message

	< option > denotes REQUIRED argument
	[ option ] denotes OPTIONAL argument

    EXAMPLE

perl saveDb.pl --db agua --outputdir E:/0.4/bin/sql/04

=cut

#### TIME
my $time = time();

#### USE LIBS
use FindBin qw($Bin);
use lib "$Bin/../../lib";
use lib "$Bin/../../lib/external";

#### HACK FOR WINDOWS
use lib "E:/agua/lib/external";

#### INTERNAL MODULES
use DBaseFactory;
use Data::Dumper;
use Timer;
use Util;
use Conf::Agua;

#### EXTERNAL MODULES
use File::Path;
use File::Copy;
use Getopt::Long;

#### GET OPTIONS
my $db;	
my $outputdir;	
my $help;
GetOptions ('db=s' => \$db,
            'outputdir=s' => \$outputdir,
            'help' => \$help) or die "No options specified. Try '--help'\n";

if ( defined $help )	{	usage();	}

#### FLUSH BUFFER
$| =1;

#### CHECK INPUTS
die "Database not defined (option --db)\n" if not defined $db;
die "Output directory not defined (option --outputdir)\n" if not defined $outputdir; 
die "File with same name as output directory already exists: $outputdir\n" if -f $outputdir;

#### CREATE OUTPUT DIRECTORY
#### 1. INCREMENT EXISTING DIR
if ( -d $outputdir )
{
	my $number = 0;
	my $archivedir = "$outputdir.$number";
	while ( -d $archivedir )	{	 $archivedir = "$outputdir." . $number++;	}
	File::Copy::move($outputdir, $archivedir);
}

File::Path::mkpath($outputdir) if not -d $outputdir;
die "Can't create output directory: $outputdir\n" if not -d $outputdir;

#### GET CONF
my $configfile = "default.conf";
if ( $^O =~ /^MSWin32$/ )
{
    $configfile = "default-win32.conf";
}
my $conf = Conf::Agua->new(inputfile=>"$Bin/../../conf/$configfile", 0);

#### GET DATABASE INFO
my $dbtype = $conf->getKeyValue("database", 'DBTYPE');
my $database = $conf->getKeyValue("database", 'DATABASE');
my $user = $conf->getKeyValue("database", 'USER');
my $password = $conf->getKeyValue("database", 'PASSWORD');


#### CREATE DB OBJECT USING DBASE FACTORY
my $dbobject;
if ( $dbtype eq "SQLite" )
{
	my $db = $conf->getKeyValue("database", 'DBFILE');
	$dbobject = DBaseFactory->new( "SQLite", { 'DBFILE' => "$db" } );
}
elsif ( $dbtype eq "MySQL" )
{
	#### CREATE DB OBJECT USING DBASE FACTORY
	$dbobject = DBaseFactory->new( 'MySQL',
		{
			'DATABASE'	=>	$database,
			'USER'      =>  $user,
			'PASSWORD'  =>  $password
		}
	) or die "Can't create database object to create database: $database. $!\n";
}

my $query = qq{SHOW TABLES};
my $tables = $dbobject->queryarray($query);
print "saveDb.pl    tables: @$tables\n";

for my $table ( @$tables )
{
	my $outputfile = "$outputdir/$table.tsv";

	if ( $^O =~ /^MSWin32$/ )   {   $outputfile =~ s/\//\\/g;  }

	$query = qq{SELECT * from $table};
	my $lines = $dbobject->querytwoDarray($query);
	next if not defined $lines;

	open(OUTFILE, ">$outputfile") or die "Can't open output file: $outputfile\n";
	for my $line ( @$lines )
	{
		next if not defined $line or scalar(@$line) == 0;
		#exit;
		no warnings;
		print OUTFILE join "\t", @$line;
		use warnings;

		print OUTFILE "\n";
	}
	close(OUTFILE);
	print "saveDb.pl    Created output file: $outputfile\n";
}


#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "\n";
print "saveDb.pl    Run time: $runtime\n";
print "saveDb.pl    Completed $0\n";
print Util::datetime(), "\n";
print "saveDb.pl    ****************************************\n\n\n";
exit;

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#									SUBROUTINES
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


sub usage
{
    print `perldoc $0`;

	exit;
}
