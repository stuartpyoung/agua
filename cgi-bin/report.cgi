#!C:\\perl\\bin\\perl -w

#!/usr/bin/perl -w


#### DEBUG

=head2

	APPLICATION 	report.cgi

	PURPOSE


        USES Report.pm TO PERFORM THE FOLLOWING TASKS:

            1. RETRIEVE A STORED REPORT SET OF PARAMETERS, INCLUDING NOTES

            2. RETRIEVE THE FINAL DATA SET OF A STORED REPORT

            3. APPLY A SERIES OF FILTERS TO REPORT DATA (E.G., SNP: dbSNP, READS,

                VARIANT FREQUENCY, EXONIC, INTRONIC, INSIDE UTR, DISTANCE UPSTREAM

                OF GENE, DISTANCE DOWNSTREAM OF GENE)

            4. FILTER REPORT DATA BASED ON PROVIDED CRITERIA AND SHOW THE NUMBER

                OF RETURNED DATA ENTRIES (E.G., SNP THRESHOLD FILTER) AT EACH

                STEP/CRITERIA


	USAGE

		./report.cgi <mode> <sessionId> <

        mode : addReport | removeReport | renameReport | etc.


    VERSION

        0.2     FIXED Dynaloader.pm ERROR WITH SQLite:

            DynaLoader.pm: Name "DBD::SQLite::sqlite_version" used only once: possible typo at C:/Perl/lib/DynaLoader.pm line 225.

    EXAMPLES


cd C:\DATA\base\cgi-bin\Bioptic0.2.5 "mode=snpFilter&sessionId=1228319084.3060.776&username=admin"

localhost:8080/cgi-bin/Bioptic0.2.5/report.cgi?mode=snpFilter&sessionId=1228319084.3060.776&username=admin

localhost:8080/Bioptic2/html/dojo.1.2.2/dojox/data/demos/stores/filestore_dojotree.php?path=.%2FReport1%2Fcode.js


perl report.cgi "module=Report::SNP&mode=reports&sessionId=1228791394.7868.158&username=admin"

saveReport

perl report.cgi "module=Report::SNP&mode=saveReport&sessionId=1228791394.7868.158&username=admin"



perl report.cgi "module=Report::SNP&mode=saveReport&sessionId=1228791394.7868.158&username=admin"



=cut

use strict;

# C CODE TO GET RID OF THIS ERROR: DynaLoader.pm: Name "DBD::SQLite::sqlite_version" used only once: possible typo at C:/Perl/lib/DynaLoader.pm line 225.
#
# if ( sv = get_sv("DBD::SQLite::sqlite_version", TRUE | GV_ADDMULTI) );


#### PRINT HEADER 


my $putdata = <STDIN>;
if ( not defined $putdata )
{
	print "{    error: 'Putdata not defined'    }";
	exit;	
}

#### REMOVE BACKSLASHES FROM PUTDATA
$putdata =~ s/\\//g;
$putdata =~ s/false/"false"/g;
$putdata =~ s/true/"true"/g;


#### USE LIBS
use FindBin qw($Bin);
use lib "$Bin/lib";
use lib "$Bin/lib/external";

#### INTERNAL MODULES
use Admin;
use DBaseFactory;
use Report;
use Util;

#### EXTERNAL MODULES
use DBI;
use DBD::SQLite;
use Data::Dumper;
#use CGI;
#use CGI::Carp qw(fatalsToBrowser);

#### GET USER, PASSWORD AND SESSION ID FROM CGI OBJECT
#my $cgi = CGI->new();

##### GET mode AND module
#my $mode = $cgi->param('mode');
#my $module = $cgi->param('module');

#### GET CONF
my $configfile = "default.conf";
$configfile = "default-win32.conf" if $^O =~ /^MSWin32$/;
my $conf = Util::conf("conf/$configfile", 0);

#### GET BIN DIRECTORY
my $bindir = $conf->{BINDIR};

#### GET CLUSTER IF DEFINED
my $cluster = $conf->{CLUSTER};

#### GET SETUID SCRIPT
my $setuid = "$Bin/msubMaster.pl";

##### CREATE DB OBJECT USING DBASE FACTORY
my $dbtype = $conf->{DBTYPE};
my $dbobject = 	DBaseFactory->new( $dbtype,
	{
		'DBFILE'	=>	$conf->{DBFILE},
		'DATABASE'	=>	$conf->{DATABASE},
		'USER'      =>  $conf->{USER},
		'PASSWORD'  =>  $conf->{PASSWORD}
	}
) or print qq{ error: 'Cannot create database object $conf->{DATABASE}: $!' } and exit;


use JSON;
my $jsonParser = JSON->new();
my $json = $jsonParser->decode($putdata);
print "Content-type: text/xml\n\n{ error: 'report.cgi    json not defined' }\n" and exit if not defined $json;

my $mode = $json->{mode};
print "Content-type: text/xml\n\n{ error: 'report.cgi    mode not defined' }\n" and exit if not defined $mode;

my $module = $json->{module};
print "Content-type: text/xml\n\n{ error: 'report.cgi    module not defined' }\n" and exit if not defined $module;

#### GET mode
my $location = $module;
$location =~ s/::/\//;
$location .= ".pm";

require $location;

my $report = $module->new(
	{
		'DBOBJECT'	=>	$dbobject,							
		#'CGI'	    =>	$cgi,
        'CONF'      =>  $conf,
        'JSON'      =>  $json
	}
);

#### RUN QUERY
no strict;
my $result = $report->$mode();
use strict;

exit;
