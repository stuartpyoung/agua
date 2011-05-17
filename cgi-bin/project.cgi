#!/usr/bin/perl -U

#### DEBUG

=head2

	APPLICATION 	project.cgi

	PURPOSE

        SUBMISSION OF REQUESTS TO project.pl, WHICH USES Project.pm

		TO PERFORM THE FOLLOWING TASKS:

            1. RETURN FILE AND DIRECTORY LIST OF A GIVEN PATH

            2. MAINTAIN A DATABASE TABLE OF THE HIERARCHY OF FOLDERS:

            3. MAINTAIN THE USER/GROUP/WORLD PERMISSIONS ON THE PROJECT FOLDERS

            4. PROVIDE THE FOLLOWING FUNCTIONALITY:

    NOTES

        MUST DO THE FOLLOWING FOR THIS SCRIPT TO WORK:

            1. SET SUID BIT FOR UID AND GID OF SCRIPT ON LINUX COMMAND LINE

                AFTER EVERY SSH TRANSFER:

chmod u+s project.cgi; chmod g+s project.cgi; ll project.cgi

	USAGE

perl -U project.cgi <username> <sessionId> <mode> 


    EXAMPLES

perl -U project.cgi "mode=fileSystem&sessionId=1228791394.7868.158&username=admin&path=.%2FProject1%2Finputdir"

perl -U project.cgi "mode=fileSystem&sessionId=1228791394.7868.158&username=admin&path=.%2FProject1%2Fcode.js"

perl -u project.cgi "mode=fileSystem&sessionId=1228791394.7868.158&username=admin&path=.%2FProject1%2Finputdir"

perl -U project.cgi "mode=fileSystem&query=Project1/Workflow1-assembly/exome&sessionId=1228791394.7868.158&username=admin"

	OK - RUN project.cgi AS apache AND EXECUTES project.pl AS USER syoung (HARD CODED IN project.cgi FOR TESTING)



perl -U project.cgi "mode=projects&sessionId=1228791394.7868.158&username=admin"

=cut

use strict;

#### DEBUG HEADER
print "Content-type: text/html\n\n";

#### USE LIBS
use FindBin qw($Bin);
use lib "$Bin/lib";
use lib "$Bin/lib/external";

#### INTERNAL MODULES
use Agua::Admin;
use Agua::DBaseFactory;
use Agua::Project;

#use Util;

#### EXTERNAL MODULES
use DBI;
use DBD::SQLite;
use Data::Dumper;
use Util;

my $input = $ARGV[0];
if ( not defined $input )
{
    $input = $ENV{'QUERY_STRING'};
}

#### CREATE CGI OBJECT TO BE PASSED TO Project.pm
#### FOR USE IN ITS cgiParam METHOD
use Agua::JSON;
my $jsonConverter = Agua::JSON->new();
my $json = $jsonConverter->cgiToJson($input);
print "Content-type: text/xml\n\n{ error: 'project.cgi    json not defined' }\n" and exit if not defined $json;

#### GET mode
my $mode = $json->{mode};
print "{ error: 'project.cgi    mode not defined' }" and exit if not defined $mode;

#### GET CONF
use FindBin qw($Bin);
use lib "$Bin/lib";
use Conf::Agua;
my $conf = Conf::Agua->new(
	inputfile	=>	"$Bin/conf/default.conf",
	backup		=>	1,
	separator	=>	"\t",
	spacer		=>	"\\s\+"
);

#### CHECK mode IS DEFINED
print "{ error: 'Mode not defined' }" and exit if not defined $mode;

my $project = Agua::Project->new(
	{
		'json'	    =>	$json,
        'conf'      =>  $conf
	}
);

#### RUN QUERY
no strict;
my $result = $project->$mode();
use strict;

exit;
