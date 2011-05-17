#!/usr/bin/perl -w

#### DEBUG

=head2

	APPLICATION 	view.cgi

	PURPOSE

		THIS APPPLICATION USES THE View.pm MODULE TO PERFORM THE

		FOLLOWING TASKS:

			1. ADD/REMOVE Views

			2. ALTER Views

`   USAGE

        INPUTS ARE 'PUT' (EQUIVALENT TO STDIN) TO THIS APPLICATION 

        IN JSON FORMAT.

        THE FOLLOWING JSON HASH KEYS ARE REQUIRED:

        mode        Name of subroutine to be executed

        username    Short username of clintt

        password    Database password

        sessionId   If session has already been created


    EXAMPLES

echo '{"mode":"addView","username":"myUsername","sessionId":"0000000000.0000.000","sessionId":"w3rasdfaw34qaerfgasdfas", ... [additional_arguments] ... }'   | ./view.cgi

=cut

use strict;

#### PRINT HEADER
print "Content-type: text/html\n\n";

#### FLUSH STDOUT SO THE MIME TYPE GETS OUT BEFORE ANY ERRORS
$| = 1;


#### USE LIBS
use FindBin qw($Bin);
use lib "$Bin/lib";
use lib "$Bin/lib/external";

#### INTERNAL MODULES
use Agua::View;
use Agua::DBaseFactory;
use Util;

#### EXTERNAL MODULES
use DBI;
use DBD::SQLite;
use Data::Dumper;

#### GET PUTDATA
my $input	= <STDIN>;
print "{ error: 'view.cgi    input not defined' }"
	and exit if not defined $input
	or not $input or $input =~ /^\s*$/;

#### GET JSON
use JSON;
my $jsonParser = JSON->new();
my $json = $jsonParser->allow_nonref->jsonToObj($input);
print "{ 'error' : 'view.cgi   json not defined' }" and exit if not defined $json;

#### GET MODE
my $mode = $json->{mode};
print "{ error: 'view.cgi    mode not defined' }" and exit if not defined $mode;
$| = 1;

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

#### INTIANTIATE View OBJECT
my $view = Agua::View->new(
	{
        'conf'      =>  $conf,
        'json'      =>  $json
		#'cluster'	=>	$cluster
	}
);

#### RUN QUERY
no strict;
my $result = $view->$mode();
use strict;
exit;
