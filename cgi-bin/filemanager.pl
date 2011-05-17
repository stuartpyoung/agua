#!/usr/bin/perl -w

#!C:\\perl\\bin\\perl -w

#### DEBUG

=head2

	APPLICATION     filemanager.pl

	PURPOSE

		RESPOND TO CLIENT QUERIES:

			1. SAVE WORKFLOW JSON SENT BY CLIENT

    		2. RUN WORKFLOWS ON REQUEST

	USAGE

		./workflow.cgi <mode> <database> [additional_arguments]

    EXAMPLES

			THE Workflow OBJECT PERFORMS THE FOLLOWING TASKS:

				1. SAVE WORKFLOWS

                2. RUN WORKFLOWS

                3. PROVIDE WORKFLOW STATUS

        SEE Workflow.pm FOR EXAMPLES


perl C:\DATA\base\cgi-bin\agua\filemanager.pl < C:\DATA\base\html\agua\fileroot\admin\tempfile-admin.json
{"workflow":"Workflow3-indels","sessionId":"1228791394.7868.158","mode":"executeWorkflow","project":"Project1","start":"0","username":"admin"}


=cut

use strict;

#### PRINT HEADER

#### FLUSH STDOUT SO THE MIME TYPE GETS OUT BEFORE ANY ERRORS
$| = 1;


#### USE LIBS
use FindBin qw($Bin);
use lib "$Bin/lib";
use lib "$Bin/lib/external";

#### INTERNAL MODULES
use FileManager;

#### GET PUTDATA
my $input	= <STDIN>;

#my $mode =	$ARGV[0];
#my $filepath = $ARGV[1];

#my $input = $ARGV[0];

if ( not defined $input or not $input or $input =~ /^\s*$/ )
{
	print "{ error: 'filemanager.pl    input not defined' }";

	#### FLUSH
	$| = 1;
	exit;
}

#### GET JSONs
use JSON;
my $jsonParser = JSON->new();
my $json = $jsonParser->allow_nonref->jsonToObj($input);
print "Content-type: text/xml\n\nfilemanager.pl    json not defined\n" and exit if not defined $json;

#### GET MODE
my $mode = $json->{mode};
print "Content-type: text/xml\n\nfilemanager.pl    modenot defined\n" and exit if not defined $mode;

my $fileManager = FileManager->new(
	{
        'JSON'      =>  $json,
	}
);

#### RUN QUERY
no strict;
my $result = $fileManager->$mode();
use strict;
exit;
