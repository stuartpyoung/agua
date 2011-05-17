#!/usr/bin/perl -w

#### DEBUG

=head2

	APPLICATION     workflow.pl

	PURPOSE

		RESPOND TO CLIENT QUERIES:

			1. SAVE WORKFLOW JSON SENT BY CLIENT

    		2. RUN WORKFLOWS ON REQUEST

	USAGE

		./workflow.pl <mode> <database> [additional_arguments]

    EXAMPLES

			THE Workflow OBJECT PERFORMS THE FOLLOWING TASKS:

				1. SAVE WORKFLOWS

                2. RUN WORKFLOWS

                3. PROVIDE WORKFLOW STATUS

        SEE Workflow.pm FOR EXAMPLES


perl C:\DATA\base\cgi-bin\agua\workflow.pl < C:\DATA\base\html\agua\fileroot\admin\tempfile-admin.json
{"workflow":"Workflow3-indels","sessionId":"1228791394.7868.158","mode":"executeWorkflow","project":"Project1","start":"0","username":"admin"}


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
use Agua::Workflow;
use Agua::DBaseFactory;
use Util;

#### EXTERNAL MODULES
use DBI;
use DBD::SQLite;
use Data::Dumper;

#### GET PUTDATA
#my $input	= <STDIN>;
my $input	= $ARGV[0];

print "{ error: 'workflow.pl    input not defined' }"
	and exit if not defined $input
	or not $input or $input =~ /^\s*$/;

#### GET JSON
use JSON;
my $jsonParser = JSON->new();
my $json = $jsonParser->allow_nonref->jsonToObj($input);
print "{ 'error' : 'workflow.pl   JSON not defined' }"
	and exit if not defined $json;

#### GET MODE
my $mode = $json->{mode};
print "{ error: 'workflow.pl    mode not defined' }"
	and exit if not defined $mode;
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


my $workflow = Agua::Workflow->new(
	{
        'conf'      =>  $conf,
        'json'      =>  $json
	}
);

#### RUN QUERY
no strict;
my $result = $workflow->$mode();
use strict;

exit;
