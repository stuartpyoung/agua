#!/usr/bin/perl -w

#!C:\\perl\\bin\\perl -w


#### DEBUG

=head2

	APPLICATION 	download.cgi

	PURPOSE

		1. VERIFY USE ACCESS TO FILES OR PROJECT DATABASE ENTRIES

		2. SEND DOWNLOAD DATA STREAM TO CLIENT

	USAGE

		./download.cgi <CGI_arguments>

	EXAMPLES

perl -U download.cgi "mode=downloadFile&username=syoung&sessionId=1266310689.20187.37&filepath=Project1/Workflow1/454HCDiffs-headers-SNP.txt"


=cut

use strict;

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

#### INSTANTIATE CGI
use CGI;
my $cgi = CGI->new();

#### GET CGI PARAMS
my $mode = $cgi->param('mode');
my $filepath = $cgi->param('filepath');
my $username = $cgi->param('username');
my $sessionId = $cgi->param('sessionId');

my $json;
$json->{username} = $username;
$json->{sessionId} = $sessionId;
$json->{mode} = $mode;
$json->{filepath} = $filepath;

#### GET CONF
my $configfile = "default.conf";
my $conf = Util::conf("$Bin/conf/$configfile", 0);

my $workflow = Agua::Workflow->new(
    {
		'conf'		=>	$conf,
        'json' 		=>	$json
    }
);


#### RUN QUERY
no strict;
my $result = $workflow->$mode();
use strict;


exit;



