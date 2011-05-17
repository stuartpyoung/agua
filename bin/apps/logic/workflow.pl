#!/usr/bin/perl -w

#### DEBUG

=head2

	APPLICATION     workflow.pl

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


perl C:\DATA\base\cgi-bin\agua\workflow.pl < C:\DATA\base\html\agua\fileroot\admin\tempfile-admin.json
{"workflow":"Workflow3-indels","sessionId":"1228791394.7868.158","mode":"executeWorkflow","project":"Project1","start":"0","username":"admin"}


=cut

use strict;

#### PRINT HEADER

#### FLUSH STDOUT SO THE MIME TYPE GETS OUT BEFORE ANY ERRORS
$| = 1;


#### USE LIB FOR Net::LDAP, JSON, ETC.
use lib "E:/agua/lib/external";


#### USE LIBS
use FindBin qw($Bin);
use lib "$Bin/lib";
use lib "$Bin/lib/external";

#### INTERNAL MODULES
use Agua::Admin;
use Agua::Workflow;
use Agua::DBaseFactory;
use Util;
use Conf::Agua;

#### EXTERNAL MODULES
use DBI;
use DBD::SQLite;
use Data::Dumper;

#### PRINT HEADER FOR DEBUG

#### GET PUTDATA
my $input	= <STDIN>;
#my $input = $ARGV[0];

if ( not defined $input or not $input or $input =~ /^\s*$/ )
{
	print "{ error: 'workflow.pl    input not defined' }";

	#### FLUSH
	$| = 1;
	exit;
}

#### GET JSONs
use JSON;
my $jsonParser = JSON->new();
my $json = $jsonParser->allow_nonref->jsonToObj($input);

#my $outfile = "workflow.out";


#### GET CONF
my $configfile = "default-linux.conf";
if ( $^O =~ /^MSWin32$/ )
{
    $configfile = "default-win32.conf";
}

my $conf = Conf::Agua->new(inputfile=>"$Bin/conf/$configfile", 0);

#### GET BIN DIRECTORY
my $bindir = $conf->getKeyValue("agua", 'BINDIR');

#### GET CLUSTER IF DEFINED
my $cluster = $conf->getKeyValue("agua", 'CLUSTERTYPE');

#### GET SETUID SCRIPT
my $setuid = "$Bin/msubMaster.pl";

##### CREATE DB OBJECT USING DBASE FACTORY
my $dbtype = $conf->getKeyValue("database", 'DBTYPE');
my $dbobject = 	DBaseFactory->new( $dbtype,
	{
		'DBFILE'	=>	$conf->getKeyValue("database", 'DBFILE'),
		'DATABASE'	=>	$conf->getKeyValue("database", 'DATABASE'),
		'USER'      =>  $conf->getKeyValue("database", 'USER'),
		'PASSWORD'  =>  $conf->getKeyValue("database", 'PASSWORD')
	}
) or print qq{ error: 'Cannot create database object $conf->getKeyValue("database", 'DATABASE'): $!' } and exit;


#### CHECK JSON
if ( not defined $json )
{
	print "{ 'error' : 'workflow.pl   JSON not defined' }";
	exit;	
}

#### GET MODE
my $mode = $json->{mode};

#### CHECK MODE
print "{ error: 'workflow.pl    mode not defined' }" and exit if not defined $mode;

$| = 1;

my $workflow = Workflow->new(
	{
		'DBOBJECT'	=>	$dbobject,							
		#'CGI'	    =>	$cgi,
        'CONF'      =>  $conf,
        'JSON'      =>  $json,
		'CLUSTER'	=>	$cluster,
		'SETUID' 	=>  $setuid,
		'QSUB'		=>	'/usr/local/bin/qsub',
		'QSTAT'		=>	'/usr/local/bin/qstat'
	}
);

#### RUN QUERY
no strict;
my $result = $workflow->$mode();
use strict;

exit;



#sub _getppid() {
#    
#    my $ppid;
#
#    if ($^O =~ /^MSWin/) {
#        my $pid = $$;
#        my $machine = "\\\\.";
#        
#        require Win32::OLE;
#        require Win32::OLE::Variant;
#    
#        # WMI Win32_Process class
#        my $class =
#"winmgmts:{impersonationLevel=impersonate}$machine\\Root\\cimv2";
#        if (my $wmi = Win32::OLE-> GetObject($class)) {
#            if(my $proc=$wmi-> Get(qq{Win32_Process.Handle="$pid"})) {
#                $ppid = $proc-> {ParentProcessId} if
#($proc-> {ParentProcessId}!=0);
#            }
#        }
#    }
#    else {
#        $ppid = getppid();
#    }
#    
#    return $ppid;
#}

#### PRINT OWN PROCESS ID
#my $process_id = $$;
#
#

#### GET THE PARENT PROCESS ID AND KILL IT
#my $parent_pid = _getppid();
#`kill -9 $parent_pid`;
