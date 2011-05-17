#!/usr/bin/perl -w

#### DEBUG

=head2

	APPLICATION 	init.pl

	PURPOSE

		EXECUTE INITIALISATION TASKS FOR AGUA:

			LOAD USER DATA INTO AGUA DATABASE

			PRINT X.509 KEY FILES FOR HTTPS AND AWS

			MOUNT AGUA DATA VOLUMES

			MOUNT USER VOLUMES	

=cut

use strict;

#### DISABLE BUFFERING OF STDOUT
$| = 1;

#### DEBUG HEADER
print "Content-type: text/html\n\n";

use CGI::Carp qw(fatalsToBrowser);

#### GET CONF
use FindBin qw($Bin);
use lib "$Bin/lib/";
use Conf::Agua;
my $conf = Conf::Agua->new(
	inputfile	=>	"$Bin/conf/default.conf",
	backup		=>	1,
	separator	=>	"\t",
	spacer		=>	"\\s\+"
);

#### PRINT TO HTML FILE
my $htmldir = $conf->{HTMLDIR};
my $filename = "initlog.html";
my $htmlfile = "$htmldir/$filename";
my $url = "./$filename";
print qq{<input type="button" onClick="window.open('$url', '_blank', 'toolbar=1,location=0,directories=0,status=0,menubar=1,scrollbars=1,resizable=1,navigation=0')" value="Click to view Progress Log"><br>
RUNNING AS: };
print `whoami`;

open(STDOUT, ">$htmlfile") or die "Can't open htmlfile: $htmlfile\n" if defined $htmlfile;

print qq{<html>
<head>
	<title>Init Progress Log</title>
</head>
<style type="text/css">
	body {
		font-size: 12px;
	}
</style>

<script>

// TIME FORMAT: minutes:seconds
var waitSeconds = 5;
function beginrefresh(){

	if (waitSeconds==1)
	{
		window.location.reload()
	}
	else
	{ 
		waitSeconds-=1
		cursec=waitSeconds;
		curtime=cursec+" seconds left until page refresh!"
	}

	window.status=curtime
	setTimeout("beginrefresh()",1000)
}

window.onload=beginrefresh
</script>

<body>

};
print "<PRE>\n";

print "whoami: ";
print `whoami`;

#### USE LIBS
use lib "lib";

#### INTERNAL MODULES
use Agua::Init;
use Util;

#### EXTERNAL MODULES
use Data::Dumper;

#### GET INPUT
my $input = $ARGV[0];
print "No JSON input<br>Exiting<br>\n" and exit if not defined $input or not $input;
$input =~ s/\s+$//;

#### CONVERT JSON TEXT TO OBJECT
use JSON;
my $jsonParser = JSON->new();
my $json = $jsonParser->allow_nonref->jsonToObj($input);

#### ADD CONF TO JSON OBJECT
$json->{data}->{conf} = $conf;

#### CREATE ADMIN OBJECT
my $aws = Init->new($json->{data});

#### EXECUTE INITIALISATION
$aws->init();

