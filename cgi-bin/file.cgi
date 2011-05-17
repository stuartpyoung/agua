#!/usr/bin/perl -w

#### DEBUG
$DEBUG = 1;

=head2

    APPLICATION     file

    PURPOSE

        PRINT DATA FROM A SPECIFIED FILE:

            1. FROM A TSV FILE:

                -   RETURN TOP OR TAIL N LINES 

                -   FOR EACH LINE, RETURN AT MOST A SPECIFIED NUMBER OF CHARACTERS

            2. FROM A NON-TSV TEXT FILE (E.G., SEQUENCE FILE)

                -   RETURN N BYTES STARTING FROM A BYTE OFFSET

                -   RETURN LAST N BYTES IF BYTE OFFSET EXCEEDS FILE LENGTH

                -   RETURN FILE CONTENTS IF N BYTES EXCEEDS FILE LENGTH

    USAGE

perl -w file.pl username=syoung&file=/usr/bin/

=cut


#### PRINT HEADER
print "Content-type: text/html\n\n";

#### USE LIB
use FindBin qw($Bin);
use lib "$Bin/lib";

#### USE LIB FOR BUNDLED EXTERNAL MODULES, Net::LDAP, JSON, ETC.
use lib "$Bin/lib/external";

#### INTERNAL MODULES
use Admin;
use Util;
use DBaseFactory;

#### EXTERNAL MODULES
use DBI;
#use DBD::SQLite;
use DBD::mysql;
use Data::Dumper;

#### GET PUTDATA
my $input	= <STDIN>;

#### CHECK POSTDATA
#### FLUSH
$| = 1;
if ( not defined $input or not $input or $input =~ /^\s*$/ )
{
	print "{ error: 'admin.pl: input not defined' }";

	#### FLUSH
	$| = 1;
	exit;
}

#### CHECK JSON
use JSON;
my $jsonObject = JSON->new();
my $json = $jsonObject->jsonToObj($input);
print "Content-type: text/xml\n\nJSON not defined\n" and exit if not defined $json;

#### GET MODE
my $mode = $json->{mode};
print "Content-type: text/xml\n\nmodenot defined\n" and exit if not defined $mode;

#### GET CONF
my $configfile = "default.conf";
if ( $^O =~ /^MSWin32$/ )
{
    $configfile = "default-win32.conf";
}
my $conf = Util::conf("$Bin/conf/$configfile", 0);

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



my $admin = Admin->new(
	{
		'DBOBJECT'	=>	$dbobject,
		#'CGI'	    =>	$cgi,
        'CONF'      =>  $conf,
        'JSON'      =>  $json
	}
);

#### RUN QUERY
no strict;
my $result = $admin->$mode();
use strict;

exit;
