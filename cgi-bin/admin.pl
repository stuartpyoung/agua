#!/usr/bin/perl -w

#### DEBUG

=head2

	APPLICATION 	admin.pl

	PURPOSE

		RESPOND TO XMLHTTPRequest LOGIN AUTHENTICATION QUERIES

        TO DO:

            - Login USING USERNAME AND PASSWORD

            - Validation USING USERNAME AND SESSION ID

            - CREATE NEW DATABASE

=cut

use strict;

#### PRINT HEADER WITH TEXT MIME TYPE
print "Content-type: text/html\n\n";


#### FLUSH STDOUT SO THE MIME TYPE GETS OUT BEFORE ANY ERRORS
$| = 1;


#### USE LIBS
use FindBin qw($Bin);
use lib "$Bin/lib";

#### USE LIB FOR BUNDLED EXTERNAL MODULES, Net::LDAP, JSON, ETC.
use lib "$Bin/lib/external";

#### INTERNAL MODULES
use Agua::Admin;
use Agua::DBaseFactory;
use Util;

#### EXTERNAL MODULES
use DBI;
use DBD::mysql;
use Data::Dumper;

#### GET PUTDATA
#my $input	= <STDIN>;
my $input	= $ARGV[0];

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
print "Content-type: text/xml\n\n{ error: 'admin.pl    json not defined' }\n" and exit if not defined $json;

#### GET MODE
my $mode = $json->{mode};
print "Content-type: text/xml\n\n{ error: 'admin.pl    mode not defined' }\n" and exit if not defined $mode;

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

my $admin = Agua::Admin->new(
	{
        'conf'      =>  $conf,
        'json'      =>  $json
	}
);

#### RUN QUERY
no strict;
my $result = $admin->$mode();
use strict;

exit;
