#!/usr/bin/perl -w
use strict;

use Net::LDAP;
use CGI;

my $cgi = CGI->new();
my $username = $cgi->param('username');
my $password = $cgi->param('password');

print "Content-type: text/html\n\n";

my $ldap = Net::LDAP->new( 'ldap.ccs.miami.edu' );

#Server: ldap.ccs.miami.edu
#Binddn: uid=USERNAME,ou=Users,dc=ccs,dc=miami,dc=edu
#Bindpw: USERPASS

# bind to a directory with dn and password
my $message = $ldap->bind(
    "uid=$username,ou=Users,dc=ccs,dc=miami,dc=edu",
    "password" => "$password"
);

print "\$message->code(): ", $message->code(), "\n";
