#!/usr/bin/perl -w
use strict;

$DEBUG = 1;

#### TIME
my $time = time();

=head2

	APPLICATION     etcgroup

	PURPOSE

		PARSE USERS AND THEIR GROUPS FROM THE /etc/groups INFILE AND

		OUTPUT THEM TO FILE IN etcgroup TABLE FORMAT

	INFILEPUT:

		/etc/group INFILE 

	OUTFILEPUT

	    1. BSV OR TSV etcusers TABLE FILE FORMAT

    USAGE

		./etcgroup.pl <--inputfile String> <--outputfile String> [--home] [-h] 

    --inputfile		:   /full/path/to/inputfile (i.e., /etc/group)
    --outputfile	:   /full/path/to/outputfile
    --format		:   Two formats current supported: (tsv|bsv) (DEFAULT=tsv)
						bsv: bar-separated values, for loading into SQLite
						tsv: tab-separated values, for loading into MySQL
    --home			:   Comma-separated list of directories containing user home dirs (DEFAULT=/home)

	EXAMPLE


perl /data/agua/0.4/bin/scripts/etcgroup.pl \
--inputfile /nethome/syoung/base/pipeline/etc/group \
--outputfile /nethome/syoung/base/pipeline/etc/users.tsv \
--format tsv \
--home "/mihg/users,/home"

perl /data/agua/0.4/bin/scripts/etcgroup.pl \
--inputfile /nethome/syoung/base/pipeline/etc/group \
--outputfile /data/agua/0.4/bin/sql/tsvfiles/users.tsv \
--format tsv \
--home "/mihg/users,/home"



perl etcgroup.pl --inputfile E:/base/pipeline/etc/group --outputfile E:/base/pipeline/etc/etcgroup.bsv

THEN RUN loader.pl

perl loader.pl --dbfile e:/0.4/bin/db/data.dbl --sqlfile e:/0.4/bin/sql/etcgroup.sql --bsvfile e:/base/pipeline/etc/etcgroup.bsv  --nodups



=cut

#### FLUSH BUFFER
$| = 1;

#### EXTERNAL MODULES
use FindBin qw($Bin);
use Data::Dumper;
use Getopt::Long;

#### USE LIB
use lib "$Bin/../../lib";

#### INFILETERNAL MODULES
use Timer;


#### GET OPTIONS
my $inputfile;	
my $outputfile;
my $home = "/home";
my $format = "tsv";
my $help;
GetOptions (
	'inputfile=s' => \$inputfile,
	'outputfile=s' => \$outputfile,
	'home=s' => \$home,
	'format=s' => \$format,
	'help' => \$help) or die "No options specified. Try '--help'\n";

if ( defined $help )	{	usage();	}

#### FLUSH BUFFER
$| =1;

#### CHECK INFILEPUTS
die "Input file not defined (option --inputfile)\n" if not defined $inputfile; 
die "Output file not defined (option --outputfile)\n" if not defined $outputfile;
die "Could not find input file (option --inputfile)\n" if not -f $inputfile;
die "format must be tsv or bsv (option --format)\n" if not $format =~ /^(tsv|bsv)$/;

open(INFILE, $inputfile) or die "Can't open input file: $inputfile\n";
open(OUTFILE, ">$outputfile") or die "Can't open output file: $outputfile\n";

my $exceptions = exceptions();
#exit;

#### GET A HASH OF UNIQUE username,gid,groupname ENTRIES
my $already_seen;
while ( <INFILE> )
{
	next if $_ =~ /^\s*$/ or $_ =~ /^\#/;
#next if $_ =~ /mihg/;

	$_ =~ /^([^:]+):[^:]+:([^:]+):([^:]+)/;
	my $groupname = $1;
	next if exists $exceptions->{$groupname};
	my $gid = $2;
	my $users = $3;
	chomp($users);

	my @usernames = split ",", $users;
	if ( defined $groupname and $groupname )
	{
		foreach my $username ( @usernames )
		{
			next if exists $exceptions->{$username};

			my $line = "$username|$gid|$groupname\n";
			$already_seen->{$line} = 1 if not exists $already_seen->{$line};
		}
	}
}

my @keys = keys( %$already_seen );
@keys = sort @keys;


#### CREATE homedirs (username:homedir) HASH
my $homedirs;
my @homes = split ",", $home;
foreach my $key ( @keys )
{
	my ($username) = $key =~ /^([^\|]+)/;
	for my $basedir ( @homes )
	{


		my $homedir = "$basedir/$username";
		if ( $^O =~ /^MSWin32$/ )   {   $homedir =~ s/\//\\/g;  }
		$homedirs->{$username} = $homedir if -d $homedir;
	}
}
#exit;



####FORMAT:
####CREATE TABLE IF NOT EXISTS etcgroup
####(
####	username			VARCHAR(20),
####	gid                 INT(9),
####	groupname			VARCHAR(20),
####    email               VARCHAR(50),
####    lastname            VARCHAR(20),
####    firstname           VARCHAR(20),
####    description         TEXT,
####    homedir		        TEXT,
####	
####	PRIMARY KEY (username, gid, groupname)
####);

foreach my $key ( @keys )
{
	my ($username) = $key =~ /^([^\|]+)/;

	#### IGNORE ANY ENTRIES FOR USERS THAT DO NOT HAVE
	#### AN EXISTING HOME DIRECTORY IN THE FILE SYSTEM
	my $homedir = $homedirs->{$username};
	next if not defined $homedir;

	$key =~ s/\n$//;
	my $email = "$username\@med.miami.edu";
	my ($firstname, $lastname) = $username =~ /^(.)(.+)$/;
	$lastname =~ s/\d+$//;
	my $description = '';

	if ( $format =~ /^bsv$/ )
	{
		print OUTFILE "$key|$email|$firstname|$lastname|$description|$homedir\n";
	}
	else
	{
		$key =~ s/\|/\t/g;
		$key =~ s/\t$//;
		print OUTFILE "$key\t$email\t$firstname\t$lastname\t$description\t$homedir\n";
	}
}
close(INFILE);
close(OUTFILE);

#### PRINFILET RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "\nRun time: $runtime\n";
print "Completed $0\n";
print Timer::datetime(), "\n";
print "****************************************\n\n\n";
exit;


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### ####             SUBROUTFILEINFILEES                 #### #### #### ####  
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 

sub exceptions
{
	my $exceptions;
	%$exceptions = map { $_ => 1 } qw{
mihgadm
mihganlstad
mihganlstaut
mihganlstms
mihganlstspg
mihganlst
mihgappd
mihgclin
mihglabtech
mihglab
mihgspg
mihgstds
mihg
stu01
stu02
stu03
stu04
stu05
stu06
stu07
stu08
stu09
stu10
stu11
stu12
stu13
stu14
stu15
stu16
stu17
stu18
stu19
stu20
stu21
stu22
stu23
stu24
stu25
stu26
stu27
stu28
stu29
stu30
stu31
stu32
stu33
stu34
stu35
stu36
stu37
stu38
stu39
stu40
testbdn
tester10
tester1
tester2
tester3
tester4
tester5
tester6
tester7
tester8
tester9

accoustics
ant
ccsofficeadmin
ccstemp
root
adm
daemon
lp
news
uucp
nobody
nogroup
wheel
daemon
kmem
sys
tty
operator
mail
bin
procview
procmod
owner
everyone
group
staff
smmsp
_lp
_postfix
_postdrop
certusers
_keytabusers
utmp
authedusers
interactusers
netusers
consoleusers
_mcxalr
_pcastagent
_pcastserver
_serialnumberd
_devdocs
_sandbox
localaccounts
netaccounts
_mdnsresponder
_uucp
_ard
dialer
network
_www
_cvs
_svn
_mysql
_sshd
_qtss
_mailman
_appserverusr
admin
_appserveradm
_clamav
_amavisd
_jabber
_xgridcontroller
_xgridagent
_appowner
_windowserver
_spotlight
accessibility
_tokend
_securityagent
_calendar
_teamsserver
_update_sharing
_installer
_atsserver
_lpadmin
_unknown
};

	return $exceptions;
}


sub usage
{
    print `perldoc $0`;

	exit;
}
