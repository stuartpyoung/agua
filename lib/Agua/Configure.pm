use MooseX::Declare;

=head2

PACKAGE		Configure

PURPOSE

    1. CONFIGURE THE Agua DATABASE

    2. CONFIGURE DATA AND APPLICATION PATHS AND SETTINGS

        E.G., PATHS TO BASIC EXECUTABLES IN CONF FILE:

        [applications]
        STARCLUSTERDEFAULT      /data/apps/starcluster/110202bal/bin/starcluster
        BOWTIE                  /data/apps/bowtie/0.12.2
        CASAVA                  /data/apps/casava/1.6.0/bin
        CROSSMATCH              /data/apps/crossmatch/0.990329/cross_match
        CUFFLINKS               /data/apps/cufflinks/0.8.2
        ...

=cut


use strict;
use warnings;
use Carp;

#### USE LIB FOR INHERITANCE
use FindBin qw($Bin);
use lib "$Bin/../";
use lib "$Bin/../external";
use Data::Dumper;

class Agua::Configure {

#### USE LIB
use FindBin qw($Bin);
use lib "$Bin/../../lib";

#### EXTERNAL MODULES
use Data::Dumper;
use Term::ReadKey;

#### INTERNAL MODULES
use Agua::DBaseFactory;
use Conf::Agua;

# STRINGS
has 'configfile'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'logfile'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'database'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'user'			=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'password'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'rootuser'		=> ( isa => 'Str|Undef', is => 'rw', default => 'root' );
has 'rootpassword'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
# OBJECTS
has 'dbobject'		=> ( isa => 'Agua::DBase::MySQL', is => 'rw', required => 0 );
has 'conf' 	=> (
	is =>	'rw',
	isa => 'Conf::Agua',
	default	=>	sub { Conf::Agua->new(	backup	=>	1, separator => "\t"	);	}
);

####/////}

method BUILD ($hash) {

	#### SET DEFAULT CONFIG FILE IF NOT SPECIFIED
	my $configfile = $Bin;
	$configfile =~ s/\/[^\/]+$//;
	$configfile =~ s/\/[^\/]+$//;
	$configfile .= "/conf/default.conf";
	$self->configfile($configfile) if not $self->configfile();
	print "Agua::Configure::BUILD    configfile: ", $self->configfile, "\n";

	#### LOAD CONFIG FILE INTO PARSER
	$self->conf()->inputfile($self->configfile());

	#### SET Agua DATABASE USER'S NAME
	$self->user($self->conf()->getKeyValue("database", "USER"));
	print "Agua::Configure::BUILD    user: ", $self->user, "\n";

	#### SET Agua DATABASE USER'S PASSWORD
	$self->password($self->conf()->getKeyValue("database", "PASSWORD"));
	print "Agua::Configure::BUILD    password: ", $self->password, "\n";

	#### SET DATABASE NAME
	$self->database($self->conf()->getKeyValue("database", "DATABASE"))
		if not $self->database();
	print "Agua::Configure::BUILD    database: ", $self->database, "\n";

	$self->logfile("$Bin/log/config.log") if not $self->logfile();
	$self->openLogfile($self->logfile());
	print "Agua::Configure::BUILD    logfile: ", $self->logfile(), "\n";
}

=head2

	SUBROUTINE		cron

	PURPOSE

		INSERT INTO /etc/crontab COMMANDS TO BE RUN AUTOMATICALLY

		ACCORDING TO THE DESIRED TIMETABLE:

			1. checkBalancers.pl	-	VERIFY LOAD BALANCERS ARE

				RUNNING AND RESTART THEM IF THEY HAVE STOPPED. THIS

				TASK RUNS ONCE A MINUTE

	NOTES

		cat /etc/crontab
		# /etc/crontab: system-wide crontab
		# Unlike any other crontab you don't have to run the `crontab'
		# command to install the new version when you edit this file
		# and files in /etc/cron.d. These files also have username fields,
		# that none of the other crontabs do.

		SHELL=/bin/sh
		PATH=/usr/local/sbin:/usr/local/bin:/sbin:/bin:/usr/sbin:/usr/bin

		# m h dom mon dow user	command
		17 *	* * *	root    cd / && run-parts --report /etc/cron.hourly
		25 6	* * *	root	test -x /usr/sbin/anacron || ( cd / && run-parts --report /etc/cron.daily )
		47 6	* * 7	root	test -x /usr/sbin/anacron || ( cd / && run-parts --report /etc/cron.weekly )
		52 6	1 * *	root	test -x /usr/sbin/anacron || ( cd / && run-parts --report /etc/cron.monthly )

=cut

method cron {
	print "Running cron configuration\n";

	my $installdir = $self->conf()->getKeyValue("agua", 'INSTALLDIR');	
	my $inserts	= [
		'* *     * * *   root    $installdir/bin/scripts/checkBalancers.pl > /tmp/loadbalancers.out'
];
	#my $crontext = crontextFromFile();
	my $crontext = `crontab -1`;

	#### REMOVE INSERTS IF ALREADY PRESENT	
	foreach my $insert ( @$inserts )
	{
		my $temp = $insert;
		$temp =~ s/\*/\\*/g;
		$temp =~ s/\$/\$/g;
		$temp =~ s/\-/\\-/g;
		$temp =~ s/\//\\\//g;
		$crontext =~ s/$temp//msg;
	}

	#### ADD INSERTS TO BOTTOM OF CRON LIST	
	$crontext =~ s/\s+$//;
	foreach my $insert ( @$inserts ) {	$crontext .= "\n$insert";	}
	$crontext .= "\n";

	`echo '$crontext' | crontab -`;
}

method crontextFromFile ($crontab) {
	$crontab = "/etc/crontab" if not defined $crontab;
	$self->backupFile($crontab);
	my $temp = $/;
	$/ = undef;
	open(FILE, $crontab) or die "Can't open crontab: $crontab\n";
	my $crontext = <FILE>;
	close(FILE) or die "Can't close crontab: $crontab\n";
	$/ = $temp;

	return $crontext;
}

method crontextToFile ($crontab, $crontext) {
	open(OUT, ">$crontab") or die "Can't open crontab: $crontab\n";
	print OUT $crontext;
	close(OUT) or die "Can't close crontab: $crontab\n";
}

method mysql {
	#### PRINT INFO
    print qq{\n
Welcome to the agua dabase configuration utility (config.pl)
You can use this application to:
	1. Set the root MySQL password (if not already set)
	2. Initialize a new Agua database
	3. Reinitialize an existing 'agua' database (deletes existing data)
	4. Set the user name and password for the Agua database
\n};

	#### CHECK FOR QUIT
    #return if not $self->yes("Type 'Y' to continue or 'N' to exit");

    print qq{\n
Please provide values for the following items (Hit RETURN to accept '[default]' value)
};

    my $database        =   $self->database();
    my $user       		=   $self->user();	
    my $password   		=   $self->password();
    my $rootuser       	=   $self->rootuser();
    my $rootpassword   	=   $self->rootpassword();

##### debug
#$rootpassword = "open4root";

    #### DATABASE NAME
	$database = $self->inputValue("Database name", $database);	
    print "Database name : $database\n";

    #### ROOT PASSWORD
    #### MASK TYPING FOR PASSWORD INPUT
    ReadMode 2;
	$rootpassword = $self->inputValue("Root password (will not appear on screen)", $rootpassword);	
    print "\nRoot password: $rootpassword\n";
    #### UNMASK TYPING
    ReadMode 0;

    $self->rootuser($rootuser);
    $self->rootpassword($rootpassword);

	#### PROMPT TO SET ROOT MYSQL PASSWORD
	my $resetroot = "Do you want to reset the root password?\nType Y to reset or N to cancel";
	$self->setRootPassword($database) if $self->yes($resetroot);

    #### Agua USER NAME
	$user = $self->inputValue("Agua username", $user);	
    print "Agua username : $user\n";

    #### Agua USER PASSWORD
	$password = $self->inputValue("Agua password", $password);	
    print "Agua password: $password\n";

    #### SAVE CONFIGURATIONS
    $self->database($database);
    $self->user($user);
    $self->password($password);

	###### RESTART MYSQL

    #### CREATE DB OBJECT USING DBASE FACTORY
    my $dbobject = Agua::DBaseFactory->new( 'MySQL',
        {
			database	=>	"mysql",
            user      	=>  "root",
            password  	=>  $rootpassword
        }
    ) or die "Can't create database object to create database: $database. $!\n";

	#### PROMPT TO DELETE IF DATABASE EXISTS
	$self->checkExists($dbobject, $database);

	#### CREATE DATABASE 
	$self->createDatabase($database);

	#### SET Agua USER AND PASSWORD
	$self->createAguaUser($database);

	#### POPULATE DATABASE WITH DUMP FILE
	$self->loadDumpfile($database);

    #### PRINT CONFIRMATION THAT THE agua DATABASE HAS BEEN CREATED
	my $configfile = $self->configfile();
	my $logfile = $self->logfile();
	my $timestamp = $dbobject->timestamp();
	print qq{
*******************************************************
         Successfully created Agua database.

Please make a note of your MySQL access credentials:

\tDatabase:\t$database
\tUsername:\t$user
\tPassword:\t$password

You can test this on the command line:

\tmysql -u $user -p$password
\tUSE $database
\tSHOW TABLES

Your updated configuration file is here:

\t$configfile

This transaction has been recorded in the log:

\t$logfile

Timestamp: $timestamp
*******************************************************\n};

}

method openLogfile ($logfile) {
	print "Agua::Configure::openLogfile    Agua::Configure::openLogfile(logfile)\n";
	print "Agua::Configure::openLogfile    logfile: $logfile\n";
	$logfile .= ".1";

	if ( -f $logfile )
	{
		my ($stub, $index) = $logfile =~ /^(.+?)\.(\d+)$/;
		$index++;
		$logfile = $stub . "." . $index;
	}
	print "Agua::Configure::openLogfile    logfile: $logfile\n";
	$self->logfile($logfile);

	return $logfile;	
}


method loadDumpfile ($database) {
    my $rootuser       	=   $self->rootuser();
    my $rootpassword   	=   $self->rootpassword();

	my $dumpfile = $self->conf()->getKeyValue("database", "DUMPFILE");
	my $installdir = $self->conf()->getKeyValue('agua', "INSTALLDIR");
	$dumpfile = $installdir . "/" . $dumpfile;
    print "Agua::Configure::loadDumpfile    dumpfile: $dumpfile\n";

	print `mysql -u $rootuser -p$rootpassword $database < $dumpfile`;
}

method createDatabase ($database) {
    my $rootuser       	=   $self->rootuser();
    my $rootpassword   	=   $self->rootpassword();
    my $user       		=   $self->user();
    my $password   		=   $self->password();

	#### CREATE DATABASE AND Agua USER AND PASSWORD
	my $sqlfile = "$Bin/../sql/createDatabase.sql";
    print "Agua::Configure::createDatabase    sqlfile: $sqlfile\n";
	my $create = qq{CREATE DATABASE $database;};
	$self->printFile($sqlfile, $create);
	my $command = "mysql -u $rootuser -p$rootpassword < $sqlfile";
	print "$command\n";
	print `$command`;
	`rm -fr $sqlfile`;
}

method createAguaUser ($database) {
    my $rootuser       	=   $self->rootuser();
    my $rootpassword   	=   $self->rootpassword();
    my $user       		=   $self->user();
    my $password   		=   $self->password();

	#### CREATE DATABASE AND Agua USER AND PASSWORD
	my $sqlfile = "$Bin/../sql/createAguaUser.sql";
    print "Agua::Configure::createAguaUser    sqlfile: $sqlfile\n";
	my $create = qq{
USE mysql;
GRANT ALL PRIVILEGES ON $database.* TO $user\@localhost IDENTIFIED BY '$password';	
FLUSH PRIVILEGES;};
	$self->printFile($sqlfile, $create);
	my $command = "mysql -u $rootuser -p$rootpassword < $sqlfile";
	print "$command\n";
	print `$command`;
	`rm -fr $sqlfile`;
}

method setRootPassword ($database) {
    my $rootpassword   	=   $self->rootpassword();
    my $user       		=   $self->user();
    my $password   		=   $self->password();

	#### RESTART MYSQL WITH  --skip-grant-tables
	my $stop = "/etc/init.d/mysql stop";

	print "$stop\n";
	print `$stop`;
	my $start = "sudo mysqld --skip-grant-tables &";
	print "$start\n";
	print `$start`;
	sleep(1);

	#### SET root USER PASSWORD
	my $sqlfile = "$Bin/../sql/setRootPassword.sql";
    print "Root::Configure::setRootPassword    sqlfile: $sqlfile\n";
	my $create = qq{
UPDATE mysql.user SET Password=PASSWORD('$rootpassword') WHERE User='root'; 
FLUSH PRIVILEGES;
};
	$self->printFile($sqlfile, $create);
	my $command = "mysql -u root < $sqlfile";
	print "Setting root user password\n";
	print "$command\n";
	print `$command`;
	`rm -fr $sqlfile`;

	##### KILL mysql	
	##print `/etc/init.d/mysql restart`;	 #### HANGS
	#my $kill = qq{kill `ps aux | grep mysql | awk -F" " '{ if (\$1 == "mysql") print \$2}'`};
	#`$kill`;


	#sleep(1);
	#my @pidfiles = </var/lib/mysql/*pid>;
	#my $pid = `cat $pidfiles[0]` if @pidfiles;
	#my $oldpid = $pid;
	#while ( defined $pid and $pid and ($oldpid != $pid) )
	#{
	#	`kill -9 $pid`;
	#	#`rm -fr /var/lib/mysql/*pid &> /dev/null`;
	#	sleep(1);
	#
	#	$oldpid = $pid;
	#	undef $pid;
	#	@pidfiles = </var/lib/mysql/*pid>;
	#	$pid = `cat $pidfiles[0]` if @pidfiles;
	#}

}

method checkExists ($dbobject, $database) {

    my $rootuser       	=   $self->rootuser();
    my $rootpassword   	=   $self->rootpassword();

    #### CHECK IF DATABASE EXISTS
    my $is_database = $dbobject->is_database($database);
    print "Agua::Configure::mysql    Agua is_database: $is_database\n";
	return if not $is_database;

	my $warning = "Database '$database' exists already.\nPress 'Y' to overwrite it or 'N' to exit";
	exit if not $self->yes($warning);

	#### PRINT DUMP COMMAND FILE
	my $timestamp = $dbobject->timestamp();
	$timestamp =~ s/\s+/-/g;
	my $dumpfile = "$Bin/../sql/dump/agua.$timestamp.dump";
	my $cmdfile = "$Bin/../sql/dump/agua.$timestamp.cmd";
	$self->printFile($cmdfile, "#!/bin/sh\n\nmysqldump -u $rootuser -p$rootpassword $database > $dumpfile\n");

	#### DUMP CONTENTS OF Agua DATABASE TO FILE
	`chmod 755 $cmdfile`;

	print "Printing dump file backup of current database: $dumpfile\n";
	print `$cmdfile`;
	`rm -fr $cmdfile`;

	#### PRINT HEAD OF DUMP FILE
	print "dumpfile head:\n\n";
	print `head $dumpfile`;
	print "\n";

	#### DELETE DATABASE
	$dbobject->drop_database($database) or die "Can't drop database: $database. $!\n";

	my $logfile = $self->logfile();

	print qq{
*******************************************************
         Successfully deleted existing Agua database.

The tables and data in the existing database have been

saved here:

$dumpfile

This transaction has been recorded in the log:

\t$logfile

Timestamp: $timestamp
*******************************************************\n};

}

method printFile ($outfile, $contents) {
	print "Agua::Configure::printFile    outfile not defined\n"
		and exit if not defined $outfile;

	open(OUT, ">$outfile") or die "Can't open outfile: $outfile\n";
	print OUT $contents;
	close(OUT);
	print `cat $outfile`;	
}

method inputValue ($message, $default) {
	print "Agua::Configure::inputValue    message is not defined\n"
		and exit if not defined $message;
	$default = '' if not defined $default;
	print "$message [$default]: ";

	my $input = '';
    while ( $input =~ /^\s*$/ )
    {
        $input = <STDIN>;
        $input =~ s/\s+//g;

		$default = $input if $input;
		return $default if $default;

        print "$message [$default]: ";
    }
}


method getKey () {
	my $key		=	$self->key();
	print "Agua::Configure::_getKey    key not defined (--key)\n"
		and exit if not $key;
	my $value = $self->_getKey();
	print "$value\n";
}

method _getKey () {
	my $key		=	$self->key();
	print "Agua::Configure::_getKey    key: $key\n";

	return $self->conf()->getValue($key);
}

=head2

    SUBROUTINE      apps

    PURPOSE

        CONFIGURE APPLICATIONS INTERACTIVELY:

            1. SET THE PATHS TO EXISTING APPLICATION ENTRIES

            2. ADD NEW APPLICATIONS

			3. DELETE APPLICATIONS

method apps {

	#### PRINT APPLICATION PATHS        
    my $apps = $self->conf()->get("apps");
    print "Type Y to use the default application paths or type N to set application paths\n";
    return if yes("Type Y to use the default application paths or type N to set custom application paths");

    my $continue;
    while ( not defined $continue )
    {
        my $counter = 1;
        foreach my $app ( @$apps )
        {
            my @keys = keys %$app;
            my $key = $keys[0];
            print "$counter $key\t:\t$app->{$key}\n";
            $counter++;
        }

        $continue = input("Type application number or 'Q' for quit");

        if ( $continue !~ /^\d+$/ and $continue !~ /^q$/i )
        {
            print "\n**** Invalid input: $continue\n\n";
            $continue = undef;
            next;
        }

        #### CONTINUE MUST BE A NUMBER OR 'Q'
        if ( $continue =~ /^q$/i )
        {
            last;
        }
        else
        {
            my $app = $$apps[$continue - 1];
            if ( not defined $app )
            {
                print "No application entry for index: $continue\n";
            }
            else
            {
                my @keys = keys %$app;
                my $key = $keys[0];

                my $value = input("Enter path to $key [$app->{$key}]");
                $$apps[$continue - 1] = { $key => $value };    
            }

            $continue = undef;
        }
    }    
}

=cut

=head2

    SUBROUTINE      yes

    PURPOSE

        PROMPT THE USER TO ENTER 'Y' OR 'N'

=cut
method yes ($message) {
	return if not defined $message;
	print "$message: ";
    my $max_times = 10;

	$/ = "\n";
	my $input = <STDIN>;
	my $counter = 0;
	while ( $input !~ /^Y$/i and $input !~ /^N$/i )
	{
		if ( $counter > $max_times ) { print "Exceeded 10 tries. Exiting...\n"; }
		print "$message: ";
		$input = <STDIN>;
		$counter++;
	}	

	if ( $input =~ /^N$/i )	{	return 0;	}
	else {	return 1;	}
}

method backupFile ($filename) {
	my $counter = 1;
	my $backupfile = "$filename.$counter";
	while ( -f $backupfile )
	{
		$counter++;
		$backupfile = "$filename.$counter";
	}
	`cp $filename $backupfile`;

	print "Could not create backupfile: $backupfile\n" and exit if not -f $backupfile;
	print "backupfile created: $backupfile\n";
}


}


