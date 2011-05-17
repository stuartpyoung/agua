package Agua::Installer;

=head2

	PACKAGE		Agua::Installer

    PURPOSE

		1. INSTALL THE DEPENDENCIES FOR Agua

		2. CREATE THE REQUIRED DIRECTORY STRUCTURE

        <INSTALLDIR>/bin
                     cgi-bin  --> /var/www/cgi-bin/agua/<VERSION>
                        conf --> <INSTALLDIR>/conf
                        lib --> <INSTALLDIR>/lib
                        sql --> <INSTALLDIR>/bin/sql
                     conf
                     html --> /var/www/html/agua/<VERSION>
                     lib
                     t
=cut

use strict;
use warnings;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our $AUTOLOAD;

#### EXTERNAL MODULES
use File::Copy::Recursive;
use File::Path;
use Term::ReadKey;
use Data::Dumper;
use FindBin qw($Bin);

#### SET SLOTS
our @DATA = qw(

INSTALLDIR
WWWDIR
LOGFILE

);
our $DATAHASH;
foreach my $key ( @DATA )	{	$DATAHASH->{lc($key)} = 1;	}

sub install {
	my $self	=	shift;
	print "Agua::Installer::install    Agua::Installer::install\n";

	#### OPEN LOG
	my $logfile = $self->get_logfile() || "$Bin/log/install.log";
	$self->openLogfile($logfile);

	#### SET VERSION AND BUILD
	$self->setVersion("$Bin/resources/agua/VERSION");

	#### COPY .bash_profile
	$self->copyBashProfile();

	#### FIX getwcd (sh: getcwd() failed: No such file or directory)
	$self->fixGetcwd();

	#### INSTALL IPTABLES
	$self->installPackage("iptables");

	#### INSTALL EC2-API-TOOLS
	$self->installEc2("http://s3.amazonaws.com/ec2-downloads/ec2-api-tools.zip");

	### INSTALL MYSQL
	$self->installMysql();

	#### INSTALL EXPAT FOR XML::SIMPLE
	$self->installPackage("expat");

	#### INSTALL CURL
	$self->installPackage("curl");

	#### INSTALL cpanminus
	$self->runCommands(["curl -L http://cpanmin.us | perl - App::cpanminus"]);

	#### INSTALL PERL MODS
	$self->installPerlMods("$Bin/resources/agua/perlmods.txt");

	#### INSTALL APACHE2 IF NOT INSTALLED
	$self->installApache();

	#### LINK INSTALLATION DIRECTORIES TO WEB DIRECTORIES
	$self->linkDirectories();

	#### SET PERMISSIONS TO ALLOW ACCESS BY www-data USER
	$self->setPermissions("$Bin/resources/agua/permissions.txt");

	###### DO SETUID ON CGI-DIR SCRIPTS
	##$self->setSuid(['agua', "init.cgi", "workflow.cgi"]);

	#### CONFIRM INSTALLATION AND PRINT INFO
	$self->installConfirmation();	
}

sub fixGetcwd {
	my $self		=   shift;

	if ( not -d "/usr/bin/getcwd" )
	{
		$self->runCommands([
			"ln -s /bin/pwd /usr/bin/getcwd",
			"ln -s /bin/pwd /bin/getcwd"
		]);

	}
}
sub openLogfile {
	my $self		=   shift;
	my $logfile 	= 	shift;

	print "Agua::Installer::openLogfile    Install::openLogfile(logfile)\n";
	print "Agua::Installer::openLogfile    logfile: $logfile\n";
	$logfile .= ".1";

	if ( -f $logfile )
	{
		my ($stub, $index) = $logfile =~ /^(.+?)\.(\d+)$/;
		$index++;
		$logfile = $stub . "." . $index;
	}

	open (STDOUT, "| tee -ai $logfile") or die "Can't split STDOUT to logfile: $logfile\n";
	print "Agua::Installer::openLogfile    Writing to logfile: $logfile\n";

	return $logfile;	
}

sub closeLogfile {
	close (STDOUT);
}

sub getArch {	
	my $self		=	shift;
	my $command = "uname -a";
	print "Agua::Installer::getArch    command: $command\n";
	my $output = `$command`;

	#### Linux ip-10-126-30-178 2.6.32-305-ec2 #9-Ubuntu SMP Thu Apr 15 08:05:38 UTC 2010 x86_64 GNU/Linux
	return "ubuntu" if $output =~ /ubuntu/i;
	#### Linux ip-10-127-158-202 2.6.21.7-2.fc8xen #1 SMP Fri Feb 15 12:34:28 EST 2008 x86_64 x86_64 x86_64 GNU/Linux
	return "centos" if $output =~ /fc\d+/;
	return "linux";
}

sub setVersion {
	my $self		=	shift;
	my $versionfile	=	shift;

	print "Agua::Installer::setVersion    versionfile is not defined\n"
		and exit if not defined $versionfile;
	print "Agua::Installer::setVersion    Can't find versionfile: $versionfile\n"
		and exit if not -f $versionfile;

	my $contents = $self->contents($versionfile);
	print "Version versionfile is empty: $versionfile\n" 
		and exit if not defined $contents or not $contents;

	my ($version, $build) = $contents =~ /^([^-]+)-(.+)$/;
	print "Agua::Installer::setVersion    version: $version\n";
	print "Agua::Installer::setVersion    build: $build\n";

	$self->set_version($version);
	$self->set_build($build);	

	print "Agua::Installer::setVersion    self->get_version(): ", $self->get_version(), "\n";

}



sub copyBashProfile {
	my $self		=	shift;

	$self->replaceFile("~/.bash_profile", "$Bin/resources/starcluster/.bash_profile");
}

sub installEc2 {
	my $self		=	shift;
	my $url			=	shift;

	print "Agua::Installer::installEc2    Agua::Installer::installEc2()\n";

	#### INSTALL ec2-api-tools
	print "Agua::Installer::installEc2    Installing EC2-API-TOOLS\n";
	my $tempdir = "/tmp/ec2-api-tools";
	File::Path::mkpath($tempdir);
	chdir($tempdir);
	print `wget $url`;

	$self->runCommands(
		[
			"mkdir -p $tempdir",
			"cd $tempdir; wget $url",
			"unzip $tempdir/*zip",
			"ls $tempdir/ec2-api-tools-*/bin/* | grep -v \\\"\.cmd\$\\\" | xargs -n 1 -iFILE cp FILE /usr/bin;",
			"rm -fr $tempdir"
		]
	)	
}

sub installMysql {
	my $self		=	shift;

	print "Agua::Installer::installMysql    Agua::Installer::installMysql()\n";

	#### INSTALL MYSQL
	my $mysql = "/etc/init.d/mysql";
	#if ( not -f $mysql )
	#{
		print "Agua::Installer::installMysql    mysql not found: $mysql. Installing.\n";

		$self->runCommands([
			"rm -fr /var/lib/dpkg/lock",
			"dpkg --configure -a",
			"rm -fr /var/cache/apt/archives/lock",
			"sleep 3",

			"export DEBIAN_FRONTEND=noninteractive; sudo apt-get -y remove mysql-server libmysqlclient-dev mysql-server-5.1 mysql-cluster-server-5.1 libmysqlclient16-dev mysql-cluster-server libaprutil1-dev apache2-prefork-dev",
			"sleep 3",
			"rm -fr /var/lib/dpkg/lock",
			"dpkg --configure -a",
			"rm -fr /var/cache/apt/archives/lock",
			"sleep 3",

			"export DEBIAN_FRONTEND=noninteractive; sudo apt-get -y remove libmysqlclient-dev mysql-cluster-client mysql-cluster-client-5.1 mysql-cluster-server-5.1  mysql-common",

			"sleep 3",
			"rm -fr /var/lib/dpkg/lock",
			"dpkg --configure -a",
			"rm -fr /var/cache/apt/archives/lock",
			"sleep 3",

			"export DEBIAN_FRONTEND=noninteractive; sudo apt-get -y remove mysql-server libmysqlclient-dev mysql-server-5.1 mysql-cluster-server-5.1 libmysqlclient16-dev mysql-cluster-server libaprutil1-dev apache2-prefork-dev",

			#"sleep 10",
			#"rm -fr /var/lib/dpkg/lock",
			#"dpkg --configure -a",
			#"rm -fr /var/cache/apt/archives/lock",
			#"sleep 10",
			#
			#"export DEBIAN_FRONTEND=noninteractive; sudo apt-get -y install libmysqlclient16-dev libmysqlclient-dev mysql-server mysql-server-5.1"

		]);
	#}

	#### INSTALL PREREQUISITES FOR DBD::mysql PERL MODULE
	$self->runCommands([
		"rm -fr /var/lib/dpkg/lock",
		"dpkg --configure -a",
		"rm -fr /var/cache/apt/archives/lock",
		"sleep 5",
		"export DEBIAN_FRONTEND=noninteractive; sudo apt-get -y install libmysqlclient16-dev libmysqlclient-dev mysql-server mysql-server-5.1"
	]);

	##### START MYSQL AUTOMATICALLY AT BOOT
	#$self->installPackage("chkconfig");
	#$self->runCommands(["chkconfig --level 345 mysql on"]);

	#### RESTART MYSQL
	$self->runCommands(["/etc/init.d/mysql restart"]);
}


sub runCommands {
	my $self		=	shift;
	my $commands 	=	shift;
	print "Agua::Installer::runCommands    Agua::Installer::runCommands(commands)\n";
	foreach my $command ( @$commands )
	{
		print "Agua::Installer::runCommands    command: $command\n";		
		print `$command` or die("Error with command: $command\n$! , stopped");
	}
}

sub installPerlMods {
	my $self		=	shift;
	my $perlmodsfile=	shift;
	print "Agua::Installer::installPerlMods    Agua::Installer::installPerlMods(perlmodsfile)\n";
	print "Agua::Installer::installPerlMods    perlmodsfile: $perlmodsfile\n";
	my $contents = $self->contents($perlmodsfile);
	print "Agua::Installer::installPerlMods    perlmodsfile is empty: $perlmodsfile\n"
		and exit if not defined $contents or not $contents;

	#### INSTALL OTHER DEPENDENCIES
	my @modules = split "\n", $contents;
	foreach my $module ( @modules )
	{
		print "Agua::Installer::installPerlMods    Problem installing module $module\n"
			if not $self->cpanminusInstall($module);
	}

	#### BACKUP CPAN CONFIG
	#### my $configfile = "~/.cpan/CPAN/MyConfig.pm";
	#### $self->backupCpanConfig($configfile);
	#### my $replacefile = "$Bin/resources/cpan-MyConfig.pm";
	###$ self->replaceCpanConfig($configfile, $replacefile);
	####
	#### INSTALL USING CPAN
	#### foreach my $module ( @modules )
	#### {
	####	print "Agua::Installer::installPerlMods    Problem installing module $module\n"
	####		if not cpanInstall($module);
	#### }
	####
	#### RESTORE CPAN CONFIG FILE
	#### my $restorefile = "$Bin/resources/cpan-MyConfig.pm";
	#### $self->restoreCpanConfig($configfile, $restorefile);
}

sub installApache {
	my $self		=	shift;

	#### INSTALL APACHE2
	$self->installPackage("apache2");
	print "Agua::Installer::installApache    Finished installing apache2\n";	

	#### REPLACE mpm-worker WITH mpm-prefork (Non-threaded) CGI DAEMON
	#### TO AVOID THIS ERROR: 'unable to connect to cgi daemon'
	$self->removePackage("apache2-mpm-prefork");
	$self->installPackage("apache2-prefork-dev");

	#### REPLACE /etc/apache2/apache2.conf
	$self->backupFile("/etc/apache2/apache2.conf", "/etc/apache2/apache2.conf.ORIGINAL");
	$self->replaceFile("/etc/apache2/apache2.conf", "$Bin/resources/apache2/apache2.conf");

	#### REPLACE /etc/apache2/sites-available/default TO:
	#### 1. SET HTML ROOT
	#### 2. ENABLE CGI-BIN
	#### 3. ALLOW FOLLOW SYMLINKS IN CGI-BIN (AVOID ERROR: 'method PUT not allowed')
	$self->backupFile("/etc/apache2/sites-available/default", "/etc/apache2/sites-available/default.ORIGINAL");
	$self->replaceFile("/etc/apache2/sites-available/default", "$Bin/resources/apache2/sites-available/default");

	#### REPLACE /etc/apache2/envvars TO:
	#### 1. SET UMASK (DEFAULT 775/664 FOR NEW FILES/DIRS)
	#### 2. SET SGE ENVIRONMENT VARIABLES
	$self->backupFile("/etc/apache2/envvars", "/etc/apache2/envvars.ORIGINAL");
	$self->replaceFile("/etc/apache2/envvars", "$Bin/resources/apache2/envvars");

	#### CREATE AGUA DIRS
	my $wwwdir = $self->get_wwwdir();
	$self->createDir("$wwwdir/html/agua");
	$self->createDir("$wwwdir/cgi-bin/agua");

	#### START APACHE
	print "Agua::Installer::installApache    Starting apache\n";
	$self->runCommands(["/etc/init.d/apache2 restart"]);
}

sub createDir {
	my $self		=	shift;
	my $directory	=	shift;
	print "Agua::Installer::createDir    Agua::Installer::createDir(directory)\n";


	print "Agua::Installer::createDir    directory not defined\n"
		and exit if not defined $directory;
	print "Agua::Installer::createDir    directory is a file: $directory\n"
		and exit if -f $directory;
	return if -d $directory;

	print `mkdir -p $directory`;
	print "Agua::Installer::createDir    Can't create directory: $directory\n"
		and exit if not -d $directory;
}

sub replaceFile {
	my $self			=	shift;
	my $originalfile 	=	shift;
	my $replacementfile	=	shift;
	my $force			=	shift;

	$self->backupFile($replacementfile, $originalfile, 1);

	#
	#	and exit if not defined $originalfile;
	#	and exit if not defined $replacementfile;
	#	and exit if not -f $replacementfile;
	#
	#return if not defined $originalfile or not $originalfile;
	#return if not defined $replacementfile or not $replacementfile;
	#return if not -f $replacementfile;
	#
	#my ($originaldir) = $originalfile =~ /^(.+?)\/[^\/]+$/;
	#if ( not -d $originaldir )
	#{
	#	my $command = "mkdir -p $originaldir";
	#}
	#
	#my $command = "cp $replacementfile $originalfile";
	#`$command`;
}


sub backupFile {
	my $self			=	shift;
	my $originalfile 	=	shift;
	my $backupfile 		=	shift;
	my $force			=	shift;

	print "Agua::Installer::backupFile    originalfile not defined\n"
		and exit if not defined $originalfile;
	print "Agua::Installer::backupFile    backupfile not defined\n"
		and exit if not defined $backupfile;
	print "Agua::Installer::backupFile    Can't find originalfile: $originalfile\n"
		and exit if not -f $originalfile;
	print "Agua::Installer::backupFile    Skipping backup as backupfile already exists: : $backupfile\n" and return if -f $backupfile and not defined $force;

	my ($backupdir) = $backupfile =~ /^(.+?)\/[^\/]+$/;
	print "Agua::Installer::replaceFile    Creating backupdir: $backupdir\n";
	if ( not -d $backupdir )
	{
		my $command = "mkdir -p $backupdir";
		print "Agua::Installer::backupFile    command: $command\n";
		print `$command`;
		print "Can't create backupdir: $backupdir\n" if not -d $backupdir;
	}

	my $command = "cp $originalfile $backupfile";
	print "Agua::Installer::backupFile    command: $command\n";
	print `$command`;
}


sub installPackage  {
	my $self		=	shift;
    my $package     =   shift;

	print "Agua::Installer::installPackage    Agua::Installer::installPackage(package)\n";
    return 0 if not defined $package or not $package;
	print "Agua::Installer::installPackage    package: $package\n";

    if ( -f "/usr/bin/apt-get" )
    {
		$self->runCommands([
		"rm -fr /var/lib/dpkg/lock",
		"dpkg --configure -a",
		"rm -fr /var/cache/apt/archives/lock"
		]);

		$ENV{'DEBIAN_FRONTEND'} = "noninteractive";
		my $command = "/usr/bin/apt-get -q -y install $package";
		print "Agua::Installer::installPackage    command: $command\n";
		system($command);
		#die("Problem with command: $command\n$!\n") if $!;
    }
	elsif ( -f "/usr/bin/yum" )
	{
		my $command = "/usr/bin/yum -y install $package";
		print "Agua::Installer::installPackage    command: $command\n";
		system($command);
		#die("Problem with command: $command\n$!\n") if $!;
	}    
}

sub removePackage  {
	my $self		=	shift;
    my $package     =   shift;

	print "Agua::Installer::removePackage    Agua::Installer::removePackage(package)\n";
    return 0 if not defined $package or not $package;
	print "Agua::Installer::removePackage    package: $package\n";

    if ( -f "/usr/bin/apt-get" )
    {
		$self->runCommands([
		"rm -fr /var/lib/dpkg/lock",
		"dpkg --configure -a",
		"rm -fr /var/cache/apt/archives/lock"
		]);

		$ENV{'DEBIAN_FRONTEND'} = "noninteractive";
		my $command = "/usr/bin/apt-get -q -y --purge remove $package";
		print "Agua::Installer::removePackage    command: $command\n";
		system($command);
		#die("Problem with command: $command\n$!\n") if $!;
    }
	elsif ( -f "/usr/bin/yum" )
	{
		my $command = "/usr/bin/yum -y remove $package";
		print "Agua::Installer::removePackage    command: $command\n";
		system($command);
		#die("Problem with command: $command\n$!\n") if $!;
	}    
}



sub cpanInstall {
	my $self		=	shift;
	my $module =    shift;
    my $logfile =    shift;

	print "Agua::Installer::cpanInstall    Agua::Installer::cpanInstall(module)\n";
    print "Agua::Installer::cpanInstall    module: $module\n";
	print "Agua::Installer::cpanInstall    logfile: $logfile\n" if defined $logfile;

	return 0 if not defined $module or not $module;

	my $command = "/usr/local/bin/cpanm $module";

	#my $command = "dh-make-pl $module; dpkg -i $module.deb";
    #my $command = "PERL_MM_USE_DEFAULT=1 /usr/bin/perl -MCPAN -e 'install $module'";
    #$command .= " &>> $logfile"  if defined $logfile;
    print "Agua::Installer::cpanInstall    command: $command\n";
    print `$command`;
}

sub cpanminusInstall {
	my $self		=	shift;
	my $module 		=    shift;
	return 0 if not defined $module or not $module;

	my $command = "/usr/local/bin/cpanm $module";
    print `$command`;
}

sub setSuid {
	my $self		=	shift;
	my $executables	=	shift;

    my $cgidir      =   $self->wwwdir() . "/cgi-bin";

	foreach my $executable ( @$executables )
	{
		#### SET chown
		my $command = "chown root $cgidir/$executable";
		print "$command\n";
		print `$command`;

		#### SETUID AND SETGID
		$command = "chmod u+s $cgidir/$executable";
		$command = "chmod g+s $cgidir/$executable";
	}
}    


sub linkDirectories {
	my $self		=	shift;

    my $installdir  =   $self->get_installdir();
    my $wwwdir      =   $self->get_wwwdir();
	my $version		=	$self->get_version();
	my $build		=	$self->get_build();

	print "Agua::Installer::linkDirectories    installdir not defined or empty\n"
		and exit if not defined $installdir or not $installdir;
	print "Agua::Installer::linkDirectories    wwwdir not defined or empty\n"
		and exit if not defined $wwwdir or not $wwwdir;
	print "Agua::Installer::linkDirectories    version not defined or empty\n"
		and exit if not defined $version or not $version;

	$self->removeDir("$installdir/$version");
    $self->removeLink("$installdir/$version");
    $self->addLink("$installdir/$version-$build", "$installdir/$version");

	#### REMOVE EXISTING LINKS
    print "Agua::Installer::linkDirectories    Removing any existing links\n";
    $self->removeLink("$wwwdir/html/agua/$version");
    $self->removeLink("$wwwdir/cgi-bin/agua/$version");
    $self->removeLink("$installdir/$version/cgi-bin/lib");
    $self->removeLink("$installdir/$version/cgi-bin/sql");
    $self->removeLink("$installdir/$version/cgi-bin/conf");

    #### LINK WEB DIR AND CGI DIR
    print "Agua::Installer::linkDirectories    Creating links\n";
    $self->addLink("$installdir/$version/html", "$wwwdir/html/agua/$version");
    $self->addLink("$installdir/$version/cgi-bin", "$wwwdir/cgi-bin/agua/$version");
    $self->addLink("$installdir/$version/lib", "$installdir/$version/cgi-bin/lib");
    $self->addLink("$installdir/$version/bin/sql", "$installdir/$version/cgi-bin/sql");
    $self->addLink("$installdir/$version/conf", "$installdir/$version/cgi-bin/conf");

	print "Agua::Installer::linkDirectories    Completed\n"
}

sub setPermissions {
	my $self		=	shift;
	my $permissionsfile=	shift;
	print "Agua::Installer::setPermissions    Agua::Installer::setPermissions(permissionsfile)\n";
	print "Agua::Installer::setPermissions    permissionsfile: $permissionsfile\n";
	my $contents = $self->contents($permissionsfile);
	print "Agua::Installer::setPermissions    permissionsfile is empty: $permissionsfile\n"
		and exit if not defined $contents or not $contents;

	my $installdir = $self->get_installdir();
	my $version = $self->get_version();
	my $base = "$installdir/$version";

	#### INSTALL OTHER DEPENDENCIES
	my @commands = split "\n", $contents;
	foreach my $command ( @commands )
	{
		$command =~ s/BASE/$base/;
		print "$command\n";
		next if $command =~ /^#/;
		print `$command`;
	}
}

sub installConfirmation {
	my $self		=	shift;

    my $installdir  =   $self->get_installdir();
    my $version  	=   $self->get_version();

    print "\n";
    print "Agua has been installed to this directory:\n";
    print "\n";
    print "    $installdir/$version\n";
    print "\n";
    print "You can browse to agua here:\n";
    print "\n";
    print "    http://localhost/agua/$version\n";
    print "		(If you are on the installation host.)\n";
    print "\n";
    print "    http://<hostname>/agua/$version\n";
    print "		(Where <hostname> is the URL of the installation host.)\n";
    print "\n";
    print "*******************************************************\n";
    print "*******************************************************\n";
    print "\n";
    exit;    
}


sub input {
	my $self		=	shift;
	my $message		=	shift;
	return if ( not defined $message );
    print "$message\n";	

	$/ = "\n";
	my $input = <STDIN>;
	while ( $input =~ /^\s*$/ )
	{
		print "$message\n";
		$input = <STDIN>;
	}	

    chop($input);
	return $input;
}

sub addLink {
	my $self		=	shift;
    my $source      =   shift;
    my $target      =   shift;
	print "Agua::Installer::addLink    source: $source\n";
	print "Agua::Installer::addLink    target: $target\n";    
	print "Agua::Installer::addLink    symlink($source, $target)\n";    
	print `ln -s $source $target`;
    print "Agua::Installer::addLink    Could not create link: $target\n"
		and exit if not -l $target;
}

sub removeLink {
	my $self		=	shift;
    my $target      =   shift;
	print "Agua::Installer::removeLink    target: $target\n";

	return if not -l $target;    
	print "Agua::Installer::removeLink    unlink($target)\n";    
    unlink($target);
    print "Agua::Installer::removeLink    Could not unlink: $target\n"
		and exit if -l $target;
}

sub removeDir {
	my $self		=	shift;
    my $target      =   shift;

	return if -l $target || -f $target;    
    print `$target`;
    print "Agua::Installer::removeLink    Could not remove target: $target\n"
		and exit if -l $target;
}


sub yes {
	my $self		=	shift;
	my $message		=	shift;
	return if ( not defined $message );
    my $max_times = 10;

	$/ = "\n";
	my $input = <STDIN>;
	my $counter = 0;
	while ( $input !~ /^Y$/i and $input !~ /^N$/i )
	{
		if ( $counter > $max_times )	{	print "Exceeded 10 tries. Exiting...\n";	}

		print "$message\n";
		$input = <STDIN>;
		$counter++;
	}	

	if ( $input =~ /^N$/i )	{	return 0;	}
	else {	return 1;	}
}



sub contents {
	my $self		=	shift;
	my $file		=	shift;

	print "Agua::Installer::contents    Agua::Installer::contents(file)\n";
	print "Agua::Installer::contents    file: $file\n";

	die("Agua::Installer::contents    file not defined\n") if not defined $file;
	die("Agua::Installer::contents    Can't find file: $file\n$!") if not -f $file;


	my $temp = $/;
	$/ = undef;
	open(FILE, $file) or die("Can't open file: $file\n$!");
	my $contents = <FILE>;
	close(FILE);
	$/ = $temp;

	return $contents;
}

sub backupCpanConfig {
	my $self		=	shift;
	my $configfile 		=	shift;
	print "Agua::Installer::backupCpanConfig    configfile: $configfile\n";

	return if not defined $configfile or not $configfile;
	my $backupfile = "$configfile.original"; 
	print "Agua::Installer::backupCpanConfig    backupfile: $backupfile\n";
	if ( not -f $backupfile and -f $configfile )
	{
		my $command = "cp $configfile $backupfile";
		print "Agua::Installer::backupCpanConfig    command: $command\n";
		`$command`;
	}
}

sub restoreCpanConfig {
	my $self		=	shift;
	my $configfile 		=	shift;
	print "Agua::Installer::restoreCpanConfig    configfile: $configfile\n";

	return if not defined $configfile or not $configfile;
	if ( -f $configfile )
	{
		print "Agua::Installer::restoreCpanConfig    configfile: $configfile\n";
		my $command = "cp $configfile.original $configfile";
		print "Agua::Installer::restoreCpanConfig    command: $command\n";
		`$command`;
	}
}

sub replaceCpanConfig {
	my $self		=	shift;
	my $configfile 		=	shift;
	my $replacement 	=	shift;
	print "Agua::Installer::replaceCpanConfig    configfile: $configfile\n";
	print "Agua::Installer::replaceCpanConfig    replacement: $replacement\n";
	return if not defined $configfile or not $configfile;
	return if not defined $replacement or not $replacement;
	return if not -f $replacement;

	my ($configdir) = $configfile =~ /^(.+?)\/[^\/]+$/;
	print "Agua::Installer::replaceCpanConfig    Creating configdir: $configdir\n";
	if ( not -d $configdir )
	{
		my $command = "mkdir -p $configdir";
		print "Agua::Installer::replaceCpanConfig    command: $command\n";
		print `$command`;
		print "Can't create configdir: $configdir\n" if not -d $configdir;
	}

	my $command = "cp $replacement $configfile";
	print "Agua::Installer::replaceCpanConfig    command: $command\n";
	`$command`;
}




sub getfiles {
	my $self		=	shift;
    my $directory   =   shift;
    my $suffix      =   shift;


    opendir(DIR, $directory) or die "Can't open directory: $directory. $!";
    my @files = readdir(DIR);
    closedir(DIR) or die "Can't close directory: $directory. $!";

    return \@files if not defined $suffix;

    for ( my $i = 0; $i < $#files + 1; $i++ )
    {
        use re 'eval';
        if ( $files[$i] !~ /$suffix$/ or not -f "$directory/$files[$i]" )
        {
            splice(@files, $i, 1);
            $i--;
        }
        no re 'eval';
    }

    return \@files;
}

sub copyDirectories {
	my $self		=	shift;
    my $installdir     =   shift;
    my $directories =   shift;

    #### COPY ALL FILES TO BASE DIR
    print "\nCopying folders to base directory: @$directories\n";
    foreach my $directory ( @$directories )
    {
        my $sourcedir = "$Bin/../../$directory";
        print "Copying $sourcedir TO $installdir/$directory\n";
        my $targetdir = "$installdir/$directory";

        #### PRINT '.' AS PROGRESS COUNTER
        print ".";

        if ( -d $sourcedir )
        {
            my $success = File::Copy::Recursive::rcopy("$sourcedir", "$targetdir");
            if ( not $success )
            {
                die "Could not copy directory '$sourcedir' to '$targetdir': $!\n";
            }
        }
        else
        {
            die "Directory is missing from agua distribution: $sourcedir\n";
        }
    }
}


################################################################################
##################			HOUSEKEEPING SUBROUTINES			################
################################################################################
=head2

	SUBROUTINE		initialise

	PURPOSE

		INITIALISE THE self OBJECT WITH USER-INPUT ARGUMENT VALUES

=cut

sub initialise {
    my $self		=	shift;
	my $arguments	=	shift;

    #### VALIDATE USER-PROVIDED ARGUMENTS
	($arguments) = $self->validate_arguments($arguments, $DATAHASH);	

	#### LOAD THE USER-PROVIDED ARGUMENTS
	foreach my $key ( keys %$arguments )
	{		
		#### LOAD THE KEY-VALUE PAIR
		$self->value($key, $arguments->{$key});
	}
}

=head2

	SUBROUTINE		new

	PURPOSE

		CREATE THE NEW self OBJECT AND INITIALISE IT, FIRST WITH DEFAULT 

		ARGUMENTS, THEN WITH PROVIDED ARGUMENTS

        THE PASSWORD FOR THE DEFAULT ROOT USER CAN BE SET ON THE MYSQL

        COMMAND LINE:

        use myEST;
        insert into users values('admin', 'mypassword', now());

=cut

sub new {
 	my $class 		=	shift;
	my $arguments 	=	shift;

	my $self = {};
    bless $self, $class;

	#### INITIALISE THE OBJECT'S ELEMENTS
	$self->initialise($arguments);

    return $self;
}


=head2

	SUBROUTINE		value

	PURPOSE

		SET A PARAMETER OF THE self OBJECT TO A GIVEN value

    INPUT

        1. parameter TO BE SET

		2. value TO BE SET TO

    OUTPUT

        1. THE SET parameter INSIDE THE BioSVG OBJECT

=cut
sub value {
    my $self		=	shift;
	my $parameter	=	shift;
	my $value		=	shift;

	$parameter = lc($parameter);

    if ( defined $value)
	{	
		$self->{"_$parameter"} = $value;
	}
}

=head2

	SUBROUTINE		validate_arguments

	PURPOSE

		VALIDATE USER-INPUT ARGUMENTS BASED ON

		THE HARD-CODED LIST OF VALID ARGUMENTS

		IN THE data ARRAY
=cut

sub validate_arguments {
	my $self		=	shift;
	my $arguments	=	shift;
	my $DATAHASH	=	shift;

	my $hash;
	foreach my $argument ( keys %$arguments )
	{
		if ( $self->is_valid($argument, $DATAHASH) )
		{
			$hash->{$argument} = $arguments->{$argument};
		}
		else
		{
			warn "'$argument' is not a known parameter\n";
		}
	}

	return $hash;
}

=head2

	SUBROUTINE		is_valid

	PURPOSE

		VERIFY THAT AN ARGUMENT IS AMONGST THE LIST OF

		ELEMENTS IN THE GLOBAL '$DATAHASH' HASH REF

=cut

sub is_valid {
	my $self		=	shift;
	my $argument	=	shift;
	my $DATAHASH	=	shift;

	#### REMOVE LEADING UNDERLINE, IF PRESENT
	$argument =~ s/^_//;

	#### CHECK IF ARGUMENT FOUND IN '$DATAHASH'
	if ( exists $DATAHASH->{lc($argument)} )
	{
		return 1;
	}

	return 0;
}

=head2

	SUBROUTINE		AUTOLOAD

	PURPOSE

		AUTOMATICALLY DO 'set_' OR 'get_' FUNCTIONS IF THE

		SUBROUTINES ARE NOT DEFINED.

=cut

sub AUTOLOAD {
    my ($self, $newvalue) = @_;
    my ($operation, $attribute) = ($AUTOLOAD =~ /(get|set)(_\w+)$/);

    # Is this a legal method name?
    unless ( defined $operation && $operation && defined $attribute && $attribute ) {
        print "Method name $AUTOLOAD is not in the recognized form (get|set)_attribute\n"
			and exit;
    }

    unless( exists $self->{$attribute} or $self->is_valid($attribute) )
	{
		#if ( not defined $operation )
		#{
	        #die "No such attribute '$attribute' exists in the class ", ref($self);
			#return;
		#}
    }

    # Turn off strict references to enable "magic" AUTOLOAD speedup
    no strict 'refs';

    # AUTOLOAD accessors
    if($operation eq 'get') {
        # define subroutine
        *{$AUTOLOAD} = sub { shift->{$attribute} };

    # AUTOLOAD mutators
    }elsif($operation eq 'set') {
        # define subroutine4

        *{$AUTOLOAD} = sub { shift->{$attribute} = shift; };

        # set the new attribute value
        $self->{$attribute} = $newvalue;
    }

    # Turn strict references back on
    use strict 'refs';

    # return the attribute value
    return $self->{$attribute};
}


=head 2

	SUBROUTINE		DESTROY

	PURPOSE

		When an object is no longer being used, this will be automatically called

		and will adjust the count of existing objects

=cut
sub DESTROY {
    my($self) = @_;

	#### TIDY UP, DISCONNECT DATABASE HANDLES, ETC.
	$self->closeLogfile();
}



=head1 LICENCE

This code is released under the GPL, a copy of which should
be provided with the code.

=end pod

=cut

1;
