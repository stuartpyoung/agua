package AWS;


=head2

	PACKAGE		Init

	PURPOSE

		CARRY OUT INITIAL Agua SETUP:

		1. CREATE, STOP/START/TERMINATE AND OTHERWISE MANAGE

			AWS AMIs, MANAGE VOLUMES AND OTHER RESOURCES

		2. STORE AWS CREDENTIALS INFO IN FILES

		3. STORE AWS CREDENTIALS INFO IN aws DATABASE TABLE

=cut

use strict;
use warnings;
use Carp;

#### USE LIB (FOR INHERITANCE)
use lib "../";
use lib "../external";

#### INHERIT FROM App CLASS
require Exporter;
our $AUTOLOAD;

#### EXTERNAL MODULES
use Data::Dumper;
use File::Copy::Recursive;
use File::Path;

#### INTERNAL MODULES
use Util;


#### SET SLOTS
our @DATA;

if ( 1 ) {

@DATA = qw(

USERNAME
AMAZONUSERID
AWSACCESSKEYID
SECRETACCESSKEY
EC2PUBLICCERT
EC2PRIVATEKEY
DATAVOLUME
DATAVOLUMESIZE
USERVOLUME
USERVOLUMESIZE
DATAVOLUMECHECKBOX
USERVOLUMECHECKBOX
CONF

);

}
our $DATAHASH;
foreach my $key ( @DATA )	{	$DATAHASH->{lc($key)} = 1;	}

=head2

	SUBROUTINE		init

	PURPOSE

		1. LOAD USER DATA INTO AGUA DATABASE

=cut

sub init {
	my $self	=	shift;
	print "Agua::Init::init     AWS::init()\n";

	#### 	1. SET ENVIRONMENT VARIABLES
	####
	print "Agua::Init::init    Doing setEnvironment()\n";
	$self->setEnvironment();

	####	2. PRINT X.509 PUBLIC CERT AND PRIVATE KEY 
	#### 		FILES FOR HTTPS AND AWS
	####	
	print "Agua::Init::init    Doing printKeyfiles()\n";
	$self->printKeyfiles();

	####	3. MOUNT /data VOLUME CONTAINING AGUA DATA
	####
	print "Agua::Init::init    Doing mountData()\n";
	$self->mountData();

	####	4. MOUNT /nethome USER DATA 
	####	
	print "Agua::Init::init    Doing mountNethome()\n";
	$self->mountNethome();

	##### 	5. MOUNT MYSQL AND START
	print "Agua::Init::init    Doing mountMysql()\n";
	$self->mountMysql();

	####	6. ADD nfs INFO TO /etc/fstab FOR LATER EXPORTS
	$self->addNfsToFstab();

	#### 	7. COPY KEYS FROM INSTALL DIR TO USER DIR
	print "Agua::Init::init    Doing copyKeyfiles()\n";
	$self->copyKeyfiles();

	####	8. GENERATE NEW PUBLIC CERTIFICATE (HTTPS)
	print "Agua::Init::init    Doing generateCACert()\n";
	$self->generateCACert();

	####    9. RESTART HTTPD SERVER
	print "Agua::Init::init    Doing restartHttpd()\n";
	$self->restartHttpd();

	###    10. SET STARCLUSTER KEYPAIR AND CONFIG FILE
	print "Agua::Init::init    Doing configStarCluster()\n";
	$self->configStarcluster();	
}

=head2

	SUBROUTINE		addNfsToFstab

	PURPOSE

		ADD nfs MOUNT INFO TO /etc/fstab. THE MOUNTPOINTS WILL 

		BE EXPORTED TO THE MASTER AND EXEC NODES OF EACH NEW

		STARCLUSTER.

=cut

sub addNfsToFstab {
	my $self	=	shift;
	my $conf 	= 	$self->get_conf();

	my $datadir = $conf->getKeyValue('agua', "DATADIR");
	my $userdir = $conf->getKeyValue('agua', "USERDIR");
	my $datadevice = $conf->getKeyValue("aws", "DATADEVICE");
	my $userdevice = $conf->getKeyValue('agua', "USERDEVICE");

	#### ADD TO FSTAB
	$self->addToFstab("$datadevice	$datadir	nfs     rw,vers=3,rsize=32768,wsize=32768,hard,proto=tcp 0 0");
	$self->addToFstab("$userdevice	$userdir	nfs     rw,vers=3,rsize=32768,wsize=32768,hard,proto=tcp 0 0");

	print "Agua::AWS::addNfsToFstab    Completed\n";
}


=head2

	SUBROUTINE 		configStarcluster

	PURPOSE

		1. COPY DEFAULT CONFIG FILE TO GIVEN USER'S CONFIG FILE

		2. CREATE KEYPAIR FILE FROM PRIVATE AND PUBLIC KEYS AND

		3. ADD KEYPAIR FILE, AWS ACCESS IDs AND OTHER CREDENTIAL

			INFO TO USER'S CLUSTER CONFIG

=cut

sub configStarcluster {
	my $self	=	shift;

	#### 1. COPY CONFIG FILE
	#### 	AND TRANSFER TO CLUSTER USING Conf::StarCluster
	$self->copyClusterConfig();

	#### 2. CREATE KEYPAIR FILE FROM PRIVATE KEY AND PUBLIC KEYS
	#### 	AND TRANSFER TO CLUSTER USING Conf::StarCluster
	$self->generateClusterKeypair();

	#### 3. ADD ACCESS KEYS AND KEYPAIR INFO
	$self->initialiseClusterConfig();
}


=head2

	SUBROUTINE 		copyClusterConfig

	PURPOSE

		CREATE KEYPAIR FILE FROM PRIVATE KEY AND PUBLIC KEYS

		AND TRANSFER TO CLUSTER USING Conf::StarCluster

=cut

sub copyClusterConfig {
	my $self	=	shift;

	#### PASSED VARIABLES
	my $conf 			= 	$self->get_conf();
	my $username 		=	$self->get_username();

	#### SET FILES
	my $installdir 		= 	$conf->getKeyValue('agua', "INSTALLDIR");
	my $userdir 		= 	$conf->getKeyValue('agua', "USERDIR");

	#### RUN starcluster.pl TO GENERATE KEYPAIR FILE IN .starcluster DIR
	my $starcluster = "$installdir/bin/apps/cluster/starcluster.pl";
	my $command = qq{$starcluster copyConfig \\
--username $username
};
	print "$command\n";
	print `$command`;
}

=head2

	SUBROUTINE 		generateClusterKeypair

	PURPOSE

		CREATE KEYPAIR FILE FROM PRIVATE KEY AND PUBLIC KEYS

		AND TRANSFER TO CLUSTER USING Conf::StarCluster

=cut

sub generateClusterKeypair {
	my $self	=	shift;

	#### PASSED VARIABLES
	my $conf 			= 	$self->get_conf();
	my $username 		=	$self->get_username();
	my $keyname			=	$self->get_keyname();

	#### SET DEFAULT keyname
	$keyname = "$username-key" if not defined $keyname;

	#### SET FILES
	my $installdir 		= 	$conf->getKeyValue('agua', "INSTALLDIR");
	my $userdir 		= 	$conf->getKeyValue('agua', "USERDIR");
	my $filedir			=	"$userdir/$username/.keypairs";
	my $privatekey		=	"$filedir/private.pem";
	my $publiccert		=	"$filedir/public.pem";

	#### RUN starcluster.pl TO GENERATE KEYPAIR FILE IN .starcluster DIR
	my $starcluster = "$installdir/bin/apps/cluster/starcluster.pl";
	my $command = qq{$starcluster generateKeypair \\
--privatekey $privatekey \\
--publiccert $publiccert \\
--keyname $keyname \\
--username $username
};
	print "$command\n";
	print `$command`;
}


=head2

	SUBROUTINE 		initialiseClusterConfig

	PURPOSE

		CREATE KEYPAIR FILE FROM PRIVATE KEY AND PUBLIC KEYS

		AND TRANSFER TO CLUSTER USING Conf::StarCluster

=cut

sub initialiseClusterConfig {
	my $self	=	shift;

	#### PASSED VARIABLES
	my $conf 			= 	$self->get_conf();
	my $username 		=	$self->get_username();		
	my $amazonuserid 	=	$self->get_amazonuserid();
	my $accesskeyid 	=	$self->get_awsaccesskeyid();
	my $secretaccesskey =	$self->get_secretaccesskey();
	my $keyname			=	$self->get_keyname();

	#### SET DEFAULT keyname
	$keyname = "$username-key" if not defined $keyname;

	my $userdir 		= 	$conf->getKeyValue('agua', "USERDIR");
	my $filedir			=	"$userdir/$username/.keypairs";
	my $privatekey		=	"$filedir/private.pem";
	my $publiccert		=	"$filedir/public.pem";

print "Agua::Init::initialiseClusterConfig    privatekey: $privatekey\n" if defined $privatekey;
print "Agua::Init::initialiseClusterConfig    publiccert: $publiccert\n" if defined $publiccert;
print "Agua::Init::initialiseClusterConfig    username: $username\n" if defined $username;
print "Agua::Init::initialiseClusterConfig    amazonuserid: $amazonuserid\n" if defined $amazonuserid;
print "Agua::Init::initialiseClusterConfig    accesskeyid: $accesskeyid\n" if defined $accesskeyid;
print "Agua::Init::initialiseClusterConfig    secretaccesskey: $secretaccesskey\n" if defined $secretaccesskey;
print "Agua::Init::initialiseClusterConfig    keyname: $keyname\n" if defined $keyname;

	#### SET FILES
	my $installdir 		= 	$conf->getKeyValue('agua', "INSTALLDIR");

	#### RUN starcluster.pl TO GENERATE KEYPAIR FILE IN .starcluster DIR
	my $starcluster = "$installdir/bin/apps/cluster/starcluster.pl";
	my $command = qq{$starcluster writeConfigfile \\
--privatekey $privatekey \\
--publiccert $publiccert \\
--amazonuserid $amazonuserid \\
--accesskeyid $accesskeyid \\
--secretaccesskey $secretaccesskey \\
--keyname $keyname \\
--username $username
};
	print "$command\n";
	print `$command`;
}


=head2

	SUBROUTINE 		copyKeyfiles

	PURPOSE

		COPY PRIVATE KEY AND PUBLIC CERT FROM INSTALL DIR TO USER DIR

=cut
sub copyKeyfiles {
	my $self	=	shift;

	#### GET FILES AND DIRECTORIES
	my $conf 			= 	$self->get_conf();
	my $username 		=	$self->get_username();
	my $privatefile 	=	$self->get_privatefile();
	my $publicfile 		=	$self->get_publicfile();
	my $userdir 		= 	$conf->getKeyValue('agua', "USERDIR");
	print "Agua::Init::copyKeyfiles    privatefile: $privatefile\n";

	#### CREATE TARGET DIRECTORY IF NOT EXISTS
	my $targetdir = "$userdir/$username/.keypairs";
	print "Agua::Init::copyKeyfiles    targetdir: $targetdir\n";
	File::Path::mkpath($targetdir) if not -d $targetdir;
	print "Agua::Init::copyKeyfiles    Can't create targetdir: $targetdir\n" and exit if not -d $targetdir;	

	#### COPY PRIVATE KEY AND PUBLIC CERT FROM INSTALL DIR
	#### TO USER DIR
	my $copykey = "cp -f $privatefile $targetdir";
	print "Agua::Init::copyKeyfiles    copykey: $copykey\n";
	print `$copykey`;

	my $copycert = "cp -f $publicfile $targetdir";
	print "Agua::Init::copyKeyfiles    copycert: $copycert\n";
	print `$copycert`;
	print "Agua::Init::copyKeyfiles    Completed\n";

	my $chmod = "chmod 600 $targetdir/*.pem";
	print "$chmod\n";
	print `$chmod`;
}

=head2

	SUBROUTINE 		generatePrivateKey

	PURPOSE

		GENERATE PRIVATE KEY FILE

=cut
sub generatePrivateKey {
	my $self	=	shift;
	my $keyname	=	shift;

	#### SET FILE NAMES
	my $conf 			= 	$self->get_conf();
	my $username 		=	$self->get_username();
	my $userdir 		= 	$conf->getKeyValue('agua', "USERDIR");
	my $filedir			=	"$userdir/$username/.starcluster";
	my $privatekey		=	"$filedir/privatekey.pem";
	my $publiccert		=	"$filedir/publiccert.pem";
	print "StarCluster::generatePrivateKey    filedir: $filedir\n";
	print "StarCluster::generatePrivateKey    privatekey: $privatekey\n";

	#### SET DEFAULT keyname
	$keyname = "$username-key" if not defined $keyname;

	#### 1. GENERATE PRIVATE KEY
	chdir($filedir);
	my $command = qq{openssl genrsa -out $privatekey-pass 1024};
	print "$command\n";
	print `$command`;

	#### 2. REMOVE PASS-PHRASE FROM KEY:
	my $remove = "openssl rsa -in $privatekey-pass -out $privatekey";
	print "$remove\n";
	print `$remove`;
    #### writing RSA key

	#### 3. SET PERMISSIONS TO 600
	my $chmod = "chmod 600 $privatekey";
	print "$chmod\n";
	print `$chmod`;

	#### 4. CREATE A PUBLIC CERTIFICATE FOR THE PRIVATE KEY
	my $public = qq{openssl rsa -in $privatekey \\
-pubout \\
-out $publiccert
};
	print "$public\n";
	print `$public`;

	#### 4. IMPORT PUBLIC KEY TO AMAZON
	my $import = qq{ec2-import-keypair \\
--debug \\
$keyname \\
--public-key-file $publiccert \\
-U https://ec2.amazonaws.com 
};
	print "$import\n";
	print `$import`;
}



=head2

	SUBROUTINE 		generateCACert

	PURPOSE

		GENERATE AUTHENTICATED CERTIFICATE FILE USING GIVEN PRIVATE KEY

		AND COPY TO APACHE conf DIR

		(NB: MUST RESTART APACHE TO USE NEW PUBLIC CERT)


	NOTES

		1. GET DOMAIN NAME

		2. CREATE CONFIG FILE

		3. GENERATE CERTIFICATE REQUEST

			openssl req

		4. GENERATE PUBLIC CERTIFICATE

			openssl x509 -req

		5. COPY TO APACHE AND RESTART APACHE

		NB: APACHE ssl.conf FILE SHOULD LOOK LIKE THIS:

			...		
			SSL Virtual Hosts
			<IfDefine SSL>

			<VirtualHost _default_:443>
			ServerAdmin webmaster@domain.com
			...
			SSLEngine on
			SSLCertificateFile /etc/httpd/conf/ssl.crt/server.crt
			SSLCertificateKeyFile /etc/httpd/conf/ssl.key/server.key
			SetEnvIf User-Agent ".*MSIE.*" nokeepalive ssl-unclean-shutdown
			CustomLog /var/log/httpd/ssl_request_log \
				"%t %h %{SSL_PROTOCOL}x %{SSL_CIPHER}x \"%r\" %b"
			</VirtualHost>

			</IfDefine>
			...

=cut
sub generateCACert {
	my $self	=	shift;

	#### SET FILE NAMES
	my $conf 				=	$self->get_conf();
	my $username 			=	$self->get_username();
	my $distinguished_name 	= 	"agua_" . $username . "_DN";
	my $userdir 			=	$conf->getKeyValue('agua', "USERDIR");
	my $apachedir 			= 	$conf->getKeyValue('agua', "APACHECONF");
	$apachedir =~ s/\/[^\/]+$//;

	my $filedir				=	"$userdir/$username/.starcluster";
	my $privatekey			=	"$filedir/privatekey.pem";
	my $pipefile			=	"$filedir/pipeline.pem";
	my $CA_certfile			=	"$filedir/CA-cert.pem";
	my $configfile			=	"$filedir/config.txt";

	#### MAKE DIRECTORY
	File::Path::mkpath($filedir) if not -d $filedir;

	#### COPY PRIVATE KEY FROM INSTALL DIR TO USER DIR
	my $privatefile = $self->get_privatefile();
	my $copy = "cp -f $privatefile $privatekey-pass";
	print "$copy\n";
	print `$copy`;	

	#### REMOVE PASS-PHRASE FROM PRIVATE KEY
	print "Conf::StarCluster::generateCACert    Removing pass-phrase from key\n";
	my $remove = "openssl rsa -in $privatekey-pass -out $privatekey";
	print "$remove\n";
	print `$remove`;
    #### writing RSA key

	#### 1. GET DOMAIN NAME
	my $instanceid = `curl -s http://169.254.169.254/latest/meta-data/instance-id`;
	$instanceid = $self->untaint($instanceid);
	my $domainname = `ec2-describe-instances | grep $instanceid | cut -f 4`;
	$domainname = $self->untaint($domainname);
	print "Conf::StarCluster::generateCACert    instanceid: $instanceid\n";
	print "Conf::StarCluster::generateCACert    domainname: $domainname\n";

	#### 2. CREATE CONFIG FILE
	open(OUT, ">$configfile") or die "Can't open configfile: $configfile\n"; 	
	print OUT qq{# SSL server cert/key parms
# Cert extensions
subjectKeyIdentifier    = hash
authorityKeyIdentifier  = keyid:always,issuer:always
basicConstraints        = CA:false
nsCertType              = server
# openssl req
[req]
default_bits            = 1024
prompt                  = no
distinguished_name      = $distinguished_name
# DN fields for SSL Server cert
[$distinguished_name]
C                       = US
ST                      = Maryland
O                       = UMCP/OIT/TSS/EIS
CN                      = $domainname
emailAddress            = trash\@trash.com
};
	close(OUT) or die "Can't close configfile: $configfile\n";

	#### 3. GENERATE CERTIFICATE REQUEST
	chdir($filedir);
	my $request = qq{openssl \\
req \\
-config $configfile \\
-newkey rsa:1024 \\
-key $privatekey \\
-out $pipefile
};
	print "Conf::StarCluster::generateCACert    request: $request\\n";
	`$request`;
	print "Conf::StarCluster::generateCACert    Can't find pipefile: $pipefile\n" and exit if not -f $pipefile;

	#### 4. GENERATE PUBLIC CERTIFICATE
	chdir($filedir);
	my $certify = qq{openssl \\
x509 -req \\
-extfile $configfile \\
-days 730 \\
-signkey $privatekey \\
-in $pipefile \\k
-out $CA_certfile
};
	print "Conf::StarCluster::generateCACert    certify: $certify\n";
	`$certify`;
	print "Conf::StarCluster::generateCACert    Can't find CA_certfile: $CA_certfile\n" and exit if not -f $CA_certfile;


	#### COPY THE PRIVATE KEY AND CERTIFICATE TO APACHE
	my $keydir = "$apachedir/ssl.key";
	print "Conf::StarCluster::generateCACert    keydir: $keydir\n";
	File::Path::mkpath($keydir) if not -d $keydir;
	print "Conf::StarCluster::generateCACert    Can't create keydir: $keydir\n"
		if not -d $keydir;
	my $certdir = "$apachedir/ssl.key";
	File::Path::mkpath($certdir) if not -d $certdir;
	print "Conf::StarCluster::generateCACert    Can't create certdir: $certdir\n"
		if not -d $certdir;

	my $copyprivate = "cp -f $privatekey $keydir/server.key";
	print "Conf::StarCluster::generateCACert    copyprivate: $copyprivate\n";
	`$copyprivate`;

	my $copypublic = "cp -f $CA_certfile $certdir/server.crt";
	print "Conf::StarCluster::generateCACert    copypublic: $copypublic\n";
	`$copypublic`;
}

=head2

	SUBROUTINE		setEnvironment

	PURPOSE

		1. SET EC2 ENVIRONMENT VARIABLES

=cut

sub setEnvironment {
	my $self	=	shift;

	my $userdata		=	$self->get_userdata();
	my $dbobject		=	$self->get_dbobject();
	my $conf 			= 	$self->get_conf();

	#### GET KEYFILES
    my $publicfile = $self->get_publicfile();
	my $privatefile = $self->get_privatefile();

	$ENV{'EC2_PRIVATE_KEY'} = $privatefile;
	$ENV{'EC2_CERT'} = $publicfile;
	$ENV{'EC2_HOME'} = $conf->getKeyValue("aws", "EC2HOME");
	$ENV{'EC2_AMITOOL_HOME'} = $conf->getKeyValue("aws", "EC2HOME");
	$ENV{'EC2_APITOOL_HOME'} = $conf->getKeyValue("aws", "EC2HOME");
	$ENV{'JAVA_HOME'} = $conf->getKeyValue("aws", "JAVAHOME");
	my $path = $ENV{'PATH'};
	$path = $conf->getKeyValue("aws", "EC2HOME") . "/bin:$path";
	$ENV{'PATH'} = $path;



}


=head2

	SUBROUTINE		mountNethome

	PURPOSE

		1. LOAD USER DATA INTO AGUA DATABASE

=cut

sub mountNethome {
	my $self		=	shift;
	my $username	=	shift;

	my $availzone = $self->get_conf()->getKeyValue(("aws", 'AVAILZONE'));
	my $mountpoint = $self->get_conf()->getKeyValue("agua", 'USERDIR');
	my $size = $self->get_conf()->getKeyValue(("aws", 'USERVOLUMESIZE'));
	my $device = $self->get_conf()->getKeyValue(("aws", 'USERDEVICE'));
	print "Agua::Init::mountNethome    availzone: $availzone\n";
	print "Agua::Init::mountNethome    mountpoint: $mountpoint\n";
	print "Agua::Init::mountNethome    size: $size\n";
	print "Agua::Init::mountNethome    device: $device\n";

	File::Path::mkpath($mountpoint) if not -d $mountpoint;
	print "Agua::Init::mountNethome    Can't create mountpoint\n"
		and exit if not -d $mountpoint;

	#### GET INPUT VOLUME SIZE IF DEFINED
	my $userdata = $self->get_userdata();
	if ( defined $self->get_volumesize() and $self->get_volumesize() )
	{
		$size = $self->get_volumesize();
	}

	#### SANITY CHECK
	my $mysqldir = "$mountpoint/mysql";
	print "Agua::Init:mountNethome    Skipping mountNethome. mysqldir present: $mysqldir\n" if -d $mysqldir;
	return if -d $mysqldir;

	#### GET KEYFILES
    my $publicfile = $self->get_publicfile();
	my $privatefile = $self->get_privatefile();

	#### GET INSTANCE ID
	my $instanceid = $self->get_instanceid();
	print "Agua::Init::mountNethome    instanceid: $instanceid\n";

	#### MAKE MOUNT POINT
	print "Agua::Init::mountNethome    Can't find mountpoint: $mountpoint\n" and exit if not -d $mountpoint;

	#### CREATE VOLUME
    my $create_command = "ec2-create-volume -v --debug -s $size -z $availzone | grep VOLUME | cut -f2";
	$create_command = $self->untaint($create_command);
	print "Agua::Init::mountNethome    create_command: $create_command\n";
	my $volumeid = `$create_command`;

	print "Agua::Init::mountNethome    volumeid: $volumeid\n";
	return if not defined $volumeid or not $volumeid;
	$volumeid =~ s/\s*\n$//;

	#### DETACH ANY EXISTING VOLUMES ATTACHED TO THIS DEVICE
	$self->detachAttached($instanceid, $device);

	#### ATTACH THE NEW VOLUME TO THE DEVICE
	my $attach_status = $self->attachVolume($instanceid, $volumeid, $device, $mountpoint);
	print "Agua::Init::mountNethome    attach_status: $attach_status\n";
	#if ( $attach_status eq "attaching" )
	#{
	#	$self->detachVolume($instanceid, $device);
	#
	#	$attach_status = $self->attachVolume($instanceid, $volumeid, $device, $mountpoint);
	#}

	#### FORMAT VOLUME
	my $format_command = "mkfs.ext3 -F $device";
	$format_command = $self->untaint($format_command);	
	print "Agua::Init::formatNethome    format_command: $format_command\n";
	my $format_success = `$format_command`;
	$format_success = $self->untaint($format_success);	
	print "Agua::Init::formatNethome    format_success: $format_success\n" if defined $format_success;

	#### MOUNT DEVICE TO MOUNTPOINT
	my $mount_command = "mount -t ext3 $device $mountpoint";
	$mount_command = $self->untaint($mount_command);
	print "Agua::Init::formatNethome    mount_command: $mount_command\n";
	my $mount_success = `$mount_command`;
	$mount_success = $self->untaint($mount_success);	
	print "Agua::Init::mountNethome    mount_success: $mount_success\n" if defined $mount_success;

	#### ADD TO FSTAB SO VOLUME DEVICE IS MOUNTED AUTOMATICALLY AFTER REBOOT
	$self->addToFstab("$device   $mountpoint      ext3    defaults,nobootwait	0	0\n");
}


=head2

	SUBROUTINE		mountMysql

	PURPOSE

		COPY MYSQL FROM /data DIR, RELOCATE POINTERS TO MYSQL

		TO NEW MYSQL LOCATION AND RESTART MYSQL

	INPUTS

		USES AWS INFORMATION IN AGUA CONF FILE, E.G.:

		PROJECTS_SNAPSHOT       snap-55fe4a3f
		USERDIR     /data
		USERVOLUMESIZE           40
		USERDEVICE         /dev/sdh

=cut

sub mountMysql {
	my $self	=	shift;

	#### COPY MYSQL FROM /data/mysql TO /nethome/mysql
	my $data_mountpoint = $self->get_conf()->getKeyValue("agua", 'DATADIR');
	my $mountpoint = $self->get_conf()->getKeyValue("agua", 'USERDIR');
	my $device = $self->get_conf()->getKeyValue(("aws", 'USERDEVICE'));
	print "Agua::Init::mountMysql    mountpoint: $mountpoint\n";
	print "Agua::Init::mountMysql    device: $device\n";

    #### STOP MYSQL
	my $mysql = "/etc/init.d/mysqld";
	$mysql = "/etc/init.d/mysql" if not -f $mysql;

	#### START MYSQL ONLY IF FOLDER ALREADY PRESENT
	my $mysqldir = "$mountpoint/mysql";
	print "Agua::Init::mountMysql    Looking for mysqldir: $mysqldir\n";
	if ( -d $mysqldir )
	{
		#### START MYSQL
		my $start_command = "$mysql restart";
		print "Agua::Init::mountMysql    start_command: $start_command\n";
		print `$start_command`;
		return;
	}

	#### OTHERWISE, COPY MYSQL FOLDER FROM /data
	print "Common::copyWorkflow    copying $data_mountpoint/mysql to $mysqldir\n";
	my $copy_success = File::Copy::Recursive::rcopy("$data_mountpoint/mysql", "$mountpoint/mysql");
	print "Common::copyWorkflow    copy result: $copy_success\n";
	if ( not defined $copy_success or not $copy_success )
	{
		print "Agua::Init::mountMysql    copy of mysql folders failed. Returning\n";
		return;
	}

    my $stop_command = "$mysql stop";
	print "Agua::Init::mountMysql    stop_command: $stop_command\n";
	print `$stop_command`;

	my $etcdir = "$mountpoint/mysql/etc/mysql";
	my $logdir = "$mountpoint/mysql/log/mysql";
	my $libdir = "$mountpoint/mysql/lib/mysql";

	#### ADD TO /etc/fstab
	$self->addToFstab("$device   $mountpoint      ext3    defaults        0 0\n");
    $self->addToFstab("$etcdir /etc/mysql     none bind\n");
    $self->addToFstab("$libdir /var/lib/mysql none bind\n");
    $self->addToFstab("$logdir /var/log/mysql none bind\n");

    #### CREATE BINDINGS TO LINK TO MYSQL ON EBS VOLUME
    ####   'apps' table is already populated
    ####   'admin' user password is set to default 'open4admin'
    ####   	Mysql root password is set to default 'open4root'

    # USE bind AND mount TO MIRROR NEW LOCATIONS IN OLD FOLDERS
    # Point MySQL to the correct database files on the EBS volume using mount bind.

    # USING FSTAB ENTRIES, MOUNT MYSQL DIRECTORIES FROM AGUA MYSQL DIRECTORY
	# MOUNT /etc/mysql
	my $etc_mysqldir = "/etc/mysql";
    my $command = "mount $etc_mysqldir";
	$command = "mkdir $etc_mysqldir;\n" . $command if not -d $etc_mysqldir; 
	print "Agua::Init::mountMysql    $command\n";
	print `$command`;

    # MOUNT /var/lib/mysql
	my $varlib_mysqldir = "/var/lib/mysql";
    $command = "mount $varlib_mysqldir";
	$command = "mkdir $varlib_mysqldir;\n" . $command if not -d $varlib_mysqldir;
	print "Agua::Init::mountMysql    $command\n";
	print `$command`;

    # MOUNT /var/log/mysql
	my $varlog_mysqldir = "/var/log/mysql";
    $command = "mount $varlog_mysqldir";
	$command = "mkdir $varlog_mysqldir;\n" . $command if not -d $varlog_mysqldir;
	print "Agua::Init::mountMysql    $command\n";
	print `$command`;

	# SET PERMISSIONS SO mysql USER HAS ACCESS TO FOLDERS
	`chown -R mysql:mysql $mountpoint/mysql/etc/mysql`;
    `chown -R mysql:mysql $mountpoint/mysql/lib/mysql`;
    `chown -R mysql:mysql $mountpoint/mysql/log/mysql`;

    #### START MYSQL
    my $restart_command = "$mysql restart";
	print "Agua::Init::mountMysql    restart_command: $restart_command\n";
	print `$restart_command`;
}


=head2

	SUBROUTINE		addToFstab

	PURPOSE

		START AN INSTANCE FOR THE GIVEN USER

=cut

sub addToFstab {
	my $self		=	shift;
	my $line 		= 	shift;


	my ($device, $mountpoint, $filesystem) = $line =~ /^\s*(\S+)\s+(\S+)\s+(\S+)/;

	#### CHECK INPUTS
	print "Agua::Init::addToFstab    device not defined\n" and exit if not defined $device;
	print "Agua::Init::addToFstab    mountpoint not defined\n" and exit if not defined $mountpoint;
	print "Agua::Init::addToFstab    filesystem not defined\n" and exit if not defined $filesystem;

	#### REMOVE NEW LINE
	$line =~ s/\s+$//;

	#### LOAD OLD FSTAB INFO
	my $fstabfile = "/etc/fstab";
	open(FSTABBKP, $fstabfile) or die "Can't open backupfile: $fstabfile";
	my $fstab;
	@$fstab = <FSTABBKP>;
	close(FSTABBKP) or die "Can't close backupfile: $fstabfile\n";

	#### REMOVE EMPTY LINES AND NEW LINE
	for ( my $i = 0; $i < @$fstab; $i++ )
	{
		if ( $$fstab[$i] =~ /^\s*$/ )
		{
			splice @$fstab, $i, 1;
			$i--;
		}
		$$fstab[$i] =~ s/\s+$//;
	}

	#### BACKUP FSTAB FILE
	my $counter = 1;
	my $backupfile = "$fstabfile.$counter";
	while ( -f $backupfile )
	{
		$counter++;
		$backupfile = "$fstabfile.$counter";
	}
	print "Agua::Init::addToFstab    backupfile: $backupfile\n";
	rename $fstabfile, $backupfile;

	#### REMOVE EXISTING FSTAB LINE FOR THIS MOUNT POINT
	for ( my $i = 0; $i < @$fstab; $i++ )
	{
		my ($current_device, $current_mountpoint, $current_filesystem) = $$fstab[$i] =~ /^\s*(\S+)\s+(\S+)\s+(\S+)/;
		next if not defined $current_device or not defined $current_mountpoint; 
		next if $current_device =~ /^#/;

		if ( $current_device eq $device
			and $current_mountpoint eq $mountpoint
			and $current_filesystem eq $filesystem )
		{
			print "Agua::Init::addToFstab    splicing line $i: $$fstab[$i]\n";
			splice @$fstab, $i, 1;
			$i--;
		}
	}

	#### ADD TO FSTAB SO VOLUME DEVICE IS MOUNTED AUTOMATICALLY AFTER REBOOT
	push @$fstab, $line;

	#### ADD NEW LINE
	for ( my $i = 0; $i < @$fstab; $i++ )
	{
		$$fstab[$i] .= "\n";
	}

    #### PRINT TO FSTAB FILE
	open(FSTAB, ">$fstabfile") or die "AWS::addToFstab    Can't open fstabfile: $fstabfile\n";
	print FSTAB @$fstab;
	close(FSTAB) or die "AWS::addToFstab    Can't close fstabfile: $fstabfile\n";
}


=head2

	SUBROUTINE		restartHttpd

	PURPOSE

		START THE HTTPD SERVER

=cut

sub restartHttpd {
	my $self	=	shift;

    #### START MYSQL
    my $restart_command = "/etc/init.d/httpd restart";
	print "Agua::Init::restartHttpd    restart_command: $restart_command\n";
	print `$restart_command`;
}

=head2

	SUBROUTINE		mountData

	PURPOSE

		START AN INSTANCE FOR THE GIVEN USER

	INPUTS

		USES AWS INFORMATION IN AGUA CONF FILE, E.G.:

		#### AWS
		DATASNAPSHOT       snap-55fe4a3f
		DATADIR     /data
		DATAVOLUMESIZE           40
		DATADEVICE         /dev/sdh

=cut

sub mountData {
	my $self		=	shift;

	my $username	=	$self->get_conf()->getKeyValue("agua", 'ADMINUSER');	
	my $snapshot	=	$self->get_conf()->getKeyValue(("aws", 'DATASNAPSHOT'));
	my $availzone 	= 	$self->get_conf()->getKeyValue(("aws", 'AVAILZONE'));
	my $mountpoint	= 	$self->get_conf()->getKeyValue("agua", 'DATADIR');
	my $size 		= 	$self->get_conf()->getKeyValue(("aws", 'DATAVOLUMESIZE'));
	my $device 		= 	$self->get_conf()->getKeyValue("aws", 'DATADEVICE');

	#### USE CUSTOM DATA SNAPSHOT IF DEFINED
	if ( defined $self->get_datavolumecheckbox() and $self->get_datavolumecheckbox() )
	{
		print "Admin::AWS::mountData    using custom datavolume: \t" , $self->get_datavolume(), "\n";
		$snapshot = $self->get_datavolume();
	}
	print "Admin::AWS::mountData    snapshot: $snapshot\t";

	#### USE DATA VOLUME SIZE IF DEFINED
	if ( defined $self->get_datavolumesize() and $self->get_datavolumesize() )
	{
		print "Admin::AWS::mountData    using custom datavolumesize: \t" , $self->get_datavolumesize(), "\n";
		$size = $self->get_datavolumesize();
	}
	print "Admin::AWS::mountData    size: $size\t";

	print "Admin::AWS::mountData    printenv:\n";
	print `printenv`;

	#### SANITY CHECK
	my $appsdir = "$mountpoint/apps";
	print "Agua::Init:mountData    Skipping mountData. appsdir present: $appsdir\n" if -d $appsdir;
	return if -d $appsdir;

	#### GET KEYFILES
    my $publicfile = $self->get_publicfile();
	my $privatefile = $self->get_privatefile();


	#### CREATE VOLUME
    my $create_command = "ec2-create-volume --snapshot $snapshot -s $size -z $availzone -K $privatefile -C $publicfile  | grep VOLUME | cut -f2";
	print "Agua::Init::mountData    create_command: $create_command\n";
	`$create_command` =~ /^(.+)$/;
	my $volumeid = $1;
	print "Agua::Init::mountData    volumeid: $volumeid\n";

	return if not defined $volumeid or not $volumeid;
	$volumeid =~ s/\s*\n$//;

	#### MAKE MOUNT POINT
	File::Path::mkpath($mountpoint) if not -d $mountpoint;
	print "Agua::Init::mountData    Can't find mountpoint: $mountpoint\n" and exit if not -d $mountpoint;

	#### GET INSTANCE ID
	my $instanceid = $self->get_instanceid();
	print "Agua::Init::mountData    instanceid: $instanceid\n";

	#### DETACH ANY EXISTING VOLUME ATTACHED TO THIS DEVICE
	$self->detachAttached($instanceid, $device);

	#### ATTACH THE NEW VOLUME TO THE DEVICE
	$self->attachVolume($instanceid, $volumeid, $device, $mountpoint);

	#### ADD LINE TO FSTAB
	my $line = "$device	$mountpoint	ext3	defaults,nobootwait	0	0\n";
	$self->addToFstab($device, $mountpoint, $line);

	#### MOUNT DEVICE TO MOUNTPOINT
	#### PAUSE TO ENSURE ATTACH IS SECURE
	while ( not -e $device )
	{
		print "Agua::Init::mountData    device not ready: $device\n";
	    sleep(3); 
	}
	my $mount_command = "mount -t ext3 $device $mountpoint";
	$mount_command = $self->untaint($mount_command);
	print "Agua::Init::mountData    mount_command: $mount_command\n";
	my $mount_success = `$mount_command`;
	print "Agua::Init::mountData    BEFORE untaint mount_success: $mount_success\n"  if defined $mount_success;
	$mount_success = $self->untaint($mount_success);
	print "Agua::Init::mountData    AFTER untaint mount_success: $mount_success\n" if defined $mount_success;

	return $mount_success;
}



=head2

	SUBROUTINE		createVolume

	PURPOSE

		CREATE AN EBS VOLUME

=cut

sub createVolume {
	my $self		=	shift;

	my $snapshot	=	$self->get_conf()->getKeyValue(("aws", 'DATASNAPSHOT'));
	my $availzone 	= 	$self->get_conf()->getKeyValue(("aws", 'AVAILZONE'));
	my $mountpoint	= 	$self->get_conf()->getKeyValue("agua", 'DATADIR');
	my $size 		= 	$self->get_conf()->getKeyValue(("aws", 'DATAVOLUMESIZE'));

	#### GET KEYFILES
    my $publicfile = $self->get_publicfile();
	my $privatefile = $self->get_privatefile();

	#### CREATE VOLUME
    my $create_command = "ec2-create-volume --snapshot $snapshot -s $size -z $availzone -K $privatefile -C $publicfile  | grep VOLUME | cut -f2";
	print "Agua::Init::createVolume    create_command: $create_command\n";
	`$create_command` =~ /^(.+)$/;
	my $volumeid = $1;
	print "Agua::Init::createVolume    volumeid: $volumeid\n";

	return $volumeid;
}


=head2

	SUBROUTINE		attachVolume

	PURPOSE

		ATTACH A VOLUME TO THIS DEVICE

=cut

sub attachVolume {
	my $self		=	shift;
	my $instanceid	=	shift;
	my $volumeid	=	shift;
	my $device		=	shift;
	my $mountpoint	=	shift;

	#### ec2-attach-volume "$VOLUME" -K "$EC2_PRIVATE_KEY" -C "$EC2_CERT" -i "$INSTANCE_ID" -d "$DEVICE"
	$instanceid =~ /^(.+)$/;
	$instanceid = $1;
	$device =~ /^(.+)$/;
	$device = $1;
	my $attach_command = "ec2-attach-volume $volumeid -i $instanceid -d $device";
	$attach_command = $self->untaint($attach_command);
	print "Agua::Init::attachVolume    attach_command: $attach_command\n";
	`$attach_command | cut -f 5` =~ /^(.+)$/;
	my $attach_success = $1;
	$attach_success =~ s/\s+$//;
	print "Agua::Init::attachVolume    attach_success: $attach_success\n";
	my $tries = 8;
	my $sleep = 3;
	my $counter = 0;
	while ( $counter < $tries and $attach_success ne "attached" )
	{
		#### FORMAT:
		#### ec2dvol 		
		####		
		#### VOLUME  vol-85f401ed    40      snap-55fe4a3f   us-east-1a      in-use  2010-11-18T19:40:43+0000
		#### ATTACHMENT      vol-85f401ed    i-b6147adb      /dev/sdh        attached        2010-11-18T19:40:49+0000

		sleep($sleep);		
		my $success_command = "ec2-describe-volumes $volumeid | grep ATTACHMENT | cut -f 5";
		$success_command = $self->untaint($success_command);
		print "$success_command\n";
		`$success_command` =~ /^(.+)$/;
		$attach_success = $1;
		$attach_success =~ s/\s+$//;
		print "Agua::Init::attachVolume    counter $counter attach_success: $attach_success\n";
		$counter = $tries if $attach_success eq "attached";
		$counter++;
	}

	return $attach_success;
}

=head2

	SUBROUTINE		detachAttached

	PURPOSE

		DETACH AN EXISTING VOLUME ATTACHED TO A DEVICE (/dev/*** ) ON THE INSTANCE 

		OUTPUT FORMAT:

			VOLUME vol-7318ed1b    10      snap-f1453d9b   us-east-1a      in-use  2010-11-19T02:14:37+0000
			ATTACHMENT      vol-7318ed1b    i-10acc57d      /dev/sda1       attached        2010-11-19T04:02:54+0000
			VOLUME  vol-751eeb1d    40      snap-55fe4a3f   us-east-1a      available       2010-11-19T02:05:52+0000
			...
			VOLUME  vol-ed05f085    40      snap-55fe4a3f   us-east-1a      in-use  2010-11-19T02:27:52+0000

			ATTACHMENT      vol-ed05f085    i-10acc57d      /dev/sdh        attached        2010-11-19T04:02:54+0000

=cut

sub detachAttached {
	my $self		=	shift;
	my $instanceid	=	shift;
	my $device		=	shift;

	my $attachedvolumeid;
	my @volumes = `ec2-describe-volumes`;
	foreach my $volume ( @volumes )
	{
		next if not ($volume =~ /$instanceid/ and $volume =~ /$device/);
		($attachedvolumeid) = $volume =~ /ATTACHMENT\s+(\S+)/;
		last;
	}

	if ( defined $attachedvolumeid and $attachedvolumeid )
	{
		print "Agua::Init::detachAttached    Doing self->detachVolume($instanceid, $device, $attachedvolumeid)\n";
		$self->detachVolume($instanceid, $device, $attachedvolumeid);
	}
}



=head2

	SUBROUTINE		detachVolume

	PURPOSE

		DETACH AN EXISTING VOLUME

=cut

sub detachVolume {
	my $self		=	shift;
	my $instanceid	=	shift;
	my $device		=	shift;
	my $volumeid	=	shift;


	#### FIRST, UNMOUNT DEVICE
	my $unmount = "umount $device";
	my $unmount_success = `$unmount`;

	my $detach_command = "ec2-detach-volume $volumeid";
	my $detach_output = `$detach_command`;
	my $detach_success = `ec2-describe-volumes $volumeid | cut -f 6`;
	return 1 if $detach_success eq "available";

	#### KEEP WAITING FOR DETACHED VOLUME TO BECOME 'available'
	my $tries = 8;
	my $sleep = 3;
	$detach_success = $self->wait_detach($volumeid, $tries, $sleep);
	return 1 if $detach_success eq "available";

	#### DID NOT DETACH SO GO NUCLEAR WITH --force
	print "Agua::Init::detachVolume    Simple detach failed.\n";	
	return 1 if $detach_success eq "available";
	$detach_command = "ec2-detach-volume --force $volumeid ";
	print "Agua::Init::detachVolume    Rerunning detach, this time with --force:\n\n$detach_command\n\n";
	$detach_output = `$detach_command`;
	$detach_success = `ec2-describe-volumes $volumeid | cut -f 6`;
	return 1 if $detach_success eq "available";
	$detach_success = $self->wait_detach($volumeid, $tries, $sleep);

	#### CHECK FILESYSTEM FOR ERRORS AFTER USING --force 
	my $check = "fsck -fy $device";

	return 1 if $detach_success eq "available";
	return 0;
}

sub wait_detach {
	my $volumeid	=	shift;
	my $tries		=	shift;
	my $sleep		=	shift;

	#### WAIT UNTIL THE DETACHING VOLUME IS 'available' TO MAKE SURE IT'S
	#### DETACHING HAS COMPLETED
	my $counter = 0;
	my $detach_success = '';
	while ( $counter < $tries )
	{
		#### FORMAT:
		#### ec2dvol 		
		####		
		#### VOLUME  vol-85f401ed    40      snap-55fe4a3f   us-east-1a      in-use  2010-11-18T19:40:43+0000
		#### ATTACHMENT      vol-85f401ed    i-b6147adb      /dev/sdh        attached        2010-11-18T19:40:49+0000

		sleep($sleep);
		$detach_success = `ec2-describe-volumes $volumeid | cut -f 6`;
		$detach_success =~ s/\s+//g;
		$counter = $tries if $detach_success eq "available";
		$counter++;
	}

	return $detach_success;
}


=head2

	SUBROUTINE		get_instanceid

	PURPOSE

		RETURN _instanceid IF DEFINED, OTHERWISE DISCOVER

		AWS INSTANCE ID FROM EXTERNAL SITE

=cut

sub get_instanceid {
	my $self	=	shift;

	return $self->{_instanceid} if defined $self->{_instanceid};

	my $instanceid =`curl -s http://169.254.169.254/latest/meta-data/instance-id`;
	$instanceid = $self->untaint($instanceid);
	$self->{_instanceid} = $instanceid if defined $instanceid;

	return $instanceid;
}

=head2

	SUBROUTINE		get_domainname

	PURPOSE

		RETURN _domainname IF DEFINED, OTHERWISE DISCOVER

		DOMAIN NAME FROM EXTERNAL SITE

=cut

sub get_domainname {
	my $self	=	shift;

	return $self->{_domainname} if defined $self->{_domainname};

	my $instanceid = $self->get_instanceid();
	my $domainname = `ec2-describe-instances | grep $instanceid | cut -f 4`;
	$domainname = $self->untaint($domainname);
 	print "Agua::Init::get_domainname    domainname: $domainname\n";
	$self->{_domainname} = $domainname if defined $domainname;

	return $domainname;
}

=head2

	SUBROUTINE		printKeyfiles

	PURPOSE

		PRINT THE PRIVATE KEY AND PUBLIC CERTIFICATE TO FILES

=cut

sub printKeyfiles {
	my $self		=	shift;


	#### GET KEYS
	my $publiccert	=	$self->get_ec2publiccert();
	my $privatekey	=	$self->get_ec2privatekey();

	#### PUBLIC CERT HAS 65-LETTER LINES
	my ($publicstring) = $publiccert =~ /^\s*\-{1,5}BEGIN CERTIFICATE\-{1,5}\s*(.+?)\s*\-{1,5}END CERTIFICATE\-{1,5}$/ms;
	$publicstring =~ s/\s+//g;
	my $public = $self->formatLines($publicstring, 64);
	$public = "-----BEGIN CERTIFICATE-----\n"
				. $public
				. "-----END CERTIFICATE-----";

	#### PUBLIC CERT HAS 65-LETTER LINES
	my ($privatestring) = $privatekey =~ /^\s*\-{1,5}BEGIN PRIVATE KEY\-{1,5}\s*(.+?)\s*\-{1,5}END PRIVATE KEY\-{1,5}$/ms;
	$privatestring =~ s/\s+//g;
	my $private = $self->formatLines($privatestring, 76);
	$private = "-----BEGIN PRIVATE KEY-----\n"
				. $private
				. "-----END PRIVATE KEY-----";

	#### GET KEYFILES
    my $publicfile = $self->get_publicfile();
	my $privatefile = $self->get_privatefile();

	#### CREATE KEYDIR
	print "Agua::Init:printKeyfiles       Can't create keydir\n" and exit if not $self->createKeydir();

	open(PUBLICFILE, ">$publicfile") or die "AWS:printKeyfiles    Can't open publicfile: $publicfile\n";
	print PUBLICFILE $public;
	close(PUBLICFILE) or die "AWS:printKeyfiles    Can't close publicfile: $publicfile\n";
	#`chmod 600 $publicfile`;

	open(PRIVATEFILE, ">$privatefile") or die "AWS:printKeyfiles    Can't open privatefile: $privatefile\n";
	print PRIVATEFILE $private;
	close(PRIVATEFILE) or die "AWS:printKeyfiles    Can't close privatefile: $privatefile\n";
	`chmod 600 $privatefile`;


}
=head2

	SUBROUTINE		get_publicfile

	PURPOSE

		SET LOCATION OF PUBLIC FILE

=cut

sub get_publicfile {
	my $self	=	shift;

	print "Agua::Init::get_publicfile    AWS::get_publicfile()\n";
	return $self->{_publicfile} if defined $self->{_publicfile};

	my $keydir = $self->get_keydir();
	my $publicfile = "$keydir/public.pem";
	$self->set_publicfile($publicfile);

	return $publicfile;
}

=head2

	SUBROUTINE		get_privatefile

	PURPOSE

		SET LOCATION OF PRIVATE FILE

=cut

sub get_privatefile {
	my $self	=	shift;

	return $self->{_privatefile} if defined $self->{_privatefile};

	my $keydir = $self->get_keydir();
	my $privatefile = "$keydir/private.pem";
	$self->set_privatefile($privatefile);

	return $privatefile;
}

sub get_keydir {
	my $self	=	shift;

	my $installdir = $self->get_conf()->getKeyValue("agua", 'INSTALLDIR');
	my $username = $self->get_username();
	print "Agua::Init::get_privatefile    username not defined. Returning\n" and return if not defined $username;

	return "$installdir/conf/$username/.keypairs";
}

sub formatLines {
	my $self	=	shift;
	my $string	=	shift;
	my $length	=	shift;


	my $lines = '';
	my $offset = 0;
	while ( $offset < length($string) )
	{
		$lines .= substr($string, $offset, $length);
		$lines .= "\n";
		$offset += $length;
	}


	return $lines;	
}


=head2

	SUBROUTINE		createKeydir

	PURPOSE

		SET LOCATION OF PUBLIC FILE

=cut

sub createKeydir {
	my $self	=	shift;

	print "Agua::Init::createKeydir    AWS::createKeydir()\n";

	return $self->{_publicfile} if defined $self->{_publicfile};

	my $keydir = $self->get_keydir();
	print "Agua::Init::createKeydir    keydir: $keydir\n";
	File::Path::mkpath($keydir) if not -d $keydir;
	print "Agua::Init::createKeydir    Can't create keydir: $keydir\n" and return 0 if not -d $keydir;

	return 1;
}


=head2

	SUBROUTINE		getKeys

	PURPOSE

		GET THE AWS KEYS FOR THE GIVEN USER

=cut

sub getKeys {
	my $self		=	shift;

	my $json		=	$self->get_json();
	my $dbobject	=	$self->get_dbobject();
	print "Agua::Init::getKeys    json: $json\n";
	print "Agua::Init::getKeys    dbobject: $dbobject\n";

	my $username	=	$json->{username};
	my $query = "SELECT * FROM aws WHERE username='$username'";
	my $keys = $dbobject->queryhash($query);
	print "keys:\n";
	print Dumper $keys;

	return $keys;
}

sub untaint {
	my $self	=	shift;
	my $string 	=	shift;

	$string =~ /^(.+)$/;
	return $1;
}


################################################################################
##################			HOUSEKEEPING SUBROUTINES			################
################################################################################


=head2

	SUBROUTINE		initialise

	PURPOSE

		INITIALISE THE self OBJECT:

			1. LOAD THE DATABASE, USER AND PASSWORD FROM THE ARGUMENTS

			2. FILL OUT %VARIABLES% IN XML AND LOAD XML

			3. LOAD THE ARGUMENTS

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

	if ( defined $arguments->{userdata} )
	{
		#### LOAD THE USER-PROVIDED ARGUMENTS
		foreach my $key ( keys %{$arguments->{userdata}} )
		{		
			#### LOAD THE KEY-VALUE PAIR
			$self->value($key, $arguments->{userdata}->{$key});
		}
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
    unless($operation && $attribute) {
        croak "Method name $AUTOLOAD is not in the recognized form (get|set)_attribute\n";
    }
    unless( exists $self->{$attribute} or $self->is_valid($attribute) )
	{
        #croak "No such attribute '$attribute' exists in the class ", ref($self);
		return;
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


# When an object is no longer being used, this will be automatically called
# and will adjust the count of existing objects
sub DESTROY {
    my($self) = @_;

	#if ( defined $self->{_databasehandle} )
	#{
	#	my $dbh =  $self->{_databasehandle};
	#	$dbh->disconnect();
	#}

#    my($self) = @_;
}



1;

