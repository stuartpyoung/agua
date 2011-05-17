use MooseX::Declare;

=head2

	PACKAGE		StarCluster

    PURPOSE

        1. START A CLUSTER, MOUNT AQUARIUS /data AND /nethome

			ON ITS NODES

        2. USERS CAN RUN small, medium, large OR CUSTOM CLUSTERS
            (ALL USERS USE admin USER'S CONFIG FILE)

        3. EACH WORKFLOW USES A SINGLE CLUSTER TO RUN ALL ITS STAGES

		4. MORE THAN ONE WORKFLOW CAN USE A CLUSTER AT THE SAME TIME

        5. THE CLUSTER SHUTS DOWN WHEN THERE ARE NO RUNNING WORKFLOWS


	NOTES

	    1. A NEW SGE CELL IS CREATED ON THE AGUA HEAD NODE FOR EACH

			CLUSTER THAT IS STARTED E.G., /opt/sge6/syoung-smallcluster

	    2. WITHIN EACH CELL, A NEW QUEUE IS CREATED FOR EACH WORKFLOW

	        E.G., syoung-project1-workflow1

		3. ADD THE threaded PARALLEL ENVIRONMENT TO EACH QUEUE SO THAT

			JOBS CAN USE MULTIPLE CPUS (LATER: CHECK IF NECESSARY)


    TO DO

        1. ALLOW USERS TO SPECIFY MAX NUMBER OF NODES, AMIs, ETC.

            BY EDITING USER-SPECIFIC CONFIG FILES CONTAINING

			[cluster myClusterName] SECTIONS	

            USERNAME    CONFIGFILE
            admin       /home/admin/.starcluster/config
            jgilbert    /home/jgilbert/.starcluster/config

        2. ALLOW USER TO SPECIFY A CLUSTER FOR EACH WORKFLOW

            (ALL STAGES IN THE SAME WORKFLOW USE THE SAME CLUSTER)

        3. IMPLEMENT AUTO-LEVELLING TO AUTOMATICALLY 

            INCREASE/DECREASE THE NUMBER OF NODES BASED ON

            THE NUMBER OF QUEUED JOBS	

=cut

class Agua::StarCluster with (Agua::Common::Base, Agua::Common::Cluster, Agua::Common::SGE, Agua::Common::Util) {


use FindBin qw($Bin);
use lib "$Bin/..";

#### INTERNAL MODULES
use Agua::DBaseFactory;
use Conf::Agua;
use Conf::StarCluster;

#### EXTERNAL MODULES
use Data::Dumper;
use File::Path;
use Getopt::Simple;
#use Moose::Util::TypeConstraints;

#### AWS INFO

#### Boolean
has 'help'			=> ( is  => 'rw', 'isa' => 'Bool', required	=>	0, documentation => "Print help message"	);
#### Int
has 'sleep'			=> ( is  => 'rw', 'isa' => 'Int', default	=>	600	);
has 'workflowpid'	=> ( is  => 'rw', 'isa' => 'Int', required	=>	0	);
#has 'clustersize'	=> ( is  => 'rw', 'isa' => 'Int', required	=>	0, default => 2	);
has 'nodes'			=> ( is  => 'rw', 'isa' => 'Str', required	=>	0	);
has 'maxnodes'		=> ( is  => 'rw', 'isa' => 'Str', required	=>	0	);
has 'minnodes'		=> ( is  => 'rw', 'isa' => 'Str', required	=>	0	);
has 'sleepinterval'	=> ( is  => 'rw', 'isa' => 'Int', required	=>	0, default => 30	);
has 'slots'			=> ( is  => 'rw', 'isa' => 'Int', required	=>	0, default => 1	);
#has 'qmasterport'=> ( is  => 'rw', 'isa' => 'Str', required	=>	0, default => 36231	);
#has 'execdport'=> ( is  => 'rw', 'isa' => 'Str', required	=>	0, default => 36232	);
#### String
has 'sgeroot'		=> ( is  => 'rw', 'isa' => 'Str', default => "/opt/sge6"	);
has 'sgecell'		=> ( is  => 'rw', 'isa' => 'Str', required	=>	0	);
has 'fileroot'		=> ( is  => 'rw', 'isa' => 'Str', required	=>	0	);
has 'username'		=> ( is  => 'rw', 'isa' => 'Str', required	=>	0	);
has 'clustertype'	=> ( is  => 'rw', 'isa' => 'Str', required	=>	0	);
has 'keyname'		=> ( is  => 'rw', 'isa' => 'Str', required	=>	0	);
has 'cluster'		=> ( is  => 'rw', 'isa' => 'Str', required	=>	0	);
has 'clusteruser'	=> ( is  => 'rw', 'isa' => 'Str', default	=> "sgeadmin"	);
has 'availzone'		=> ( is  => 'rw', 'isa' => 'Str', required	=>	0	);
has 'instancetype'	=> ( is  => 'rw', 'isa' => 'Str', required	=>	0	);
has 'nodeimage'		=> ( is  => 'rw', 'isa' => 'Str', required	=>	0	);
has 'starcluster'	=> ( is  => 'rw', 'isa' => 'Str', required	=>	0	);
has 'plugins'		=> ( is  => 'rw', 'isa' => 'Str', default	=> "automount,sge"	);
has 'amazonuserid'	=> ( is  => 'rw', 'isa' => 'Str', required	=>	0	);
has 'accesskeyid'	=> ( is  => 'rw', 'isa' => 'Str', required	=>	0	);
has 'secretaccesskey'=> ( is  => 'rw', 'isa' => 'Str', required	=>	0	);
has 'privatekey'	=> ( is  => 'rw', 'isa' => 'Str', required	=>	0	);
has 'publiccert'	=> ( is  => 'rw', 'isa' => 'Str', required	=>	0	);
has 'keypairfile'	=> ( is  => 'rw', 'isa' => 'Str|Undef', required	=>	0	);
has 'outputdir'		=> ( is  => 'rw', 'isa' => 'Str', default 	=>	''	);
has 'sources'		=> ( is  => 'rw', 'isa' => 'Str', default 	=>	''	);
has 'mounts'		=> ( is  => 'rw', 'isa' => 'Str', default 	=>	''	);
has 'devs'			=> ( is  => 'rw', 'isa' => 'Str', default 	=>	''	);
has 'configfile'	=> ( is  => 'rw', 'isa' => 'Str'	);
has 'outputfile'	=> ( is  => 'rw', 'isa' => 'Str'	);
#### Object
has 'dbobject'		=> ( isa => 'Agua::DBase::MySQL', is => 'rw', required => 0 );
has 'monitor'		=> ( is  => 'rw', 'isa' => 'Maybe' );
has 'sourcedirs'	=> ( is  => 'rw', 'isa' => 'ArrayRef[Str]' );
has 'mountpoints'	=> ( is  => 'rw', 'isa' => 'ArrayRef[Str]' );
has 'devices'		=> ( is  => 'rw', 'isa' => 'ArrayRef[Str]' );
has 'fields'    	=> ( isa => 'ArrayRef[Str|Undef]', is => 'ro', default => sub { ['help', 'amazonuserid', 'accesskeyid', 'secretaccesskey', 'privatekey', 'publiccert', 'username', 'keyname', 'cluster', 'nodes', 'outputdir', 'maxnodes', 'minnodes', 'sources', 'mounts', 'devs', 'instancetype', 'availzone', 'nodeimage', 'configfile'] } );
has 'instancetypeslots' =>  ( isa => 'HashRef', is => 'ro', default => sub {
	{
		"t1.micro"	=>	1,
		"m1.small"	=>  1,
		"m1.large"	=>	4,
		"m1.xlarge"	=>	8,
		"m2.xlarge"	=>	7,
		"m2.2xlarge"=>	13,
		"m2.4xlarge"=>	26,
		"c1.medium"	=>	5,
		"c1.xlarge"	=>	20,
		"cc1.4xlarge"=>	34,
		"cg1.4xlarge"=>	34
	}
});
#### http://aws.amazon.com/ec2/instance-types

has 'conf' 	=> (
	is =>	'rw',
	'isa' => 'Conf::Agua',
	default	=>	sub { Conf::Agua->new(	backup	=>	1, separator => "\t"	);	}
);
has 'config'=> (
	is =>	'rw',
	'isa' => 'Conf::StarCluster',
	default	=>	sub { Conf::StarCluster->new(	backup	=>	1, separator => "="	);	}
);

####/////}}

=head2

	SUBROUTINE		BUILD

	PURPOSE

		GET AND VALIDATE INPUTS, AND INITIALISE OBJECT

=cut

method BUILD ($hash) {

	$self->getopts();
	$self->initialise();
}

method initialise () {

	my $conf 		= 	$self->conf();
	my $username 	= 	$self->username();
	my $cluster 	= 	$self->cluster();

	#### ADD CONFIG FILE TO config OBJECT
	my $configfile = $self->getConfigfile($username, $cluster);
	$self->config()->inputfile($configfile);

	my $userdir 	= 	$conf->getKeyValue('agua', 'USERDIR');
	my $aguadir 	= 	$conf->getKeyValue('agua', 'AGUADIR');
	print "Agua::StarCluster::initialise    username not defined\n" and exit if not defined $username;
	print "Agua::StarCluster::initialise    userdir not defined\n" and exit if not defined $userdir;

	#### SET DATABASE HANDLE
	$self->setDbh();

	#### SET CLUSTER TYPE
	$self->clustertype($conf->getKeyValue('agua', "CLUSTERTYPE"));

	#### SET KEYNAME IF NOT DEFINED
	$self->keyname("$username-key") if not $self->keyname();

	#### SET KEYPAIR FILE IF NOT DEFINED
	$self->keypairfile("$userdir/$username/.starcluster/id_rsa-" . $self->keyname()) if not $self->keypairfile();

	#### SET OUTPUT DIR IF NOT DEFINED
	$self->outputdir("$userdir/$aguadir/$username/.cluster") if not $self->outputdir();

	#### SET SOURCEDIRS AND MOUNTPOINTS
	my @sourcedirs = 	split ",", 	$self->sources() 	|| () if $self->sources();
	my @mountpoints = 	split ",", 	$self->mounts() 	|| () if $self->mounts();
	my @devices = 		split ",", 	$self->devs() 		|| () if $self->devs();
	$self->sourcedirs(\@sourcedirs) if $self->sources();
	$self->mountpoints(\@mountpoints) if $self->mounts();
	$self->devices(\@devices) if $self->devs();
}

method getopts () {
	my @temp = @ARGV;
	my $args = $self->args();

	my $olderr;
	open $olderr, ">&STDERR";	
	open(STDERR, ">/dev/null") or die "Can't redirect STDERR to /dev/null\n";
	my $options = Getopt::Simple->new();
	$options->getOptions($args, "Usage: blah blah blah"); 
	open STDERR, ">&", $olderr;
	my $switch = $options->{switch};
	foreach my $key ( keys %$switch )
	{
		$self->$key($switch->{$key}) if defined $switch->{$key};
	}

	@ARGV = @temp;

}

method args() {
	my $meta = $self->meta();

	my %option_type_map = (
		'Bool'     => '!',
		'Str'      => '=s',
		'Int'      => '=i',
		'Num'      => '=f',
		'ArrayRef' => '=s@',
		'HashRef'  => '=s%',
		'Maybe'    => ''
	);

	my $attributes = $self->fields();
	my $args = {};
	foreach my $attribute_name ( @$attributes )
	{
		my $attr = $meta->get_attribute($attribute_name);
		my $attribute_type  = $attr->{isa};

		$attribute_type =~ s/\|.+$//;
		$args -> {$attribute_name} = {  type => $option_type_map{$attribute_type}  };
	}

	return $args;
}



=head2

	SUBROUTINE 		generateKeypair

	PURPOSE

		GENERATE A KEYPAIR USING PRIVATE AND PUBLIC KEYS

=cut

method generateKeypair () {
	#### SET DEFAULT KEYNAME
	$self->keyname() = "$self->username()-key" if not defined $self->keyname();

	print "Agua::StarCluster::generateKeypair    privatekey: ", $self->privatekey(), "\n";
	print "Agua::StarCluster::generateKeypair    publiccert: ", $self->publiccert(), "\n";
	print "Agua::StarCluster::generateKeypair    username: ", $self->username(), "\n";
	print "Agua::StarCluster::generateKeypair    keyname: ", $self->keyname(), "\n";

	#### PRINT HELP
	if ( defined $self->help() )
	{
		print qq{

		$0 start <--privatekey String> <--publiccert String> [--help]

		--privatekey 	Location of private key
		--privatecert 	Location of public key
		--username		Name of user 
		--keyname		Name of key
		--help			Print help info
};
	}

	#### SET TARGET KEYPAIR FILE
	my $conf = $self->conf();
	my $userdir 		= 	$conf->getKeyValue('agua', "USERDIR");
	my $keypairdir 		=	"$userdir/" . $self->username() . "/.starcluster";
	my $keypairfile		=	"$keypairdir/id_rsa-" . $self->keyname();

	`mkdir -p $keypairdir` if not -d $keypairdir;
	print "Agua::StarCluster::generatedKeypair    can't create keypairdir: $keypairdir\n" and exit if not -d $keypairdir;

	my $delete = qq{ec2-delete-keypair \\
-K } . $self->privatekey() . qq{ \\
-C } . $self->publiccert()  . qq{ \\
id_rsa-} . $self->keyname();
	print "$delete\n";
	print `$delete`;

	my $add = qq{ec2-add-keypair \\
-K } . $self->privatekey() . qq{ \\
-C } . $self->publiccert() . qq{ \\
id_rsa-} . $self->keyname() . qq{ \\
> $keypairfile
};
	print "$add\n";
	print `$add`;
}

=head2

	SUBROUTINE 		writeConfigfile

	PURPOSE

		ADD KEYPAIR FILE, AWS ACCESS IDs AND OTHER CREDENTIAL

		INFO TO USER'S CLUSTER CONFIG

=cut

method writeConfigfile () {


	#### PRINT HELP
	if ( defined $self->help() )
	{
		print qq{

	$0 writeConfigfile <--amazonuserid String> <--accesskeyid String> [--help]

	--privatekey 		Location of private key
	--privatecert 		Location of public key
	--amazonuserid 		AWS user ID
	--accesskeyid		AWS access key ID
	--secretaccesskey 	AWS secret access key
	--keyname			Name of keypair
	--username			Change this user's config file
	--help				Print help info
	--cluster 			Name of cluster (e.g., microcluster) 
	--nodes 			Number of nodes to begin with
	--maxnodes 			Max. number of nodes (using load balancer)
	--minnodes 			Min. number of nodes (using load balancer)
	--outputdir 		Print load balancer STDOUT/STDERR to file in this directory
	--sources 			Comma-separated list of source dirs (e.g., /data,/nethome)
	--mounts 			Comma-separated list of mount points (e.g., /data,/nethome)
	--devs 				Comma-separated list of devices (e.g., /dev/sdh,/dev/sdi)
	--configfile 		Print config to this file (e.g., /nethome/admin/.starcluster/config)
};
	}

	print "Agua::StarCluster::writeConfigfile    cluster not defined\n"
		and exit if not $self->cluster();
	print "Agua::StarCluster::writeConfigfile    username not defined\n"
		and exit if not $self->username();
	print "Agua::StarCluster::writeConfigfile    configfile not defined\n"
		and exit if not $self->configfile();

	#### SET DEFAULT KEYNAME
	$self->keyname() = $self->username(). "-key" if not defined $self->keyname();

	#### SET TARGET KEYPAIR FILE
	my $keypairfile		=	$self->keypairfile();

	#### SET INPUTFILE
	my $configfile	=	$self->configfile();
	my $conf 		= 	$self->conf();
	my $username 	= 	$self->username();
	my $cluster 	= 	$self->cluster();
	my $userdir 	= 	$conf->getKeyValue('agua', "USERDIR");

	my $config = Conf::StarCluster->new(
		separator	=>	"=",
		inputfile	=>	$configfile
	);

	print "Agua::StarCluster::writeConfigfile    privatekey: ", $self->privatekey(), "\n";
	print "Agua::StarCluster::writeConfigfile    publiccert: ", $self->publiccert(), "\n";
	print "Agua::StarCluster::writeConfigfile    amazonuserid: ", $self->amazonuserid(), "\n";
	print "Agua::StarCluster::writeConfigfile    accesskeyid: ", $self->accesskeyid(), "\n";
	print "Agua::StarCluster::writeConfigfile    secretaccesskey: ", $self->secretaccesskey(), "\n";
	print "Agua::StarCluster::writeConfigfile    keyname: ", $self->keyname(), "\n";
	print "Agua::StarCluster::writeConfigfile    username: ", $self->username(), "\n";
	print "Agua::StarCluster::writeConfigfile    sources: ", $self->sources(), "\n";
	print "Agua::StarCluster::writeConfigfile    mounts: ", $self->mounts(), "\n";
	print "Agua::StarCluster::writeConfigfile    devs: ", $self->devs(), "\n";
	print "Agua::StarCluster::writeConfigfile    instancetype: ", $self->instancetype(), "\n";
	print "Agua::StarCluster::writeConfigfile    nodes: ", $self->nodes(), "\n";
	print "Agua::StarCluster::writeConfigfile    availzone: ", $self->availzone(), "\n";
	print "Agua::StarCluster::writeConfigfile    nodeimage: ", $self->nodeimage(), "\n";
	print "Agua::StarCluster::writeConfigfile    cluster: ", $self->cluster(), "\n";
	print "Agua::StarCluster::writeConfigfile    configfile: ", $self->configfile(), "\n";

	#### SET [global]
	$config->addKeyValue("global", "", "DEFAULT_TEMPLATE", $cluster);
	$config->addKeyValue("global", "", "ENABLE_EXPERIMENTAL", "True");

	#### SET [cluster <clusterName>]
	$config->addKeyValue("cluster", $cluster, "KEYNAME", "id_rsa-" . $self->keyname());
	$config->addKeyValue("cluster", $cluster, "AVAILABILITY_ZONE", $self->availzone());
	$config->addKeyValue("cluster", $cluster, "CLUSTER_SIZE", $self->nodes());
	$config->addKeyValue("cluster", $cluster, "CLUSTER_USER", $self->clusteruser());
	$config->addKeyValue("cluster", $cluster, "NODE_IMAGE_ID", $self->nodeimage());
	$config->addKeyValue("cluster", $cluster, "NODE_INSTANCE_TYPE", $self->instancetype());
	$config->addKeyValue("cluster", $cluster, "CLUSTER_USER", $self->clusteruser());
	$config->addKeyValue("cluster", $cluster, "PLUGINS", $self->plugins());

	#### SET [aws info]
	$config->addKeyValue("aws", "info", "AWS_USER_ID", $self->amazonuserid());
	$config->addKeyValue("aws", "info", "AWS_ACCESS_KEY_ID", $self->accesskeyid());
	$config->addKeyValue("aws", "info", "AWS_SECRET_ACCESS_KEY", $self->secretaccesskey());	
	#### SET [key <keyname>]
	$config->addKeyValue("key", "id_rsa-" . $self->keyname(), "KEYNAME", "id_rsa-". $self->keyname());
	$config->addKeyValue("key", "id_rsa-" . $self->keyname(), "KEY_LOCATION", $keypairfile);

	#### SET [plugin automount]
	my ($internalip, $internallongip) = $self->getHeadInternalIps();
	my $portmapport = $self->conf()->getKeyValue("starcluster:nfs", "PORTMAPPORT");
	my $nfsport 	= $self->conf()->getKeyValue("starcluster:nfs", "NFSPORT");
	my $mountdport 	= $self->conf()->getKeyValue("starcluster:nfs", "MOUNTDPORT");
	$config->addKeyValue("plugin", "automount", "setup_class", "automount.NfsShares");
	$config->addKeyValue("plugin", "automount", "head_ip", $internalip);
	$config->addKeyValue("plugin", "automount", "cluster", $cluster);
	$config->addKeyValue("plugin", "automount", "portmapport", $portmapport);
	$config->addKeyValue("plugin", "automount", "nfsport", $nfsport);
	$config->addKeyValue("plugin", "automount", "mountdport", $mountdport);
	$config->addKeyValue("plugin", "automount", "interval", $self->sleepinterval())
		if $self->sleepinterval();
	$config->addKeyValue("plugin", "automount", "devices", $self->devs())
		if $self->devs();
	$config->addKeyValue("plugin", "automount", "mountpoints", $self->mounts())
		if $self->mounts();
	$config->addKeyValue("plugin", "automount", "sourcedirs", $self->sources())
		if $self->sources();
	print "Agua::StarCluster::writeConfigfile    devs: ", $self->devs(), "\n";
	print "Agua::StarCluster::writeConfigfile    mounts: ", $self->mounts(), "\n";

	#### SET sge PLUGIN
	my $slots = $self->setSlots($self->instancetype());
	my ($qmasterport, $execdport) = $self->getPorts();
	print "Agua::StarCluster::writeConfigfile    slots: $slots\n"; 
	print "Agua::StarCluster::writeConfigfile    qmasterport: $qmasterport\n"; 
	print "Agua::StarCluster::writeConfigfile    execdport: $execdport\n"; 
	$config->addKeyValue("plugin", "sge", "setup_class", "sge.CreateCell");
	$config->addKeyValue("plugin", "sge", "head_ip", $internalip);
	$config->addKeyValue("plugin", "sge", "head_long_ip", $internallongip);
	$config->addKeyValue("plugin", "sge", "root", $self->sgeroot());
	$config->addKeyValue("plugin", "sge", "cell", $cluster);
	$config->addKeyValue("plugin", "sge", "qmasterport", $qmasterport);
	$config->addKeyValue("plugin", "sge", "execdport", $execdport);
	$config->addKeyValue("plugin", "sge", "slots", $slots);

	print "Configfile printed: $configfile\n";
}

=head2

	SUBROUTINE		_createVolume

	PURPOSE

		CREATE AN EBS VOLUME

=cut

method _createVolume ($privatekey, $publiccert, $snapshot, $availzone, $size) {

	#### CREATE VOLUME
    my $create_command = "ec2-create-volume --snapshot $snapshot -s $size -z $availzone -K ". $self->privatekey() . " -C " . $self->publiccert() . " | grep VOLUME | cut -f2";
	print "Agua::StarCluster::_createVolume    create_command: $create_command\n";
	my $volumeid = `$create_command`;
	print "Agua::StarCluster::_createVolume    volumeid: $volumeid\n";

	return $volumeid;
}


=head2

	SUBROUTINE 		start

	PURPOSE

		START UP A STARCLUSTER cluster FOR THE GIVEN USER

=cut
method start () {
	print "Agua::StarCluster::start    StarCluster::start(privatekey, publiccert, username, cluster, keyname)\n";

	$self->getopts();

	print "Agua::StarCluster::start    privatekey not defined\n" and exit if not $self->privatekey();
	print "Agua::StarCluster::start    publiccert not defined\n" and exit if not $self->publiccert();
	print "Agua::StarCluster::start    keyname not defined\n" and exit if not $self->keyname();
	print "Agua::StarCluster::start    username not defined\n" and exit if not $self->username();
	print "Agua::StarCluster::start    cluster not defined\n" and exit if not $self->cluster();

	#### PRINT HELP
	if ( defined $self->help() )
	{
		print qq{

	$0 start <--username String> <--cluster String> [--help]

	--privatekey 	Location of private key file
	--publiccert	Location of public certificate file
	--keyname 		Name of the key (e.g., 'admin-key')
	--username 		User account name
	--cluster		Name of cluster to run (e.g., 'smallcluster')
	--outputdir		Location of public certificate file
	--nodes			Number of nodes to run in cluster
	--maxnodes		Maximum number of nodes to run (optional load balancing)
	--minnodes		Minimum number of nodes to run (optional load balancing)
	--help			Print help info
};
	}

#	#### LAUNCH THE CLUSTER
#	$self->launchCluster();

	#### WAIT FOR CLUSTER TO COME UP
	sleep($self->sleep());

	#### LAUNCH LOAD BALANCER
 	print "Agua::StarCluster::start    Doing self->launchBalancer()\n";
	$self->launchBalancer() if $self->maxnodes() or $self->minnodes();

	#### **** DEPRECATED: automount.py STARCLUSTER PLUGIN TAKES CARE OF THIS **
	#### **** THIS IS TO MAKE SURE THAT THE NODES ARE FULLY OPERATIONAL
	#### **** BEFORE THE MOUNTS ARE ATTEMPTED
	####
	#### OPEN NFS AND SGE PORTS IN SECURITY GROUP
	#### my $cluster 	= 	$self->cluster();
	#### print "Agua::StarCluster::start    cluster: $cluster\n";
	#### print "Agua::StarCluster::start    Doing self->openPorts()\n";
	#### $self->openPorts("\@sc-$cluster");
	####
	#### SETUP SHARES FROM HEAD
	#### $self->addNfsMounts();
	####
	#### MOUNT SHARES ON MASTER AND ALL NEW NODES
	#### $self->mountShares($self->privatekey(), $self->publiccert(), $self->username(), $self->keyname(), $self->cluster(), $nodes);
	####	

	#### **** DEPRECATED: sge.py STARCLUSTER PLUGIN TAKES CARE OF THIS **
	##### SET THE DEFAULT QUEUE ON MASTER
	#####	$self->setQueue("default");
	####	
	##### SET threaded PARALLEL ENVIRONMENT ON MASTER
	##### 	print "Agua::StarCluster::start    Doing self->setPE('threaded')\n";
	#####	$self->setPE("threaded", "default");	


}


method getPorts () {

	#### SELECT UNUSED PORTS FOR SGE_QMASTER_PORT AND SGE_EXECD_PORT
	#### NB: MAX 100 CLUSTERS, QMASTER PORT RANGE: 36241-37231
	my $minmaster = 36231;
	my $maxmaster = 37231;
	my $qmasterport;
	my $execdport;

	#### GET THE PORTS FOR THE CLUSTER IF WE ARE SIMPLY
	#### REINITIALISING ITS 
	my $username = $self->username();
	my $cluster = $self->cluster();
	my $root = $self->sgeroot();

	my $query = qq{SELECT qmasterport
FROM clustervars
WHERE username = '$username'
AND cluster = '$cluster'};

	$qmasterport = $self->dbobject()->query($query);
	$execdport = $qmasterport + 1 if defined $qmasterport;

	if ( not defined $qmasterport )
	{
		$query = qq{SELECT MAX(qmasterport)
FROM clustervars
WHERE qmasterport > $minmaster};

		$qmasterport = $self->dbobject()->query($query);
		$qmasterport += 10 if defined $qmasterport;
		$execdport = $qmasterport + 1 if defined $qmasterport;

		#### IF WE HAVE REACHED MAX PORT NUMBER, qmasterport IS UNDEFINED.
		#### IN THIS CASE, RECYCLE UNUSED PORTS IN ALLOWED RANGE
		if ( not defined $qmasterport )
		{
			my $query = qq{SELECT COUNT(qmasterport) FROM clustervars};
			my $count = $self->dbobject()->query($query);
			if ( $count == 0 )
			{
				$qmasterport = $minmaster + 10;
			}
			else
			{
				my $query = qq{SELECT qmasterport FROM clustervars};
				print "$query\n";
				my $ports = $self->dbobject()->queryarray($query);
				for ( my $i = 0; $i < @$ports - 1; $i++ )
				{
					if ( ($$ports[$i + 1] - $$ports[$i]) > 10 )
					{
						$qmasterport = $$ports[$i] + 10;
					}
				}
			}

			$execdport = $qmasterport + 1 if defined $qmasterport;
		}

		#### ADD NEW PORTS TO DATABASE TABLE
		$query = qq{INSERT INTO clustervars
VALUES ('$username', '$cluster', '$qmasterport', '$execdport', '$root', '$username-$cluster')};
		print "$query\n";
		my $success = $self->dbobject()->do($query);
	}


	return $qmasterport, $execdport;	
}

method setSlots ($nodetype) {
	#### GET THE SPECIFIED NODE TYPE OR, IF MISSING, USE DEFAULT
	my $cluster = $self->cluster();
	print "Agua::StarCluster::setSlots    cluster: $cluster\n";
	print "Agua::StarCluster::setSlots    nodetype: $nodetype\n" if defined $nodetype;

	$nodetype = $self->instancetype();
	print "Agua::StarCluster::setSlots    nodetype: $nodetype\n";

	my $slots = $self->instancetypeslots()->{$nodetype};
	print "Agua::StarCluster::setSlots    slots: $slots\n";

	$self->slots($slots);
}

method setQueue ($queue) {
	my $slots = $self->setSlots(undef);
	print "Agua::StarCluster::setQueue    slots: $slots\n"; 

	my $queuefile = $self->getQueuefile("queue-$queue");
	print "Agua::StarCluster::setQueue    queuefile: $queuefile\n"; 
	$self->addQueue($queue, $queuefile, $slots);
}

method setPE ($pe, $queue) {
 	print "Agua::StarCluster::setPE    pe: $pe\n";
 	print "Agua::StarCluster::setPE    queue: $queue\n";

	my $slots = $self->setPE(undef);
	print "Agua::StarCluster::setPE    slots: $slots\n"; 

	my $pefile = $self->getQueuefile("pe-$pe");
	print "Agua::StarCluster::setPE    pefile: $pefile\n"; 
	my $queuefile = $self->getQueuefile("queue-$queue");
	print "Agua::StarCluster::setPE    queuefile: $queuefile\n"; 

	$self->addPE($pe, $pefile, $slots);

	#$self->addPEToQueue($pe, $queue, $queuefile);

	print "Agua::StarCluster::setPE    Completed\n"; 
}

=head2

	SUBROUTINE		getMonitor

	PURPOSE

		OVERRIDES Agua::Common::Cluster::getMonitor

=cut

method getMonitor {
	print "Agua::Workflow::getMonitor    Agua::Workflow::getMonitor()\n";

	return $self->monitor() if $self->monitor();
	my $classfile = "Agua/Monitor/" . uc($self->clustertype()) . ".pm";
	my $module = "Agua::Monitor::" . $self->clustertype();
	require $classfile;

	my $monitor = $module->new(
		{
			'conf' 		=>	$self->conf()
		}
	);
	$self->monitor($monitor);

	return $monitor;
}

method getBalancerOutputfile ($cluster) {

#### DEBUG
#return "/tmp/loadbal";

	return $self->outputdir() . "/$cluster-loadbalancer.out";
}

method balancerOutput ($cluster, $lines){
	my $outputfile = $self->getBalancerOutputfile($cluster);
	$lines	=	20 if not defined $lines;

	my $tail = `tail -$lines $outputfile`;

	#### GET THE LAST RECORD
	my ($record) = $tail =~ /^.+(>>> \*\* \d+ \*\* Oldest job is.+)$/s;

	return $record;
}

method launchBalancer () {

	#### SET STATUS, PID, ETC. IN clusterstatus TABLE
	my $username = $self->username();
	my $cluster = $self->cluster();
	my $configfile = $self->configfile();

	#### START THE LOAD BALANCER 
	my $command = '';
	$command = "#starcluster -c $configfile bal $cluster";
	$command .= " -m ". $self->maxnodes();
	$command .= " -n ". $self->minnodes();

	my $outputfile = $self->getBalancerOutputfile($cluster);
	$command .= " &> $outputfile";

	File::Path::mkpath($self->outputdir()) if not -d $self->outputdir();
	print "Agua::StarCluster::launchBalancer    Can't create outputdir: ", $self->outputdir(), "\n" and exit if not -d $self->outputdir();

	#### NB: MUST CHANGE TO CONFIGDIR FOR PYTHON PLUGINS TO BE FOUND
	my ($configdir) = $self->configfile() =~ /^(.+?)\/[^\/]+$/;

	my $pid = fork();
	if ( $pid == 0 )
	{
		print "Running child\n";
		chdir("$configdir/plugins");
		exec("cd $configdir/plugins; $command");
		exit(0);
	}

	my $query = qq{SELECT 1 FROM clusterstatus
WHERE username='$username'
AND cluster='$cluster'};
	my $exists = $self->dbobject()->query($query);
	print "Agua::StarCluster::launchBalancer    exists: $exists\n" if defined $exists;

	my $table = "clusterstatus";	
	my $now = $self->dbobject->now();
	if ( defined $exists )
	{
		$query = qq{UPDATE $table
SET pid='$pid',
polled=$now
WHERE username='$username'
AND cluster='$cluster'};
		print "$query\n";
		$self->dbobject()->do($query);
	}
	else
	{
		$query = qq{SELECT *
FROM cluster
WHERE username='$username'
AND cluster='$cluster'};
		print "$query\n";
		my $object = $self->dbobject()->queryhash($query);
		$object->{pid} = $pid;
		$object->{started} = $now;
		$object->{polled} = $now;

		my $required = ["username","cluster"];
		my $required_fields = ["username", "cluster"];

		#### CHECK REQUIRED FIELDS ARE DEFINED
		my $not_defined = $self->dbobject()->notDefined($object, $required_fields);
		print "{ error: 'Admin::addToGroup    undefined values: @$not_defined' }" and exit if @$not_defined;

		#### DO THE ADD
		my $inserted_fields = $self->dbobject()->fields($table);
		my $success = $self->_addToTable($table, $object, $required_fields, $inserted_fields);
		print "Agua::StarCluster::launchBalancer    success: $success\n" if defined $success;
	}
}


method launchCluster () {
	print "Agua::StarCluster::launchCluster    StarCluster::launchCluster()\n";

	#### GET STARCLUSTER EXECUTABLE
	my $starcluster = $self->starcluster();
	if ( not $starcluster )
	{
		my $conf = $self->conf();
		print "Agua::StarCluster::launchCluster    self->conf not defined\n" and exit if not defined $conf;
		$starcluster = $conf->getKeyValue("applications", "STARCLUSTERDEFAULT");
	}
	print "Agua::StarCluster::launchCluster    starcluster: $starcluster\n";

	### START THE CLUSTER	
	my $command = "$starcluster -c ". $self->configfile() . " start ". $self->cluster();
	$command .= " -s " . $self->nodes() if $self->nodes() and not $self->minnodes();

	my $outputfile = $self->outputdir() . "/" . $self->cluster() . "-starcluster.out";
	$command .= " &> $outputfile";
	print "Agua::StarCluster::launchCluster    $command\n";

	my ($configdir) = $self->configfile() =~ /^(.+?)\/[^\/]+$/;
	print "Agua::StarCluster::launchCluster    configdir: $configdir\n";

	File::Path::mkpath($self->outputdir()) if not -d $self->outputdir();
	print "Agua::StarCluster::launchCluster    outputdir: ", $self->outputdir(), "\n";
	print "Agua::StarCluster::launchCluster    Can't create outputdir: ", $self->outputdir(), "\n" and exit if not -d $self->outputdir();

	my $pid = fork();
	if ( $pid == 0 )
	{
		print "Running start child\n";
		chdir("$configdir/plugins");
		exec("cd $configdir/plugins; $command");
		exit(0);
	}	
#exit;

}

=head2

	SUBROUTINE 		addNfsMounts

	PURPOSE

		CONFIGURE NFS SHARES ON HEAD NODE (I.E., NOT MASTER)

=cut

method addNfsMounts () {
	#### CHECK INPUTS
	print "Agua::StarCluster::addNfsMounts    privatekey not defined\n"and exit if not $self->privatekey();
	print "Agua::StarCluster::addNfsMounts    publiccert not defined\n"and exit if not $self->publiccert();
	print "Agua::StarCluster::addNfsMounts    username not defined\n"and exit if not $self->username();
	print "Agua::StarCluster::addNfsMounts    cluster not defined\n"and exit if not $self->cluster();
	print "Agua::StarCluster::addNfsMounts    keyname not defined\n"and exit if not $self->keyname();
	print "Agua::StarCluster::addNfsMounts    sourcedirs not defined\n"and exit if not $self->sourcedirs();
	print "Agua::StarCluster::addNfsMounts    devices not defined\n"and exit if not $self->devices();

	##### SET FIXED PORT FOR HEAD NODE NFS
	$self->setMountdPort(undef, undef);

	#### ADD ENTRIES TO /etc/fstab
	my $inserts = [];
	my $removex = [];
	for ( my $i = 0; $i < @{$self->sourcedirs()}; $i++ )
	{
		my $sourcedir = ${$self->sourcedirs()}[$i];
		my $device = ${$self->devices()}[$i];

		push @$inserts, "$device  $sourcedir    nfs     rw,vers=3,rsize=32768,wsize=32768,hard,proto=tcp 0 0";
		push @$removex, "$device\\s+$sourcedir\\s+nfs";
	}

	#### RESTART PORTMAP AND NFS DAEMONS
	$self->restartDaemons(undef, undef);

	#### SET /etc/fstab ENTRIES FOR SHARES
	$self->addFstab($removex, $inserts, undef, undef);
}


method setNfsExports ($volumes, $internalips) {

	$self->addToExports($volumes, $internalips);

	#### RESTART PORTMAP AND NFS DAEMONS
	$self->restartDaemons();
}			




=head2

	#### *** automount.py STARCLUSTER PLUGIN TAKES CARE OF THIS ***
	#### KEEP JUST IN CASE THERE ARE PROBLEMS LATER

	SUBROUTINE		monitorAddedNodes

	PURPOSE

		MONITOR ADDED STARCLUSTER NODES BY TRACKING LOAD

		BALANCER OUTPUT, E.G:

			...
			>>> *** ADDING 1 NODES at 2011-02-09 02:42:52.871015.
			>>> Launching node(s): node002
			...

	NOTES

		1. PERIODICALLY CHECK FOR '*** ADDING 1 NODES ' IN 
		   <cluster>-starcluster.out FILE

		2. FOR EACH FOUND, GET NODE NAME (E.G., 'node002') AND
		   USE starcluster listclusters <cluster> TO RETRIEVE
		   INTERNAL IP

		3. CHECK AGAINST LIST OF COMPLETED NODES

		4. IF NOT COMPLETED, SET NFS EXPORTS ON THE HEAD NODE

			FOR THIS NEWLY ADDED NODE

		5. ADD NODE TO LIST OF COMPLETED NODES

=cut
method monitorAddedNodes  ($outputfile) {
	#### SET KEYPAIR FILE
	my $keypairfile = $self->keypairfile();
	print "Agua::StarCluster::monitorAddedNodes   keypairfile: $keypairfile\n";

	#### SET CONFIG FILE
	my $configfile = $self->configfile();
	print "Agua::StarCluster::monitorAddedNodes    StarCluster::monitorAddedNodes(username, keyname, cluster)\n";

	my $sleepinterval = $self->sleepinterval();
	print "sleepinterval: $sleepinterval\n";

	####	1. PERIODICALLY CHECK FOR '*** ADDING 1 NODES ' IN 
	####	   <cluster>-starcluster.out FILE
	my $completed_nodes = [];
	while ( 1 )
	{
		print "in loop\n";
		open(FILE, $outputfile) or die "Can't open outputfile: $outputfile\n";
		$/ = "*** ADDING";
		####	2. FOR EACH FOUND, GET NODE NAME (E.G., 'node002') AND
		####	   USE starcluster listclusters <cluster> TO RETRIEVE
		####	   INTERNAL IP
		my @nodes = <FILE>;
		close(FILE);
		shift @nodes;

		####	3. CHECK AGAINST LIST OF COMPLETED NODES
		foreach my $node ( @nodes )
		{
			my $completed = 0;
			foreach my $completed_node ( @$completed_nodes )
			{
				$completed = 1 if $completed_node->{internalip} eq $node->{internalip};
				last if $completed;
			}
			####	4. IF NOT COMPLETED, SET NFS EXPORTS FOR THIS NODE ON HEAD
			$self->setNfsExports($self->sourcedirs(), [$node->{internalip}]);			

			####	5. ADD NODE TO LIST OF COMPLETED NODES
			print "StarCluslter::monitorAddedNodes    completed node: $node->{name}\n" and next if $completed;
			push @$completed_nodes, $node;
		}

		sleep($self->sleepinterval());
	}
}

method getLaunchedNodes () {
	my $command = "starcluster -c " . $self->configfile . " listclusters ". $self->cluster();
	print "Agua::StarCluster::getMasterIp    command: $command\n";
	####	    Cluster nodes:
	####         master running i-9b3e5ff7 ec2-50-17-20-70.compute-1.amazonaws.com 
	####        node001 running i-953e5ff9 ec2-67-202-10-108.compute-1.amazonaws.com 
	####    Total nodes: 2
	my $list = `$command`;
	print "Agua::StarCluster::getMasterIp    list result not defined or empty\n"
		and return if not defined $list or not $list;
	print "Agua::StarCluster::getMasterIp    list: $list\n";

	my ($nodelist) = $list =~ /Cluster nodes:\s*\n(.+)$/ms;
	my @lines = split "\n", $nodelist;
	my $nodes = [];
	foreach my $line ( @lines )
	{
		next if $line =~ /^\s*$/;
		my ($name, $status, $instanceid, $internalip) =
			$line =~ /^\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;
		push @$nodes, {
			name=> $name,
			status=>$status,
			instanceid=> $instanceid,
			internalip=> $internalip
		};
	}

	return $nodes;
}

method mountShares () {
	#### SET KEYPAIR FILE
	my $keypairfile = $self->keypairfile();
	print "Agua::StarCluster::mountShares    keypairfile: $keypairfile\n";

	#### SET CONFIG FILE
	my $configfile = $self->configfile();

	#### WHILE CLUSTER IS STARTING UP, SEARCH FOR IPS OF NODES
	#### AND ADD THEM TO /etc/exports
	my $completed_nodes = [];
	while ( scalar(@$completed_nodes) < $self->nodes() )
	{
		#### GET INTERNAL IPS, ETC. OF ALL NODES IN CLUSTER
		my $launched_nodes = $self->getLaunchedNodes();

		#### IGNORE ALREADY COMPLETED NODES
		foreach my $launched_node ( @$launched_nodes )
		{
			my $completed = 0;
			foreach my $completed_node ( @$completed_nodes )
			{
				$completed = 1 if $completed_node->{internalip} eq $launched_node->{internalip};
				last if $completed;
			}
			next if $completed;

			push @$completed_nodes, $launched_node;

			if ( $launched_node->{name} eq "master" )
			{
				##### AUTOMATE RESTART sge MASTER AND EXECD AFTER REBOOT MASTER
				my $masterip = $launched_node->{internalip};
				$self->setSgeStartup($self->username(), $masterip, $keypairfile);

				#### EXCLUDE MASTER NODE FROM EXEC NODES LIST
				$self->excludeMasterHost($masterip);
			}

			#### ADD ENTRIES TO /etc/exports
			$self->setNfsExports($self->sourcedirs, $launched_node->{internalip});			}
	}
}

=head

	SUBROUTINE		autoMount

	PURPOSE

		MAKE REMOTE CALLS TO NODES TO SET UP NFS MOUNTS

	#############################################################
	####
	####	DEPRECATED		DEPRECATED		DEPRECATED
	####
	#### automount.py STARCLUSTER PLUGIN TAKES CARE OF THIS
	####
	####	KEEP JUST IN CASE THERE ARE PROBLEMS LATER
	####
	#############################################################
=cut

method autoMount () {
	#### SET KEYPAIR FILE
	my $keypairfile = $self->keypairfile();
	print "Agua::StarCluster::mountShares    keypairfile: $keypairfile\n";

	### GET INTERNAL IP OF HEAD NODE
	my ($externalip, $headip) = $self->getLocalIps();
	print "Agua::StarCluster::addNfsMounts    headip: $headip\n";

	#### GET INTERNAL IPS OF ALL NODES
	my $nodeips = $self->getInternalIps();
	print "Agua::StarCluster::addNfsMounts    nodeips: @$nodeips\n";

	#### MOUNT ON MASTER AND EXEC NODES
	foreach my $nodeip ( @$nodeips )
	{
		my $inserts = [];
		my $removex = [];
		for ( my $i = 0; $i < @{$self->sourcedirs()}; $i++ )
		{
			my $sourcedir = ${$self->sourcedirs()}[$i];
			my $mountpoint = ${$self->mountpoints()}[$i];

			$self->mountNfs($sourcedir, $headip, $mountpoint, $keypairfile, $nodeip);

			push @$inserts, "$headip:$sourcedir  $mountpoint nfs nfsvers=3,defaults 0 0";
			push @$removex, "$headip:$sourcedir";
		}

		$self->addFstab($removex, $inserts, $keypairfile, $nodeip);
	}
}


=head2

	SUBROUTINE 		excludeMasterHost

	PURPOSE

		1. REMOVE MASTER FROM EXEC HOST LIST

		2. REMOVE MASTER FROM CONFIGURATION LIST:

		3. REMOVE MASTER FROM HOST LIST

		4. RESTART sgemaster DAEMON

	NOTES

		MUST BE RUN AS root

=cut

method excludeMasterHost ($targetip) {
	#### targetip IS THE HOST MOUNTING THE SHARE
	print "Agua::StarCluster::excludeMasterHost    StarCluster::excludeMasterHost(targetip)\n";
	print "Agua::StarCluster::excludeMasterHost    targetip: $targetip\n";

	#### SET CONFIG FILE
	my $configfile 		= 	$self->configfile();
	my $keypairfile 		= 	$self->keypairfile();

	#### SET REMOTE SSH 
	my $remotessh = $self->setRemoteSsh($targetip);

	#### SET QCONF
	my $conf 			= 	$self->conf();
	my $sge_root		= 	$conf->getKeyValue("cluster", "SGE_ROOT");
	my $sge_qmaster_port= 	$conf->getKeyValue("cluster", "SGE_QMASTER_PORT");
	my $qconf			= 	$conf->getKeyValue("cluster", "QCONF");

	#### 1. REMOVE MASTER FROM EXEC HOST LIST
	#	Host object "ip-10-124-245-118.ec2.internal" is still referenced in cluster queue "all.q"
	my $exec = qq{$remotessh "export SGE_ROOT=$sge_root; export SGE_QMASTER_PORT=$sge_qmaster_port; $qconf -de $targetip"};
	print "\n$exec\n\n";
	my $exec_result = `$exec`;
	print "Agua::StarCluster::excludeMasterHost    BEFORE exec_result: $exec_result\n";


	#### 2. REMOVE MASTER FROM CONFIGURATION LIST:
	#    root@ip-10-124-245-118.ec2.internal removed "ip-10-124-245-118.ec2.internal" from configuration list
	my $config = qq{$remotessh "export SGE_ROOT=$sge_root; export SGE_QMASTER_PORT=$sge_qmaster_port; $qconf -dconf $targetip"};
	print "\n$config\n\n";
	my $config_result = `$config`;
	print "Agua::StarCluster::excludeMasterHost    BEFORE config_result: $config_result\n";
#
	#### 3. REMOVE MASTER FROM HOST LIST
	#### GET CURRENT GROUP CONFIG
	#    group_name @allhosts
	#    hostlist ip-10-124-245-118.ec2.internal ip-10-124-247-224.ec2.internal
	my $hostlist = `$remotessh qconf -shgrp \@allhosts`;
	$hostlist =~ s/\s+$//;
	print "Agua::StarCluster::excludeMasterHost    BEFORE hostlist: $hostlist\n";

	#### WRITE TEMP FILE TO USER-OWNED /tmp DIRECTORY
	#### IN PREPARATION FOR COPY AS ROOT TO REMOTE HOST
	my $tempdir = "/tmp/" . $self->username();
	File::Path::mkpath($tempdir) if not -d $tempdir;
	my $groupfile = "allhosts.group";
	open(OUT, ">/$tempdir/$groupfile") or die "Can't open groupfile: $groupfile\n";
	print OUT $hostlist;
	close(OUT) or die "Can't close groupfile: /$tempdir/$groupfile\n";

	#### COPY GROUP CONFIG FILE TO '~' ON TARGET HOST
	my $copy;
	$copy = qq{scp -i $keypairfile /$tempdir/$groupfile root\@$targetip:$groupfile} if defined $keypairfile;
	$copy = qq{cp /$tempdir/$groupfile ~/$groupfile} if not defined $keypairfile;
	print "Agua::StarCluster::excludeMasterHost    copy: $copy\n";
	my $result = `$copy`;	
	print "Agua::StarCluster::excludeMasterHost    result: $result\n";

	#### SET GROUP CONFIG FROM FILE ON TARGET HOST
	my $group = qq{$remotessh "export SGE_ROOT=$sge_root; export SGE_QMASTER_PORT=$sge_qmaster_port; $qconf -Mhgrp ~/$groupfile"};
	print "\n$group\n\n";
	my $qconf_result = `$group`;
	#    root@ip-10-124-245-118.ec2.internal modified "@allhosts" in host group list
	print "Agua::StarCluster::excludeMasterHost    qconf_result: $qconf_result\n";

	#### RESTART sgemaster
	my $restart = "$remotessh /etc/init.d/sgemaster.starcluster stop; /etc/init.d/sgemaster.starcluster start";
	print "\n$restart\n\n";
	my $restart_result = `$restart`;
	print "Agua::StarCluster::excludeMasterHost    restart_result: $restart_result\n";
}


=head2

	SUBROUTINE		setMountdPort

	PURPOSE

		ADD FIXED PORT FOR mountd IN /etc/default/nfs-kernel-server

			- 	ON REMOTE HOST IF KEYPAIRFILE AND IP PROVIDED

			-	OTHERWISE, ON LOCAL HOST

=cut

method setMountdPort ($keypairfile, $targetip) {	
	print "Agua::StarCluster::setMountdPort    StarCluster::setMountdPort(keypairfile, targetip)\n";

	#### SET REMOTE SSH 
	my $remotessh = $self->setRemoteSsh($keypairfile, $targetip);

	my $uname = `$remotessh uname -a`;
	print "Agua::StarCluster::setMountdPort    uname: $uname\n";

	my $conf = $self->conf();
	my $mountdport = $conf->getKeyValue("starcluster:nfs", "MOUNTDPORT");	
	print "Agua::StarCluster::setMountdPort    mountdport: $mountdport\n";

	my $insert = qq{MOUNTD_PORT=32767};
	my $filter = qq{MOUNTD_PORT};
	my $nfsconfigfile = "/etc/sysconfig/nfs";
	if ( $uname =~ /ubuntu/i )
	{
		$insert = qq{RPCMOUNTDOPTS="--port $mountdport --manage-gids"};
		$filter = qq{RPCMOUNTDOPTS};
		$nfsconfigfile = "/etc/default/nfs-kernel-server";
	}
	print "Agua::StarCluster::setMountdPort    insert: $insert\n";
	print "Agua::StarCluster::setMountdPort    filter: $filter\n";
	print "Agua::StarCluster::setMountdPort    nfsconfigfile: $nfsconfigfile\n";

	#### GET NFS CONFIG FILE CONTENTS
	my $cat = qq{$remotessh cat $nfsconfigfile}; 
	print "Agua::StarCluster::setMountdPort    cat: $cat\n";
	my $nfsconfig = `$cat`;	
	print "Agua::StarCluster::setMountdPort    BEFORE nfsconfig: $nfsconfig\n";

	#### BACKUP NFS CONFIG FILE
	my $move = qq{$remotessh mv -f $nfsconfigfile $nfsconfigfile.bkp}; 
	print "Agua::StarCluster::setMountdPort    move: $move\n";
	print `$move`;	

	#### COMMENT OUT EXISTING RPCMOUNTOPTS LINES
	my @lines = split "\n", $nfsconfig;
	for ( my $i = 0; $i < $#lines + 1; $i++ )
	{
		if ( $lines[$i] =~ /^$filter/ )
		{
			splice @lines, $i, 1;
			$i--;
		}
	}

	#### ADD NEW RPCMOUNTOPTS LINE
	push @lines, "$insert\n";

	#### WRITE TEMP FILE TO USER-OWNED /tmp DIRECTORY
	#### IN PREPARATION FOR COPY AS ROOT TO REMOTE HOST
	my $tempdir = "/tmp/" . $self->username();
	File::Path::mkpath($tempdir) if not -d $tempdir;
	my $tempfile = "$tempdir/nfs-kernel-server";
	open(OUT, ">$tempfile") or die "Can't open tempfile: $tempfile\n";
	my $content = join "\n", @lines;
	print OUT $content;
	close(OUT) or die "Can't close tempfile: $tempfile\n";

	my $copy = qq{cp $tempfile $nfsconfigfile};
	if ( defined $keypairfile )
	{
		$copy = qq{scp -i $keypairfile $tempfile root\@$targetip:$nfsconfigfile};
	}
	print "Agua::StarCluster::setMountdPort    copy: $copy\n";
	my $result = `$copy`;	
	print "Agua::StarCluster::setMountdPort    result: $result\n";	

	##### DEBUG
	#my $catnew = qq{$remotessh cat $nfsconfigfile}; 
}

=head2

	SUBROUTINE		addToExports

	PURPOSE

		ON MASTER, SET UP EXPORT TO HEAD INSTANCE IN /etc/exports:

		/home ip-10-124-247-224.ec2.internal(async,no_root_squash,no_subtree_check,rw)
		/opt/sge6 ip-10-124-247-224.ec2.internal(async,no_root_squash,no_subtree_check,rw)
		/data ip-10-124-247-224.ec2.internal(async,no_root_squash,no_subtree_check,rw)
	*** /data ip-10-127-158-202.ec2.internal(async,no_root_squash,no_subtree_check,rw)

=cut

method addToExports ($volumes, $recipientips, $targetip, $keypairfile) {
	#### sourceip IS THE HOST DOING THE SHARING
	#### targetip IS THE HOST MOUNTING THE SHARE
	print "Agua::StarCluster::addToExports    StarCluster::addToExports(username, keypairfile, targetip, sourceip, volume)\n";
	print "Agua::StarCluster::addToExports    username: $self->username()\n";
	print "Agua::StarCluster::addToExports    recipientips: @$recipientips\n";
	print "Agua::StarCluster::addToExports    volume: @$volumes\n";
	print "Agua::StarCluster::addToExports    keypairfile: $keypairfile\n" if defined $keypairfile;
	print "Agua::StarCluster::addToExports    targetip: $targetip\n" if defined $targetip;

	#### SET CONFIG FILE
	my $configfile = $self->configfile();

	#### SET REMOTE SSH 
	my $remotessh = $self->setRemoteSsh($targetip);

	#### GET CONTENTS OF /etc/exports
	my $exportsfile = "/etc/exports";
	my $cat = qq{$remotessh cat $exportsfile}; 
	my $exports = `$cat`;
	$exports =~ s/\s+$//;

	#### REMOVE EXISTING ENTRY FOR THESE VOLUMES
	my @lines = split "\n", $exports;
	foreach my $volume ( @$volumes )
	{
		foreach my $recipientip ( @$recipientips )
		{
			for ( my $i = 0; $i < $#lines + 1; $i++ )
			{
				if ( $lines[$i] =~ /^$volume\s+$recipientip/ )
				{
					splice @lines, $i, 1;
					$i--;
				}
			}
		}
	}

	foreach my $volume ( @$volumes )
	{
		foreach my $recipientip ( @$recipientips )
		{
			push @lines, "$volume $recipientip(async,no_root_squash,no_subtree_check,rw)";
		}
	}
	my $output = join "\n", @lines;

	my $move = qq{$remotessh mv -f $exportsfile $exportsfile.bkp}; 
	`$move`;	

	#### WRITE TEMP FILE TO USER-OWNED /tmp DIRECTORY
	#### IN PREPARATION FOR COPY AS ROOT TO REMOTE HOST
	my $tempdir = "/tmp/" . $self->username();
	File::Path::mkpath($tempdir) if not -d $tempdir;
	my $tempfile = "/$tempdir/exports";
	open(OUT, ">$tempfile") or die "Can't open tempfile: $tempfile\n";
	print OUT $output;
	close(OUT) or die "Can't close tempfile: $tempfile\n";

	my $copy;
	$copy = qq{scp -i $keypairfile $tempfile root\@$targetip:$exportsfile} if defined $keypairfile;
	$copy = qq{cp $tempfile $exportsfile} if not defined $keypairfile;

	my $result = `$copy`;	

	##### DEBUG
	#my $catnew = qq{$remotessh cat $exportsfile}; 
}



=head2

	SUBROUTINE		addFstab

	PURPOSE

		ADD ENTRIES TO /etc/fstab TO AUTOMOUNT NFS IMPORTS, I.E.:

			/dev/sdh  /data      nfs     rw,vers=3,rsize=32768,wsize=32768,hard,proto=tcp 0 0
			/dev/sdi  /nethome      nfs     rw,vers=3,rsize=32768,wsize=32768,hard,proto=tcp 0 0

=cut

method addFstab ($removex, $inserts, $keypairfile, $targetip) {
	#### targetip IS THE HOST MOUNTING THE SHARE
	print "Agua::StarCluster::addFstab    StarCluster::addFstab(keypairfile, targetip, sourceip, volume)\n";
	print "Agua::StarCluster::addFstab    removex: $removex\n";
	print "Agua::StarCluster::addFstab    inserts: @$inserts\n";
	print "Agua::StarCluster::addFstab    targetip: $targetip\n" if defined $targetip;

	#### SET CONFIG FILE
	my $configfile = $self->configfile();

	#### SET REMOTE SSH 
	my $remotessh = $self->setRemoteSsh($keypairfile, $targetip);

	#### GET CONTENTS OF /etc/exports
	my $fstabfile = "/etc/fstab";
	my $cat = qq{$remotessh cat $fstabfile}; 
	my $exports = `$cat`;
	$exports =~ s/\s+$//;

	#### REMOVE EXISTING ENTRY FOR THESE VOLUMES
	my @lines = split "\n", $exports;
	for ( my $i = 0; $i < $#lines + 1; $i++ )
	{
		foreach my $remove ( @$removex )
		{
			if ( $lines[$i] =~ /^$remove/ )
			{
				splice @lines, $i, 1;
				$i--;
			}
		}
	}

	foreach my $insert ( @$inserts )
	{
		push @lines, $insert;
	}
	my $output = join "\n", @lines;

	my $move = qq{$remotessh mv -f $fstabfile $fstabfile.bkp}; 
	print "Agua::StarCluster::addFstab    move: $move\n";
	`$move`;	

	#### WRITE TEMP FILE TO USER-OWNED /tmp DIRECTORY
	#### IN PREPARATION FOR COPY AS ROOT TO REMOTE HOST
	my $tempdir = "/tmp/" . $self->username();
	File::Path::mkpath($tempdir) if not -d $tempdir;
	my $tempfile = "$tempdir/exports";
	open(OUT, ">$tempfile") or die "Can't open tempfile: $tempfile\n";
	print OUT $output;
	close(OUT) or die "Can't close tempfile: $tempfile\n";

	my $copy;
	$copy = qq{scp -i " . $self->keypairfile() . " $tempfile root\@$targetip:$fstabfile} if defined $targetip;
	$copy = qq{cp $tempfile $fstabfile} if not defined $targetip;

	print "Agua::StarCluster::addFstab    copy: $copy\n";
	my $result = `$copy`;	

	#### DEBUG
	my $catnew = qq{$remotessh cat $fstabfile}; 
	print "Agua::StarCluster::addFstab    cat: $catnew\n";
	print `$catnew`;
	print "\n";
}


=head2

	SUBROUTINE 		setSgeStartup

	PURPOSE

		ADD TO /etc/init.d/rc.local

=cut

method setSgeStartup ($targetip) {
	my $removex = [
	"",
	"/etc/init.d/sgemaster.starcluster start",
	"echo 'sgemaster.starcluster started'",
	"/etc/init.d/sgeexecd.starcluster start",
	"echo 'sgeexecd.starcluster started'"
];
	my $inserts;
	@$inserts = @$removex;

	#### targetip IS THE HOST MOUNTING THE SHARE
	print "Agua::StarCluster::setSgeStartup    StarCluster::setSgeStartup(username, keypairfile, targetip, sourceip, volume)\n";
	print "Agua::StarCluster::setSgeStartup    username: $self->username()\n";
	print "Agua::StarCluster::setSgeStartup    removex: $removex\n";
	print "Agua::StarCluster::setSgeStartup    inserts: @$inserts\n";
	print "Agua::StarCluster::setSgeStartup    targetip: $targetip\n";

	#### SET REMOTE SSH 
	my $remotessh = $self->setRemoteSsh($targetip);

	#### GET CONTENTS OF /etc/exports
	my $startupfile = "/etc/init.d/rc.local";
	my $cat = qq{$remotessh cat $startupfile}; 
	my $contents = `$cat`;
	$contents =~ s/\s+$//;

	#### REMOVE EXISTING ENTRY FOR THESE VOLUMES
	my @lines = split "\n", $contents;
	for ( my $i = 0; $i < $#lines + 1; $i++ )
	{
		foreach my $remove ( @$removex )
		{
			if ( $lines[$i] =~ /^$remove/ )
			{
				splice @lines, $i, 1;
				$i--;
			}
		}
	}

	foreach my $insert ( @$inserts )
	{
		push @lines, $insert;
	}
	my $output = join "\n", @lines;

	my $move = qq{$remotessh mv -f $startupfile $startupfile.bkp}; 
	print "Agua::StarCluster::setSgeStartup    move: $move\n";
	`$move`;	

	#### WRITE TEMP FILE TO USER-OWNED /tmp DIRECTORY
	#### IN PREPARATION FOR COPY AS ROOT TO REMOTE HOST
	my $tempdir = "/tmp/" . $self->username();
	File::Path::mkpath($tempdir) if not -d $tempdir;
	my $tempfile = "$tempdir/startup";
	open(OUT, ">$tempfile") or die "Can't open tempfile: $tempfile\n";
	print OUT $output;
	close(OUT) or die "Can't close tempfile: $tempfile\n";

	my $copy;
	$copy = qq{scp -i ". $self->keypairfile() . " $tempfile root\@$targetip:$startupfile} if defined $targetip;
	$copy = qq{cp $tempfile $startupfile} if not defined $targetip;

	print "Agua::StarCluster::setSgeStartup    copy: $copy\n";
	my $result = `$copy`;	

	#### DEBUG
	my $catnew = qq{$remotessh cat $startupfile}; 
	print "Agua::StarCluster::setSgeStartup    cat: $catnew\n";
	print `$catnew`;	
}





=head2

	SUBROUTINE		mountNfs

	PURPOSE

		MOUNT EXPORTED NFS volume FROM MASTER ON AQ-7

	NOTES

		CHECK IF MASTER'S MOUNT IS SEEN BY AQ-7:

		showmount -e ip-10-124-245-118

			Export list for ip-10-124-245-118:
			/data     ip-10-127-158-202.ec2.internal,ip-10-124-247-224.ec2.internal
			/opt/sge6 ip-10-124-247-224.ec2.internal
			/home     ip-10-124-247-224.ec2.internal

=cut

method mountNfs ($source, $sourceip, $mountpoint, $keypairfile, $targetip) {	
	#### targetip IS THE HOST MOUNTING THE SHARE
	print "Agua::StarCluster::mountNfs    StarCluster::mountNfs(source, sourceip, mountpoint, keypairfile, targetip)\n";
	print "Agua::StarCluster::mountNfs    source: $source\n";
	print "Agua::StarCluster::mountNfs    sourceip: $sourceip\n";
	print "Agua::StarCluster::mountNfs    mountpoint: $mountpoint\n";
	print "Agua::StarCluster::mountNfs    keypairfile: $keypairfile\n";
	print "Agua::StarCluster::mountNfs    targetip: $targetip\n";

	#### SET REMOTE SSH 
	my $remotessh = $self->setRemoteSsh($targetip);

	#### CREATE MOUNTPOINT DIRECTORY
	my $mkdir = "$remotessh mkdir -p $mountpoint";
	print "$mkdir\n";
	print `$mkdir`;

	#### MOUNT NFS SHARE TO MOUNTPOINT
	my $mount = qq{$remotessh mount -t nfs $sourceip:$source $mountpoint}; 
	print "$mount\n";
	print `$mount`;
}


=head2

	SUBROUTINE		setKeypairfile

	PURPOSE

		SET SSH COMMAND IF KEYPAIRFILE, ETC. ARE DEFINED
=cut

method setKeypairfile () {
	print "Agua::StarCluster::setKeypairfile    username not defined\n" and exit if not defined $self->username();
	$self->keyname() = $self->username() . "-key" if not defined $self->keyname();

	print "Agua::StarCluster:::setKeypairfile     username not defined: $self->username()\n" and exit if not defined $self->username();	

	my $conf 		= 	$self->conf();
	my $userdir 	= 	$conf->getKeyValue('agua', "USERDIR");
	print "Agua::StarCluster::setKeypairfile    userdir not defined\n"and exit if not defined $userdir;

	return "$userdir/" . $self->username() . "/.starcluster/id_rsa-". $self->keyname();
}



=head2

	SUBROUTINE		setRemoteSsh

	PURPOSE

		SET SSH COMMAND IF KEYPAIRFILE, ETC. ARE DEFINED
=cut

method setRemoteSsh ($keypairfile, $targetip) {
	my $remotessh = '';
	$remotessh = "ssh -o StrictHostKeyChecking=no -i " . $self->keypairfile() . " root\@$targetip" if defined $targetip;

	return $remotessh;
}


=head2

	SUBROUTINE		restartDaemons

	PURPOSE

		ON MASTER, RESTART NFS

=cut

method restartDaemons ($keypairfile, $targetip) {
	#### SET REMOTE SSH 
	my $remotessh = $self->setRemoteSsh($keypairfile, $targetip);

	#### CHECK LINUX FLAVOUR
	my $uname = `$remotessh uname -a`;

	#### SET BINARIES ACCORDING TO FLAVOUR
	my $portmapbinary = "service portmap";
	my $nfsbinary = "service nfs";
	if ( $uname =~ /ubuntu/i )
	{
		$portmapbinary = "/etc/init.d/portmap";
		$nfsbinary = "/etc/init.d/nfs";
	}

	#### RESTART SERVICES
	my $portmap = qq{$remotessh $portmapbinary restart};
	print `$portmap`;

	my $nfs = qq{$remotessh $nfsbinary restart};
	print `$nfs`;
}

=head2

	SUBROUTINE 		getInstanceInfo

	PURPOSE

		RETURN A HASH OF INSTANCE VARIABLES

	NOTES

		RESERVATION     r-e30f5d89      558277860346    default
		INSTANCE        i-b42f3fd9      ami-90af5ef9    ec2-75-101-214-196.compute-1.amazonaws.com      ip-10-127-158-202.ec2.internal running aquarius        0               t1.micro        2010-12-24T09:51:37+0000        us-east-1a     aki-b51cf9dc    ari-b31cf9da            monitoring-disabled     75.101.214.196  10.127.158.202        ebs
		BLOCKDEVICE     /dev/sda1       vol-c6e346ae    2010-12-24T09:51:40.000Z
		BLOCKDEVICE     /dev/sdh        vol-266dc84e    2010-12-24T23:03:04.000Z
		BLOCKDEVICE     /dev/sdi        vol-fa6dc892    2010-12-24T23:05:50.000Z

=cut

method getInstanceInfo () {

	my $instanceid = `curl -s http://169.254.169.254/latest/meta-data/instance-id`;
	my $command = qq{ec2-describe-instances \\
-K } . $self->privatekey() . qq{ \\
-C } . $self->publiccert();	
	print "$command\n";
	my $result = `$command`;

	my $instance;
	$instance->{blockdevices} = [];
	my @reservations = split "RESERVATION", $result;
	foreach my $reservation ( @reservations )
	{
		next if $reservation !~ /\s+$instanceid\s+/;
		my @lines = split "\n", $reservation;
		foreach my $line ( @lines )
		{
			push @{$instance->{blockdevices}}, $line if $line =~ /^BLOCKDEVICE/;
			$instance->{instance} = $line if $line =~ /^INSTANCE/;
			$instance->{reservation} = $line if $line =~ /^\s+/;
		}
	}

	return $instance;
}

method getLocalIps () {
	my $instance = $self->getInstanceInfo();
	my @elements = split " ", $instance->{instance};

	return $elements[3], $elements[4];
}


=head2

	SUBROUTINE		openPorts

	PURPOSE

		OPEN NFS PORTS FOR THE GIVEN GROUP ON EC2

=cut

method openPorts ($group) {
	print "Agua::StarCluster::openPorts    StarCluster::openPorts(privatekey, publiccert, group)\n";
	print "Agua::StarCluster::openPorts    group not defined or empty\n"
		and exit if not defined $group or not $group;

	my $conf 	= 	$self->conf();
	my $portmap	= 	$conf->getKeyValue("starcluster:nfs", "PORTMAPPORT");
	my $nfs 	= 	$conf->getKeyValue("starcluster:nfs", "NFSPORT");
	my $mountd 	= 	$conf->getKeyValue("starcluster:nfs", "MOUNTDPORT");
	my $sge 	= 	$conf->getKeyValue("cluster", "SGEQMASTERPORT");
	my $ec2 	= 	$conf->getKeyValue("applications", "EC2");	
	my $java_home= 	$conf->getKeyValue("aws", "JAVAHOME");	
	print "Agua::StarCluster::openPorts    ec2 not defined\n" and exit if not defined $ec2;

	#### SET EC2_HOME ENVIRONMENT VARIABLE
	my $ec2_home = $ec2;
	$ec2_home =~ s/\/bin$//;

	$ENV{'EC2_HOME'} = $ec2_home;
	$ENV{'JAVA_HOME'} = $java_home;

	##### CREATE SECURITY GROUP
	#my $creategroup = "$ec2/ec2-add-group $group -d 'StarCluster group'";

	my $tasks = [
		#### PORTMAP
		"$group -p $portmap -P tcp",
		"$group -p $portmap -P udp",

		#### NFS
		"$group -p $nfs -P udp",
		"$group -p $nfs -P tcp",

		#### MOUNTD
		"$group -p $mountd -P udp",
		"$group -p $mountd -P tcp",

		#### SGE_QMASTER_PORT
		"$group -p $sge -P udp",
		"$group -p $sge -P tcp"
	];

	#### RUN COMMANDS
	foreach my $task ( @$tasks )
	{
		my $command = qq{$ec2/ec2-authorize \\
-K } . $self->privatekey() . qq{ \\
-C } . $self->publiccert() . qq{ \\
$task\n};
		print $command;
		print `$command`;
	}
}

method getMasterIp () {
	#### SET CONFIG FILE
	my $configfile = $self->configfile();

	print "Agua::StarCluster::getMasterIp    StarCluster::getMasterIp(username, keyname)\n";

    #### StarCluster - (http://web.mit.edu/starcluster)
    #### Software Tools for Academics and Researchers (STAR)
    #### Please submit bug reports to starcluster@mit.edu
	####        
    #### -----------------------------------------------
    #### smallcluster (security group: @sc-smallcluster)
    #### -----------------------------------------------
    #### Launch time: 2010-10-15T02:53:52.000Z
    #### Zone: us-east-1a
    #### Keypair: gsg-keypair
    #### Cluster nodes:
    ####      master running i-9edc26f3 ec2-174-129-54-141.compute-1.amazonaws.com 
    ####     node001 running i-98dc26f5 ec2-75-101-225-244.compute-1.amazonaws.com 
	my $command = "starcluster -c $configfile listclusters";
	print "Agua::StarCluster::getMasterIp    command: $command\n";
	my $list = `$command`;
	print "Agua::StarCluster::getMasterIp    list: $list\n";
	my ($entry) = $list =~ /" . $self->cluster()\s+\(security group:.+?\n\-+\n(.+)$/ms;
	print "Agua::StarCluster::getMasterIp    entry: $entry\n";
	my ($masterip) = $entry =~ /^.+?master.+?(ec2\S+)/ms;
	print "Agua::StarCluster::getMasterIp    masterip: $masterip\n";

	return $masterip;
}


=head

	SUBROUTINE		getInstanceIds

	PURPOSE

		RETRIEVE THE INTERNAL IPS FOR ALL EXECUTION NODES

		IN THE GIVEN CLUSTER OWNED BY THIS USER

	NOTES

		1. EXTRACT EXTERNAL IPS FOR STARCLUSTER EXEC NODES

		starcluster -c /path/to/configfile listclusters

			-----------------------------------------------
			smallcluster (security group: @sc-smallcluster)
			-----------------------------------------------
			Launch time: 2011-01-06T14:11:09.000Z
			Zone: us-east-1a
			Keypair: id_rsa-admin-key
			EBS volumes:
				vol-fc5af194 on master:/dev/sdj (status: attached)
			Cluster nodes:
				 master running i-6f21d203 ec2-72-44-59-38.compute-1.amazonaws.com 
				node001 running i-6921d205 ec2-67-202-9-15.compute-1.amazonaws.com 


		2. LOOK UP INTERNAL IPS USING EXTERNAL IPS IN RESERVATION

			RECORDS OUTPUT BY ec2-describe-instances

=cut

method getInstanceIds () {
	print "Agua::StarCluster::getInstanceIds    StarCluster::getInstanceIds()\n";
	#### 1. GET LIST OF NODES FOR THIS CLUSTER	
    #### -----------------------------------------------
    #### smallcluster (security group: @sc-smallcluster)
    #### -----------------------------------------------
    #### Launch time: 2010-10-15T02:53:52.000Z
    #### Zone: us-east-1a
    #### Keypair: gsg-keypair
    #### Cluster nodes:
    ####      master running i-9edc26f3 ec2-174-129-54-141.compute-1.amazonaws.com 
    ####     node001 running i-98dc26f5 ec2-75-101-225-244.compute-1.amazonaws.com 

	my $command = "starcluster -c ". $self->configfile() ." listclusters " . $self->cluster();
	print "Agua::StarCluster::start    command: $command\n";
	my $list = `$command`;
	print "Agua::StarCluster::start    list: $list\n";
	my ($entry) = $list =~ /" . $self->cluster()\s+\(security group:.+?\n[-]+\n(.+)([-]+)*/ms;
	print "Agua::StarCluster::getInstanceIds    entry: $entry\n";
	my ($lines) = $entry =~ /Cluster nodes:\s*\n(.+)$/ms;
	print "Agua::StarCluster::getInstanceIds    lines: $lines\n";
	my $instanceids = [];
	my @array = split "\n", $lines;
	foreach my $line ( @array )
	{
		next if $line =~ /^\s*$/;
		push @$instanceids, $2 if $line =~ /^\s+(\S+)\s+\S+\s+(i-\S+)/;
	}
	print "Agua::StarCluster::getInstanceIds    instanceids: @$instanceids\n";

	return $instanceids;
}


=head

	SUBROUTINE		getInternalIps

	PURPOSE

		RETRIEVE THE INTERNAL IPS FOR ALL EXECUTION NODES

		IN THE GIVEN CLUSTER OWNED BY THIS USER

	NOTES

		1. EXTRACT EXTERNAL IPS FOR STARCLUSTER EXEC NODES

		starcluster -c /path/to/configfile listclusters

			-----------------------------------------------
			smallcluster (security group: @sc-smallcluster)
			-----------------------------------------------
			Launch time: 2011-01-06T14:11:09.000Z
			Zone: us-east-1a
			Keypair: id_rsa-admin-key
			EBS volumes:
				vol-fc5af194 on master:/dev/sdj (status: attached)
			Cluster nodes:
				 master running i-6f21d203 ec2-72-44-59-38.compute-1.amazonaws.com 
				node001 running i-6921d205 ec2-67-202-9-15.compute-1.amazonaws.com 


		2. LOOK UP INTERNAL IPS USING EXTERNAL IPS IN RESERVATION

			RECORDS OUTPUT BY ec2-describe-instances

=cut

method getInternalIps () {
	#### SET CONFIG FILE
	my $configfile = $self->configfile();

	print "Agua::StarCluster::getInternalIps    StarCluster::getInternalIps()\n";
	print "Agua::StarCluster::getInternalIps    privatekey: ", $self->privatekey(), "\n";
	print "Agua::StarCluster::getInternalIps    publiccert: ", $self->publiccert(), "\n";

	#### 1. GET LIST OF NODES FOR THIS CLUSTER	
    #### -----------------------------------------------
    #### smallcluster (security group: @sc-smallcluster)
    #### -----------------------------------------------
    #### Launch time: 2010-10-15T02:53:52.000Z
    #### Zone: us-east-1a
    #### Keypair: gsg-keypair
    #### Cluster nodes:
    ####      master running i-9edc26f3 ec2-174-129-54-141.compute-1.amazonaws.com 
    ####     node001 running i-98dc26f5 ec2-75-101-225-244.compute-1.amazonaws.com 

	my $instanceids = $self->getInstanceIds();
	print "Agua::StarCluster::getInternalIps    instanceids: @$instanceids\n";

	my $clusternodes = $self->getClusterNodes();

	my $internalips = [];
	foreach my $instanceid ( @$instanceids )
	{
		foreach my $clusternode ( @$clusternodes )
		{
			if ( $instanceid eq $clusternode->{instanceid} )
			{
				my @elements = split "\t", $clusternode->{instance};
				push @$internalips,  $elements[4];
			}
		}
	}
	print "Agua::StarCluster::getInternalIps    internalips: @$internalips\n";

	return $internalips;
}



method getHeadInternalIps () {
	my $instanceInfo = $self->getInstanceInfo();
	print "Agua::StarCluster::getHeadInternalIP    Could not retrieve instanceInfo\n" and exit if not defined $instanceInfo;
	my @instanceRow = split "\t", $instanceInfo->{instance};
print "instanceRow: @instanceRow\n";

	return  $instanceRow[17], $instanceRow[4];
}

=head

	SUBROUTINE		getClusterNodes

	PURPOSE

		RETRIEVE THE INSTANCE RECORDS FOR ALL EXECUTION NODES

		IN THE GIVEN CLUSTER OWNED BY THIS USER

	NOTES

			ec2din OUTPUT EXAMPLE:

			RESERVATION     r-e30f5d89      558277860346    default
			INSTANCE        i-b42f3fd9      ami-90af5ef9    ec2-75-101-214-196.compute-1.amazonaws.com      ip-10-127-158-202.ec2.internal  running aquarius        0   t1.micro 2010-12-24T09:51:37+0000        us-east-1a      aki-b51cf9dc    ari-b31cf9da            monitoring-disabled     75.101.214.196  10.127.158.202      ebs
			BLOCKDEVICE     /dev/sda1       vol-c6e346ae    2010-12-24T09:51:40.000Z
			BLOCKDEVICE     /dev/sdh        vol-266dc84e    2010-12-24T23:03:04.000Z
			BLOCKDEVICE     /dev/sdi        vol-fa6dc892    2010-12-24T23:05:50.000Z
			RESERVATION     r-6dfecb07      558277860346    @sc-masters,@sc-smallcluster
			INSTANCE        i-6f21d203      ami-a5c42dcc    ec2-72-44-59-38.compute-1.amazonaws.com ip-10-124-245-118.ec2.internal  running id_rsa-admin-key        0   m1.large 2011-01-06T14:11:09+0000        us-east-1a      aki-fd15f694    ari-7b739e12            monitoring-disabled     72.44.59.38     10.124.245.118      instance-store
			BLOCKDEVICE     /dev/sdj        vol-fc5af194    2011-01-06T14:14:11.000Z
			RESERVATION     r-63fecb09      558277860346    @sc-smallcluster
			INSTANCE        i-6921d205      ami-a5c42dcc    ec2-67-202-9-15.compute-1.amazonaws.com ip-10-124-247-224.ec2.internal  running id_rsa-admin-key        0   m1.large 2011-01-06T14:11:09+0000        us-east-1a      aki-fd15f694    ari-7b739e12            monitoring-disabled     67.202.9.15     10.124.247.224      instance-store

=cut

method getClusterNodes () {
	#### SET CONFIG FILE
	my $configfile = $self->configfile();

	print "Agua::StarCluster::getClusterNodes    StarCluster::getClusterNodes(username, cluster, privatekey, publiccert)\n";
	print "Agua::StarCluster::getClusterNodes    username: ", $self->username(). "\n";
	print "Agua::StarCluster::getClusterNodes    cluster: ", $self->cluster(). "\n";
	print "Agua::StarCluster::getClusterNodes    privatekey: ", $self->privatekey(). "\n";
	print "Agua::StarCluster::getClusterNodes    publiccert: ", $self->publiccert(). "\n";

	my $ec2din = qq{ec2-describe-instances \\
-K } . $self->privatekey() . qq{\\
-C } . $self->publiccert();	
	my $return = `$ec2din`;

	print "Agua::StarCluster::getClusterNodes    return from ec2din is undefined or empty\n" and exit if not defined $return or not $return;

	my $nodes = [];
	my @reservations = split "RESERVATION", $return;
	shift @reservations;
	foreach my $reservation ( @reservations )
	{
		my @instances = split "INSTANCE", $reservation;
		shift @instances; #### DISCARD RESERVATION LINE
		foreach my $instance ( @instances )
		{
			my @lines = split "\n", $instance;
			my $node;
			$node->{instance} = shift @lines;
			($node->{instanceid}) = $node->{instance} =~ /^\s+(\S+)/;
			$node->{blockdevices} = [];
			foreach my $line ( @lines )
			{
				push @{$node->{blockdevices}}, $line if $line =~ /^BLOCKDEVICE/;
			}
			push @$nodes, $node;
		}
	}

	#### DEBUG
	#foreach my $node ( @$nodes )
	#{
	#}

	return $nodes;
}


method addConfig ($section_name, $section_value, $key, $value, $comment) {
	print "Agua::StarCluster::addConfig    section_name not defined\n" and exit if not defined $section_name;
	print "Agua::StarCluster::addConfig    key not defined\n" and exit if not defined $key;
	print "Agua::StarCluster::addConfig    value not defined\n" and exit if not defined $value;


	my $config = $self->config();

	$config->addConfig($section_name, $section_value, $key, $value, $comment);
}


method removeConfig ($section_name, $section_value, $key, $value) {
	print "Agua::StarCluster::removeConfig    section_name not defined\n" and exit if not defined $section_name;
	print "Agua::StarCluster::removeConfig    key not defined\n" and exit if not defined $key;


	my $config = $self->config();

	$config->removeConfig($section_name, $section_value, $key, $value);
}



method getKeyValue ($name, $value, $key) {

	my $config = $self->config();

	return $config->getKeyValue($name, $value, $key);
}


=head2

	SUBROUTINE 		copyConfig

	PURPOSE

		GENERATE A USER-SPECIFIC STARCLUSTER CONFIG FILE:

			1. READ DEFAULT CONFIG FILE

			2. ADD USER AWS KEYS, KEYPAIR FILES

			3. PRINT TO NEW CONFIG FILE


			AWS_ACCESS_KEY_ID=AKIAIZXZ6S7ARZ44TTHQ
			AWS_SECRET_ACCESS_KEY=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
			AWS_USER_ID=728213020069

			[key starcluster-1]
			# Section name should match KEYNAME
			KEY_LOCATION=/agua/home/admin/.keypairs/id_rsa-starcluster-1


=cut
method copyConfig () {
	print "Agua::StarCluster::copyConfig    StarCluster::copyConfig(username)\n";

	#### GET config INSTANCE	
	my $conf = $self->conf();
	my $config = $self->config();

	#### SET OUTPUT CONFIG FILE
	my $userdir = $conf->getKeyValue('agua', "USERDIR");
	my $datadir = $conf->getKeyValue('agua', "DATADIR");
	my $configfile = "$datadir/starcluster/config";
	my $outputfile = "$userdir/" . $self->username() . "/.starcluster/config";
	print "Agua::StarCluster::copyConfig    configfile: $configfile\n";
	print "Agua::StarCluster::copyConfig    outputfile: $outputfile\n";

	#### COPY
	$config->inputfile($configfile);
	$config->copy($outputfile);
}


=head2

	SUBROUTINE 		stop

	PURPOSE

		stop UP A STARCLUSTER cluster FOR THE GIVEN USER

=cut
method stop () {
	print "Agua::StarCluster::stop    StarCluster::stop(privatekey, publiccert, username, cluster, keyname)\n";

	print "Agua::StarCluster::stop    privatekey: ", $self->privatekey() , "\n";
	print "Agua::StarCluster::stop    publiccert: ", $self->publiccert() , "\n";
	print "Agua::StarCluster::stop    username: ", $self->username() , "\n";
	print "Agua::StarCluster::stop    cluster: ", $self->cluster() , "\n";
	print "Agua::StarCluster::stop    keyname: ", $self->keyname() , "\n";

	#### PRINT HELP
	if ( defined $self->help() )
	{
		print qq{

	$0 stop <--username String> <--cluster String> [--help]

	--username 		User account name
	--cluster		Name of cluster to run (e.g., 'smallcluster')
	--privatekey	Location of private key file
	--publiccert	Location of public certificate file
	--keyname		Name of key (e.g., 'admin-key')
	--help			Print help info
};
	}

	#### MOUNT SHARES FROM AQUARIUS ON MASTER AND NODES
	$self->removeNfsMounts($self->privatekey(), $self->publiccert(), $self->username(), $self->cluster(), $self->keyname());

	### stop THE CLUSTER
	my $command = "starcluster -c ". $self->configfile() . " stop ". $self->cluster();
	print "Agua::StarCluster::stop    command: $command\n";
	my $success = `$command`;
	print "Agua::StarCluster::stop    success: $success\n";		
}


=head2

	SUBROUTINE 		removeNfsMounts

	PURPOSE

		REMOVE NFS MOUNT INFORMATION FROM SYSTEM FILES

		AND RESTART NFS DAEMONS

=cut

method removeNfsMounts () {
	print "Agua::StarCluster::removeNfsMounts    StarCluster::removeNfsMounts(username, cluster, keyname)\n";

	print "Agua::StarCluster::removeNfsMounts    privatekey: " , $self->privatekey() , "\n";
	print "Agua::StarCluster::removeNfsMounts    publiccert: " , $self->publiccert() , "\n";
	print "Agua::StarCluster::removeNfsMounts    username: " , $self->username() , "\n";
	print "Agua::StarCluster::removeNfsMounts    cluster: " , $self->cluster() , "\n";
	print "Agua::StarCluster::removeNfsMounts    keyname: " , $self->keyname() , "\n";

	#### CHECK INPUTS
	print "Agua::StarCluster::removeNfsMounts    privatekey not defined\n"and exit if not defined $self->privatekey();
	print "Agua::StarCluster::removeNfsMounts    publiccert not defined\n"and exit if not defined $self->publiccert();
	print "Agua::StarCluster::removeNfsMounts    username not defined\n"and exit if not defined $self->username();
	print "Agua::StarCluster::removeNfsMounts    cluster not defined\n"and exit if not defined $self->cluster();
	print "Agua::StarCluster::removeNfsMounts    keyname not defined\n"and exit if not defined $self->keyname();

	#### SET DEFAULT KEYNAME
	my $keypairfile = $self->keypairfile();
	print "Agua::StarCluster::removeNfsMounts    keypairfile: $keypairfile\n";

	#### GET INTERNAL IPS OF ALL NODES IN CLUSTER
	my $nodeips = $self->getInternalIps($self->username(), $self->cluster(), $self->privatekey(), $self->publiccert());
	#my $nodeips = [ "ip-10-124-241-66.ec2.internal" ];
	print "Agua::StarCluster::removeNfsMounts    nodeips: @$nodeips\n";

	#### REMOVE ENTRIES FROM /etc/exports
	my $volumes = [ "/agua", "/data", "/nethome" ];
	$self->removeExports($self->username(), $volumes, $nodeips);
}



=head2

	SUBROUTINE		removeExports

	PURPOSE

		ON MASTER, SET UP EXPORT TO HEAD INSTANCE IN /etc/exports:

		/home ip-10-124-247-224.ec2.internal(async,no_root_squash,no_subtree_check,rw)
		/opt/sge6 ip-10-124-247-224.ec2.internal(async,no_root_squash,no_subtree_check,rw)
		/data ip-10-124-247-224.ec2.internal(async,no_root_squash,no_subtree_check,rw)
	*** /data ip-10-127-158-202.ec2.internal(async,no_root_squash,no_subtree_check,rw)

=cut

method removeExports ($volumes, $recipientips, $targetip) {
	#### sourceip IS THE HOST DOING THE SHARING
	#### targetip IS THE HOST MOUNTING THE SHARE
	print "Agua::StarCluster::removeExports    StarCluster::removeExports(volumes, recipientips, targetip)\n";
	print "Agua::StarCluster::removeExports    recipientips: @$recipientips\n";
	print "Agua::StarCluster::removeExports    volume: @$volumes\n";
	print "Agua::StarCluster::removeExports    targetip: $targetip\n" if defined $targetip;

	#### SET CONFIG FILE
	my $configfile = $self->configfile();
	my $keypairfile = $self->keypairfile();

	#### SET REMOTE SSH 
	my $remotessh = $self->setRemoteSsh($targetip);

	#### GET CONTENTS OF /etc/exports
	my $exportsfile = "/etc/exports";
	my $cat = qq{$remotessh cat $exportsfile}; 
	my $exports = `$cat`;
	$exports =~ s/\s+$//;

	#### REMOVE EXISTING ENTRY FOR THESE VOLUMES
	my @lines = split "\n", $exports;
	foreach my $volume ( @$volumes )
	{
		foreach my $recipientip ( @$recipientips )
		{
			for ( my $i = 0; $i < $#lines + 1; $i++ )
			{
				if ( $lines[$i] =~ /^$volume\s+$recipientip/ )
				{
					splice @lines, $i, 1;
					$i--;
				}
			}
		}
	}
	my $output = join "\n", @lines;

	my $move = qq{$remotessh mv -f $exportsfile $exportsfile.bkp}; 
	`$move`;	

	#### WRITE TEMP FILE TO USER-OWNED /tmp DIRECTORY
	#### IN PREPARATION FOR COPY AS ROOT TO REMOTE HOST
	my $tempdir = "/tmp/" . $self->username();
	File::Path::mkpath($tempdir) if not -d $tempdir;
	my $tempfile = "/$tempdir/exports";
	open(OUT, ">$tempfile") or die "Can't open tempfile: $tempfile\n";
	print OUT $output;
	close(OUT) or die "Can't close tempfile: $tempfile\n";

	my $copy;
	$copy = qq{scp -i $keypairfile $tempfile root\@$targetip:$exportsfile} if defined $keypairfile;
	$copy = qq{cp $tempfile $exportsfile} if not defined $keypairfile;

	my $result = `$copy`;	

	##### DEBUG
	#my $catnew = qq{$remotessh cat $exportsfile}; 
}



#### DUMPER
method dump () { 
	print "Doing StarCluster::dump(@_)\n";
    require Data::Dumper;
    $Data::Dumper::Maxdepth = shift if @_;
    print Data::Dumper::Dumper $self;
}





}	#### class Agua::StarCluster

=head2


	SUBROUTINE		getConfigfile

	PURPOSE

		SET SSH COMMAND IF KEYPAIRFILE, ETC. ARE DEFINED

method getConfigfile () {
	print "Agua::StarCluster:::getConfigfile     username not defined: $self->username()\n" and exit if not defined $self->username();	

	print "Agua::StarCluster:::getConfigfile     Agua::StarCluster:::getConfigfile()\n";	
	return $self->configfile() if $self->configfile();

	my $conf 			= 	$self->conf();
	my $userdir 		= 	$conf->getKeyValue('agua', "USERDIR");
	print "Agua::StarCluster:::getConfigfile     userdir not defined: $userdir\n" and exit if not defined $userdir;	

	return "$userdir/" . $self->username() . "/.starcluster/config";
}

=head2

	SUBROUTINE		addVolume

	PURPOSE

		ADD A VOLUME AND ITS MOUNT POINT TO THE aws TABLE

#=cut

method addVolume () {
	my $json		=	$self->json();
	print "Agua::StarCluster::addVolume    json: $json\n";
	print Dumper $json;	

	my $self->username()	=	$json->{username};
	my $mountpoint	=	$json->{mountpoint};
	my $device		=	$json->{device};
	my $volume		=	$json->{volume};
	my $snapshot	=	$json->{snapshot};

	#### GET EXISTING KEYS IF AVAILABLE
	my $keys = $self->getKeys($self->username());	
	print "Agua::StarCluster::addVolume    keys:\n";
	print Dumper $keys;

	my $self->accesskeyid() = $keys->{awsaccesskeyid};
	my $self->secretaccesskey() = $keys->{secretaccesskey};
	print "Agua::StarCluster::addVolume    accesskeyid: ". $self->accesskeyid()\n";
	print "Agua::StarCluster::addVolume    secret: ". $self->accesskeyid()\n";

	my $ec2 = Net::Amazon::EC2->new(
		AWSAccessKeyId	=>	$self->accesskeyid(), 
		SecretAccessKey	=>	$self->secretaccesskey(),
		debug => 1
	);




}





=head2

	SUBROUTINE		runInstance

	PURPOSE

		START AN INSTANCE FOR THE GIVEN USER

#=cut

method runInstance () {
	my $self->accesskeyid()		=	shift;	
	my $self->secretaccesskey()	=	shift;	


	my $ec2 = Net::Amazon::EC2->new(
		AWSAccessKeyId	=>	$self->accesskeyid(), 
		SecretAccessKey	=>	$self->secretaccesskey()
	);


	# Start 1 new instance from AMI: ami-XXXXXXXX
	my $instance = $ec2->runInstances(
		InstanceType => 'm1.large',
		ImageId => 'ami-0a59bb63',
		MinCount => 1,
		MaxCount => 1
	);

	#### 32-bit
	#ImageId => 'ami-0859bb61',

	print "instance: \n";
	print Dumper $instance;


	my $running_instances = $ec2->describe_instances;

	foreach my $reservation (@$running_instances) {
	   foreach my $instance ($reservation->instances_set) {
		   print "instance->instance_id: " , $instance->instance_id . "\n";
	   }
	}

	my $instance_id = $instance->instances_set->[0]->instance_id;
	print "instance_id: $instance_id\n";

	return $instance_id;
}	


method terminateInstance () {
	my $self->accesskeyid()		=	shift;	
	my $self->secretaccesskey()	=	shift;	

	my $ec2 = Net::Amazon::EC2->new(
		AWSAccessKeyId => $self->accesskeyid(), 
		SecretAccessKey => $self->secretaccesskey()
	);

	##### Terminate instance
	#my $result = $ec2->terminateInstances(InstanceId => $instance_id);

	#return $result;
}



TAKEN FROM writeConfigfile

	######## MASTER TO AQUARIUS MOUNT
	######## CREATE VOLUME
	####my $snapshot = $conf->getKeyValue("starcluster:data", "SNAPSHOT");
	####my $size = $conf->getKeyValue("starcluster:data", "VOLUMESIZE");
	####my $availzone = $conf->getKeyValue("aws", "AVAILZONE");
	####my $volumeid = $self->_createVolume($self->privatekey(), $self->publiccert(), $snapshot, $availzone, $size);
	####
	######## SET [vol XXX] INFORMATION
	####my $datadir = $conf->getKeyValue('agua', "DATADIR");
	####my $datadevice = $conf->getKeyValue("aws", "STARCLUSTERDATADEVICE");
	####my $datapartition = "1";
	####$config->addKeyValue("volume", "data", "MOUNT_PATH", $datadir);
	####$config->addKeyValue("volume", "data", "DEVICE", $datadevice);
	####$config->addKeyValue("volume", "data", "PARTITION", $datapartition);
	####$config->addKeyValue("volume", "data", "VOLUME_ID", $volumeid);
	####$config->addKeyValue("cluster", "smallcluster", "VOLUMES", "data");

exit;

