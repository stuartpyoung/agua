package Agua::Common::Cluster;
use Moose::Role;
#use Moose::Util::TypeConstraints;

=head2

	PACKAGE		Agua::Common::Cluster

	PURPOSE

		CLUSTER METHODS FOR Agua::Common

=cut


use Data::Dumper;
use File::Path;
use Conf::StarCluster;

sub clusterStatus {
	my $self		=	shift;
	my $username 	=	shift;
	my $cluster 	=	shift;

 	print "Agua::Common::Cluster::getClusterStatus    username not defined\n" and exit if not defined $username;
 	print "Agua::Common::Cluster::getClusterStatus    cluster not defined\n" and exit if not defined $cluster;

	my $query = qq{SELECT *
FROM clusterstatus
WHERE username='$username'
AND cluster='$cluster'};
	print "$query\n";
	return $self->dbobject()->queryhash($query);
}

sub getAdminKey {
	my $self		=	shift;
	my $username	=	shift;

 	print "Agua::Common::Cluster::getAdminKey    username not defined\n" and exit if not defined $username;

	return $self->adminkey() if $self->can('adminkey') and defined $self->adminkey();

	my $adminkey_names = $self->conf()->getKeyValue('aws', 'ADMINKEY');
	$adminkey_names = '' if not defined $adminkey_names;
	my @names = split ",", $adminkey_names;
	my $adminkey = 0;
	foreach my $name ( @names )
	{
		if ( $name eq $username )	{	return $adminkey = 1;	}
	}

	$self->adminkey($adminkey) if $self->can('adminkey');

	return $adminkey;
}

sub getConfigfile {
	my $self		=	shift;
	my $username	=	shift;
	my $cluster		=	shift;
	my $lender		=	shift;

	my $userdir = $self->conf()->getKeyValue('agua', 'USERDIR');
	my $configfile;
	$configfile = "$userdir/$username/.starcluster/$cluster.config" if not defined $lender;
	$configfile = "$userdir/$lender/.starcluster/$cluster.config" if defined $lender;

	return $configfile;
}

=head2

	SUBROUTINE		updateCluster

	PURPOSE

		SAVE A CLUSTER TO THE cluster TABLE

=cut

sub updateCluster {
	my $self		=	shift;


    my $dbobject        =	$self->dbobject();
    my $json 			=	$self->json();


	my $query = qq{SELECT * FROM clusters
WHERE username='$json->{username}'
AND cluster='$json->{cluster}'
};
	my $cluster = $dbobject->query($query);
	print Dumper $cluster;

	$cluster->{minnodes} = $json->{minnodes};
	$cluster->{maxnodes} = $json->{maxnodes};

	my $success = $self->_removeCluster();	
 	print "{ error: 'Agua::Common::Cluster::updateCluster    Could not delete cluster $json->{cluster} from cluster table'}" and exit if not defined $success;

	#### SET TABLE AND REQUIRED FIELDS	
	my $table = "cluster";
	my $required_fields = ["username", "cluster", "minnodes", "maxnodes", "instancetype", "amiid"];
	#### CHECK REQUIRED FIELDS ARE DEFINED
	my $not_defined = $dbobject->notDefined($cluster, $required_fields);
    print "{ error: 'Agua::Common::Cluster::_addCluster    not defined: @$not_defined' }" and exit if @$not_defined;

	#### DO ADD
	$success = $self->_addToTable($table, $cluster, $required_fields);	
	print "{ status: 'Agua::Common::Cluster::updateCluster    Could not add cluster $json->{cluster} into cluster table' }" if not $success;
	print "{ status: 'Agua::Common::Cluster::updateCluster    Successful insert of cluster $json->{cluster} into cluster table' }" if $success;
	exit;
}

=head2

	SUBROUTINE		addCluster

	PURPOSE

		ADD A CLUSTER TO THE cluster TABLE

		IF THE USE HAS NO AWS CREDENTIALS INFORMATION,

		USE THE 'admin' USER'S AWS CREDENTIALS AND STORE

		THE CONFIGFILE IN THE ADMIN USER'S .starcluster

		DIRECTORY
=cut

sub addCluster {
	my $self		=	shift;


    my $dbobject    =	$self->dbobject();
    my $json 		=	$self->json();
    my $conf 		=	$self->conf();
	my $username 	= 	$json->{username};	
	my $cluster 	= 	$json->{cluster};	
	my $userdir 	= 	$conf->getKeyValue('agua', 'USERDIR');


	#### GET AWS
	my $aws = $self->getAws($username);	

	$json->{amazonuserid} = $aws->{amazonuserid};
	$json->{accesskeyid} = $aws->{awsaccesskeyid};
	$json->{secretaccesskey} = $aws->{secretaccesskey};

	##### ADD TO cluster TABLE
	#$self->_removeCluster();	
	#my $success = $self->_addCluster();

	#### SET PRIVATE KEY AND PUBLIC CERT FILE LOCATIONS	
	my $adminkey = $self->getAdminKey($username);
	$json->{privatekey} = $self->getPrivatefile() if not $adminkey;
	$json->{publiccert} = $self->getPublicfile() if not $adminkey;
	$json->{privatekey} = $self->getPrivatefile('admin') if $adminkey;
	$json->{publiccert} = $self->getPublicfile('admin') if $adminkey;

	#### SET KEYPAIR FILE IF NOT DEFINED
	$json->{keypairfile} = "$userdir/$username/.starcluster/id_rsa-$username-key" if not $adminkey;
	$json->{keypairfile} = "$userdir/admin/.starcluster/id_rsa-admin-key" if $adminkey;

	#### SET KEYNAME
	$json->{keyname} = "$username-key";
	$json->{keyname} = "admin-key" if $adminkey;

	#### SET NFS MOUNT INFO
	$json->{sources} = $conf->getKeyValue('starcluster:mounts', 'SOURCEDIRS')
		if not defined $json->{sources};
	$json->{mounts} = $conf->getKeyValue('starcluster:mounts', 'MOUNTPOINTS')
		if not defined $json->{mounts};
	$json->{devs} = $conf->getKeyValue('starcluster:mounts', 'DEVICES')
		if not defined $json->{devs};

	#### SET CONFIGFILE
	$json->{configfile} = $self->getConfigfile($username, $cluster) if not $adminkey;
	$json->{configfile} = $self->getConfigfile($username, $cluster, 'admin') if $adminkey;

	#### SET NODES TO GREATER OF 1 OR minnodes
	#### IF maxnodes IS GREATER THAN minnodes, THE LOAD BALANCER WILL 
	#### MAINTAIN THE NUMBER OF NODES BETWEEN minnodes AND maxnodes
	#### DEPENDING ON THE JOB LOAD ONCE THE CLUSTER HAS STARTED
	$json->{nodes} = 1;
	$json->{nodes} = $json->{minnodes} if $json->{minnodes} > $json->{nodes};

	#### SET NODE IMAGE ID
	$json->{nodeimage} = $json->{amiid};

	#### ADD CONF
	$json->{conf} = $self->conf();

	#### DEBUG

	#### CREATE STARCLUSTER config FILE
 	print "Agua::Common::Cluster::addCluster    Doing Agua::StarCluster->new(json)\n";
	my $starcluster = Agua::StarCluster->new($json);
	$starcluster->writeConfigfile();
	exit;
}


=head2

	SUBROUTINE		saveCluster

	PURPOSE

		SAVE A CLUSTER TO THE cluster TABLE

=cut

sub saveCluster {
	my $self		=	shift;

    my $dbobject    =	$self->dbobject();
    my $json 		=	$self->json();


	my $success = $self->_removeCluster();	
	$success = $self->_addCluster();
	print "{ status: 'Agua::Common::Cluster::saveCluster    Could not add cluster $json->{cluster} into cluster table' }" if not $success;
	print "{ status: 'Agua::Common::Cluster::saveCluster    Successful insert of cluster $json->{cluster} into cluster table' }" if $success;
	exit;
}

=head2

	SUBROUTINE		_addCluster

	PURPOSE

		INTERNAL USE ONLY: ADD A CLUSTER TO THE cluster TABLE

=cut

sub _addCluster {
	my $self		=	shift;


    my $dbobject        =	$self->dbobject();
    my $json 			=	$self->json();


    #### VALIDATE
    print "{ error: 'Agua::Common::Cluster::_addCluster    User session not validated' }" and exit unless $self->validate();

	#### SET TABLE AND REQUIRED FIELDS	
	my $table = "cluster";
	my $required_fields = ["username", "cluster", "minnodes", "maxnodes", "instancetype", "amiid"];

	#### CHECK REQUIRED FIELDS ARE DEFINED
	my $not_defined = $dbobject->notDefined($json, $required_fields);
    print "{ error: 'Agua::Common::Cluster::_addCluster    not defined: @$not_defined' }" and exit if @$not_defined;

	#### DO ADD
	my $success = $self->_addToTable($table, $json, $required_fields);	

#### DEBUG COMMENTED OUT

	#### ADD THE NEW CELL DIRECTORY TO $SGE_ROOT
	my $celldir = $self->getCelldir();
	my $defaultdir = $self->getDefaultCelldir();
	if ( not -d $celldir )
	{
		my $copy = "cp -pr $defaultdir $celldir";
		print `$copy`;

		my $chown = "chown -R sgeadmin:sgeadmin $celldir";
		print `$chown`;
	}
	$success = -d $celldir;

	my $cluster = $json->{cluster};
	print "{ error: 'Agua::Common::Cluster::_addCluster    Could not copy cluster $cluster directory' }" and exit if not defined $success or not $success;

	print "{ error: 'Agua::Common::Cluster::_addCluster    Successfully added cluster $cluster' }" if not defined $success or not $success;
}

sub getCelldir {
	my $self 		=	shift;

	my $json		=	$self->json();	
	my $cluster = $json->{cluster};
	my $username = $json->{username};
	my $sgeroot = $self->conf()->getKeyValue('cluster', 'SGEROOT');

	return "$sgeroot/$username-$cluster";	
}

sub getDefaultCelldir {
	my $self 		=	shift;

	my $json		=	$self->json();	
	my $sgeroot = $self->conf()->getKeyValue('cluster', 'SGEROOT');

	return "$sgeroot/default";
}

=head2

	SUBROUTINE		removeCluster

	PURPOSE

		REMOVE A CLUSTER FROM cluster TABLE

=cut

sub removeCluster  {
	my $self		=	shift;


    my $dbobject        =	$self->dbobject();
    my $json 			=	$self->json();

    #### VALIDATE
    print "{ error: 'Agua::Common::Cluster::removeCluster    User session not validated' }" and exit unless $self->validate();

	#### REMOVE FROM cluster
	my $success = $self->_removeCluster();
 	print "{ error: 'Agua::Common::Cluster::removeCluster    Could not remove cluster $json->{cluster}'}" and exit if not defined $success;
	$success = $self->_removeCelldir();
 	print "{ error: 'Agua::Common::Cluster::removeCluster    Could not remove celldir for cluster $json->{cluster}'}" and exit if not defined $success;


	#### REMOVE CLUSTER CONFIG FILE
	my $userdir = $self->conf()->getKeyValue('agua', 'USERDIR');
	my $configfile = "$userdir/" . $self->username() . "/.starcluster/config";
	`rm -fr $configfile`;
	print "{ error: 'Agua::Common::Cluster::_removeCluster    Can't remove cluster $json->{cluster} configfile: $configfile' }" and exit if -f $configfile;

 	print "{ error: 'Agua::Common::Cluster::removeCluster    Successfully removed cluster $json->{cluster}'}" if $success;
	exit;

}	#### removeCluster



=head2

	SUBROUTINE		_removeCluster

	PURPOSE

		REMOVE A CLUSTER FROM cluster TABLE

=cut

sub _removeCluster {
	my $self		=	shift;


    my $dbobject        =	$self->dbobject();
    my $json 			=	$self->json();

    #### VALIDATE
    print "{ error: 'Agua::Common::Cluster::_removeCluster    User session not validated' }" and exit unless $self->validate();

	#### SET CLUSTER, TABLE AND REQUIRED FIELDS	
	my $cluster = $json->{cluster};
	my $table = "cluster";
	my $required_fields = ["username", "cluster"];

	#### CHECK REQUIRED FIELDS ARE DEFINED
	my $not_defined = $dbobject->notDefined($json, $required_fields);
    print "{ error: 'Agua::Common::Cluster::_removeCluster    undefined values: @$not_defined' }" and exit if @$not_defined;

	#### REMOVE FROM cluster
	my $success = $self->_removeFromTable($table, $json, $required_fields);
 	print "{ error: 'Agua::Common::Cluster::_removeCluster    Could not remove cluster $cluster from cluster table'}" and exit if not defined $success;
}

sub _removeCelldir {
	my $self		=	shift;

	#### REMOVE THE CELL DIRECTORY
	my $celldir = $self->getCelldir();
	if ( -d $celldir )
	{
		my $remove = "rm -fr $celldir";
		print `$remove`;
	}
	my $success = not (-d $celldir);

	print "{ error: 'Agua::Common::Cluster::_removeCelldir    Could not delete celldir' }" and exit if not defined $success or not $success;

	return $success;
}	#### _removeCelldir

=head2

    SUBROUTINE:     getClusters

    PURPOSE:

		RETURN AN ARRAY OF cluster HASHES

=cut

sub getClusters {
	my $self		=	shift;


    my $dbobject		=	$self->dbobject();
	my $json			=	$self->json();


    #### VALIDATE
    my $username = $json->{'username'};
    print "{ error: 'User session not validated' }" and exit unless $self->validate();

	#### GET ALL SOURCES
	my $query = qq{SELECT * FROM cluster
WHERE username='$username'
ORDER BY maxnodes ASC, cluster};
	my $clusters = $dbobject->queryhasharray($query);

	#### SET TO EMPTY ARRAY IF NO RESULTS
	$clusters = [] if not defined $clusters;

	return $clusters;
}

sub activeClusters {
	my $self		=	shift;


    #### GET DATABASE OBJECT
    my $dbobject     =	$self->dbobject();

	my $query = qq{SELECT * FROM clusterstatus
WHERE status
='running'};
	my $clusters = $dbobject->queryhasharray($query);

	#### SET TO EMPTY ARRAY IF NO RESULTS
	$clusters = [] if not defined $clusters;

	return $clusters;	
}

sub clusterInputs {	
	my $self		=	shift;
	my $object		=	shift;


	my $username 	=	$object->{username};
	my $cluster 	=	$object->{cluster};

	#### ADD CONF OBJECT	
	$object->{conf} = $self->conf();

	#### ADD min/maxnodes FROM cluster TABLE
	my $clusternodes = $self->getCluster($username, $cluster);
	$object->{minnodes} = $clusternodes->{minnodes};
	$object->{maxnodes} = $clusternodes->{maxnodes};

	#### DETERMINE WHETHER TO USE ADMIN KEY FILES
	my $adminkey = $self->getAdminKey($username);
	if ( $adminkey )
	{
		#### SET KEYNAME
		$object->{keyname} = 	"admin-key";

		#### IF USER HAS AWS CREDENTIALS, USE THEIR CONFIGFILE.
		$object->{configfile} =  $self->getConfigfile($username, $cluster, 'admin')
			if $adminkey;

		#### SET PRIVATE KEY AND PUBLIC CERT FILE LOCATIONS	
		$object->{privatekey} =  $self->getPrivatefile('admin');
		$object->{publiccert} =  $self->getPublicfile('admin');
	}
	else
	{
		#### SET KEYNAME
		$object->{keyname} =  "$username-key";

		#### QUIT IF USER HAS NO AWS CREDENTIALS
		my $aws = $self->getAws($username);
		print "{ error: 'Agua::Workflow::clusterInputs    No AWS credentials for user $username' }" and exit if not defined $aws;

		#### IF USER HAS AWS CREDENTIALS, USE THEIR CONFIGFILE.
		$object->{configfile} =  $self->getConfigfile($username, $cluster);

		#### SET PRIVATE KEY AND PUBLIC CERT FILE LOCATIONS	
		$object->{privatekey} =  $self->getPrivatefile($username);
		$object->{publiccert} =  $self->getPublicfile($username);
	}

	print "Agua::Common::Cluster::clusterInputs    configfile: $object->{configfile}\n";

	return $object;
}



=head2

    SUBROUTINE:     getCluster

    PURPOSE:

		THE cluster TABLE ENTRY FOR THE SPECIFIED USER'S TABLE

=cut

sub getCluster {
	my $self		=	shift;
	my $username	=	shift;
	my $cluster		=	shift;


	#### GET ALL SOURCES
	my $query = qq{SELECT * FROM cluster
WHERE username='$username'
AND cluster='$cluster'};

	return $self->dbobject()->queryhash($query);
}

=head2

    SUBROUTINE:     getClusterWorkflows

    PURPOSE:

		RETURN AN ARRAY OF cluster HASHES

=cut

sub getClusterWorkflows {
	my $self		=	shift;


    my $dbobject	=	$self->dbobject();
	my $username 	=	$self->username();

	#### GET ALL SOURCES
	my $query = qq{SELECT * FROM clusterworkflow
WHERE username='$username'
ORDER BY project, workflow};
	my $clusters = $dbobject->queryhasharray($query);

	#### SET TO EMPTY ARRAY IF NO RESULTS
	$clusters = [] if not defined $clusters;

	return $clusters;
}

=head2

    SUBROUTINE:     getClusterByWorkflow

    PURPOSE:

		RETURN THE CLUSTER FOR A GIVEN PROJECT WORKFLOW

=cut

sub getClusterByWorkflow {
	my $self		=	shift;
	my $username	=	shift;
	my $project		=	shift;
	my $workflow	=	shift;


	print "Agua::Common::Cluster::getClusterByWorkflow    username not defined\n" and exit if not defined $username;
	print "Agua::Common::Cluster::getClusterByWorkflow    project not defined\n" and exit if not defined $project;
	print "Agua::Common::Cluster::getClusterByWorkflow    workflow not defined\n" and exit if not defined $workflow;

    my $dbobject		=	$self->dbobject();

	#### GET ALL SOURCES
	my $query = qq{SELECT cluster FROM clusterworkflow
WHERE username='$username'
AND project='$project'
AND workflow='$workflow'};
	my $cluster = $dbobject->query($query);

	return $cluster;
}


=head2

	SUBROUTINE		saveClusterWorkflow

	PURPOSE

		ADD AN ENTRY OR UPDATE AN EXISTING ENTRY IN clusterworkflow TABLE
=cut

sub saveClusterWorkflow {
	my $self		=	shift;

    my $dbobject        =	$self->dbobject();
    my $json 			=	$self->json();


	my $success = $self->_removeClusterWorkflow();	

	$success = $self->_removeQueue();

	$success = $self->_addClusterWorkflow();
	print "{ status: 'Agua::Common::Cluster::saveClusterWorkflow    Could not add cluster $json->{cluster} into cluster table' }" if not $success;

	$success = $self->_addClusterQueue();

	print "{ status: 'Agua::Common::Cluster::saveClusterWorkflow    Successful insert of cluster $json->{cluster} into cluster table' }" if $success;
	exit;
}

=head2

	SUBROUTINE		_addClusterWorkflow

	PURPOSE

		INTERNAL USE ONLY: ADD A CLUSTER TO THE clusterworkflow TABLE

=cut

sub _addClusterWorkflow {
	my $self		=	shift;


    my $dbobject        =	$self->dbobject();
    my $json 			=	$self->json();


    #### VALIDATE
    print "{ error: 'Agua::Common::Cluster::_addClusterWorkflow    User session not validated' }" and exit unless $self->validate();

	#### SET TABLE AND REQUIRED FIELDS	
	my $table = "clusterworkflow";
	my $required_fields = ["username", "project", "workflow"];

	#### CHECK REQUIRED FIELDS ARE DEFINED
	my $not_defined = $dbobject->notDefined($json, $required_fields);
    print "{ error: 'Agua::Common::Cluster::_addClusterWorkflow    not defined: @$not_defined' }" and exit if @$not_defined;

	#### DO ADD
	my $success = $self->_addToTable($table, $json, $required_fields);	
	print "{ status: 'Agua::Common::Cluster::_addClusterWorkflow    Failed to insert cluster $json->{cluster} into cluster table' }" if not defined $success or not $success;

	print "{ status: 'Agua::Common::Cluster::_addClusterWorkflow    Successful insert of cluster $json->{cluster} into cluster table' }" if $success;
}

=head2

	SUBROUTINE		_removeClusterWorkflow

	PURPOSE

		REMOVE A CLUSTER FROM clusterworkflow TABLE

=cut

sub _removeClusterWorkflow {
	my $self		=	shift;

$DEBUG = 1;

    my $dbobject        =	$self->dbobject();
    my $json 			=	$self->json();

    #### VALIDATE
    print "{ error: 'Agua::Common::Cluster::_removeClusterWorkflow    User session not validated' }" and exit unless $self->validate();

	#### SET TABLE AND REQUIRED FIELDS	
	my $table = "clusterworkflow";
	my $required_fields = ["username", "project", "workflow"];

	#### CHECK REQUIRED FIELDS ARE DEFINED
	my $not_defined = $dbobject->notDefined($json, $required_fields);
    print "{ error: 'Agua::Common::Cluster::_removeClusterWorkflow    undefined values: @$not_defined' }" and exit if @$not_defined;

	#### REMOVE FROM cluster
	return $self->_removeFromTable($table, $json, $required_fields);

}	#### _removeClusterWorkflow


=head2

	SUBROUTINE		_addQueue

	PURPOSE

		INTERNAL USE ONLY: ADD A CLUSTER TO THE clusterworkflow TABLE

=cut

sub _addQueue {
	my $self		=	shift;

    my $dbobject        =	$self->dbobject();
    my $json 			=	$self->json();


    #### VALIDATE
    print "{ error: 'Agua::Common::Cluster::_addQueue    User session not validated' }" and exit unless $self->validate();

	my $queue = $json->{username} . "-" . $json->{project} . "-" . $json->{workflow};	
 	print "Agua::Common::Cluster::_addQueue    queue: $queue\n";

	my $queuefile = $self->getQueuefile($queue);
	print "Agua::Common::Cluster::_addQueuefile    queuefile: $queuefile\n";

	return $self->addQueue($queue, $queuefile);
}

=head2

	SUBROUTINE		_removeQueue

	PURPOSE

		REMOVE A CLUSTER FROM clusterworkflow TABLE

=cut

sub _removeQueue {
	my $self		=	shift;

$DEBUG = 1;

    my $dbobject        =	$self->dbobject();
    my $json 			=	$self->json();

    #### VALIDATE
    print "{ error: 'Agua::Common::Cluster::_removeQueue    User session not validated' }" and exit unless $self->validate();

	my $queue = $json->{username} . "-" . $json->{project} . "-" . $json->{workflow};	
 	print "Agua::Common::Cluster::_removeQueue    queue: $queue\n";

	my $monitor = $self->getMonitor();
 	print "Agua::Common::Cluster::_removeQueue    monitor: $monitor\n";

	my $queuefile = $self->getQueuefile($queue);
	print "Agua::Common::Cluster::_removeQueue    Can't remove queuefile: $queuefile\n" if 

	return $monitor->removeQueue($queue, $queuefile);

}	#### _removeQueue



sub getMonitor {
	my $self		=	shift;
	my $clustertype	=	shift;
	print "Agua::Common::Cluster::getMonitor    (clustertype)\n";
	print "Agua::Common::Cluster::getMonitor    clustertype: $clustertype\n";

	return $self->monitor() if $self->monitor();

	$clustertype =  $self->conf()->getKeyValue('agua', 'CLUSTERTYPE') if not defined $clustertype;
	my $classfile = "Agua/Monitor/" . uc($clustertype) . ".pm";
	my $module = "Agua::Monitor::$clustertype";
	require $classfile;

	my $monitor = $module->new(
		{
			'pid'		=>	$self->workflowpid(),
			'conf' 		=>	$self->conf(),
			'dbobject'	=>	$self->dbobject()
		}
	);
	$self->monitor($monitor);

	return $monitor;
}



1;