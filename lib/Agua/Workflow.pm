#use Moose;
use MooseX::Declare;


=head2

	PACKAGE		Workflow

	PURPOSE

		THE Workflow OBJECT PERFORMS THE FOLLOWING TASKS:

			1. SAVE WORKFLOWS

			2. RUN WORKFLOWS

			3. PROVIDE WORKFLOW STATUS

	NOTES

		Workflow::executeWorkflow
			|
			|
			|
			|
		Workflow::runStages
				|
				|
				|
				-> 	my $stage = Agua::Stage->new()
					...
					|
					|
					-> $stage->run() 
						|
						|
						? DEFINED 'CLUSTER' AND 'SUBMIT'
						|
						|
						YES ->  Agua::Stage::clusterSubmit() 
						|
						|
						NO ->  Agua::Stage::execute()

	EXAMPLES

runWorkflow

perl workflow.cgi < workflow-runWorkflow.json
{
    'sessionId' : '1228791394.7868.158',
    'username' : 'admin',
    'project' : 'Project1',
    'workflow' : 'Workflow3-indels',
    'mode' : 'runWorkflow',
    'start' : 0
}

executeWorkflow

perl ./workflow.cgi < workflow-executeWorkflow.json

"sessionId=1228791394.7868.158&username=admin&project=Project1&workflow=Workflow3-indels&mode=executeWorkflow&start=0"


saveWorkflow

perl workflow.cgi < workflow-saveWorkflow.json

"sessionId=1228791394.7868.158&username=admin&mode=saveWorkflow&project=Project1&workflow=Workflow3-indels"


runStatus

perl workflow.cgi < workflow-runStatus.json
print "\n";
"sessionId=1228791394.7868.158&username=admin&mode=runStatus&project=Project1&workflow=Workflow3-indels"


workflowApplications

perl workflow.cgi < workflow-workflowApplications.json

"sessionId=1228791394.7868.158&username=admin&mode=workflowApplications&project=Project1&workflow=Workflow3-indels"

"sessionId=1228791394.7868.158&username=admin&mode=workflowApplications&workflow=Workflow3-indels&project=Project1"


=cut

use strict;
use warnings;
use Carp;

##### USE LIBS
#use FindBin qw($Bin);
#use lib "$Bin/../";
#use lib "$Bin/../external";

class Agua::Workflow extends Agua::Common {

#### EXTERNAL MODULES
use Data::Dumper;

use FindBin qw($Bin);
use lib "$Bin/..";

#### INTERNAL MODULES	
use Agua::DBaseFactory;
use Conf::Agua;
use Agua::Stage;

#use Agua::Monitor;

# INTS
has 'workflowpid'	=> ( isa => 'Int|Undef', is => 'rw', required => 0 );
has 'workflownumber'=>  ( isa => 'Str', is => 'rw' );
has 'start'     	=>  ( isa => 'Int', is => 'rw' );
has 'submit'  		=>  ( isa => 'Int', is => 'rw' );

# STRINGS
has 'fileroot'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'qstat'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );

# STRINGS
has 'queue'			=>  ( isa => 'Str|Undef', is => 'rw', default => 'default' );
has 'cluster'		=>  ( isa => 'Str|Undef', is => 'rw' );
has 'username'  	=>  ( isa => 'Str', is => 'rw' );
has 'workflow'  	=>  ( isa => 'Str', is => 'rw' );
has 'project'   	=>  ( isa => 'Str', is => 'rw' );
#has 'queue_options'=>  ( isa => 'Str|Undef', is => 'rw', default => '' );
#has 'outputdir'	=>  ( isa => 'Str', is => 'rw' );
#has 'setuid'		=>  ( isa => 'Str|Undef', is => 'rw', default => '' );
#has 'scriptfile'	=>  ( isa => 'Str', is => 'rw' );
#has 'installdir'   =>  ( isa => 'Str', is => 'rw' );
#has 'qsub'			=>  ( isa => 'Str', is => 'rw' );
#has 'qstat'		=>  ( isa => 'Str', is => 'rw' );

# OBJECTS
has 'json'		=> ( isa => 'HashRef', is => 'rw', required => 0 );
has 'dbobject'	=> ( isa => 'Agua::DBase::MySQL', is => 'rw', required => 0 );
has 'stages'		=> 	( isa => 'ArrayRef', is => 'rw', required => 0 );
has 'stageobjects'	=> 	( isa => 'ArrayRef', is => 'rw', required => 0 );
has 'monitor'		=> 	( isa => 'Maybe|Undef', is => 'rw', required => 0 );

has 'conf' 	=> (
	is =>	'rw',
	'isa' => 'Conf::Agua',
	default	=>	sub { Conf::Agua->new(	backup	=>	1, separator => "\t"	);	}
);

####////}	

method BUILD ($hash) {
#exit;

	#### SET DATABASE HANDLE
	$self->setDbh();

    #### VALIDATE
    print "{ error: 'Agua::Workflow::BUILD    User session not validated' }" and exit unless $self->validate();

	#### SET WORKFLOW PROCESS ID
	$self->workflowpid($$);	

	#### SET CLUSTER IF DEFINED
	print "Agua::Workflow::BUILD    conf->getKeyValue(agua, CLUSTERTYPE) not defined\n"
		and exit if not defined $self->conf()->getKeyValue('agua', 'CLUSTERTYPE');    
	print "Agua::Workflow::BUILD    conf->getKeyValue(cluster, QSUB) not defined\n"
		and exit if not defined $self->conf()->getKeyValue('cluster', 'QSUB');
	print "Agua::Workflow::BUILD    conf->getKeyValue(cluster, QSTAT) not defined\n"
		and exit if not defined $self->conf()->getKeyValue('cluster', 'QSTAT');

	#### IF JSON IS DEFINED, ADD VALUES TO SLOTS
	if ( $self->json() )
	{
		foreach my $key ( keys %{$self->{json}} )
		{
			$self->$key($self->{json}->{$key}) if $self->can($key);
		}
	}
}


=head2

	SUBROUTINE		checkBalancers

	PURPOSE

		MONITOR THE LOAD BALANCERS RUNNING ON EACH CLUSTER (I.E., SGE CELL)

		AND RESTART THEM IF THEY BECOME INACTIVE (E.G., DUE TO A RECURRENT

		ERROR IN STARCLUSTER WHEN PARSING XML OUTPUT AND/OR XML ERROR MESSAGES

		CAUSED BY A BUG IN SGE)

	NOTES    

		1. CRON JOB RUNS THIS SCRIPT EVERY MINUTE

		* * * * * /agua/0.6/bin/scripts/checkBalancers.pl --logfile /tmp/checkbalancers.out

		2. CHECKS clusterstatus DATABASE TABLE FOR RUNNING CLUSTERS

		3. FOR RUNNING CLUSTERS, CHECKS IF LOCK FILE IS OLDER THAN 1 MINUTE

		4. IF LOCK FILE STALE, RESTART bal FOR THE CLUSTER:

			-   KILL PREVIOUS PROCESS USING PID IN PID FILE

			-   START bal: STORE PID IN PID FILE, CONCAT OUTPUT TO OUTPUT FILE

=cut

method checkBalancers {

    #### GET DATABASE OBJECT
    my $dbobject     =	$self->dbobject();
    my $json         =	$self->json();

	#### GET ALL status='running' CLUSTERS
	my $clusterobjects = $self->activeClusters();

	#### KEEP ONLY CLUSTERS WHERE THE LOAD BALANCER HAS FAILED	
	for ( my $i = 0; $i < @$clusterobjects; $i++ )
	{
		my $clusterobject = $$clusterobjects[$i];
		my $pid = $clusterobject->{balancerpid};
		print "Agua::Workflow::checkBalancers    No pid for clusterobject $clusterobject->{clusterobject}\n" and next if not defined $pid or not $pid;
		my $running = kill 0, $pid;
		print "Agua::Workflow::checkBalancers    running: $running\n";
		if ( $running )
		{
			splice @$clusterobjects, $i, 1;
			$i--;
		}
	}

print "HERE\n";
exit;

	#### RESTART LOAD BALANCERS FOR REMAINING CLUSTERS
	foreach my $clusterobject (@$clusterobjects)
	{
		$self->restartLoadBalancer($clusterobject);
	}
}

method restartLoadBalancer ($clusterobject) {

    #### GET DATABASE OBJECT
    my $dbobject     =	$self->dbobject();

	my $username 	=	$clusterobject->{username};
	my $cluster 	=	$clusterobject->{cluster};

	$clusterobject = $self->getStarClusterInputs($clusterobject);

	#### RUN starcluster.pl TO GENERATE KEYPAIR FILE IN .starcluster DIR
	#### NB: StarCluster WILL SET STATUS IN clusterstatus TABLE
	my $starcluster = Agua::StarCluster->new($clusterobject);
	print "Agua::Workflow::startCluster    Doing StarCluster->start(clusterobject)\n";
	$starcluster->launchBalancer();
}

method getStarClusterInputs ($object) {

	my $username 	=	$object->{username};
	my $cluster 	=	$object->{cluster};

	#### DETERMINE WHETHER TO USE ADMIN KEY FILES
	my $adminkey = $self->getAdminKey($username);

	#### SET KEYNAME
	$object->{keyname} =  "$username-key";
	$object->{keyname} = 	"admin-key" if $adminkey;

	#### QUIT IF USER HAS NO AWS CREDENTIALS
	my $aws = $self->getAws($username);
	print "{ error: 'Agua::Workflow::getStarClusterInputs    No AWS credentials for user $username' }" and exit if not defined $aws;

	#### IF USER HAS AWS CREDENTIALS, USE THEIR CONFIGFILE.
	print "Agua::Workflow::getStarClusterInputs    adminkey: $adminkey\n";
	$object->{configfile} =  $self->getConfigfile($username, $cluster);
	$object->{configfile} =  $self->getConfigfile($username, $cluster, 'admin')
		if $self->getAdminKey($username);
	print "Agua::Workflow::getStarClusterInputs    configfile: $object->{configfile}\n";

	#### SET PRIVATE KEY AND PUBLIC CERT FILE LOCATIONS	
	$object->{privatekey} =  $self->getPrivatefile();
	$object->{publiccert} =  $self->getPublicfile();
	$object->{privatekey} =  $self->getPrivatefile('admin') if $adminkey;
	$object->{publiccert} =  $self->getPublicfile('admin') if $adminkey;

	#### GET CLUSTER NODES INFO FROM cluster TABLE
	my $clusternodes = $self->getCluster($username, $cluster);
	$object->{minnodes} = $clusternodes->{minnodes};
	$object->{maxnodes} = $clusternodes->{maxnodes};

	#### ADD CONF OBJECT	
	$object->{conf} = $self->conf();

	#### RENAME CLUSTER TO REFLECT USER
	$object->{cluster} = "$username-$cluster";


	return $object;
}


=head2

	SUBROUTINE		getStatus

	PURPOSE

		GET WORKFLOW STATUS

=cut

method getStatus {

    my $username	= 	$self->username();
    my $project 	= 	$self->project();
    my $workflow 	= 	$self->workflow();
    my $json  		=	$self->json();
    my $start 		= 	$json->{start};

	#### GET RUN STATUS FROM stage TABLE
    my $query = qq{SELECT *, NOW() AS now
FROM stage
WHERE username ='$username'
AND project = '$project'
AND workflow = '$workflow'
ORDER BY number
};
    my $stages = $self->dbobject()->queryhasharray($query);

	print "{ error: 'No stages with run status for username: $username, project: $project, workflow: $workflow from start: $start' }" and exit if not defined $stages;

	foreach my $stage ( @$stages )
	{
	}

	#### SET CLUSTER - RETRIEVE FROM clusterworkflow TABLE
	my $cluster = $self->getClusterByWorkflow($username, $project, $workflow);
	$self->cluster($cluster) if defined $cluster;

    #### IF CLUSTER IS NOT DEFINED THEN JOB WAS SUBMITTED LOCALLY 
	#### AND THE stage TABLE SHOULD BE UPDATED BY THE PROCESS ON EXIT.
	#### SO JUST PRINT STATUS AND QUIT
	if ( not defined $cluster)
	{
		my $status;
		$status->{stages} = $stages;
		require JSON;
		my $jsonObject = JSON->new();
		my $statusJson = $jsonObject->objToJson($status);
		$statusJson =~ s/"/'/g;    
		print "$statusJson\n";
		exit;    
	}

	#### GET CLUSTER STATUS (INCLUDES THIS USER'S PROJECT WORKFLOW)
	my $clusterstatus = $self->clusterStatus();

	#### GET QSTAT SUMMARY FOR THIS USER'S PROJECT WORKFLOW
	my $monitor = $self->getMonitor();
	my $queuestatus = $monitor->queueStatus();

	#### UPDATE stage TABLE WITH JOB STATUS FROM QSTAT
	#### BEFORE RETURNING THE STATUS OF STAGES 
	$self->updateStageStatus($monitor, $stages);

	my $status;
	$status->{stages} 	= $stages;
	$status->{cluster} 	= $clusterstatus;
	$status->{queue} 	= $queuestatus;

    require JSON;
    my $jsonObject = JSON->new();
    my $statusJson = $jsonObject->objToJson($status);
    $statusJson =~ s/"/'/g;    
    print "$statusJson\n";
    exit;    
}

#### UPDATE stage TABLE WITH JOB STATUS FROM QSTAT
method updateStageStatus($monitor, $stages) {
	my $statusHash = $monitor->statusHash();
	#exit;

	foreach my $stage ( @$stages )
	{
		my $stagejobid = $stage->{stagejobid};
		next if not defined $stagejobid or not $stagejobid;

		#### GET STATUS
		my $status;
		if ( defined $statusHash )
		{
			$status = $statusHash->{$stagejobid};
			next if not defined $status;

			#### SET TIME ENTRY TO BE UPDATED
			my $datetime = "queued";
			$datetime = "started" if defined $status and $status eq "running";

			$datetime = "completed" if not defined $status;
			$status = "completed" if not defined $status;

			#### UPDATE THE STAGE ENTRY IF THE STATUS HAS CHANGED
			if ( $status ne $stage->{status} )
			{
				my $query = qq{UPDATE stage
SET status='$status',
$datetime=NOW()
WHERE username ='$stage->{username}'
AND project = '$stage->{project}'
AND workflow = '$stage->{workflow}'
AND number='$stage->{number}'};
				my $result = $self->dbobject()->do($query);
			}
		}
	}

}

method clusterStatus {
	#### INSTANTIATE STARCLUSTER OBJECT
	my $starcluster = Agua::StarCluster->new(
		{
			username => $self->username(),
			project => $self->project(),
			workflow => $self->workflow(),
			cluster => $self->cluster(),
			conf	=>	$self->conf()
		}
	);

	return $starcluster->balancerOutput($self->cluster());
}

=head2

	SUBROUTINE		executeWorkflow

	PURPOSE

		EXECUTE THE WORKFLOW

		Agua::Workflow.pm -> executeWorkflow()
		|
		|
		-> Agua::Workflow.pm -> runStages()
			|
			| has many Agua::Stages
			|
			-> Agua::Stage -> run()
				| 
				|
				-> Agua::Stage -> execute() LOCAL JOB
				|
				OR
				|
				-> Agua::Stage -> clusterSubmit()  CLUSTER JOB

=cut

method executeWorkflow {


    my $dbobject     =	$self->dbobject();
    my $json         =	$self->json();
	my $username 	=	$json->{username};
	my $cluster 	=	$json->{cluster};
	my $project 	=	$json->{project};
	my $workflow 	=	$json->{workflow};

	#### GET THE CLUSTER
	print "Agua::Workflow::executeWorkflow    json:\n";
	print Dumper $json;

	my $status = $self->clusterStatus($username, $cluster);
	print "Agua::Workflow::executeWorkflow    status:\n";
	print Dumper $status;

	#### CREATE STARCLUSTER IF NOT STARTED
	print "{ status: 'Workflow::executeWorkflow    starting cluster: $cluster' } /**/\n\n";
	$self->startCluster() if $status->{status} ne "running";


	#### CREATE A UNIQUE QUEUE FOR THIS WORKFLOW
	my $envars = $self->getEnvars($username, $cluster);
	$self->createQueue($project, $workflow, $username, $cluster, $envars);

	#### SET STAGES
	my $stages = $self->setStages($json);
	print "Agua::Workflow::executeWorkflow    stages:\n";
	print Dumper $stages;

	#### CHECK PREVIOUS STAGE COMPLETED CORRECTLY
	$self->checkPrevious($stages, $json);

	#### GET STAGE PARAMETERS FOR THESE STAGES
	$stages = $self->setStageParameters($stages, $json);

	#### SET START
	$self->setStart($stages, $json);

	#### GET FILEROOT
	my $fileroot = $self->getFileroot($json->{username});	

	#### SET SCRIPT FILE
	my $scriptsdir = "$fileroot/$project/$workflow/scripts";
	File::Path::mkpath($scriptsdir) if not -d $scriptsdir;
	print "Can't create scripts directory: $scriptsdir\n"
		and exit if not -d $scriptsdir;

	#### SET STDOUT FILE 
	my $stdoutdir = "$fileroot/$project/$workflow/stdout";
	File::Path::mkpath($stdoutdir) if not -d $stdoutdir;
	print "Can't create stdout directory: $stdoutdir\n"
		and exit if not -d $stdoutdir;

	#### SET STDERR FILE
	my $stderrdir = "$fileroot/$project/$workflow/stderr";
	File::Path::mkpath($stderrdir) if not -d $stderrdir;
	print "Can't create stderr directory: $stderrdir\n"
		and exit if not -d $stderrdir;

	#### GET SETUID
	my $setuid = $self->json()->{setuid};
	print "Agua::Workflow::executeWorkflow    setuid: $setuid\n"

	#### WORKFLOW PROCESS ID
	my $workflowpid = $self->workflowpid();

	#### CLUSTER, QUEUE AND QUEUE OPTIONS
	my $queue = "$project-$workflow";

	my $queue_options = $self->json()->{queue_options};
	my $submit = $self->json()->{submit};

	#### SET OUTPUT DIR
	my $outputdir =  "$fileroot/$project/$workflow/";

	my $conf = $self->conf();
	print "{ status: 'Workflow::executeWorkflow    running workflow: $project.$workflow' } /**/\n\n";

	#my $workflow_number = $self->workflownumber();
	#if ( not defined $workflow_number )	{	$workflow_number = 1;	}

	#### GET MONITOR
	my $monitor = $self->getMonitor();

	my $stageobjects = [];
	my $number = $json->{number};
    my $start = $number - 1;
	for ( my $counter = $start; $counter < @$stages; $counter++ )
	{
		my $stage = $$stages[$counter];

		#### QUIT IF NO STAGE PARAMETERS
		print "Agua::Workflow::executeWorkflow    stageparameters not defined for stage $stage->{name}\n" and exit if not defined $stage->{stageparameters};

		my $stage_number = $counter + 1;

		#### SET MONITOR
		$stage->{monitor} = $monitor;

		#### SET SGE ENVIRONMENT VARIABLES
		$stage->{envars} = $envars;

        #### SET SCRIPT, STDOUT AND STDERR FILES
		$stage->{scriptfile} =  "$scriptsdir/$stage->{number}-$stage->{name}.sh";
        $stage->{stdoutfile} = "$stdoutdir/$stage->{number}-$stage->{name}.stdout";
        $stage->{stderrfile} = "$stderrdir/$stage->{number}-$stage->{name}.stderr";
#exit;

		$stage->{cluster}	=  $cluster;
		$stage->{submit}	=  $submit ;

		$stage->{workflowpid}=	$workflowpid;
		$stage->{dbobject}	=	$dbobject;
		$stage->{conf}		=  $conf;
		$stage->{fileroot}	=  $fileroot;
		$stage->{queue}		=  $queue;
		$stage->{queue_options}	=  $queue_options;
		$stage->{outputdir}	=  $outputdir;
		$stage->{setuid}	=  $setuid;
		$stage->{installdir}=  $conf->getKeyValue('agua', "INSTALLDIR");
		$stage->{qsub}		=  $conf->getKeyValue("cluster", "QSUB");
		$stage->{qstat}		=  $conf->getKeyValue("cluster", "QSTAT");

		my $stageobject = Agua::Stage->new($stage);
		$stageobject->toString();

		push @$stageobjects, $stageobject;
	}

	$self->stages($stageobjects);

	$self->runStages();
}


method startCluster {

    my $dbobject     =	$self->dbobject();
    my $json         =	$self->json();
	my $username 	=	$json->{username};
	my $cluster 	=	$json->{cluster};

	$json = $self->clusterInputs($json);



exit;


	#### USE StarCluster OBJECT TO START THE CLUSTER
	my $starcluster = Agua::StarCluster->new($json);
	print "Agua::Workflow::startCluster    Doing StarCluster->start()\n";
	$starcluster->start();

	#### SET STATUS IN clusterstatus TABLE

}

method createQueue ($project, $workflow, $username, $cluster, $envars) {
	print "Agua::Workflow::createQueue    project not defined\n" and exit if not defined $project;
	print "Agua::Workflow::createQueue    workflow not defined\n" and exit if not defined $workflow;

	#### GET ENVIRONMENT VARIABLES FOR THIS CLUSTER/CELL
	my $json = $self->json();
	$json->{qmasterport} = $envars->{qmasterport}->{value};
	$json->{execdport} = $envars->{execdport}->{value};
	$json->{sgecell} = $envars->{sgecell}->{value};
	$json->{sgeroot} = $envars->{sgeroot}->{value};

	#### DETERMINE WHETHER TO USE ADMIN KEY FILES
	my $adminkey = $self->getAdminKey($username);
	print "Agua::Workflow::createQueue    adminkey: $adminkey\n";
	$json->{configfile} =  $self->getConfigfile($username, $cluster);
	$json->{configfile} =  $self->getConfigfile($username, $cluster, 'admin')
		if $self->getAdminKey($username);
	print "Agua::Workflow::createQueue    configfile: $json->{configfile}\n";

	#### GET CLUSTER NODES INFO FROM cluster TABLE
	my $clusternodes = $self->getCluster($username, $cluster);
	$json->{nodetype} = $clusternodes->{instancetype};

	#### ADD CONF OBJECT	
	$json->{conf} = $self->conf();

	#### RENAME nodetype TO instancetype
	$json->{instancetype} = $json->{nodetype};


	#### RUN starcluster.pl TO GENERATE KEYPAIR FILE IN .starcluster DIR
	my $queue = "$project-$workflow";
	my $starcluster = Agua::StarCluster->new($json);
	print "Agua::Workflow::createQueue    Doing StarCluster->setQueue($queue)\n";
	$starcluster->setQueue($queue);
}

method getMonitor {
	return $self->monitor() if $self->monitor();

	my $clustertype =  $self->conf()->getKeyValue('agua', 'CLUSTERTYPE');
	my $classfile = "Agua/Monitor/" . uc($clustertype) . ".pm";
	my $module = "Agua::Monitor::$clustertype";
	require $classfile;

	my $monitor = $module->new(
		{
			pid			=>	$self->workflowpid(),
			conf 		=>	$self->conf(),
			dbobject	=>	$self->dbobject(),
			username	=>	$self->username(),
			project		=>	$self->project(),
			workflow	=>	$self->workflow(),
			cluster		=>	$self->cluster()
		}
	);

	$self->monitor($monitor);

	return $monitor;
}

method setStages ($json) {


	my $username = $json->{username};
    my $project = $json->{project};
    my $workflow = $json->{workflow};
    my $workflownumber = $json->{workflownumber};
	my $number = $json->{number};

	#### GET ALL STAGES FOR THIS WORKFLOW
    my $query = qq{SELECT * FROM stage
    WHERE username ='$username'
    AND project = '$project'
    AND workflownumber = '$workflownumber'
    AND workflow = '$workflow'
    ORDER BY number};
    my $stages = $self->dbobject()->queryhasharray($query);
	print "{ error: 'Agua::Workflow::setStages    stages not defined for username: $username' }" and exit if not defined $stages;	

	return $stages;
}

method checkPrevious ($stages, $json) {
	#### IF NOT STARTING AT BEGINNING, CHECK IF PREVIOUS STAGE COMPLETED SUCCESSFULLY

	my $number = $json->{number};
    my $start = $number - 1;	
	return if $start == 0;

	my $stage_number = $start - 1;
	$$stages[$stage_number]->{appname} = $$stages[$stage_number]->{name};
	$$stages[$stage_number]->{appnumber} = $$stages[$stage_number]->{number};
	my $keys = ["username", "project", "workflow", "name", "number"];
	my $where = $self->dbobject()->where($$stages[$stage_number], $keys);
	my $query = qq{SELECT status FROM stage $where};
	my $status = $self->dbobject()->query($query);

	return 1 if not defined $status or not $status;
	print "{ error: 'Agua::Workflow::checkPrevious    previous stage not completed: $stage_number' }"	and exit if $status ne "completed";
}


method setStageParameters ($stages, $json) {
	#### GET THE PARAMETERS FOR THE STAGES WE WANT TO RUN

	#### GET THE PARAMETERS FOR THE STAGES WE WANT TO RUN
	my $number = $json->{number};
    my $start = $number - 1;
	for ( my $i = $start; $i < @$stages; $i++ )
	{
		$$stages[$i]->{appname} = $$stages[$i]->{name};
		$$stages[$i]->{appnumber} = $$stages[$i]->{number};
		my $keys = ["username", "project", "workflow", "appname", "appnumber"];
		my $where = $self->dbobject()->where($$stages[$i], $keys);
		my $query = qq{SELECT * FROM stageparameter
$where AND paramtype='input'};
		my $stageparameters = $self->dbobject()->queryhasharray($query);
		$$stages[$i]->{stageparameters} = $stageparameters;

	}

	return $stages;
}

method setStart ($stages, $json) {
	print "{ error: 'Agua::Workflow::setStart    json->{number} not defined' }"
		and exit if not defined $json->{number};
	print "{ error: 'Agua::Workflow::setStart    json->{number} is non-numeric: ", $json->{number}, "}" if $json->{number} !~ /^\d+$/;
	print "error: 'Agua::Workflow::setStart    Runner starting stage $json->{number} is greater than the number of stages' }" and exit if $json->{number} > @$stages;

	$self->start($json->{number} - 1);
}

=head2

	SUBROUTINE		runStages

	PURPOSE

		RUN ALL OF THE STAGES OF THE WORKFLOW IN SUCCESSION
=cut

method runStages {

    #### RUN STAGES FROM start
	my $start 	=	$self->start();
	my $stages 	= 	$self->stages();
	for ( my $stage_counter = 0; $stage_counter < @$stages; $stage_counter++ )
	{
		my $stage = $$stages[$stage_counter];
		my ($completed) = $stage->run();
		print "Successfully completed stage $stage_counter: $completed\n" if defined $completed;

		#### STOP IF THIS STAGE DIDN'T COMPLETE SUCCESSFULLY (OUTPUTS '0' FOR SUCCESS)
		print "Agua::Workflow::runStages    Stage ", $stage->number(), ": ", $stage->name(), "did not complete successfully\n" and exit if not defined $completed;

		if ( defined $completed and $completed != 0 )
		{
			print "{ error: 'failed', name: '", $stage->name(), "', number: '", $stage->number(), "'}\n";

            $stage->setStatus('error');
			last;
		}

		#### OTHERWISE, REPORT COMPLETED
	}    

	return 1;
}

=head2

	SUBROUTINE		stopWorkflow

	PURPOSE

		KILL ANY PROCESSES FOR THE CURRENT STAGE AND STOP THE WORKFLOW

		LOCAL APPS:

			GET THE PROCESS IDS OF ANY RUNNING STAGE OF THE WORKFLOW AND

			'kill -9' THEM (INCLUDES STAGE PID, App PARENT PID AND App CHILD PID)

		CLUSTER APPS:

			1. IN THE CASE OF A SINGLE JOB, CANCEL THAT JOB ID

			2. IF AN APPLICATION IS RUNNING LOCALLY AND SUBMITTING JOBS TO THE

			CLUSTER, KILL ITS PROCESS ID AND CANCEL ANY JOBS IT HAS SUBMITTED

			(SHOULD BE REGISTERED IN THE stagejobs TABLE)

=cut

method stopWorkflow {

    my $dbobject     =	$self->dbobject();
    my $json         =	$self->json();

    #### VALIDATE
    print "{ error: 'Workflow::stopWorkflow    User session not validated' }" and exit unless $self->validate();

	#### SET EXECUTE WORKFLOW COMMAND
    my $bindir = $self->conf()->getKeyValue("agua", 'INSTALLDIR') . "/cgi-bin";

    my $username = $json->{username};
    my $project = $json->{project};
    my $workflow = $json->{workflow};
	my $number = $json->{number};
    my $start = $number - 1;

	#### GET ALL STAGES FOR THIS WORKFLOW
    my $query = qq{SELECT * FROM stage
WHERE username ='$username'
AND project = '$project'
AND workflow = '$workflow'
AND status='running'
ORDER BY number
};
	my $running = $dbobject->queryhasharray($query);

	#### EXIT IF NO PIDS
	print qq{ error: 'No running stages in workflow  $project.$workflow. Quitting' } and exit if not defined $running;

	#### WARNING IF MORE THAN ONE STAGE RETURNED (SHOULD NOT HAPPEN 
	#### AS STAGES ARE EXECUTED CONSECUTIVELY)
	print qq{ error: 'More than one stage running in workflow $project.$workflow. Continuing with stopWorkflow' } if scalar(@$running) > 1;

	my $cluster = $$running[0]->{cluster};
	my $messages;
	if ( defined $cluster and $cluster )
	{
		$messages = $self->killLocaljob($running);
	}
	else
	{
		$messages = $self->kill_clusterjob($running);
	}

	#### SET THIS STAGE STATUS AS 'waiting'
	my $update_query = qq{UPDATE stage
SET status = 'waiting'
WHERE username ='$username'
AND project = '$project'
AND workflow = '$workflow'
AND status='running'
};
	print "Agua::Workflow::$update_query\n";
	my $success = $dbobject->do($update_query);
	my $result;
	$result->{error} = qq{ error: 'Could not update 'running' stages for project $project workflow $workflow in the stage table' } if not $success;
	$result->{status} = qq{ error: 'Updated 'running' stages for project $project workflow $workflow in the stage table' } if $success;

	#### RETURN status AND messages
	$result->{messages} = $messages;

	#### PRINT RESULT	
	$self->printJSON($messages);
	exit;
}


=head2

	SUBROUTINE		killClusterjob

	PURPOSE

		1. CANCEL THE JOB IDS OF ANY RUNNING STAGE OF THE WORKFLOW:

			1. IN THE CASE OF A SINGLE JOB, CANCEL THAT JOB ID

			2. IF AN APPLICATION IS RUNNING LOCALLY AND SUBMITTING JOBS TO THE

			CLUSTER, KILL ITS PROCESS ID AND CANCEL ANY JOBS IT HAS SUBMITTED

			(SHOULD BE REGISTERED IN THE stagejobs TABLE)

=cut

method killClusterjob {
    my $running			=   shift;


    my $dbobject     =	$self->dbobject();
    my $json         =	$self->json();

    #### VALIDATE
    print "{ error: 'Workflow::killClusterJob    User session not validated' }" and exit unless $self->validate();

	my $messages = [];
	foreach my $stage ( @$running )
	{
		#### OTHERWISE, KILL ALL PIDS
		push @$messages, $self->cancelJob($stage->{childpid});
		push @$messages, $self->killPid($stage->{parentpid});
		#push @$messages, $self->killPid($stage->{stagepid});
		#push @$messages, $self->killPid($stage->{workflowpid});
	}

	return $messages;
}

=head2

	SUBROUTINE		cancelJob

	PURPOSE

		1. CANCEL A CLUSTER JOB BY JOB ID

=cut

method cancelJob {
    my $jobid			=	shift;


    #### VALIDATE
    print "{ error: 'Workflow::stopWorkflow    User session not validated' }" and exit unless $self->validate();

	my $canceljob = $self->conf()->getKeyValue(("cluster", 'CANCELJOB'));

	my $command = "$canceljob $jobid";

	return `$command`;
}

=head2

	SUBROUTINE		killLocaljob

	PURPOSE

		1. 'kill -9' THE PROCESS IDS OF ANY RUNNING STAGE OF THE WORKFLOW

		2. INCLUDES STAGE PID, App PARENT PID AND App CHILD PID)
=cut

method killLocaljob {
    my $running			=   shift;



    #### VALIDATE
    print "{ error: 'Workflow::killLocalJob    User session not validated' }" and exit unless $self->validate();

	my $messages = [];
	foreach my $stage ( @$running )
	{
		#### OTHERWISE, KILL ALL PIDS
		push @$messages, $self->killPid($stage->{childpid});
		push @$messages, $self->killPid($stage->{parentpid});
		push @$messages, $self->killPid($stage->{stagepid});
		push @$messages, $self->killPid($stage->{workflowpid});
	}

	return $messages;
}

=head2

	SUBROUTINE		killPid

	PURPOSE

		KILL THE PROCESS AND RETURN THE RESULT

=cut

method killPid {
    my $pid		=	shift;

    #### VALIDATE
    print "{ error: 'Workflow::killPid    User session not validated' }" and exit unless $self->validate();

	my $kill_command = "kill -9 $pid";
	print $kill_command, "\n";

	`$kill_command &> /dev/null`;

	return $kill_command;
}

=head2

	SUBROUTINE		workflow

	PURPOSE

		RETURN THE XML FOR ALL STAGES OF THE GIVEN WORKFLOW

=cut

method workflowStatus {
	my $workflow_number	=	shift;

    my $dbh	    =	$self->dbh();

	my $query = qq{SELECT status FROM stage
	WHERE workflownumber='$workflow_number'
	LIMIT 1};
	my $workflowStatus = Database::simple_query($dbh, $query);

	return $workflowStatus;
}



=head2

	SUBROUTINE		workflows

	PURPOSE

		RETURN A REFERENCE TO THE HASHARRAY OF WORKFLOW NUMBER AND NAMES

        IN THE collectionworkflow TABLE OF THE myEST DATABASE

=cut

method workflows {
    my $dbh     = $self->dbh();
	my $query = qq{SELECT DISTINCT project, workflownumber, workflow, number, name
    FROM stage ORDER BY workflownumber, number};
	my $workflows = Database::simple_queryhasharray($dbh, $query);

	return $workflows;
}


}

1;


