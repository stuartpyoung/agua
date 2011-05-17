use MooseX::Declare;

#### DEBUG
our $VERSION = 0.1;

=head2

	PACKAGE		Stage

	PURPOSE:

		A Stage IS ONE STEP IN A WORKFLOW. IT HAS THE FOLLOWING

		CHARACTERISTICS:

		1. EACH Stage RUNS ITSELF AND LOGS ITS STATUS TO THE stage

			DATABASE TABLE.

		2. A Stage WILL RUN LOCALLY BY DEFAULT. IF THE submit VARIABLE IS NOT

			ZERO AND cluster VARIABLE IS NOT EMPTY, IT WILL RUN ON A

			CLUSTER.

		3. EACH Stage DYNAMICALLY SETS ITS STDOUT, STDERR, INPUT AND OUTPUT

			FILES.

=cut 

use strict;
use warnings;

#### USE LIB FOR INHERITANCE
use FindBin qw($Bin);
use lib "$Bin/../";
use lib "$Bin/../external";

class Agua::Stage with (Agua::Common::Base, Agua::Common::Util, Agua::Cluster::Jobs) {

#### INTERNAL MODULES
use Agua::Monitor::PBS;

#### EXTERNAL MODULES
use IO::Pipe;
use Data::Dumper;
use overload '==' => 'identical';
use overload 'eq' => 'equal';

# INTS
has 'workflowpid'	=>	( isa => 'Int|Undef', is => 'rw' );
has 'stagepid'		=>	( isa => 'Int|Undef', is => 'rw' );
has 'stagejobid'	=>	( isa => 'Int|Undef', is => 'rw' );
has 'number'		=>  ( isa => 'Str', is => 'rw');
has 'workflownumber'=>  ( isa => 'Str', is => 'rw');
has 'start'     	=>  ( isa => 'Int', is => 'rw' );
has 'submit'     	=>  ( isa => 'Int|Undef', is => 'rw' );

# STRINGS
has 'clustertype'	=>  ( isa => 'Str|Undef', is => 'rw', default => "SGE" );
has 'fileroot'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'executor'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'location'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'username'  	=>  ( isa => 'Str', is => 'rw', required => 1  );
has 'workflow'  	=>  ( isa => 'Str', is => 'rw', required => 1  );
has 'project'   	=>  ( isa => 'Str', is => 'rw', required => 1  );
has 'name'   		=>  ( isa => 'Str', is => 'rw', required => 1  );
has 'queue'			=>  ( isa => 'Str', is => 'rw', required => 1  );
has 'queue_options'	=>  ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'outputdir'		=>  ( isa => 'Str', is => 'rw', required => 1  );
has 'setuid'		=>  ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'scriptfile'	=>  ( isa => 'Str', is => 'rw', required => 1 );
has 'stdoutfile'	=>  ( isa => 'Str', is => 'rw' );
has 'stderrfile'	=>  ( isa => 'Str', is => 'rw' );
has 'installdir'   	=>  ( isa => 'Str', is => 'rw', required => 1  );
has 'cluster'		=>  ( isa => 'Str|Undef', is => 'rw' );
has 'qsub'			=>  ( isa => 'Str', is => 'rw' );
has 'qstat'			=>  ( isa => 'Str', is => 'rw' );
has 'resultfile'	=>  ( isa => 'Str', is => 'ro', default => sub { "/tmp/result-$$" });

# OBJECTS
has 'envars'		=> ( isa => 'HashRef', is => 'rw', required => 1 );
has 'conf'			=> ( isa => 'Conf::Agua', is => 'rw', required => 1 );
has 'dbobject'		=> ( isa => 'Agua::DBase::MySQL', is => 'rw', required => 0 );
has 'monitor'		=> 	( isa => 'Maybe', is => 'rw', required => 0 );
has 'stageparameters'=> ( isa => 'ArrayRef', is => 'rw', required => 1 );

####////}

method BUILD ($args) {

	#use Data::Dumper;
	#$Data::Dumper::Deepcopy = 1;
}

=head2

	SUBROUTINE		run

	PURPOSE

		1. RUN THE STAGE APPLICATION AND UPDATE STATUS TO 'running'

		2. UPDATE THE PROGRESS FIELD PERIODICALLY (CHECKPROGRESS OR DEFAULT = 10 SECS)

		3. UPDATE STATUS TO 'complete' WHEN EXECUTED APPLICATION HAS FINISHED RUNNING

=cut

method run {

	#### TO DO: START PROGRESS UPDATER

	#### EXECUTE APPLICATION
	my $success;


	my $submit = $self->submit();
	my $cluster = $self->cluster();

	print "Agua::Stage::run    submit: $submit\n" if defined $submit;
	print "Agua::Stage::run    submit not defined\n" if not defined $submit;
	print "Agua::Stage::run    cluster: $cluster\n" if defined $cluster;
	print "Agua::Stage::run    cluster not defined\n" if not defined $cluster;

	#### SET STATUS TO running
	$self->setStatus('running');

	#### RUN ON CLUSTER
	if ( defined $cluster and $cluster and defined $submit and $submit )
	{
		$success = $self->runOnCluster();

	}
	#### RUN LOCALLY	
	else
	{
		$success = $self->runLocally();
	}
	#exit;

	#### REGISTER PROCESS IDS SO WE CAN MONITOR THEIR PROGRESS
	#### OR kill THEM LATER IF NECESSARY
	$self->registerIds();

	#### WAIT UNTIL PROCESS HAS FINISHED
	wait;

	#### UPDATE STATUS
	#$self->setStatus('completed');

	return $success;
}


=head2

	SUBROUTINE		runLocally

	PURPOSE

		EXECUTE THE APPLICATION

{"username":"admin","sessionId":"9999999999.9999.999","project":"Project1","workflow":"Workflow1","mode":"runLocallyWorkflow","number":1}



=cut

method runLocally {

    my $stageparameters =	$self->stageparameters();
	print "Agua::Stage::runLocally    stageparemeters not defined\n" and exit if not defined $stageparameters;

	#### GET FILE ROOT
	my $username = $self->username();
	my $fileroot = $self->fileroot();

	#### CONVERT ARGUMENTS INTO AN ARRAY IF ITS A NON-EMPTY STRING
	my $arguments = $self->setArguments($stageparameters);

	#### SET ENVIRONMENT VARIABLES AND EXECUTOR
	my $executor = $self->envars()->{tostring} . "; ";
	$executor .= $self->executor() if $self->executor();

	my $application = $self->location();

	#### ADD THE INSTALLDIR IF THE LOCATION IS NOT AN ABSOLUTE PATH
	my $installdir = $self->conf()->getKeyValue("agua", 'INSTALLDIR');
	if ( $application !~ /^\//
		and $application !~ /^[A-Z]:/i )
	{
		$application = "$installdir/bin/$application";
	}
	my @system_call = ($executor, $application, @$arguments);

	#### SET STDOUT AND STDERR FILES
	my $stdoutfile = $self->stdoutfile;
	my $stderrfile = $self->stderrfile;
	push @system_call, "1> $stdoutfile" if defined $stdoutfile;
	push @system_call, "2> $stderrfile" if defined $stderrfile;


    #### SET CHANGE DIR TO APPLICATION DIRECTORY, JUST IN CASE
    my ($changedir) = $application =~ /^(.+?)\/[^\/]+$/;

	print "Agua::Stage::runLocally    \$self->dbobject()->dbh(): ",  $self->dbobject()->dbh(), "\n";

	#### NO BUFFERING
	$| = 1;

	#### COMMAND
	my $command = join " ", @system_call;

	#### RUN APP BY FORKING
	my $stagepid;
	if ( $stagepid = fork() ) #### ****** Parent ****** 
	{
		print "Agua::Stage::runLocally    stagepid: $stagepid\n";
		$self->stagepid($stagepid);
	}
	elsif ( defined $stagepid ) #### ****** Child ******
	{
		#### SET InactiveDestroy ON DATABASE HANDLE
		$self->dbobject()->dbh()->{InactiveDestroy} = 1;
		my $dbh = $self->dbobject()->dbh();
		undef $dbh;
		#$self->dbobject()->dbh()->disconnect();

		`echo '$command' > /tmp/command.txt`;

		#### CHANGE DIR
        chdir($changedir);
		#my $result = system(@system_call);
		my $result = system($command);
		my $resultfile = $self->resultfile();
		`echo $result > $resultfile`;
		exit;
	}

    #### IF NEITHER CHILD NOR PARENT THEN COULDN'T OPEN PIPE
    else
    {
        print "{ error: 'Agua::Stage::runLocally    Could not open pipe (fork failed): $!' }";
		exit;
    }


	#### PAUSE AND THEN GET RESULT.
	#### CONFIRM APPLICATION RAN PROPERLY IN INITIAL STAGE
	wait;
	sleep(3);
	my $resultfile = $self->resultfile();
	open(RESULT, $resultfile);
	my $result = <RESULT>;
	close(RESULT);
	print "Agua::Stage::runLocally    resultfile: $resultfile\n";
	print "Agua::Stage::runLocally    result: $result\n";

	#### GET MAIN PARAMS
	my $project 	= 	$self->project();
	my $workflow 	= 	$self->workflow();
	my $workflowpid	=	$self->workflowpid();

	#### SET STATUS TO 'error' IF result IS NOT ZERO

	$self->setStatus('error') if not defined $result;

	return $result;
}

=head2

	SUBROUTINE		runOnCluster

	PURPOSE

		SUBMIT THE SHELLSCRIPT FOR EXECUTION ON A PBS- OR MOAB-ENABLED CLUSTER.

		IF THE 'SETUID' ARGUMENT EXISTS, RUN THE CLUSTER SUBMISSION VIA

		THIS SCRIPT (IT USES SETUID AS THE SPECIFIED USER PROVIDED THE USER

		EXISTS IN THE SYSTEM)

=cut

method runOnCluster {

	#### CLUSTER MONITOR
	my $monitor		=	$self->monitor();	

	#### GET MAIN PARAMS
	my $username 	= $self->username();
	my $project 	= $self->project();
	my $workflownumber 	= $self->workflownumber();
	my $workflow 	= $self->workflow();
	my $number 		= $self->number();
	my $queue 		= $self->queue();
	my $cluster		= $self->cluster();
	my $qstat		= $self->qstat();
	my $qsub		= $self->qsub();
	my $workflowpid = $self->workflowpid();

	#### SET DEFAULTS
	$queue = '' if not defined $queue;

	#### GET AGUA DIRECTORY FOR CREATING STDOUTFILE LATER
	my $aguadir 	= $self->conf()->getKeyValue("agua", 'AGUADIR');

	#### GET FILE ROOT
	my $fileroot = $self->getFileroot($username);

	#### GET ARGUMENTS ARRAY
    my $stageparameters =	$self->stageparameters();
    $stageparameters =~ s/\'/"/g;
	my $arguments = $self->setArguments($stageparameters);    

	#### GET EXECUTOR AND APPLICATION
	my $executor = $self->executor();
	my $application = $self->location();	

	#### ADD THE INSTALLDIR IF THE LOCATION IS NOT AN ABSOLUTE PATH
	my $installdir = $self->conf()->getKeyValue("agua", 'INSTALLDIR');
	if ( $application !~ /^\// and $application !~ /^[A-Z]:/i )
	{
		$application = "$installdir/bin/$application";
	}

	#### SET SYSTEM CALL
	my @system_call = ($application, @$arguments);
	my $command = "$executor @system_call";

    #### GET OUTPUT DIR
    my $outputdir = $self->outputdir();

	#### SET JOB NAME AS project-workflow-number
	my $label =	$project;
	$label .= "-" . $workflownumber;
	$label .= "-" . $workflow;
	$label .= "-" . $number;

	#### SET *** BATCH *** JOB 
	my $job = $self->setJob([$command], $label, $outputdir);

	#### GET FILES
	my $commands = $job->{commands};
	my $scriptfile = $job->{scriptfile};
	my $stdoutfile = $job->{stdoutfile};
	my $stderrfile = $job->{stderrfile};
	my $lockfile = $job->{lockfile};

	#### PRINT SHELL SCRIPT	
	$self->printSgeScriptfile($scriptfile, $commands, $label, $stdoutfile, $stderrfile, $lockfile);
	print "Agua::Cluster::Jobs::execute    scriptfile: $scriptfile\n";

	#### SET QUEUE
	$job->{queue} = $self->queue();

	#### SET QSUB
	$job->{qsub} = $self->qsub();

	#### SET SGE ENVIRONMENT VARIABLES
	$job->{envars} = $self->envars() if $self->envars();


	#### SUBMIT TO CLUSTER AND GET THE JOB ID 
	my $stagepid = $monitor->submitJob($job);
	$self->stagepid($stagepid);

	if ( not defined $stagepid or $stagepid =~ /^\s*$/ )
	{
		print "Agua::Stage::runOnCluster    Cluster submission failed - stagepid is not defined or empty\n";
		exit;
	}

	#### UPDATE STAGE WITH PARENT AND CHILD PIDs
	my $query = qq{UPDATE stage
SET workflowpid = '$workflowpid',
status='queued',
queued=NOW(),
started='',
completed='',
stagepid = '$stagepid'
WHERE username = '$username'
AND project = '$project'
AND workflow = '$workflow'
};
	my $dbobject = $self->dbobject();
	my $update_success = $dbobject->do($query);
	if ( not $update_success )
	{
		print "{ error: 'Could not update stage table with parent pid ($self->workflowpid() and stagepid ('$self->stagepid()') }";
		exit;
	}
	else
	{
	}	


	#### MONITOR UNTIL JOB STATUS IS running
	my $jobstatus = $monitor->jobStatus($stagepid);

	#### GET JOB STATUS
	my $sleep = $self->conf()->getKeyValue("cluster", 'SLEEP');
	$sleep = 5 if not defined $sleep;
	while ( $jobstatus->{status} eq "Idle" ) 
	{
		sleep($sleep);
		$jobstatus = $monitor->jobStatus($stagepid);
	}

	my $startttime = $jobstatus->{starttime};
	my $submitttime = $jobstatus->{submittime};
	my $runttime = $jobstatus->{runtime};

	if ( $jobstatus->{status} eq "Error" )
	{
		my $query = qq{UPDATE stage
SET workflowpid = '$workflowpid',
status='error',
completed=NOW(),
stagepid = '$self->stagepid()'
WHERE username = '$username'
AND project = '$project'
AND workflow = '$workflow'
};
		my $dbobject = $self->dbobject();
		my $update_success = $dbobject->do($query);
		if ( not $update_success )
		{
			print "{ error: 'Could not update stage table with parent pid ($self->workflowpid() and stagepid ('$self->stagepid()') }";
			exit;
		}
		else
		{
			print "{ error: 'Error with cluster submission, job id: $stagepid' }";
		}
	}
	elsif ( $jobstatus->{status} eq "Running" )
	{
		my $query = qq{UPDATE stage
SET workflowpid = '$workflowpid',
status='running',
started=NOW(),
completed='',
stagepid = '$self->stagepid()'
WHERE username = '$username'
AND project = '$project'
AND workflow = '$workflow'
};
		my $dbobject = $self->dbobject();
		my $update_success = $dbobject->do($query);
		if ( not $update_success )
		{
			print "{ error: 'Could not update stage table with parent pid ($self->workflowpid() and stagepid ('$self->stagepid()') }";
			exit;
		}
		else
		{
			print "{ status: 'Running.Submitted to cluster with job id $stagepid' }";
		}
	}
	elsif ( $jobstatus->{status} eq "Completed" )
	{
		my $query = qq{UPDATE stage
SET workflowpid = '$workflowpid',
status='completed',
completed=NOW(),
stagepid = '$self->stagepid()'
WHERE username = '$username'
AND project = '$project'
AND workflow = '$workflow'
};
		my $dbobject = $self->dbobject();
		my $update_success = $dbobject->do($query);
		if ( not $update_success )
		{
			print "{ error: 'Could not update stage table with parent pid ($self->workflowpid() and stagepid ('$self->stagepid()') }";
			exit;
		}
		else
		{
			print "{ status: 'Running. Submitted to cluster with job id $stagepid' }";
		}
	}	

}	#	runOnCluster

method updateStatus ($set, $username, $project, $workflow) {

	my $dbobject 	=	$self->dbobject();

	my $query = qq{UPDATE stage
SET $set
WHERE username = '$username'
AND project = '$project'
AND workflow = '$workflow'
};
	my $update_success = $dbobject->do($query);
	if ( not $update_success )
	{
		print "{ error: 'Could not update stage table for username $username, project $project, workflow $workflow with set clause: $set' }";
		exit;
	}
}

=head2

	SUBROUTINE		clusterArguments

	PURPOSE

		SET ARGUMENTS AND GET VALUES FOR ALL PARAMETERS

=cut

method setArguments ($stageparameters) {


	#### SANITY CHECK
	return if not defined $stageparameters;
	return if ref($stageparameters) eq '';

	#### GET FILEROOT
	my $username 	= $self->username();
	my $fileroot 	= $self->getFileroot($username);

	#### GENERATE ARGUMENTS ARRAY
	my $arguments = [];
	foreach my $stageparameter (@$stageparameters)
	{

		my $name	 	=	$stageparameter->{name};
		my $argument 	=	$stageparameter->{argument};
		my $value 		=	$stageparameter->{value};
		my $valuetype 	=	$stageparameter->{valuetype};
		my $discretion 	=	$stageparameter->{discretion};


		#### SKIP EMPTY FLAG OR ADD 'checked' FLAG
		if ( $valuetype eq "flag" )
		{
			if (not defined $value or not $value)
			{
				next;
			}

			push @$arguments, $argument;
			next;
		}

		if ( $value =~ /^\s*$/ and $discretion ne "required" )
		{
			next;
		}

		if ( defined $value )
		{

			#### ADD THE FILE ROOT FOR THIS USER TO FILE/DIRECTORY PATHS
			#### IF IT DOES NOT BEGIN WITH A '/', I.E., AN ABSOLUTE PATH
			if ( $valuetype =~ /^(file|directory)$/ and $value =~ /^[^\/]/ )
			{	
				print "Agua::Stage::setArguments    Adding fileroot to $valuetype: $value\n";
				$value =~ s/^\///;
				$value = "$fileroot/$value";
			}


			#### ADD THE FILE ROOT FOR THIS USER TO FILE/DIRECTORY PATHS
			#### IF IT DOES NOT BEGIN WITH A '/', I.E., AN ABSOLUTE PATH
			if ( $valuetype =~ /^(files|directories)$/ and $value =~ /^[^\/]/ )
			{	
				print "Agua::Stage::setArguments    Adding fileroot to $valuetype: $value\n";
				my @subvalues = split ",", $value;
				foreach my $subvalue ( @subvalues )
				{
					$subvalue =~ s/^\///;
					$subvalue = "$fileroot/$subvalue";
				}

				$value = join ",", @subvalues;
			}

			#### SINGLE '-' OPTIONS (E.G., -i)
			if ( $argument =~ /^\-[^\-]/ )
			{
				push @$arguments, qq{$argument $value};
			}

			#### DOUBLE '-' OPTIONS (E.G., --inputfile)
			else
			{
				push @$arguments, $argument;
				push @$arguments, $value;
			}

		}
	}



	return $arguments;
}


=head2

	SUBROUTINE		registerIds

	PURPOSE

		SET THE PROCESS IDS FOR:

			- THE STAGE ITSELF

			- THE PARENT OF THE STAGE'S APPLICATION (SAME AS STAGE)

			- THE CHILD OF THE STAGE'S APPLICATION

=cut

method registerIds {
	print "Agua::Stage::registerIds    Agua::Stage::registerIds()\n";

    my $dbobject 	= $self->dbobject();
	my $workflowpid = $self->workflowpid();
	my $stagepid 	= $self->stagepid() || '';
	my $stagejobid = $self->stagejobid() || '';
	my $username 	= $self->username();
	my $project 	= $self->project();
	my $workflow 	= $self->workflow();
	my $workflownumber = $self->workflownumber();
	my $number		 = $self->number();


	#### UPDATE status TO waiting IN TABLE stage
    my $query = qq{UPDATE stage
    SET workflowpid='$workflowpid',
	stagepid='$stagepid',
	stagejobid='$stagejobid'	
	WHERE username = '$username'
    AND project = '$project'
    AND workflow = '$workflow'
    AND workflownumber = '$workflownumber'
    AND number = '$number'};
    my $success = $dbobject->do($query);
	if ( not $success )
	{
		#warn "Stage::register    Could not insert entry for stage $self->stagenumber() into 'stage' table\n";
        return 0;
    }

	return 1;
}


=head2

	SUBROUTINE		register

	PURPOSE

		UPDATE THE status TO waiting FOR A STAGE IN THE stage TABLE

=cut

method register {
	my $status	=	shift;


    my $dbobject = $self->dbobject();

	#### SET SELF _status TO waiting
	$self->status('waiting');

	my $username = $self->username();
	my $project = $self->project();
	my $workflow = $self->workflow();
	my $workflownumber = $self->workflownumber();
	my $number = $self->number();

	#### UPDATE status TO waiting IN TABLE stage
    my $query = qq{UPDATE stage
    SET status='waiting'
	WHERE username = '$username'
    AND project = '$project'
    AND workflow = '$workflow'
    AND workflownumber = '$workflownumber'
    AND number = '$number'};
    my $success = $dbobject->do($query);
	if ( not $success )
	{
		warn "Stage::register    Could not insert entry for stage $self->stagenumber() into 'stage' table\n";
        return 0;
    }

	return 1;
}

=head2

	SUBROUTINE		isComplete

	PURPOSE

		CHECK IF THIS STAGE HAS STATUS 'complete' IN THE stage

	INPUT

		WORKFLOW NAME (workflow) AND STAGE NAME (name)

	OUTPUT

		RETURNS 1 IF COMPLETE, 0 IF NOT COMPLETE

=cut

method isComplete {
    my $dbobject = $self->dbobject();

	my $project = $self->workflownumber();
	my $workflow = $self->workflownumber();
	my $workflow_number = $self->workflownumber();
	my $number = $self->number();

	my $query = qq{SELECT status
	FROM stage
	WHERE project='$project'
	AND workflownumber = '$workflow_number'
	AND workflow = '$workflow'
	AND number = '$number'};
	my $complete = $dbobject->query($query);

	return 0 if not defined $complete or not $complete;
	return 1;
}


=head2

	SUBROUTINE		setStatus

	PURPOSE

		SET THE status FIELD IN THE stage TABLE FOR THIS STAGE

=cut

method setStatus ($status) {	

    my $dbobject = $self->dbobject();

	#### GET KEY INFO FOR MYSQL TABLE
	my $project = $self->project();
	my $workflow = $self->workflow();
	my $workflownumber = $self->workflownumber();
	my $number = $self->number();
	my $now = $dbobject->now();

	#### SET STATUS
	my $update_query = qq{UPDATE stage
	SET
	status = '$status',
	started = $now
	WHERE project = '$project'
	AND workflow = '$workflow'
	AND workflownumber = '$workflownumber'
	AND number = '$number'};
	my $success = $dbobject->do($update_query);
	if ( not $success )
	{
		print "Agua::Stage::setStatus    Could not update stage (project: $project, workflow: $workflow, workflownumber: $workflownumber, number: $number) status : '$status'\n";
		exit;
	}
}


method toString () {
	print $self->_toString();
}

method _toString () {
	my @keys = qw[ username project workflownumber workflow name number start executor location fileroot queue queue_options outputdir scriptfile stdoutfile stderrfile workflowpid stagepid stagejobid submit setuid installdir cluster qsub qstat resultfile];
	my $string = '';
	foreach my $key ( @keys )
	{
		my $filler = " " x (20 - length($key));
		$string .= "$key$filler:\t";
		$string .= $self->$key() || '';
		$string .= "\n";
	}
	$string .= "\n\n";
}


} #### Agua::Stage

