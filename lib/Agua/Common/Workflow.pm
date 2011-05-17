package Agua::Common::Workflow;
use Moose::Role;
use Moose::Util::TypeConstraints;

=head2

	PACKAGE		Agua::Common::Workflow

	PURPOSE

		WORKFLOW METHODS FOR Agua::Common

=cut


use Data::Dumper;

=head2

	SUBROUTINE		downloadHistory

	PURPOSE

		SEND A FILE DOWNLOAD OF WORKFLOW APPLICATIONS

echo {"username":"syoung","sessionId":"1228791394.7868.158","mode":"downloadHistory"} | perl -U download.cgi

=cut


sub downloadHistory {
	my $self		=	shift;	


    #### GET DATABASE OBJECT
    my $dbobject     =	$self->dbobject();
    my $json         =	$self->json();

    #### VALIDATE    
    my $username = $json->{username};
    my $session_id = $json->{sessionId};
	my $validated = $self->validate();
    print "{ error: 'User session not validated' }" and exit unless $validated;

	#### GET PROJECTS AND WORKFLOWS
	my $query = qq{select DISTINCT project,workflow FROM stage WHERE username='$username' order by started};
	my $project_workflows = $dbobject->queryhasharray($query);
	print "{}" and exit if not defined $project_workflows;

	#### CONVERT TO PROJECT:WORKFLOWS ARRAY HASH
	my $projects;

	my $current_project = $$project_workflows[0]->{project};
	my $workflow_array;
	foreach my $row ( @$project_workflows )
	{
		my $project = $row->{project};
		if ( $project eq $current_project )
		{
			push @$workflow_array, $row->{workflow};
		}
		else
		{
			push @$projects, { $current_project => $workflow_array };
			$current_project = $project;
			$workflow_array = undef;
			push @$workflow_array, $row->{workflow};
		}		
	}
	push @$projects, { $current_project => $workflow_array } if defined $workflow_array;

	my $output;
	foreach my $projecthash ( @$projects)
	{
		my ($project) = keys %$projecthash;
		my $workflows = $projecthash->{$project};

		foreach my $workflow ( @$workflows )
		{
			my $query = qq{select * FROM stage WHERE username='$username' AND project='$project' and workflow='$workflow' ORDER BY number DESC};
			my $applications = $dbobject->queryhasharray($query);

			$output .= qq{Project: $project\tWorkflow: $workflow\n};
			$output .= qq{\nApplication\tNumber\tStatus\tStarted\tCompleted\n};

			#### PRINT APPLICATION NAME, NUMBER, STATUS, STARTED AND COMPLETED
			foreach my $application ( @$applications )
			{
				#username            VARCHAR(30),
				#project             VARCHAR(30),
				#workflow            VARCHAR(60),
				#workflownumber      INT(12),
				#name                VARCHAR(60),
				#number              VARCHAR(10),
				#arguments           TEXT,
				#outputs             TEXT,
				#inputs              TEXT,
				#inputfiles          TEXT,
				#started             DATETIME,
				#completed           DATETIME,
				#workflowpid         INT(12),
				#parentpid           INT(12),
				#childpid            INT(12),
				#status              VARCHAR(20),
				#runnumber           INT(20),

				$output .= qq{$application->{name}\t$application->{number}\t$application->{status}\t$application->{started}\t$application->{completed}\n};

				my $inputshash = $application->{arguments};
				if ( defined $inputshash )
				{
					my $array = $self->parseHash($inputshash);
					my $lines = "\nInputs\n";
					$lines .= "name\tvalue\n";
					$lines .= join "", @$array;
					$output .= $lines;
				}

				my $outputshash = $application->{outputs};
				if ( defined $outputshash )
				{
					my $array = $self->parseHash($outputshash);
					my $lines = "\nOutputs\n";
					$lines .= "name\tvalue\n";
					$lines .= join "", @$array;
					$output .= $lines;
				}
			}

			$output .= "\n";
		}
	}

	#### PRINT HISTORY
	my $filesize = length($output);
	my $datetime = $self->datetime();	
	my $filename = "$username-agua-history-$datetime.txt";

	print qq{Content-type: application/x-download\n};
	print qq{Content-Disposition: attachment;filename="$filename"\n};
	print qq{Content-Length: $filesize\n\n};
	print $output;
	exit;
}

=head2

	SUBROUTINE		getHistory

	PURPOSE

		RETURN THE LIST OF WORKFLOW APPLICATIONS INCLUDING

		STARTED, COMPLETED, ETC. INFORMATON

echo {"username":"syoung","sessionId":"1228791394.7868.158","mode":"getHistory"} | perl -U workflow.cgi

=cut

sub getHistory {
	my $self		=	shift;	



    #### GET DATABASE OBJECT
    my $dbobject     =	$self->dbobject();
    my $json         =	$self->json();

    #### VALIDATE    
    my $username = $json->{username};
    my $session_id = $json->{sessionId};
	my $validated = $self->validate();
    print "{ error: 'Agua::Workflow::getHistory    User session not validated' }" and exit unless $validated;

	#### GET PROJECT WORKFLOWS THAT HAVE BEEN UPDATED MOST RECENTLY
	my $query = qq{select DISTINCT project,workflow
FROM stage
WHERE username='$username'
ORDER BY started,project,workflow};
	my $project_workflows = $dbobject->queryhasharray($query);
	print "{}" and exit if not defined $project_workflows;

	#### GET HISTORY FOR EACH PROJECT WORKFLOW
	my $outputs;
	foreach my $projecthash ( @$project_workflows )
	{
		my $project = $projecthash->{project};
		my $workflow = $projecthash->{workflow};
		my $query = qq{select project, workflow, number, name, status, started, completed FROM stage
WHERE username='$username'
AND project='$project'
AND workflow='$workflow'
ORDER BY number DESC};

		my $stages = $dbobject->queryhasharray($query);
		push @$outputs, $stages if defined $stages;
	}

	#### PRINT HISTORY
	use JSON -support_by_pp; 
	my $jsonParser = JSON->new();
	my $history = $jsonParser->allow_singlequote(1)->allow_nonref->loose(1)->pretty->encode($outputs);
	print $history;
	exit;
}






=head2

	SUBROUTINE		addWorkflow

	PURPOSE

		ADD A WORKFLOW TO THE workflow TABLE

=cut

sub addWorkflow {
	my $self		=	shift;


    my $dbobject        =	$self->dbobject();
    my $json 			=	$self->json();


	my $success = $self->_addWorkflow();
	print "{ status: 'Successful insert of workflow $json->{name} into project $json->{project} in workflow table' }" if $success;
	exit;
}




=head2

	SUBROUTINE		_addWorkflow

	PURPOSE

		INTERNAL USE ONLY: ADD A WORKFLOW TO THE workflow TABLE

=cut

sub _addWorkflow {
	my $self		=	shift;


    my $dbobject        =	$self->dbobject();
    my $json 			=	$self->json();


    #### VALIDATE
    print "{ error: 'Agua::Common::Workflow::_addWorkflow    User session not validated' }" and exit unless $self->validate();

	#### SET TABLE AND REQUIRED FIELDS	
	my $table = "workflow";
	my $required_fields = ["username", "project", "name"];

	#### CHECK REQUIRED FIELDS ARE DEFINED
	my $not_defined = $dbobject->notDefined($json, $required_fields);
    print "{ error: 'Agua::Common::Workflow::_addWorkflow    not defined: @$not_defined' }" and exit if @$not_defined;

	#### REMOVE IF EXISTS ALREADY
	$self->_removeFromTable($table, $json, $required_fields);

	#### GET MAX WORKFLOW NUMBER IF NOT DEFINED
    my $username = $json->{'username'};
	my $query = qq{SELECT MAX(number)
	FROM workflow
	WHERE username = '$username'
	AND project = '$json->{project}'};
	my $number = $dbobject->query($query);
	$number = 1 if not defined $number;
	$number++ if defined $number;
	$json->{number} = $number;

	#### DO ADD
	my $success = $self->_addToTable($table, $json, $required_fields);	
 	print "{ error: 'Agua::Common::Workflow::_addWorkflow    Could not add workflow $json->{workflow} into project $json->{project} in workflow table'}" and exit if not defined $success;

	#### ADD THE PROJECT DIRECTORY TO THE USER'S agua DIRECTORY
	my $fileroot = $self->getFileroot($username);
	my $filepath = "$fileroot/$json->{project}/$json->{name}";

	if ( $^O =~ /^MSWin32$/ )   {   $filepath =~ s/\//\\/g;  }

	File::Path::mkpath($filepath);
	print "{ error: 'Agua::Common::Workflow::_addWorkflow    Could not create the fileroot directory: $fileroot' }" and exit if not -d $filepath;

	return 1;
}




=head2

	SUBROUTINE		removeWorkflow

	PURPOSE

		REMOVE A WORKFLOW FROM workflow, stage AND stageparameter

=cut

sub removeWorkflow {
	my $self		=	shift;

    my $dbobject        =	$self->dbobject();
    my $json 			=	$self->json();

    #### VALIDATE
    print "{ error: 'Agua::Common::Workflow::removeWorkflow    User session not validated' }" and exit unless $self->validate();

	#### SET TABLE AND REQUIRED FIELDS	
	my $table = "workflow";
	my $required_fields = ["username", "project", "name"];

	#### CHECK REQUIRED FIELDS ARE DEFINED
	my $not_defined = $dbobject->notDefined($json, $required_fields);
    print "{ error: 'Agua::Common::Workflow::removeWorkflow    undefined values: @$not_defined' }" and exit if @$not_defined;

	#### REMOVE FROM workflow
	my $success = $self->_removeFromTable($table, $json, $required_fields);
 	print "{ error: 'Agua::Common::Workflow::removeWorkflow    Could not delete workflow $json->{workflow} from project $json->{project} in workflow table'}" and exit if not defined $success;

	#### SET 'workflow' FIELD
	$json->{workflow} = $json->{name};

	#### REMOVE FROM stage
	$table = "stage";
	my $stage_fields = ["username", "project", "workflow"];
	$self->_removeFromTable("stage", $json, $stage_fields);

	#### REMOVE FROM clusters
	$table = "cluster";
	my $clusters_fields = ["username", "project", "workflow"];
	$self->_removeFromTable("cluster", $json, $clusters_fields);

	#### REMOVE FROM stageparameter
	$table = "stageparameter";
	my $stageparameter_fields = ["username", "project", "workflow"];
	$self->_removeFromTable($table, $json, $stageparameter_fields);

	#### REMOVE FROM view
	$table = "view";
	my $view_fields = ["username", "project"];
	$self->_removeFromTable($table, $json, $view_fields);

	#### REMOVE WORKFLOW DIRECTORY
    my $username = $json->{'username'};
	my $fileroot = $self->getFileroot($username);

	my $filepath = "$fileroot/$json->{project}/$json->{workflow}";
	if ( $^O =~ /^MSWin32$/ )   {   $filepath =~ s/\//\\/g;  }

	print "{ error: 'Agua::Common::Workflow::removeWorkflow    Cannot remove directory: $filepath' }" and exit if not File::Remove::rm(\1, $filepath);

	print "{ status: 'Agua::Common::Workflow::removeWorkflow    Successfully deleted workflow $json->{workflow} from project $json->{project} in workflow table' }";

}	#### removeWorkflow



=head2

	SUBROUTINE		renameWorkflow

	PURPOSE

		RENAME A WORKFLOW IN workflow, stage AND stageparameter

=cut

sub renameWorkflow {
	my $self		=	shift;



    my $dbobject        =	$self->dbobject();
    my $json 			=	$self->json();


	#### GET NEW NAME
	my $newname = $json->{newname};
    print "{ error: 'Agua::Common::Workflow::renameWorkflow    No newname parameter. Exiting' }" and exit if not defined $newname;

    #### VALIDATE
    print "{ error: 'Agua::Common::Workflow::renameWorkflow    User session not validated' }" and exit unless $self->validate();

	#### QUIT IF NEW NAME EXISTS ALREADY
	my $query = qq{SELECT name FROM workflow
WHERE project='$json->{project}'
AND name='$newname'};
	my $already_exists = $dbobject->query($query);
	if ( $already_exists )
	{
	    print "{ error: 'Agua::Common::Workflow::renameWorkflow    New name $newname already exists in workflow table' }";
		exit;
	}

	#### SET TABLE AND REQUIRED FIELDS	
	my $table = "workflow";
	my $required_fields = ["username", "project", "name"];

	#### CHECK REQUIRED FIELDS ARE DEFINED
	my $not_defined = $dbobject->notDefined($json, $required_fields);
    print "{ error: 'Agua::Common::Workflow::renameWorkflow    undefined values: @$not_defined' }" and exit if @$not_defined;

	#### UPDATE workflow
	my $set_hash = { name => $newname };
	my $set_fields = ["name"];
	my $success = $self->_updateTable($table, $json, $required_fields, $set_hash, $set_fields);
 	print "{ error: 'Agua::Common::Workflow::renameWorkflow    Could not rename workflow '$json->{workflow}' to '$newname' in $table table'}" and exit if not defined $success;

	#### SET 'workflow' FIELD FOR STAGE AND STAGEPARAMETER TABLES
	$json->{workflow} = $json->{name};

	#### UPDATE stage
	$table = "stage";
	$set_hash = { workflow => $newname };
	$set_fields = ["workflow"];
	$required_fields = ["username", "project", "name"];
	$self->_updateTable($table, $json, $required_fields, $set_hash, $set_fields);

	#### UPDATE stage
	$table = "stageparameter";
	$self->_updateTable($table, $json, $required_fields, $set_hash, $set_fields);

	#### UPDATE clusters
	$table = "cluster";
	$set_hash = { workflow => $newname };
	$set_fields = ["workflow"];
	$required_fields = ["username", "project", "workflow"];
	$self->_updateTable($table, $json, $required_fields, $set_hash, $set_fields);

	#### RENAME WORKFLOW DIRECTORY
	my $fileroot = $self->getFileroot();
	my $old_filepath = "$fileroot/$json->{project}/$json->{name}";
	my $new_filepath = "$fileroot/$json->{project}/$json->{newname}";
	if ( $^O =~ /^MSWin32$/ )   {   $old_filepath =~ s/\//\\/g;  }

	#### CHECK IF WORKFLOW DIRECTORY EXISTS
	print "{ error: 'Agua::Common::Workflow::renameWorkflow    Cannot find old workflow directory: $old_filepath' }" and exit if not -d $old_filepath;

	#### RENAME WORKFLOW DIRECTORY
	File::Copy::move($old_filepath, $new_filepath);

	print "{ error: 'Agua::Common::Workflow::renameWorkflow    Could not rename directory: $old_filepath to $new_filepath' }" and exit if not -d $new_filepath;

	print "{ status: 'Agua::Common::Workflow::renameWorkflow    Successfully renamed workflow $json->{workflow} to $newname in workflow table' }";

}	#### renameWorkflow




=head2

    SUBROUTINE:     getWorkflows

    PURPOSE:

		RETURN AN ARRAY OF workflow WORKFLOW HASHES

=cut


sub getWorkflows {
	my $self		=	shift;


    my $dbobject		=	$self->dbobject();
	my $json			=	$self->json();

    #### VALIDATE
    my $username = $json->{'username'};
    print "{ error: 'User session not validated' }" and exit unless $self->validate();

	#### GET ALL SOURCES
	my $query = qq{SELECT * FROM workflow
WHERE username='$username'
ORDER BY project, name};
	my $workflows = $dbobject->queryhasharray($query);

	######## IF NO RESULTS:
	####	2. INSERT DEFAULT WORFKLOW INTO workflow TABLE
	####	2. CREATE DEFAULT WORKFLOW FOLDERS
	return $self->_defaultWorkflow() if not defined $workflows;

	return $workflows;
}

=head2

	SUBROUTINE		_defaultWorkflow

	PURPOSE

		1. INSERT DEFAULT WORKFLOW INTO workflow TABLE

		2. CREATE DEFAULT PROJECT AND WORKFLOW FOLDERS

	INPUT

		1. USERNAME

		2. SESSION ID

	OUTPUT

		1. JSON HASH { project1 : { workflow}

=cut

sub _defaultWorkflow {
	my $self		=	shift;	



    #### GET DATABASE OBJECT
    my $dbobject     =	$self->dbobject();
    my $json         =	$self->json();

    #### VALIDATE    
    print "{ error: 'Agua::Common::Workflow::_defaultWorkflow    User session not validated' }" and exit unless $self->validate();

	#### SET DEFAULT WORKFLOW
	$self->json()->{project} = "Project1";
	$self->json()->{name} = "Workflow1";
	$self->json()->{number} = 1;

	#### ADD WORKFLOW1
	my $success = $self->_addWorkflow();
 	print "{ error: 'Agua::Common::Workflow::_defaultWorkflow    Could not add workflow $json->{workflow} into  workflow table' }" and exit if not defined $success;

	#### DO QUERY
    my $username = $json->{'username'};
	my $query = qq{SELECT * FROM workflow
WHERE username='$username'
ORDER BY project, name};
	my $workflows = $dbobject->queryhasharray($query);

	return $workflows;
}




=head2

	SUBROUTINE		copyWorkflow

	PURPOSE

        COPY A WORKFLOW TO ANOTHER (NON-EXISTING) WORKFLOW:

			1. UPDATE THE workflow TABLE TO ADD THE NEW WORKFLOW

			3. COPY THE WORKFLOW DIRECTORY TO THE NEW WORKFLOW IF

                copyFile IS DEFINED

                 echo '{"sourceuser":"admin","targetuser":"syoung","sourceworkflow":"Workflow0","sourceproject":"Project1","targetworkflow":"Workflow9","targetproject":"Project1","username":"syoung","sessionId":"9999999999.9999.999","mode":"copyWorkflow"}' |  ./workflow.cgi

=cut
sub copyWorkflow {
	my $self				=	shift;
    my $dbobject        =	$self->dbobject();
    my $json 			=	$self->json();


    print "{ error: 'Agua::Common::Workflow::copyWorkflow    No data provided' }"
        if not defined $json;

	my $sourceuser     = $json->{sourceuser};
	my $targetuser     = $json->{targetuser};
	my $sourceproject  = $json->{sourceproject};
	my $sourceworkflow = $json->{sourceworkflow};
	my $targetproject  = $json->{targetproject};
	my $targetworkflow = $json->{targetworkflow};
	my $copyfiles      = $json->{copyfiles};

    print "Agua::Common::Workflow::copyWorkflow    sourceuser: $sourceuser\n";
	print "Agua::Common::Workflow::copyWorkflow    targetuser: $targetuser\n";
    print "Agua::Common::Workflow::copyWorkflow    sourceworkflow: $sourceworkflow\n";
	print "Agua::Common::Workflow::copyWorkflow    targetworkflow: $targetworkflow\n";

    print "{ error: 'Agua::Common::Workflow::copyWorkflow    User not validated: $targetuser' }" and exit if not $self->validate();
    print "{ error: 'Agua::Common::Workflow::copyWorkflow    targetworkflow not defined: $targetworkflow' }" and exit if not defined $targetworkflow or not $targetworkflow;

    my $can_copy;
	$can_copy = 1 if $sourceuser eq $targetuser;
	$can_copy = $self->canCopy($sourceuser, $sourceproject, $targetuser) if $sourceuser ne $targetuser;
	print "Agua::Common::Workflow::copyWorkflow    can_copy: $can_copy\n";

    print "{ error: 'Agua::Common::Workflow::copyWorkflow    Insufficient privileges for user: $targetuser' }"
        and exit if not $can_copy;

	#### CONFIRM THAT SOURCE PROJECT EXISTS IN project TABLE
	my $query = qq{SELECT *
FROM workflow
WHERE username = '$sourceuser'
AND project = '$sourceproject'
AND name = '$sourceworkflow'};
	my $workflowObject = $dbobject->queryhash($query);
    print "{ error: 'Agua::Common::Workflow::copyProject    Source workflow does not exist: $targetworkflow' }"
        and exit if not defined $workflowObject;
    my $description = $workflowObject->{description};
    my $notes = $workflowObject->{notes};

#	#### CHECK IF TARGET PROJECT EXISTS IN project TABLE
#	$query = qq{SELECT project
#FROM workflow
#WHERE username = '$targetuser'
#AND project = '$targetproject'
#AND name = '$targetworkflow'};
#	my $project_exists = $dbobject->query($query);
#
#        and exit if defined $project_exists and $project_exists;

	#### CHECK IF WORKFLOW ALREADY EXISTS IN project TABLE
	$query = qq{SELECT 1
FROM workflow
WHERE username = '$targetuser'
AND project = '$targetproject'
AND name = '$targetworkflow'};
	my $workflow_exists = $dbobject->query($query);
	print "{ error: 'Agua::Common::Workflow::copyWorkflow    Workflow already exists: $targetworkflow' }"
        and  exit if $workflow_exists;

	#### GET WORKFLOW NUMBER FOR THE WORKFLOW TO BE INSERTED	
	$query = qq{SELECT MAX(number)
FROM workflow
WHERE username = '$targetuser'
AND project = '$targetproject'};
	my $workflownumber = $dbobject->query($query);
	$workflownumber = 0 if not defined $workflownumber;
	$workflownumber++;

	#### ADD NEW NUMBER TO WORKFLOW OBJECT
    $workflowObject->{number} = $workflownumber;

    #### COPY WORKFLOW OBJECT AND RELATED INFORMATION (STAGES, STAGE
    #### PARAMETERS, ETC.)
    $self->_copyWorkflow($sourceuser, $targetuser, $sourceproject, $targetproject, $sourceworkflow, $targetworkflow, $workflowObject);

    #### COPY FILES IF FLAGGED BY copyfiles 
    if ( defined $copyfiles and $copyfiles )
    {
        #### SET DIRECTORIES
        my $aguadir = $self->conf()->getKeyValue("agua", 'AGUADIR');
        my $userdir = $self->conf()->getKeyValue("agua", 'USERDIR');
        my $sourcedir = "$userdir/$sourceuser/$aguadir/$sourceproject/$sourceworkflow";
        my $targetdir = "$userdir/$targetuser/$aguadir/$targetproject/$targetworkflow";

        #### COPY DIRECTORY
        my $copy_result = $self->copyFilesystem($sourcedir, $targetdir);
        print "{ status: 'Agua::Common::Workflow::copyWorkflow    Copied to $targetworkflow' }"
            and exit if $copy_result;
        print "{ status: 'Agua::Common::Workflow::copyWorkflow    Could not copy to '$targetworkflow' }";
    }

    print "{ status: 'Completed copy to $targetworkflow' }";
}


#### COPY DIRECTORY
sub copyFilesystem {
    my $self    =   shift;
    my $source  =   shift;
    my $target  =   shift;

	require File::Copy::Recursive;
    my $result = File::Copy::Recursive::rcopy($source, $target);

    return $result;
}


=head2

	SUBROUTINE		copyProject

	PURPOSE

       THIS SUBROUTINE HAS TWO ROLES:

       1. COPY A PROJECT TO A (NON-EXISTING) DESTINATION PROJECT:

			1. ADD PROJECT TO project TABLE

            2. ADD ANY WORKFLOWS TO THE workflow TABLE

			2. OPTIONALLY, COPY THE PROJECT DIRECTORY

echo '{"sourceuser":"admin","targetuser":"syoung","sourceproject":"Project1","targetproject":"Project1","username":"syoung","sessionId":"9999999999.9999.999","mode":"copyProject"}' |  ./workflow.cgi

=cut

sub copyProject {
	my $self		=	shift;
    my $dbobject    =	$self->dbobject();
    my $json 		=	$self->json();

    print "{ error: 'Agua::Common::Workflow::copyProject    No data provided' }"
        if not defined $json;

	my $sourceuser = $json->{'sourceuser'};
	my $targetuser = $json->{'targetuser'};
	my $sourceproject = $json->{'sourceproject'};
	my $targetproject = $json->{'targetproject'};
	my $copyfiles = $json->{'copyfiles'};

    print "{ error: 'Agua::Common::Workflow::copyProject    User not validated: $targetuser' }"
        and exit if not $self->validate();

    my $can_copy = $self->projectPrivilege($sourceuser, $sourceproject, $targetuser, "groupcopy");
    print "{ error: 'Insufficient privileges for user: $targetuser' }" and exit if not $can_copy;

	#### CONFIRM THAT SOURCE PROJECT EXISTS IN project TABLE
	my $query = qq{SELECT description, notes
FROM project
WHERE username = '$sourceuser'
AND name = '$sourceproject'};
	my $details = $dbobject->queryhash($query);
    print "{ error: 'Agua::Common::Workflow::copyProject    Source project does not exist: $targetproject' }"
        and exit if not defined $details;
    my $description = $details->{description};
    my $notes = $details->{notes};

	#### EXIT IF TARGET PROJECT ALREADY EXISTS IN project TABLE
	$query = qq{SELECT 1
FROM project
WHERE username = '$targetuser'
AND name = '$targetproject'};
	my $exists = $dbobject->query($query);
    print "{ error: 'Agua::Common::Workflow::copyProject    Project already exists: $targetproject ' }"
        and exit if $exists;

    #### ADD PROJECT TO project TABLE
	$query = qq{INSERT INTO project
VALUES ('$targetuser', '$targetproject', '$description', '$notes')};
    my $insert_success = $dbobject->do($query);
    print "{ error: 'Agua::Common::Workflow::copyProject    Could not insert project: $targetproject' }"
        and exit if not $insert_success;

    #### GET SOURCE WORKFLOW INFORMATION
    $query = qq{SELECT * FROM workflow
WHERE username='$sourceuser'
AND project='$sourceproject'};
	my $workflowObjects = $dbobject->queryhasharray($query);

	#### COPY SOURCE WORKFLOW TO TARGET WORKFLOW
    #### COPY ALSO STAGE, STAGEPARAMETER, VIEW AND REPORT INFO
    foreach my $workflowObject ( @$workflowObjects )
    {
        $self->_copyWorkflow($sourceuser, $targetuser, $sourceproject, $targetproject, $workflowObject->{name}, $workflowObject->{name}, $workflowObject);
    }

	#### CREATE PROJECT DIRECTORY
	my $aguadir = $self->conf()->getKeyValue("agua", 'AGUADIR');
	my $userdir = $self->conf()->getKeyValue("agua", 'USERDIR');
	my $targetdir = "$userdir/$targetuser/$aguadir/$targetproject";
	File::Path::mkpath($targetdir);

    #### COPY FILES AND SUBDIRS IF FLAGGED BY copyfiles
    if ( defined $copyfiles and $copyfiles )
    {
        #### SET DIRECTORIES
        my $sourcedir = "$userdir/$sourceuser/$aguadir/$sourceproject";

        #### COPY DIRECTORY
        my $copy_result = $self->copyFilesystem($sourcedir, $targetdir);


        #### DEBUG    
        print `echo 'test' > $targetdir/apache.txt`;    


		print "{ status: 'Agua::Common::Workflow::copyProject    Copied to $targetproject' }"
            and exit if $copy_result;
        print "{ status: 'Agua::Common::Workflow::copyProject    Could not copy to '$targetproject' }";
    }

    print "{ status: 'Agua::Common::Workflow::copyProject    Copied to $targetproject' }";    
}


sub _copyWorkflow {
    my $self            =   shift;
    my $sourceuser      =   shift;
    my $targetuser      =   shift;
    my $sourceproject   =   shift;
    my $targetproject   =   shift;
    my $sourceworkflow  =   shift;
    my $targetworkflow  =   shift;
    my $workflowobject  =   shift;

    print "Agua::Common::Workflow::_copyWorkflow    sourceuser: $sourceuser\n";
    print "Agua::Common::Workflow::_copyWorkflow    targetuser: $targetuser\n";
    print "Agua::Common::Workflow::_copyWorkflow    sourceproject: $sourceproject\n";
    print "Agua::Common::Workflow::_copyWorkflow    targetproject: $targetproject\n";
    print "Agua::Common::Workflow::_copyWorkflow    workflowobject: $workflowobject\n";
    print "Agua::Common::Workflow::_copyWorkflow    targetworkflow: $targetworkflow\n";

    my $dbobject    =	$self->dbobject();

	#### CREATE PROJECT DIRECTORY
	my $aguadir = $self->conf()->getKeyValue("agua", 'AGUADIR');
	my $userdir = $self->conf()->getKeyValue("agua", 'USERDIR');
	my $targetdir = "$userdir/$targetuser/$aguadir/$targetproject/$targetworkflow";
	File::Path::mkpath($targetdir);


    #### INSERT COPY OF WORKFLOW INTO TARGET 
    $workflowobject->{username} = $targetuser;
    $workflowobject->{project} = $targetproject;
    $workflowobject->{name} = $targetworkflow;
    $self->insertWorkflow($workflowobject);


	my $query;
    #### COPY STAGES
    $query = qq{SELECT * FROM stage
WHERE username='$sourceuser'
AND project='$sourceproject'
AND workflow='$sourceworkflow'};
    my $stages = $dbobject->queryhasharray($query);
    $stages = [] if not defined $stages;
    foreach my $stage ( @$stages )
    {
        $stage->{username} = $targetuser;
        $stage->{project} = $targetproject;
        $stage->{workflow} = $targetworkflow;
    }
    $self->insertStages($stages);

    #### COPY STAGE PARAMETERS
    foreach my $stage ( @$stages )
    {
        my $query = qq{SELECT * FROM stageparameter
WHERE username='$sourceuser'
AND project='$sourceproject'
AND workflow='$sourceworkflow'
AND appnumber='$stage->{number}'};
        my $stageparams = $dbobject->queryhasharray($query);
        $stageparams = [] if not defined $stageparams;
        foreach my $stageparams ( @$stageparams )
        {
            $stageparams->{username} = $targetuser;
            $stageparams->{project} = $targetproject;
	        $stageparams->{workflow} = $targetworkflow;
        }
        $self->insertStageParameters($stageparams);
    }

    #### COPY INFORMATION IN view TABLE
    $query = qq{SELECT * FROM view
WHERE username='$sourceuser'
AND project='$sourceproject'};
    my $views = $dbobject->queryhasharray($query);
    $views = [] if not defined $views;
    foreach my $view ( @$views )
    {
        $view->{username} = $targetuser;
        $view->{project} = $targetproject;
    }
    $self->insertViews($views); 
}

sub insertViews {
    my $self        =   shift;
    my $hasharray   =   shift;

	print "Agua::Common::Workflow::insertViews    Agua::Common::Workflow::insertViews(hasharray)\n";

	#### SET TABLE AND REQUIRED FIELDS	
    my $dbobject    =	$self->dbobject();
	my $table       =   "view";
	my $required_fields = ["username", "project"];
	my $inserted_fields = $dbobject->fields($table);

    foreach my $hash ( @$hasharray )
    {    
        #### CHECK REQUIRED FIELDS ARE DEFINED
        my $not_defined = $dbobject->notDefined($hash, $required_fields);
        print "{ error: 'Agua::Common::Workflow::insertViews    undefined values: @$not_defined' }" and exit if @$not_defined;

        #### DO ADD
        my $success = $self->_addToTable($table, $hash, $required_fields, $inserted_fields);	
    }
}

sub insertReports {
    my $self        =   shift;
    my $hasharray   =   shift;


	#### SET TABLE AND REQUIRED FIELDS	
    my $dbobject    =	$self->dbobject();
	my $table       =   "report";
	my $required_fields = ["username", "project", "name", "number"];
	my $inserted_fields = $dbobject->fields($table);

    foreach my $hash ( @$hasharray )
    {    
        #### CHECK REQUIRED FIELDS ARE DEFINED
        my $not_defined = $dbobject->notDefined($hash, $required_fields);
        print "{ error: 'Agua::Common::Workflow::insertReports    undefined values: @$not_defined' }" and exit if @$not_defined;

        #### DO ADD
        my $success = $self->_addToTable($table, $hash, $required_fields, $inserted_fields);	
    }
}

sub insertStageParameters {
    my $self        =   shift;
    my $stageparameters   =   shift;


	#### SET TABLE AND REQUIRED FIELDS	
    my $dbobject    =	$self->dbobject();
	my $table       =   "stageparameter";
	my $required_fields = ["username", "project", "workflow", "appname", "appnumber", "name"];
    my $inserted_fields = $dbobject->fields($table); 
    foreach my $stageparameter ( @$stageparameters )
    {    
        #### CHECK REQUIRED FIELDS ARE DEFINED
        my $not_defined = $dbobject->notDefined($stageparameter, $required_fields);
        print "{ error: 'Agua::Common::Workflow::insertStageParameters    undefined values: @$not_defined' }" and exit if @$not_defined;

        #### DO ADD
        my $success = $self->_addToTable($table, $stageparameter, $required_fields, $inserted_fields);	
    }
}

sub insertStages {
    my $self        =   shift;
    my $stages   =   shift;


	#### SET TABLE AND REQUIRED FIELDS	
    my $dbobject    =	$self->dbobject();
	my $table       =   "stage";
	my $required_fields = ["username", "project", "workflow", "number"];
    my $inserted_fields = $dbobject->fields($table);    

    foreach my $stage ( @$stages )
    {    
        #### CHECK REQUIRED FIELDS ARE DEFINED
        my $not_defined = $dbobject->notDefined($stage, $required_fields);
        print "{ error: 'Agua::Common::Workflow::insertStages    undefined values: @$not_defined' }" and exit if @$not_defined;

        #### DO ADD
        my $success = $self->_addToTable($table, $stage, $required_fields, $inserted_fields);	
    }
}

sub insertWorkflow {
    my $self            =   shift;
    my $workflowObject  =   shift;


	#### SET TABLE AND REQUIRED FIELDS	
    my $dbobject    =	$self->dbobject();
	my $table       =   "workflow";
	my $required_fields = ["username", "project", "name", "number"];
	my $inserted_fields = $dbobject->fields($table);

    #### CHECK REQUIRED FIELDS ARE DEFINED
    my $not_defined = $dbobject->notDefined($workflowObject, $required_fields);
    print "{ error: 'Agua::Common::Workflow::insertWorkflows    undefined values: @$not_defined' }" and exit if @$not_defined;

    #### DO ADD
    my $success = $self->_addToTable($table, $workflowObject, $required_fields, $inserted_fields);	
}






1;