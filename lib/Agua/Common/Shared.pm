package Agua::Common::Shared;

=head2

	PACKAGE		Agua::Common::Shared

	PURPOSE

		PROJECT AND RESOURCE SHARING METHODS FOR Agua::Common

=cut


use Moose::Role;
use Moose::Util::TypeConstraints;


has 'totalprojects'		=> ( isa => 'ArrayRef|Undef', is => 'rw', default => undef );
has 'sharedprojects'	=> ( isa => 'ArrayRef|Undef', is => 'rw', default => undef );
has 'worldprojects'		=> ( isa => 'ArrayRef|Undef', is => 'rw', default => undef );
has 'userprojects'		=> ( isa => 'HashRef|Undef', is => 'rw', default => undef );

use Data::Dumper;

=head2

    SUBROUTINE:     getSharedParameters

    PURPOSE:

        RETURN A JSON LIST OF PARAMETERS TO ACCOMPANY THE SHARED

		APPS RETURNED BY getApps
=cut

sub getSharedParameters {

	my $self		=	shift;

    my $dbobject		=	$self->dbobject();
	my $json			=	$self->json();

	#### VALIDATE    
    my $username = $json->{username};
    print "{ error: 'Agua::Common::App::getSharedParameters    User not validated: $username ' }"
		and exit unless $self->validate();

	#### GET admin USER'S APPS	
	my $admin = $self->conf()->getKeyValue("agua", 'ADMINUSER');
    my $query = qq{SELECT * FROM parameter
WHERE owner = '$admin'
ORDER BY appname, apptype};
    my $apps = $dbobject->queryhasharray($query);
	$apps = [] if not defined $apps;

	return $apps;
}


sub getUserProjects {
    my $self        =   shift;
	my $username	=	shift;


=head2

    SUBROUTINE     sharedViewFeatures

    PURPOSE

        1. RETURN A USERNAME:VIEWS HASH ARRAY CONTAINING ALL VIEWS

			(GROUP AND WORLD) SHARED WITH THIS USER

        2. EACH 

				A. VIEW - CAN BEEN SEEN FROM SHARER'S HOME DIRECTORY

				B. COPY - CAN BE COPIED TO THE USER'S HOME DIRECTORY 

				C. READ - BY DEFINITION, A SHARED VIEW HAS READ PERMISSION

echo '{"sourceuser":"admin","targetuser":"syoung","sourceworkflow":"Workflow0","sourceproject":"Project1","targetworkflow":"Workflow9","targetproject":"Project1","username":"syoung","sessionId":"9999999999.9999.999","mode":"getData"}' |  ./workflow.cgi

=cut

sub sharedViewFeatures {
    my $self        =   shift;


    #### VALIDATE
    my $validated = $self->validate();
    print "{ error: 'Agua::Common::Shared::sharedViewFeatures    User session not validated' }" and exit unless $validated;

    #### GET DATABASE OBJECT
    my $dbobject = $self->dbobject();

    #### GET A USERNAME:PROJECTS HASH ARRAY OF PROJECTS
	#### SHARED WITH THIS USER 
    my $username = $self->json->{username};
	my $userprojects = $self->getUserProjects($username);
	return [] unless defined $userprojects and %$userprojects;

	#### GET THE VIEWS FOR THESE PROJECTS
	my $sharedviewfeatures = {};
	foreach my $username ( keys %$userprojects )
	{
		my $projects = $userprojects->{$username};
		next if scalar @$projects == 0;

		my $select = "username,project,view,feature,type";
		my $extra = "ORDER BY project, view, feature,type";
		my $shared_data = $self->_sharedProjectData($projects, "viewfeature", $select, $extra);
		$sharedviewfeatures->{$username} = $shared_data if defined $shared_data;
	}

	return $sharedviewfeatures;
}



    #### VALIDATE
    print "{ error: 'Agua::Common::Shared::getUserProjects    User session not validated' }"
		and exit unless $self->validate();

    #### GET PROJECTS SHARED WITH THE USER
    my $shared_projects = $self->getSharedProjects($username);	

    #### GET PROJECTS SHARED WITH THE USER
    my $world_projects = $self->getWorldProjects($username);	

	#### ADD GROUP SHARED TO WORLD SHARED
	my $total_projects;
	@$total_projects = (@$shared_projects, @$world_projects);

	#### CONVERT TO USERNAME:PROJECTS HASH
	my $userprojects = $self->hasharrayToHash($total_projects, "owner");

	return $userprojects;
}


=head2

    SUBROUTINE     sharedViews

    PURPOSE

        1. RETURN A USERNAME:VIEWS HASH ARRAY CONTAINING ALL VIEWS

			(GROUP AND WORLD) SHARED WITH THIS USER

        2. EACH 

				A. VIEW - CAN BEEN SEEN FROM SHARER'S HOME DIRECTORY

				B. COPY - CAN BE COPIED TO THE USER'S HOME DIRECTORY 

				C. READ - BY DEFINITION, A SHARED VIEW HAS READ PERMISSION

echo '{"sourceuser":"admin","targetuser":"syoung","sourceworkflow":"Workflow0","sourceproject":"Project1","targetworkflow":"Workflow9","targetproject":"Project1","username":"syoung","sessionId":"9999999999.9999.999","mode":"getData"}' |  ./workflow.cgi

=cut

sub sharedViews {
    my $self        =   shift;


    #### VALIDATE
    my $validated = $self->validate();
    print "{ error: 'Agua::Common::Shared::sharedViews    User session not validated' }"
		and exit unless $validated;

    #### GET DATABASE OBJECT
    my $dbobject = $self->dbobject();

    #### GET A USERNAME:PROJECTS HASH ARRAY OF PROJECTS
	#### SHARED WITH THIS USER 
    my $username = $self->json->{username};
	my $userprojects = $self->getUserProjects($username);
	return [] unless defined $userprojects and %$userprojects;

	#### GET THE VIEWS FOR THESE PROJECTS
	my $sharedviews = {};
	foreach my $username ( keys %$userprojects )
	{
		my $projects = $userprojects->{$username};
		next if scalar @$projects == 0;

		my $select = "*";
		my $extra = "ORDER BY project, view, species, build, chromosome, start, stop";
		my $shared_data = $self->_sharedProjectData($projects, "view", $select, $extra);
		$sharedviews->{$username} = $shared_data if defined $shared_data;
	}

	return $sharedviews;
}


sub addProjectRights {
	my $self		=	shift;
	my $hasharray	=	shift;
	my $project		=	shift;

	#### ADD THE PERMISSIONS FOR THIS PROJECT TO EACH VIEW

	foreach my $entry ( @$hasharray )
	{
		my $permissions = ["groupwrite", "groupcopy", "groupview", "worldwrite", "worldcopy", "worldview"];
		foreach my $permission ( @$permissions )
		{
			$entry->{$permission} = $project->{rights}->{$permission}
				if defined $project->{rights}->{$permission};
		}
	}

	return $hasharray;	
}

=head2

	SUBROUTINE		sharedStages

	PURPOSE

		RETURN ARRAY OF STAGES THAT HAVE BEEN SHARED WITH

		THE USER BY OTHER USERS

	INPUT

		1. USERNAME

		2. SESSION ID

	OUTPUT

		1. USER-KEYED HASH:

			{
				USER1 : [
					{ project: project1, workflow: workflow1, name: app1, ...},
					{ .. }
				],
				USER2 : [ ... ],
				...
			}

=cut

sub sharedStages  {
	my $self		=	shift;	



	my $jsonParser = $self->json_parser();	

	#### GET JSON
    my $json         =	$self->json();

    #### VALIDATE    
    print "{ error: 'User session not validated' }" and exit unless $self->validate();

	#### GET THE ARRAY OF SHARED PROJECTS CONVERTED TO
	#### A USERNAME:PROJECT KEYED HASH
    my $username = $json->{username};
	my $userprojects = $self->getUserProjects($username);
	return [] unless defined $userprojects and %$userprojects;

	#### GET THE STAGES FOR EACH SHARED PROJECT
	my $shared_stages = {};
	foreach my $username ( keys %$userprojects )
	{
		my $projects = $userprojects->{$username};

		my $shared_data = $self->_sharedProjectData($projects, "stage");
		$shared_stages->{$username} = $shared_data if defined $shared_data;
	}

	$shared_stages = [] if not defined $shared_stages;

	return $shared_stages;
}


=head2

	SUBROUTINE		_sharedProjectData

	PURPOSE

		RETURN ARRAY OF DATA VALUES FOR THE SHARED PROJECTS

		WITH PERMISSIONS ATTACHED TO EACH ENTRY

	INPUT

		1. PROJECTS HASH (username|owner, project, etc.)

		2. TABLE FROM WHICH TO GET THE DATA

	OUTPUT

		1. USER-KEYED HASH:

			{
				USER1 : [
					{ project: project1, workflow: workflow1, appname: app1, ...},
					{ .. }
				],
				USER2 : [ ... ],
				...
			}

=cut

sub _sharedProjectData {
	my $self		=	shift;	
	my $projects	=	shift;
	my $table		=	shift;
	my $select		=	shift;
	my $extra		=	shift;


	#### DEFAULTS
	$select = '*' if not defined $select;
	$extra = '' if not defined $extra;

    #### GET DATABASE OBJECT
    my $dbobject     =	$self->dbobject();

    #### VALIDATE    
    print "{ error: 'User session not validated' }" and exit unless $self->validate();

	#### 1. RETRIEVE USER PROJECTS (STORED EARLIER IN sharedStages)
	my $user_projects = $self->userprojects();

	#### GET THE DATA
	my $shared_project_stages = [];

	foreach my $project ( @$projects )
	{
		#### ADD username FIELD
		$project->{username} = $project->{owner};


		my $unique_keys = ["username", "project"];
		my $where = $dbobject->where($project, $unique_keys);
		my $query = qq{SELECT $select FROM $table
$where $extra};
		my $hasharray = $dbobject->queryhasharray($query);
		next if not defined $hasharray;

		#### ADD THE PERMISSIONS FOR THIS PROJECT TO EACH ENTRY
		$hasharray = $self->addProjectRights($hasharray, $project);

		(@$shared_project_stages) = (@$shared_project_stages, @$hasharray);
	}


	return $shared_project_stages;
}




=head2

	SUBROUTINE		sharedStageParameters

	PURPOSE

		RETURN ARRAY OF STAGE PARAMETERS THAT HAVE BEEN SHARED WITH

		THE USER BY OTHER USERS

	INPUT

		1. USERNAME

		2. SESSION ID

	OUTPUT

		1. USER-KEYED HASH:

			{
				USER1 : [
					{ project: project1, workflow: workflow1, appname: app1, ...},
					{ .. }
				],
				USER2 : [ ... ],
				...
			}

=cut

sub sharedStageParameters  {
	my $self		=	shift;



    #### VALIDATE    
    print "{ error: 'User session not validated' }" and exit unless $self->validate();

	#### GET JSON
    my $json         =	$self->json();

	#### STORE THE USER-KEYED SHARED PROJECTS HERE
	my $shared_stageparameters = {};

	#### GET THE ARRAY OF SHARED PROJECTS CONVERTED TO
	#### A USERNAME:PROJECT KEYED HASH
    my $username = $json->{username};
	my $userprojects = $self->getUserProjects($username);
	return [] unless defined $userprojects and %$userprojects;

	#### GET THE STAGE PARAMETERS FOR EACH SHARED PROJECT
	foreach my $username ( keys %$userprojects )
	{
		my $projects = $userprojects->{$username};

		my $shared_data = $self->_sharedProjectData($projects, "stageparameter");
		$shared_stageparameters->{$username} = $shared_data if defined $shared_data;
	}


	return $shared_stageparameters;
}

=head2

    SUBROUTINE     getSharedProjects

    PURPOSE

		RETURN THE LIST OF PROJECTS BELONGING TO GROUPS OF WHICH

		THE USER IS A MEMBER

	NOTES

		RETURNS THE FOLLOWING DATA STRUCTURE:

		$sharedprojects = [
			  {
				'rights' => {
							  'groupview' => '1',
							  'groupcopy' => '1',
							  'groupwrite' => '1'
							},
				'owner' => 'admin',
				'groupname' => 'bioinfo',
				'project' => 'Project1',
				'description' => ''
			  },
			  {
				'rights' => {
							  'groupview' => '1',
							  'groupcopy' => '1',
							  'groupwrite' => '1'
							},
				'owner' => 'admin',
				'groupname' => 'bioinfo',
				'project' => 'ProjectX',
				'description' => ''
			  },
			  ...
		]

=cut

sub getSharedProjects {
    my $self        =   shift;
    my $username	=   shift;


	return $self->sharedprojects() if $self->sharedprojects();

    #### GET USERNAME IF NOT DEFINED
    my $json = $self->json();
	$username = $json->{username} if not defined $username;

    #### GET DATABASE OBJECT
    my $dbobject = $self->dbobject();

	####  1. GET USER-ACCESSIBLE SHARED PROJECTS
	my $query = qq{SELECT DISTINCT owner, groupname, description
FROM groupmember
WHERE name='$username'
AND type='user'
ORDER BY owner};
    my $ownergroups = $dbobject->queryhasharray($query);

	my $shared_projects;
	my $already_seen;
    foreach my $ownergroup( @$ownergroups )
    {
        my $owner = $ownergroup->{owner};
        my $group = $ownergroup->{groupname};

        $query = qq{select groupname, name, description from groupmember where owner = '$owner' and groupname = '$group' and type='project'};
        my $projects = $dbobject->queryhasharray($query);
        foreach my $project ( @$projects )
        {
			#### SKIP IT IF WE'VE SEEN IT BEFORE
			next if $already_seen->{$project->{name}};
			$already_seen->{$project->{name}} = 1;

			my $hash;
            $hash->{groupname} = $project->{groupname};
            $hash->{project} = $project->{name};
            $hash->{description} = $project->{description};
            $hash->{owner} = $owner;
            push @$shared_projects, $hash;
        }
    }


    #### GET THE PERMISSIONS OF THE SHARED PROJECTS
    for my $shared_project ( @$shared_projects )
    {
        my $owner = $shared_project->{owner};
        my $groupname = $shared_project->{groupname};

        my $query = qq{select groupwrite, groupcopy, groupview from access where owner = '$owner' and groupname = '$groupname'};
        my $grouprights = $dbobject->queryhash($query);


		$grouprights = {} if not defined $grouprights;        
        $shared_project->{rights} = $grouprights;
    }
	$shared_projects = [] if not defined $shared_projects;


	#### SET _sharedprojects
	$self->sharedprojects($shared_projects);

	return $shared_projects;	
}

=head2

    SUBROUTINE     getWorldProjects

    PURPOSE

		RETURN THE LIST OF WORLD ACCESSIBLE PROJECTS ___MINUS___ THE

		PROJECTS BELONGING TO GROUPS OF WHICH THE USER IS A MEMBER

=cut

sub getWorldProjects {
    my $self        		=   shift;
    my $username			=   shift;
    my $shared_projects 	=   shift;


	$shared_projects = $self->getSharedProjects($username)
		if not defined $shared_projects;

    #### GET DATABASE OBJECT
    my $dbobject = $self->dbobject();

    #### GET GROUPS SHARED WITH WORLD
    my $query = qq{SELECT DISTINCT owner, groupname, worldwrite, worldcopy, worldview
FROM access
WHERE worldview = 1
AND owner != "$username"};
    my $world_ownergroups = $dbobject->queryhasharray($query);
	$world_ownergroups = [] if not defined $world_ownergroups;

    #### GET PROJECTS IN WORLD GROUPS
    my $world_projects = [];
    foreach my $world_ownergroup ( @$world_ownergroups )
    {
        my $owner = $world_ownergroup->{owner};
        my $group = $world_ownergroup->{groupname};
        my $query = qq{SELECT owner, groupname, name AS project, description 
FROM groupmember
WHERE owner='$owner'
AND groupname='$group'
AND type='project'};
        my $projects = $dbobject->queryhasharray($query);
        $projects = [] if not defined $projects;
		@$world_projects = (@$world_projects, @$projects) if @$projects;
    }

	#### REMOVE ANY WORLD PROJECTS THAT HAVE ALREADY BEEN SHARED WITH THIS USER
	if ( defined $world_projects and @$world_projects )
	{
		for ( my $i = 0; $i < @$world_projects; $i++ )
		{
			my $world_owner = $$world_projects[$i]->{owner};
			my $world_project = $$world_projects[$i]->{project};
			foreach my $shared_project ( @$shared_projects )
			{
				my $shared_owner = $shared_project->{owner};
				my $shared_project = $shared_project->{project};

				if ( $shared_owner eq $world_owner
					&& $shared_project eq $world_project)
				{
					splice @$world_projects, $i, 1;
					$i--;
				}
			}
		}
	}

    #### GET THE PERMISSIONS OF THE SHARED PROJECTS
    for my $world_project ( @$world_projects )
    {
        my $owner = $world_project->{owner};
        my $groupname = $world_project->{groupname};

        my $query = qq{SELECT worldwrite, worldcopy, worldview
FROM access
WHERE owner = '$owner'
AND groupname = '$groupname'};
        my $grouprights = $dbobject->queryhash($query);		

		$grouprights = {} if not defined $grouprights;        
        $world_project->{rights} = $grouprights;
    }
	$world_projects = [] if not defined $world_projects;


	#### FOR PROJECTS THAT BELONG TO MULTIPLE GROUPS, SELECT THE
	#### GROUP IN WHICH THIS USER HAS THE GREATEST PRIVILEGES    
	$world_projects = $self->highestPrivilegeProjects($world_projects);


	return $world_projects;
}


=head2

	SUBROUTINE		highestPrivilegeProjects

	PURPOSE

		FOR PROJECTS THAT BELONG TO MULTIPLE GROUPS, SELECT THE

		GROUP IN WHICH THIS USER HAS THE GREATEST PRIVILEGES    
=cut

sub highestPrivilegeProjects {
	my $self		=	shift;
	my $projects	=	shift;

	print "Agua::Common::Shared::highestPrivilegeProjects    projects not defined\n"
		and exit if not defined $projects;

	for ( my $i = 0; $i < @$projects; $i++ )
	{
		my $project = $$projects[$i];
		my $privilege = $self->sumPrivilege($project->{rights}, "world");

		for ( my $j = 0; $j < @$projects; $j++ )
		{
			my $current_project = $$projects[$j];
			my $current_privilege = $self->sumPrivilege($current_project->{rights}, "world");

			if ( $current_privilege <= $privilege
			and $project->{owner} eq $current_project->{owner}
			and $project->{project} ne $current_project->{project} )
			{

				splice @$projects, $i, 1;
				$i--;
			}
		}
	}

	return $projects;
}





1;
