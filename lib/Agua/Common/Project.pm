package Agua::Common::Project;
use Moose::Role;
use Moose::Util::TypeConstraints;

=head2

	PACKAGE		Agua::Common::Project

	PURPOSE

		PROJECT METHODS FOR Agua::Common

=cut


use Data::Dumper;

=head2

	SUBROUTINE		saveProject

	PURPOSE

		ADD A PROJECT TO THE project TABLE

=cut

sub saveProject {
	my $self		=	shift;


    my $dbobject        =	$self->dbobject();
    my $json 			=	$self->json();

	my $success = $self->_removeProject();
	print "{ status: 'Could not remove project $json->{project} from project table' }" if not $success;
	$success = $self->_addProject();	
	print "{ status: 'Successful insert of project $json->{project} into project table' }" if $success;
	exit;
}

=head2

	SUBROUTINE		addProject

	PURPOSE

		ADD A PROJECT TO THE project TABLE

=cut

sub addProject {
	my $self		=	shift;


    my $dbobject        =	$self->dbobject();
    my $json 			=	$self->json();


	my $success = $self->_addProject();	
	print "{ status: 'Successful insert of project $json->{project} into project table' }" if $success;
	exit;
}

=head2

	SUBROUTINE		_addProject

	PURPOSE

		INTERNAL USE ONLY: ADD A PROJECT TO THE project TABLE

=cut

sub _addProject {
	my $self		=	shift;


    my $dbobject        =	$self->dbobject();
    my $json 			=	$self->json();


    #### VALIDATE
    print "{ error: 'Agua::Common::Project::_addProject    User session not validated' }" and exit unless $self->validate();

	#### SET TABLE AND REQUIRED FIELDS	
	my $table = "project";
	my $required_fields = ["username", "name"];

	#### CHECK REQUIRED FIELDS ARE DEFINED
	my $not_defined = $dbobject->notDefined($json, $required_fields);
    print "{ error: 'Agua::Common::Project::_addProject    undefined values: @$not_defined' }" and exit if @$not_defined;

	#### REMOVE IF EXISTS ALREADY
	$self->_removeFromTable($table, $json, $required_fields);

	#### DO ADD
	my $success = $self->_addToTable($table, $json, $required_fields);	
 	print "{ error: 'Agua::Common::Project::_addProject    Could not add project $json->{project} into project $json->{project} in project table'}" and exit if not defined $success;

	#### ADD THE PROJECT DIRECTORY TO THE USER'S agua DIRECTORY
    my $username = $json->{'username'};
	my $fileroot = $self->getFileroot($username);
	my $filepath = "$fileroot/$json->{name}";

	if ( $^O =~ /^MSWin32$/ )   {   $filepath =~ s/\//\\/g;  }

	File::Path::mkpath($filepath);
	print "{ error: 'Agua::Common::Project::_addProject    Could not create the fileroot directory: $fileroot' }" and exit if not -d $filepath;

	return 1;
}
=head2

	SUBROUTINE		_removeProject

	PURPOSE

		REMOVE A PROJECT FROM THE project, workflow, groupmember, stage AND

		stageparameter TABLES, AND REMOVE THE PROJECT FOLDER AND DATA FILES

=cut

sub _removeProject {
	my $self		=	shift;


    my $dbobject        =	$self->dbobject();
    my $json 			=	$self->json();

    #### VALIDATE
    print "{ error: 'Agua::Common::Project::_removeProject    User session not validated' }" and exit unless $self->validate();

	#### CHECK REQUIRED FIELDS ARE DEFINED
	my $required_fields = ["username", "name"];
	my $not_defined = $dbobject->notDefined($json, $required_fields);
    print "{ error: 'Agua::Common::Project::_removeProject    undefined values: @$not_defined' }" and exit if @$not_defined;

	#### REMOVE FROM project
	my $table = "project";
	return $self->_removeFromTable($table, $json, $required_fields);
}
=head2

	SUBROUTINE		removeProject

	PURPOSE

		REMOVE A PROJECT FROM THE project, workflow, groupmember, stage AND

		stageparameter TABLES, AND REMOVE THE PROJECT FOLDER AND DATA FILES

=cut

sub removeProject {
	my $self		=	shift;


    my $dbobject        =	$self->dbobject();
    my $json 			=	$self->json();

    #### VALIDATE
    print "{ error: 'Agua::Common::Project::removeProject    User session not validated' }" and exit unless $self->validate();

	#### CHECK REQUIRED FIELDS ARE DEFINED
	my $required_fields = ["username", "name"];
	my $not_defined = $dbobject->notDefined($json, $required_fields);
    print "{ error: 'Agua::Common::Project::removeProject    undefined values: @$not_defined' }" and exit if @$not_defined;

	#### REMOVE FROM project TABLE
    my $success = $self->_removeProject();
    print "{ error: 'Agua::Common::Project::removeProject    Can't remove project' }"
        and exit if not $success;

	#### REMOVE FROM workflow
	my $table = "workflow";
	$required_fields = ["username", "project"];
	$success = $self->_removeFromTable($table, $json, $required_fields);
 	print "{ error: 'Agua::Common::Project::removeProject    Could not delete project $json->{project} from the $table table'}" and exit if not defined $success;

	#### REMOVE FROM stage
	$table = "stage";
	$success = $self->_removeFromTable($table, $json, $required_fields);
 	print "{ error: 'Agua::Common::Project::removeProject    Could not delete project $json->{project} from the $table table'}" and exit if not defined $success;

	#### REMOVE FROM stageparameter
	$table = "stageparameter";
	$success = $self->_removeFromTable($table, $json, $required_fields);
 	print "{ error: 'Agua::Common::Project::removeProject    Could not delete project $json->{project} from the $table table'}" and exit if not defined $success;

	#### REMOVE FROM groupmember
	$table = "groupmember";
	$success = $self->_removeFromTable($table, $json, $required_fields);
 	print "{ error: 'Agua::Common::Project::removeProject    Could not delete project $json->{project} from the $table table'}" and exit if not defined $success;

	#### REMOVE FROM clusters
	$table = "clusters";
	$success = $self->_removeFromTable($table, $json, $required_fields);
 	print "{ error: 'Agua::Common::Project::removeProject    Could not delete project $json->{project} from the $table table'}" and exit if not defined $success;

	#### REMOVE PROJECT DIRECTORY
    my $username = $json->{'username'};
	my $fileroot = $self->getFileroot($username);
	my $filepath = "$fileroot/$json->{project}";
	if ( $^O =~ /^MSWin32$/ )   {   $filepath =~ s/\//\\/g;  }

	print "{ error: 'Agua::Common::Project::removeProject    Cannot remove directory: $filepath' }" and exit if not File::Remove::rm(\1, $filepath);

	print "{ status: 'Agua::Common::Project::removeProject    Successfully deleted project $json->{project} from project $json->{project} in project table' }";

}	#### removeProject





=head2

    SUBROUTINE:     getProjects

    PURPOSE:

		RETURN AN ARRAY OF project HASHES

			E.G.:
			[
				{
				  'name' : 'NGS',
				  'desciption' : 'NGS analysis team',
				  'notes' : 'This project is for ...',
				},
				{
					...
			]

=cut


sub getProjects {
	my $self		=	shift;


    my $dbobject		=	$self->dbobject();
	my $json			=	$self->json();

    #### VALIDATE
    print "{ error: 'User session not validated' }" and exit unless $self->validate();

	#### GET PROJECTS FOR THIS USER
    my $username = $json->{'username'};
	my $query = qq{SELECT * FROM project
WHERE username='$username'
ORDER BY name};
	my $projects = $dbobject->queryhasharray($query);

	######## IF NO RESULTS:
	####	1. INSERT DEFAULT PROJECT INTO project TABLE
	####	2. CREATE DEFAULT PROJECT FOLDERS
	return $self->_defaultProject() if not defined $projects;

	return $projects;
}



=head2

	SUBROUTINE		_defaultProject

	PURPOSE

		1. INSERT DEFAULT PROJECT INTO project TABLE

		2. RETURN QUERY RESULT OF project TABLE

=cut

sub _defaultProject {
	my $self		=	shift;	



    #### GET DATABASE OBJECT
    my $dbobject     =	$self->dbobject();
    my $json         =	$self->json();

    #### VALIDATE    
    print "{ error: 'Agua::Common::Project::_defaultProject    User session not validated' }" and exit unless $self->validate();

	#### SET DEFAULT PROJECT
	$self->json()->{name} = "Project1";

	#### ADD PROJECT
	my $success = $self->_addProject();
	print "{ error: 'Agua::Common::Project::_defaultProject    Could not add project $json->{project} into  project table'}" and exit if not defined $success;

	#### DO QUERY
    my $username = $json->{'username'};
	my $query = qq{SELECT * FROM project
WHERE username='$username'
ORDER BY name};
	my $projects = $dbobject->queryhasharray($query);

	return $projects;
}


1;