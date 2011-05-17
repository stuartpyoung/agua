package Agua::Common::Base;
use Moose::Role;
use Moose::Util::TypeConstraints;

=head2

	PACKAGE		Agua::Common::Base

	PURPOSE

		BASE METHODS FOR Agua::Common

=cut




#### USE LIB FOR INHERITANCE
use FindBin qw($Bin);
use lib "$Bin/../../";
use lib "$Bin/../../external";

use Data::Dumper;

##############################################################################
#				DATA METHODS
=head2

	SUBROUTINE		getData

	PURPOSE

		RETURN JSON STRING OF ALL WORKFLOWS FOR EACH

		PROJECT IN THE application TABLE BELONGING TO

		THE USER

	INPUT

		1. USERNAME

		2. SESSION ID

	OUTPUT

		1. JSON HASH { "projects":[ {"project": "project1","workflow":"...}], ...}

=cut
sub getData {
	my $self		=	shift;	


#use Time::HiRes qw[gettimeofday tv_interval];
#my $time1=[gettimeofday()];	

    #### GET DATABASE OBJECT
    my $dbobject     =	$self->dbobject();
    my $json         =	$self->json();

    #### VALIDATE    
    print "{ error: 'User session not validated' }" and exit unless $self->validate();

	#### GET ALL USER PROJECTS
	#### ( TAKES 2 ms, COMPARED TO 12 ms USING MULTIPLE MYSQL QUERIES )
    my $username = $json->{username};

	my $output;

	#### GET ARRAY OF PROJECT HASHES: NAME, NUMBER AND DESCRIPTION 	
	$output->{projects} = $self->getProjects();

	#### CONVERT PROJECTS INTO ARRAY OF { projectName : [ workflows] } HASHES
	$output->{workflows} = $self->getWorkflows();

	#### GET ARRAY OF { groupName : [ members ] } HASHES
	$output->{groupmembers} = $self->getGroupMembers();

    #### GET HASH OF ADMIN MENU HEADINGS
    $output->{headings} = $self->getHeadings();

    #### GET ARRAY OF [ username: ..., email: ..., ... ] USER HASHES
    $output->{access} = $self->getAccess();

    #### GET ARRAY OF [ username: ..., email: ..., ... ] USER HASHES
    $output->{users} = $self->getUsers();

	#### GET ARRAY OF [ { groupname: ..., groupdesc: ..., ... ] GROUP HASHES
	$output->{groups} = $self->getGroups();

	#### GET ARRAY OF [ name: ..., description: ..., location: ... ] SOURCE HASHES
	$output->{sources} = $self->getSources();

	#### GET ARRAY OF [ name: ..., description: ..., location: ... ] SOURCE HASHES
	$output->{sharedsources} = $self->sharedSources();

	#### GET ARRAY OF [ name: ..., description: ..., ... ] APPLICATION HASHES
	$output->{apps} = $self->getApps();

	#### GET ARRAY OF [ name: ..., description: ..., ... ] APP PARAMETER KEY:VALUE PAIRS
	$output->{parameters} = $self->getParameters();

	#### GET ARRAY OF HASHES [ stagename1: [  ], stagename2 : [...], ...  ] 
	$output->{stages} = $self->getStages();

	#### GET ARRAY OF HASHES: [ { project: ..., workflow: ..., appname:, ..., name: ... }, ...] 
	$output->{stageparameters} = $self->getStageParameters();

	$output->{views} = $self->getViews();
	$output->{viewfeatures} = $self->getViewFeatures();
	$output->{features} = $self->getFeatures();

	#### GET AWS AUTHENTICATION INFO
	$output->{aws} = $self->getAws($username);

	#### GET (STAR)CLUSTERS
	$output->{clusters} = $self->getClusters();

	#### GET CLUSTER ALLOTED TO WORKFLOWS
	$output->{clusterworkflows} = $self->getClusterWorkflows();


	#### GET ARRAY OF SHARED PARAMETERS
	#### TO GO WITH SHARED APPS IN getApps
	$output->{sharedparameters} = $self->getSharedParameters();

	#### GET ARRAY OF **SHARED*** PROJECTS
	$output->{sharedprojects} = $self->getSharedProjects();

	#### GET ARRAY OF **SHARED*** STAGES
	$output->{sharedstages} = $self->sharedStages();

	#### GET ARRAY OF **SHARED*** STAGE PARAMETERS
	$output->{sharedstageparameters} = $self->sharedStageParameters();

	#### GET ARRAY OF **SHARED*** VIEWS
	$output->{sharedviews} = $self->sharedViews();

	#### GET ARRAY OF **SHARED*** VIEW FEATURES
	$output->{sharedviewfeatures} = $self->sharedViewFeatures();

	#### PRINT JSON AND EXIT
	use JSON -support_by_pp; 

	my $jsonParser = JSON->new();
    #my $jsonText = $jsonParser->objToJson($output, {pretty => 1, indent => 4});
    my $jsonText = $jsonParser->objToJson($output, {pretty => 0});

	#### THIS ALSO WORKS ON LINUX
	#my $jsonText = $jsonParser->encode->allow_nonref->pretty->get_utf8->($output);
	####my $apps = $jsonParser->allow_singlequote(1)->allow_nonref->loose(1)->encode($output);

	#### TO AVOID HIJACKING --- DO NOT--- PRINT AS 'json-comment-optional'
	print "{}&&$jsonText";
	exit;
}


=head2

	SUBROUTINE		getTable

	PURPOSE

		RETURN THE JSON STRING FOR ALL USER-RELATED ENTRIES IN

        THE DESIGNATED TABLE PROXY NAME (NB: ACTUAL TABLE NAMES

        DIFFER -- SOME ARE MISSING THE LAST 'S')

	INPUT

		1. USERNAME

		2. SESSION ID

        3. TABLE PROXY NAME, E.G., "stageParameters" RETURNS RELATED

            'stageparameter' TABLE ENTRIES

	OUTPUT

		1. JSON HASH { "projects":[ {"project": "project1","workflow":"...}], ...}

=cut

sub getTable {
	my $self		=	shift;	


    #### GET DATABASE OBJECT
    my $dbobject     =	$self->dbobject();
    my $json         =	$self->json();

    #### VALIDATE    
    print "{ error: 'User session not validated' }" and exit unless $self->validate();

    my $username = $json->{username};
    my $table = $json->{table};

    #### CONVERT TO get... COMMAND
    my $get = "get" . $self->cowCase($table);

    #### QUIT IF TABLE PROXY NAME IS INCORRECT
    print "{ error: 'Agua::Common::Base::getTable    method $get not defined' }" and exit unless defined $self->can($get);

    my $output = $self->$get();

    my $data = {};
    $data->{lc($table)} = $output;

	#### PRINT JSON AND EXIT
	use JSON -support_by_pp; 

	my $jsonParser = JSON->new();
    my $jsonText = $jsonParser->objToJson($data, {pretty => 1, indent => 4});
    #my $jsonText = $jsonParser->objToJson($output, {pretty => 0});

	#### THIS ALSO WORKS ON LINUX
	#my $jsonText = $jsonParser->encode->allow_nonref->pretty->get_utf8->($output);
	####my $apps = $jsonParser->allow_singlequote(1)->allow_nonref->loose(1)->encode($output);

	#### TO AVOID HIJACKING --- DO NOT--- PRINT AS 'json-comment-optional'
	print "{}&&$jsonText";
	exit;
}



##############################################################################
#				DATABASE TABLE METHODS
=head2

	SUBROUTINE		_updateTable

	PURPOSE

		UPDATE ONE OR MORE ENTRIES IN A TABLE

	INPUTS

		1. NAME OF TABLE      

		2. HASH CONTAINING OBJECT TO BE UPDATED

		3. HASH CONTAINING TABLE FIELD KEY-VALUE PAIRS

=cut
sub _updateTable {
	my $self		=	shift;
	my $table		=	shift;
	my $hash		=	shift;
	my $required_fields		=	shift;
	my $set_hash	=	shift;
	my $set_fields	=	shift;


    print "{ error: 'Agua::Common::Base::_updateTable    hash not defined' }" and exit if not defined $hash;
    print "{ error: 'Agua::Common::Base::_updateTable    required_fields not defined' }" and exit if not defined $required_fields;
    print "{ error: 'Agua::Common::Base::_updateTable    set_hash not defined' }" and exit if not defined $set_hash;
    print "{ error: 'Agua::Common::Base::_updateTable    set_fields not defined' }" and exit if not defined $set_fields;
    print "{ error: 'Agua::Common::Base::_updateTable    table not defined' }" and exit if not defined $table;

    my $dbobject        =	$self->dbobject();

	#### CHECK REQUIRED FIELDS ARE DEFINED
	my $not_defined = $dbobject->notDefined($hash, $required_fields);
    print "{ error: 'Agua::Common::Base::_updateTable    undefined values: @$not_defined' }" and exit if @$not_defined;

	#### GET WHERE
	my $where = $dbobject->where($hash, $required_fields);

	#### GET SET
	my $set = $dbobject->set($set_hash, $set_fields);
	print "{ error: 'Agua::Common::Base::_updateTable    set values not defined' }" and exit if not defined $set;

	##### UPDATE TABLE
	my $query = qq{UPDATE $table $set $where};           
	my $result = $dbobject->do($query);
}



=head2

	SUBROUTINE		_addToTable

	PURPOSE

		ADD AN ENTRY TO A TABLE

	INPUTS

		1. NAME OF TABLE      

		2. ARRAY OF KEY FIELDS THAT MUST BE DEFINED 

		3. HASH CONTAINING TABLE FIELD KEY-VALUE PAIRS

=cut

sub _addToTable {
	my $self		=	shift;
	my $table		=	shift;
	my $hash		=	shift;
	my $required_fields		=	shift;
	my $inserted_fields		=	shift;


	#### CHECK FOR ERRORS
    print "{ error: 'Agua::Common::Base::_addToTable    hash not defined' }" and exit if not defined $hash;
    print "{ error: 'Agua::Common::Base::_addToTable    required_fields not defined' }" and exit if not defined $required_fields;
    print "{ error: 'Agua::Common::Base::_addToTable    table not defined' }" and exit if not defined $table;

    my $dbobject        =	$self->dbobject();

	#### CHECK REQUIRED FIELDS ARE DEFINED
	my $not_defined = $dbobject->notDefined($hash, $required_fields);
    print "{ error: 'Agua::Common::Base::_addToTable    undefined values: @$not_defined' }" and exit if @$not_defined;

	#### GET ALL FIELDS BY DEFAULT IF INSERTED FIELDS NOT DEFINED
	$inserted_fields = $dbobject->fields($table) if not defined $inserted_fields;

	print "{ error: 'Agua::Common::Base::_addToTable    fields not defined' }" and exit if not defined $inserted_fields;
	my $fields_csv = join ",", @$inserted_fields;

	##### INSERT INTO TABLE
	my $values_csv = $dbobject->fields_csvline($inserted_fields, $hash);
	my $query = qq{INSERT INTO $table ($fields_csv)
VALUES ($values_csv)};           
	my $result = $dbobject->do($query);

	return $result;
}




=head2

	SUBROUTINE		_removeFromTable

	PURPOSE

		REMOVE AN ENTRY FROM A TABLE

	INPUTS

		1. HASH CONTAINING TABLE FIELD KEY-VALUE PAIRS

		2. ARRAY OF KEY FIELDS THAT MUST BE DEFINED 

		3. NAME OF TABLE      
=cut

sub _removeFromTable {
	my $self		=	shift;	
	my $table		=	shift;
	my $hash		=	shift;
	my $required_fields		=	shift;



    #### CHECK INPUTS
    print "{ error: 'Agua::Common::Base::_removeFromTable    hash not defined' }" and exit if not defined $hash;
    print "{ error: 'Agua::Common::Base::_removeFromTable    required_fields not defined' }" and exit if not defined $required_fields;
    print "{ error: 'Agua::Common::Base::_removeFromTable    table not defined' }" and exit if not defined $table;

    my $dbobject        =	$self->dbobject();

	#### CHECK REQUIRED FIELDS ARE DEFINED
	my $not_defined = $dbobject->notDefined($hash, $required_fields);
    print "{ error: 'Agua::Common::Base::_removeFromTable    undefined values: @$not_defined' }" and exit if @$not_defined;

	#### DO DELETE 
	my $where = $dbobject->where($hash, $required_fields);

	my $query = qq{DELETE FROM $table
$where};
	my $result = $dbobject->do($query);

	return 1 if defined $result;
	return 0;
}


=head2

    SUBROUTINE:     arrayToArrayhash

    PURPOSE:

        CONVERT AN ARRAY INTO AN ARRAYHASH, E.G.:

		{
			key1 : [ entry1, entry2 ],
			key2 : [ ... ]
			...
		}



=cut

sub arrayToArrayhash {
	my $self		=	shift;
	my $array 		=	shift;
	my $key			=	shift;


	my $arrayhash = {};
	for my $entry ( @$array )
	{
		if ( not defined $entry->{$key} )
		{
			print "Agua::Common::Base::arrayToArrayhash    entry->{$key} not defined in entry. Returning.\n";
			return;
		}
		$arrayhash->{$entry->{$key}} = [] if not exists $arrayhash->{$entry->{$key}};
		push @{$arrayhash->{$entry->{$key}}}, $entry;		
	}

	return $arrayhash;
}





##############################################################################
#				FILESYSTEM METHODS
=head2

	SUBROUTINE		getFileroot

	PURPOSE

		RETURN THE FULL PATH TO THE agua FOLDER WITHIN THE USER'S HOME DIR

=cut

sub getFileroot {
	my $self		=	shift;
	my $username	=	shift;


	#### RETURN FILE ROOT FOR THIS USER IF ALREADY DEFINED
	return $self->fileroot() if $self->fileroot() and not defined $username;

	#### OTHERWISE, GET USERNAME FROM JSON IF NOT PROVIDED
	$username = $self->json()->{username} if not defined $username;
	return if not defined $username;

	my $userdir = $self->conf()->getKeyValue('agua', 'USERDIR');
	my $aguadir = $self->conf()->getKeyValue('agua', 'AGUADIR');
	my $fileroot = "$userdir/$username/$aguadir";

	$self->fileroot($fileroot);

	return $fileroot;	
}


1;

