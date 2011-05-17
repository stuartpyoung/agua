#package Agua::Admin;
#use Moose;
use MooseX::Declare;

=head2

		PACKAGE		Admin

		PURPOSE

			THE Admin OBJECT PERFORMS THE FOLLOWING	TASKS:

				1. AUTHENTICATE USER ACCESS (PASSWORD AND SESSION ID)

				2. CREATE, MODIFY OR DELETE USERS

=cut


use strict;
use warnings;
use Carp;

#### USE LIB FOR INHERITANCE
use FindBin qw($Bin);
use lib "$Bin/../";
use lib "$Bin/../external";
use Data::Dumper;

#### INTERNAL MODULES
#use Agua::Common;

class Agua::Admin with (Agua::Common::Util,
	Agua::Common::Base,
	Agua::Common::Cluster,
	Agua::Common::Workflow,
	Agua::Common::Privileges,
	Agua::Common::User)
{

# STRINGS
has 'fileroot'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'validated'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'username'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'project'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'workflow'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
# OBJECTS
has 'json'		=> ( isa => 'HashRef|Undef', is => 'rw', default => undef );
has 'dbobject'	=> ( isa => 'Agua::DBase::MySQL', is => 'rw', required => 0 );
has 'conf' 	=> (
	is =>	'rw',
	'isa' => 'Conf::Agua',
	default	=>	sub { Conf::Agua->new(	backup	=>	1, separator => "\t"	);	}
);

### /////}

method BUILD ($hash) {

	#### SET DATABASE HANDLE
	$self->setDbh();	

    #### VALIDATE
    print "{ error: 'Admin::Project::BUILD    User session not validated' }" and exit unless $self->validate();

}

##############################################################################
#				GROUP METHODS
=head2

	SUBROUTINE		deleteProject

	PURPOSE

		ADD A GROUP OBJECT TO THE projects TABLE

=cut

method deleteProject {


    my $dbobject		=	$self->dbobject();
	my $json			=	$self->json();

	#### VALIDATE    
    my $username = $json->{'username'};
	print "{ error: 'Admin::deleteProject    User $username not validated' }" and exit unless $self->validate($username);

	#### DELETE FROM project TABLE
	my $name = $json->{data}->{name};
	my $description = $json->{data}->{description};
	my $notes = $json->{data}->{notes};
	my $query = qq{DELETE FROM project
	WHERE username='$username'
	AND name='$name'};
	my $success = $dbobject->do($query);	

	#### DELETE FROM groupmember TABLE
	$query = qq{DELETE FROM groupmember
	WHERE owner='$username'
	AND name='$name'};
	$dbobject->do($query);	

	if ( $success == 1 )
	{
		#### REMOVE PROJECT FROM workflow TABLE
		$self->removeProjectWorkflow();
	}
	else
	{
		print "{ error : 'Admin::deleteProject    Could not delete project $name from project table' }";
	}

	exit;
}


=head2

    SUBROUTINE:     removeFromGroup

    PURPOSE:

        VALIDATE THE admin USER THEN REMOVE A SOURCE, user, ETC.

		FROM THE groupmember TABLE UNDER THE GIVEN groupname

=cut

method removeFromGroup {


    my $dbobject		=	$self->dbobject();
	my $json			=	$self->json();

	#### VALIDATE USER USING SESSION ID	
    my $username = $json->{'username'};
    print "{ error: 'User $username not validated' }" and exit unless $self->validate($username);

	#### PARSE JSON INTO OBJECT
	my $jsonParser = JSON->new();
	my $data = $json->{data};

	#### PRIMARY KEYS FOR apps TABLE:
	####    name, type, location
	my $groupname = $data->{groupname};
	my $name = $data->{name};
	my $description = $data->{description};
	my $location = $data->{location};
	my $type = $data->{type};

	if ( not defined $groupname or not defined $name or not defined $type )
	{
		print "{   error: 'Either groupname, name or type not defined'   }";
		exit;
	}
	my $query = qq{DELETE FROM groupmember
WHERE owner = '$username'
AND groupname = '$groupname'
AND name = '$name'
AND type = '$type'};


	my $success = $dbobject->do($query);
	if ( $success == 1 )
	{
		print "{ status: 'Successfully deleted source $name from groupmember $groupname in groupmember table' }";
	}
	else
	{
		print "{ error: 'Could not delete member $name (type $type) from groupmember $groupname in groupmember table' }";
	}
	exit;
}




=head2

    SUBROUTINE:     addToGroup

    PURPOSE:

        VALIDATE THE admin USER THEN ADD A SOURCE TO

		THE groupmember TABLE UNDER THE GIVEN groupname

=cut

method addToGroup {


    my $dbobject		=	$self->dbobject();
	my $json			=	$self->json();

	#### VALIDATE    
	print "{ error: 'Admin::addToGroup    User not validated' }" and exit unless $self->validate();

	#### PARSE JSON INTO OBJECT
	my $jsonParser = JSON->new();
	my $data = $json->{data};
	$data->{owner} = $json->{username} if not defined $data->{owner};

	#### SET TABLE AND REQUIRED FIELDS	
	my $table = "groupmember";
	my $required_fields = ["owner", "groupname", "name", "type"];

	#### CHECK REQUIRED FIELDS ARE DEFINED
	my $not_defined = $dbobject->notDefined($data, $required_fields);
    print "{ error: 'Admin::addToGroup    undefined values: @$not_defined' }" and exit if @$not_defined;

	#### DO THE ADD
	my $inserted_fields = $dbobject->fields($table);
	my $success = $self->_addToTable($table, $data, $required_fields, $inserted_fields);
 	print "{ error: 'Admin::addToGroup    Could not insert $json->{type} $json->{name} into groupmember table' }" and exit if not defined $success;

	print "{ status: 'Admin::addToGroup    Successfully inserted $data->{type} $data->{name} into groupmember table' }";
	exit;
}






=head2

    SUBROUTINE:     deleteSource

    PURPOSE:

        VALIDATE THE admin USER THEN DELETE A SOURCE

=cut

method deleteSource {


    my $dbobject		=	$self->dbobject();
	my $json			=	$self->json();

	#### VALIDATE    
    my $username = $json->{'username'};
    my $session_id = $json->{'sessionId'};

	#### VALIDATE USER USING SESSION ID	
    if ( not $self->validate($username, $session_id) )
    {
        print "{ error: 'User $username not validated' }";
        exit;
    }

	my $data = $json->{data};

	#### PRIMARY KEYS FOR source TABLE:
	####    name, type, location
	my $name = $data->{name};

	if ( not defined $name )
	{
		print "{   error: 'Name of source not defined'   }";
		exit;
	}

	#### DELETE SOURCE
	my $query = qq{DELETE FROM source
WHERE username='$username'
AND name='$name'};
	my $success = $dbobject->do($query);
	if ( $success )
	{
		my $query = qq{DELETE FROM groupmember
WHERE owner='$username'
AND name = '$name'
AND type = 'source'};
		$success = $dbobject->do($query);
		if ( $success )
		{
			print "{	status: 'Succesfully deleted source $name from source table' }";
		}
		else
		{
			print "{ error: 'Could not delete source $name from groupmember table' }";
		}
	}
	else
	{
		print "{ error: 'Could not delete source $name from source table' }";
	}

	exit;
}


=head2

    SUBROUTINE:     saveSource

    PURPOSE:

        VALIDATE THE admin USER THEN SAVE SOURCE INFORMATION

=cut

method saveSource {


    my $dbobject		=	$self->dbobject();
	my $json			=	$self->json();

	#### VALIDATE    
    my $username = $json->{'username'};
    my $session_id = $json->{'sessionId'};

	#### VALIDATE USER USING SESSION ID	
    if ( not $self->validate($username, $session_id) )
    {
        print "{ error: 'Admin::saveSource    User $username not validated' }";
        exit;
    }

	#### PARSE JSON INTO OBJECT
	my $jsonParser = JSON->new();
	my $data = $json->{data};
	$data->{username} = $username;

	#### PRIMARY KEYS FOR apps TABLE:
	####    name, type, location
	my $name = $data->{name};
	my $description = $data->{description};
	my $type = $data->{type};
	my $location = $data->{location};

	if ( not defined $name or not defined $location )
	{
		print "{   error: 'Admin::saveSource    Either name or location not defined'   }";
		exit;
	}

	my $fields = $dbobject->fields('source');

	#### CHECK IF THIS ENTRY EXISTS ALREADY
	my $query = qq{SELECT 1 FROM source
WHERE username='$username'
AND name='$name'};
#AND type= '$type'
#AND location = '$location'};
    my $already_exists = $dbobject->query($query);

	#### UPDATE THE source TABLE ENTRY IF EXISTS ALREADY
	if ( defined $already_exists )
	{
		#### UPDATE THE source TABLE
		my $set = '';
		foreach my $field ( @$fields )
		{
			$data->{$field} =~ s/'/"/g;
			$set .= "$field = '$data->{$field}',\n";
		}
		$set =~ s/,\s*$//;

		my $query = qq{UPDATE source
SET $set
WHERE username='$username'
AND name='$name'};

		my $success = $dbobject->do($query);

		if ( $success )
		{
			#### UPDATE THE groupmember TABLE
			$query = qq(UPDATE groupmember
SET
description = '$description',
location = '$location'
WHERE owner='$username'
AND name='$name'
AND type = 'source');
			$success = $dbobject->do($query);

			if ( $success == 1 )
			{
				print "{ status: 'Successfully updated source and groupmember tables  with source $name' }";
			}
			else
			{
				print "{ error: 'Could not update groupmember table with source $name' }";
			}
			exit;
		}
		else
		{
			print "{ error: 'Could not update source table xwith source $name' }";
		}
	}

	#### OTHERWISE, INSERT THE ENTRY
	else
	{
		my $values = '';
		foreach my $field ( @$fields )
		{
			my $value = $data->{$field};
			$value = '' if not defined $value;
			$values .= "'$value',\n";
		}
		$values =~ s/,\s*$//;

		my $query = qq{INSERT INTO source
VALUES ($values) };

		my $success = $dbobject->do($query);
		if ( $success == 1 )
		{
			print "{ status: 'Successfully inserted into source table' }";
		}
		else
		{
			print "{ error: 'Could not insert into source table' }";
		}
		exit;
	}
}

##############################################################################
#				APP METHODS
=head2

    SUBROUTINE:     deleteApp

    PURPOSE:

        VALIDATE THE admin USER THEN DELETE AN APPLICATION

=cut

method deleteApp {

    my $dbobject		=	$self->dbobject();
	my $json			=	$self->json();



	#### VALIDATE    
	print "{ error: 'Admin::deleteApp    User not validated' }" and exit unless $self->validate();

	#### GET DATA 
	my $data = $json->{data};
	$data->{owner} = $json->{username};
	my $table = "app";
	my $required_fields = ["owner", "name", "type", "location"];

	#### CHECK REQUIRED FIELDS ARE DEFINED
	my $not_defined = $dbobject->notDefined($data, $required_fields);
    print "{ error: 'Admin::deleteApp    undefined values: @$not_defined' }" and exit if @$not_defined;

	#### REMOVE IF EXISTS ALREADY
	my $success = $self->_removeFromTable($table, $data, $required_fields);
	print "{ status: 'Admin::deleteApp    Deleted application $data->{name} in apps table' }" if $success;
	print "{ error: 'Admin::deleteApp    Could not delete application $data->{name} from the apps table' }" if not $success;
	exit;
}


=head2

    SUBROUTINE:     saveApp

    PURPOSE:

        VALIDATE THE admin USER THEN SAVE APPLICATION INFORMATION

=cut

method saveApp {
    my $dbobject		=	$self->dbobject();
	my $json			=	$self->json();


	#### VALIDATE    
    my $username = $json->{'username'};
	print "{ error: 'Admin::saveApp    User $username not validated' }" and exit unless $self->validate($username);

	#### GET DATA FOR PRIMARY KEYS FOR apps TABLE:
	####    name, type, location
	my $data = $json->{data};
	my $name = $data->{name};
	my $type = $data->{type};
	my $location = $data->{location};

	#### CHECK INPUTS
	print "{ error: 'Admin::saveApp    Name $name not defined or empty' }" and exit if not defined $name or $name =~ /^\s*$/;
	print "{ error: 'Admin::saveApp    Name $type not defined or empty' }" and exit if not defined $type or $type =~ /^\s*$/;
	print "{ error: 'Admin::saveApp    Name $location not defined or empty' }" and exit if not defined $location or $location =~ /^\s*$/;



	#### SET owner AS USERNAME IN data
	$data->{owner} = $username;

	#### EXIT IF ONE OR MORE PRIMARY KEYS IS MISSING	
	print "{   error: 'Admin::saveApp    Either name, type or location not defined'   }" and exit if not defined $name or not defined $type or not defined $location;

	#### DELETE EXISTING ENTRY IF PRESENT
	my $query = qq{SELECT 1 FROM app
WHERE owner='$username'
AND name='$name'};
    my $exists = $dbobject->query($query);

	#### DELETE ENTRY IF EXISTS IN app TABLE
	if ( defined $exists and $exists )
	{
		$query = qq{DELETE FROM app
WHERE owner='$username'
AND name='$name'};
		my $success = $dbobject->do($query);
	}

	#### GET VALUE FOR EACH TABLE FIELD
	my $fields = $dbobject->fields('app');
	my $values = '';
	foreach my $field ( @$fields )
	{
		my $value = $data->{$field};
		$value = '' if not defined $value;
		$values .= "'$value',\n";
	}
	$values =~ s/,\s*$//;

	#### DO INSERT
	$query = qq{INSERT INTO app
VALUES ($values) };
	my $success = $dbobject->do($query);
	if ( $success != 1 )
	{
		print "{ error: 'Admin::saveApp    Could not insert application $name into app table '}";
	}
	else
	{
		print "{ status: 'Admin::saveApp    Inserted application $name into app table' }";
	}

	exit;
}


=head2

    SUBROUTINE:     saveParameter

    PURPOSE:

        VALIDATE THE admin USER THEN SAVE APPLICATION INFORMATION

=cut

method saveParameter {


    my $dbobject		=	$self->dbobject();
	my $json			=	$self->json();

	#### VALIDATE    
    my $username = $json->{'username'};
	print "{ error: 'Admin::saveParameter    User $username not validated' }" and exit unless $self->validate($username);

	#### GET DATA FOR PRIMARY KEYS FOR parameters TABLE:
	####    name, valuetype, location
	my $data = $json->{data};
	my $appname = $data->{appname};
	my $name = $data->{name};
	my $valuetype = $data->{valuetype};

	#### CHECK INPUTS
	print "{ error: 'Admin::saveParameter    Name $appname not defined or empty' }" and exit if not defined $appname or $appname =~ /^\s*$/;
	print "{ error: 'Admin::saveParameter    Name $name not defined or empty' }" and exit if not defined $name or $name =~ /^\s*$/;
	print "{ error: 'Admin::saveParameter    valuetype $valuetype not defined or empty' }" and exit if not defined $valuetype or $valuetype =~ /^\s*$/;


	#### SET owner AS USERNAME IN data
	$data->{owner} = $username;

	#### EXIT IF ONE OR MORE PRIMARY KEYS IS MISSING	
	print "{   error: 'Admin::saveParameter    Either name, valuetype or appname not defined'   }" and exit if not defined $name or not defined $valuetype or not defined $appname;

	#### DELETE EXISTING ENTRY IF PRESENT
	my $query = qq{SELECT 1 FROM parameter
WHERE owner='$username'
AND appname='$appname'
AND name='$name'};
    my $exists = $dbobject->query($query);

	#### DELETE ENTRY IF EXISTS IN parameter TABLE
	if ( defined $exists and $exists )
	{
		$query = qq{DELETE FROM parameter
WHERE owner='$username'
AND appname='$appname'
AND name='$name'};
		my $success = $dbobject->do($query);
	}

	#### GET VALUE FOR EACH TABLE FIELD
	my $fields = $dbobject->fields('parameter');
	my $values = '';
	foreach my $field ( @$fields )
	{
		my $value = $data->{$field};
		$value = '' if not defined $value;

		#### CONVERT ' TO \'
		$value =~ s/'/\\'/g;

		$values .= "'$value',\n";
	}
	$values =~ s/,\s*$//;


	#### DO INSERT
	$query = qq{INSERT INTO parameter
VALUES ($values) };
	my $success = $dbobject->do($query);
	if ( $success != 1 )
	{
		print "{ error: 'Admin::saveParameter    Could not insert parameter $name into parameter table '}";
	}
	else
	{
		print "{ status: 'Admin::saveParameter    Inserted parameter $name into parameter table' }";
	}

	exit;
}


=head2

    SUBROUTINE:     deleteParameter

    PURPOSE:

        VALIDATE THE admin USER THEN DELETE AN APPLICATION

=cut

method deleteParameter {


    my $dbobject		=	$self->dbobject();
	my $json			=	$self->json();

	#### VALIDATE    
    my $username = $json->{'username'};
	print "{ error: 'Admin::deleteParameter    User $username not validated' }" and exit unless $self->validate($username);

	#### GET DATA FOR PRIMARY KEYS FOR parameters TABLE:
	####    name, type, appname
	my $data = $json->{data};
	my $name = $data->{name};
	my $appname = $data->{appname};

	print "{   error: 'Admin::deleteParameter    name not defined'   }" and exit if not defined $name;
	print "{   error: 'Admin::deleteParameter    type not defined'   }" and exit if not defined $appname;

	my $query = qq{SELECT 1 FROM parameter
WHERE owner='$username'
AND appname='$appname'
AND name='$name'};
    my $exists = $dbobject->query($query);

	#### EXIT IF ENTRY DOESN'T EXIST IN TABLE
	print "{ error: 'Admin::deleteParameter    Parameterlication $name not found in parameter table' }" and exit if $exists != 1;

	$query = qq{DELETE FROM parameter
WHERE owner='$username'
AND appname='$appname'
AND name='$name'};
	my $success = $dbobject->do($query);

	if ( $success )
	{
		print "{ status: 'Admin::deleteParameter    Deleted parameter $name from parameter table' }";
	}
	else
	{
		print "{ error: 'Admin::deleteParameter    Could not delete parameter $name from parameter table' }";
	}
	exit;
}








##############################################################################
#				ACCESS METHODS
=head2

	SUBROUTINE		saveAccess

	PURPOSE

		UPDATE THE 'ACCESS' TABLE ENTRIES SPECIFIED IN THE

		JSON OBJECT SENT FROM THE CLIENT:

			1. UPDATE ALL ROWS WHERE ALL ENTRIES ARE DEFINED IN THE JSON

			2. IGNORE ROWS WHERE ANY ENTRIES ARE UNDEFINED

=cut

method saveAccess {



    my $dbobject		=	$self->dbobject();
	my $json			=	$self->json();

	my $access_entries = $json->{data};

	#### VALIDATE    
    my $username = $json->{'username'};
    my $session_id = $json->{'sessionId'};
    print "{ error: 'User $username not validated' }"
		and exit if not $self->validate();

	#### GET FIELDS	
	my $fields = $dbobject->fields("access");
	my $fieldstring = join ",", @$fields;

	#### SET TABLE AND REQUIRED FIELDS	
	my $table = "access";
	my $required_fields = ['owner', 'groupname'];

	#### CHECK WHICH USERS TO ADD
	foreach my $access_entry ( @$access_entries )
	{
		#### CHECK REQUIRED FIELDS ARE DEFINED
		my $not_defined = $dbobject->notDefined($access_entry, $required_fields);
		print "{ error: 'Admin::saveAccess    undefined values: @$not_defined' }" and exit if @$not_defined;

		#### IGNORE THIS PROJECT ENTRY IF A FIELD IS UNDEFINED
		my $defined = 1;
		foreach my $field ( @$fields )
		{
			if ( not defined $access_entry->{$field} )
			{
				$defined = 0;
				last;
			}
		}
		next if not defined $defined;

		my $set = '';
		foreach my $field ( @$fields )
		{
			$set .= "$field = '$access_entry->{$field}',\n";
		}
		$set =~ s/,$//;

		my $query = qq{UPDATE access
SET $set WHERE owner = '$username'
AND groupname = '$access_entry->{groupname}'};

		my $success = $dbobject->do($query);
	}

	print "{ 'message' : 'saveAccess completed' }";
}


##############################################################################
#				GROUPMEMBER METHODS
=head2

	SUBROUTINE		saveGroups

	PURPOSE

		SAVE ALL THE USERS IN THE 'GROUP' TABLE SPECIFIED IN THE

		JSON OBJECT SENT FROM THE CLIENT:

			1. ADD ANY USERS NOT PRESENT IN THE TABLE BUT

				PRESENT IN THE JSON OBJECT

			2. DELETE ANY USERS PRESENT IN THE TABLE BUT

				NOT PRESENT IN THE JSON OBJECT

=cut

method saveGroups {


    my $dbobject		=	$self->dbobject();
	my $json			=	$self->json();

	my $jsonParser = JSON->new();
	my $users = $jsonParser->jsonToObj($json->{data});


	#### VALIDATE    
    my $username = $json->{'username'};
    my $session_id = $json->{'sessionId'};
    if ( not $self->validate($username, $session_id) )
    {
        print "{ error: 'User $username not validated' }";
        exit;
    }

    my $query = qq{SELECT DISTINCT groupname, groupdesc
	FROM groupmember
	WHERE owner='$username'};
	my $groups = $dbobject->queryhasharray($query);
	$groups = [] if not defined $groups;

	my $add_users;
	my $remove_users;

	#### CHECK WHICH USERS TO ADD
	foreach my $user ( @$users )
	{
		my $groupname = $user->{groupname};
		my $groupdesc = $user->{groupdesc};

		my $matched = 0;
		foreach my $group ( @$groups )
		{
			my $existing_name = $group->{groupdesc};
			my $existing_groupname = $group->{groupname};
			if ( $existing_name eq $groupdesc
				and $existing_groupname eq $groupname )
			{
				$matched = 1;
				$user->{groupdesc} = $group->{groupdesc};
				last;
			}
		}

		$user->{owner} = $username;
		push @$add_users, $user if not $matched;
	}

	#exit;

	#### CHECK WHICH USERS TO REMOVE
	foreach my $group ( @$groups )
	{
		my $existing_name = $group->{groupdesc};
		my $existing_groupname = $group->{groupname};


		my $matched = 0;
		foreach my $user ( @$users )
		{
			my $groupname = $user->{groupname};
			my $groupdesc = $user->{groupdesc};


			if ( $existing_name eq $groupdesc
				and $existing_groupname eq $groupname )
			{				

				$matched = 1;
			}
		}

		push @$remove_users, $group if not $matched;
	}


	my $fields = $dbobject->fields("groups");
	my $fieldstring = join ",", @$fields;
	foreach my $user ( @$add_users )
	{
		my $tsvline = $dbobject->fields_tsvline($fieldstring, $user);
		$tsvline =~ s/\t/', '/g;
		$tsvline = "'" . $tsvline . "'";

		$query = qq{INSERT INTO groups VALUES ($tsvline)};
		my $success = $dbobject->do($query);
	}

	foreach my $user ( @$remove_users )
	{
		my $groupdesc = $user->{groupdesc};
		my $groupname = $user->{groupname};
		my $query = qq{DELETE FROM groupmember WHERE groupdesc = '$groupdesc' AND groupname = '$groupname' AND owner = '$username'};
		my $success = $dbobject->do($query);
	}

	print "{ 'message' : 'saveGroups completed' }";
}



=head2

	SUBROUTINE		saveGroup

	PURPOSE

		1. ADD A GROUP OBJECT TO THE groups TABLE

		2. ADD A CORRESPONDING ENTRY IN THE access TABLE

=cut

method saveGroup {


    my $dbobject		=	$self->dbobject();
	my $json			=	$self->json();

	#### VALIDATE    
	print "{ error: 'Admin::saveGroup    User not validated' }" and exit unless $self->validate();

	#### ADD USERNAME TO NEW GROUP
	my $newgroup = $json->{data};
    my $username = $json->{'username'};
	$newgroup->{username} = $username;

	#### SET UP ADD/REMOVE VARIABLES
	my $table = "groups";
	my $required_fields = [ "username", "name" ];
	my $inserted_fields = [ "username", "name", "description", "notes" ];

	#### DO THE DELETE
	my $deleted = $self->_removeFromTable($table, $newgroup, $required_fields);
	print "{ error : 'Admin::saveGroup    Could not add group $newgroup->{name} to groups table' }" if not $deleted;	

	#### DO THE ADD
	my $added = $self->_addToTable($table, $newgroup, $required_fields, $inserted_fields);
	print "{ status : 'Group $newgroup->{name} added to groups table' }" if $added;
	print "{ error : 'Admin::saveGroup    Could not add group $newgroup->{name} to groups table' }" if not $added;


	#### ADD ENTRY IN access TABLE
	#### DO THE DELETE


}



=head2

	SUBROUTINE		deleteGroup

	PURPOSE

		DELETE A GROUP OBJECT FROM

		1. THE group TABLE

		2. THE groupmember TABLE IF PRESENT

=cut

method deleteGroup {


    my $dbobject		=	$self->dbobject();
	my $json			=	$self->json();

	#### VALIDATE    
	print "{ error: 'Admin::deleteGroup    User not validated' }" and exit unless $self->validate();

	#### ADD USERNAME TO NEW GROUP
	my $newgroup = $json->{data};
    my $username = $json->{'username'};
	$newgroup->{username} = $username;

	#### SET UP ADD/REMOVE VARIABLES
	my $table = "groups";
	my $required_fields = [ "username", "name" ];
	my $inserted_fields = [ "username", "name", "description", "notes" ];

	#### DO THE DELETE
	my $deleted = $self->_removeFromTable($table, $newgroup, $required_fields);
	print "{ error : 'Admin::deleteGroup    Could not add group $newgroup->{name} to groups table' }" if not $deleted;	
}


=head2

	SUBROUTINE		saveGroupUsers

	PURPOSE

		SAVE ALL THE USERS IN THE 'GROUP' TABLE SPECIFIED IN THE

		JSON OBJECT SENT FROM THE CLIENT:

			1. ADD ANY USERS NOT PRESENT IN THE TABLE BUT

				PRESENT IN THE JSON OBJECT

			2. DELETE ANY USERS PRESENT IN THE TABLE BUT

				NOT PRESENT IN THE JSON OBJECT

=cut

method saveGroupUsers {


    my $dbobject		=	$self->dbobject();
	my $json			=	$self->json();

	my $jsonParser = JSON->new();
	my $users = $jsonParser->jsonToObj($json->{data});


	#### VALIDATE    
    my $username = $json->{'username'};
    my $session_id = $json->{'sessionId'};
    if ( not $self->validate($username, $session_id) )
    {
        print "{ error: 'User $username not validated' }";
        exit;
    }

    #my $query = qq{SELECT * FROM groupmember WHERE owner='$owner'};
    my $query = qq{SELECT DISTINCT name, groupname, type
	FROM groupmember
	WHERE owner='$username'};
	my $groupmember = $dbobject->queryhasharray($query);
    if ( not defined $groupmember )
    {
        my $json = "{   'error' : 'No groupmember for owner: $username'    }";
        print "$json\n";
        exit;
    }

	my $add_users;
	my $remove_users;

	#### CHECK WHICH USERS TO ADD
	foreach my $user ( @$users )
	{
		my $groupname = $user->{groupname};
		my $name = $user->{name};

		my $matched = 0;
		foreach my $groupuser ( @$groupmember )
		{
			my $existing_name = $groupuser->{name};
			my $existing_groupname = $groupuser->{groupname};
			if ( $existing_name eq $name
				and $existing_groupname eq $groupname )
			{
				$matched = 1;
				$user->{groupdesc} = $groupuser->{groupdesc};
				last;
			}
		}

		$user->{owner} = $username;
		push @$add_users, $user if not $matched;
	}

	#exit;

	#### CHECK WHICH USERS TO REMOVE
	foreach my $groupuser ( @$groupmember )
	{
		my $existing_name = $groupuser->{name};
		my $existing_groupname = $groupuser->{groupname};


		my $matched = 0;
		foreach my $user ( @$users )
		{
			my $groupname = $user->{groupname};
			my $name = $user->{name};		

			if ( $existing_name eq $name
				and $existing_groupname eq $groupname )
			{


				$matched = 1;
			}
		}

		push @$remove_users, $groupuser if not $matched;
	}

	#### ADD USERS
	my $fields = $dbobject->fields("groupmember");
	my $fieldstring = join ",", @$fields;

	my $fileroot = $self->getFileroot($username);

	foreach my $user ( @$add_users )
	{
		my $groupname = $user->{groupname};
		my $query = qq{SELECT groupdesc FROM groupmember WHERE groupname = '$groupname'};
		my $groupdesc = $dbobject->query($query);

		my $location = "$fileroot/$user->{name}";

		$user->{groupdesc} = $groupdesc;
		$user->{location} = $location;

		my $tsvline = $dbobject->fields_tsvline($fieldstring, $user);
		$tsvline =~ s/\t/', '/g;
		$tsvline = "'" . $tsvline . "'";

		$query = qq{INSERT INTO groupmember VALUES ($tsvline)};
		my $success = $dbobject->do($query);

		#### UPDATE access TABLE
		#### SO THAT IS IS CONSTANTLY IN SYNC WITH groupmember TABLE
		#### NB: DEFAULT PERMISSIONS: 751
		if ( $user->{type} eq 'project' )
		{
			my $project = $user->{name};
			$query = qq{INSERT INTO access
VALUES ('$project', '$groupname', '$username', 7, 5, 1, '$username/$project')};
		my $success = $dbobject->do($query);
		}
	}


	#### REMOVE USERS
	foreach my $user ( @$remove_users )
	{
		my $name = $user->{name};
		my $groupname = $user->{groupname};
		my $query = qq{DELETE FROM groupmember
WHERE name = '$name'
AND groupname = '$groupname'
AND owner = '$username'};
		my $success = $dbobject->do($query);

		#### UPDATE access TABLE
		#### SO THAT IS IS CONSTANTLY IN SYNC WITH groupmember TABLE
		#### NB: DEFAULT PERMISSIONS: 751
		if ( $user->{type} eq 'project' )
		{
			my $project = $user->{name};
			$query = qq{DELETE FROM access
WHERE project = '$project'
AND groupname = '$groupname'
AND owner = '$username'};
		print "$query\n";
		my $success = $dbobject->do($query);
		}
	}

	print "{ 'message' : 'saveGroupUsers completed' }";
}





##############################################################################
#				LOGIN METHODS
=head2

    SUBROUTINE      ldap

    PURPOSE

        USE LDAP SERVER TO VALIDATE USER

	INPUTS

		1. JSON->USERNAME

		2. JSON->PASSWORD

		3. CONF->LDAP_SERVER

	OUTPUTS

		1. RETURN 1 ON SUCCESS, 0 ON FAILURE

=cut

method ldap {



	use Net::LDAP;

	my $dbobject	=	$self->dbobject();
	my $json 		=	$self->json();
	my $conf 		=	$self->conf();

	my $username	=	$json->{username};
	my $password 	=	$json->{password};

	#### 'ldap.ccs.miami.edu'
	#### Server: ldap.ccs.miami.edu
	#### Binddn: uid=USERNAME,ou=Users,dc=ccs,dc=miami,dc=edu
	#### Bindpw: USERPASS
	my $ldap_server = $conf->getKeyValue('agua', "LDAP_SERVER");

	####  RETURN 1 IF NO LDAP SERVER
	return 1 if not defined $ldap_server;

	#### CREATE Net::LDAP OBJECT
	my $ldap = Net::LDAP->new($ldap_server);

	#### TEST BIND TO A DIRECTORY WITH DN AND PASSWORD.
	#### WILL RETURN '0' IF AUTHENTICATED
	my $message = $ldap->bind(
		"uid=$username,ou=Users,dc=ccs,dc=miami,dc=edu",
		"password" => "$password"
	);

	my $result = $message->code();

	#### CONVERT TO 1 IF SUCCES (I.E., RESULT = 0)
	if ( $result == 0 )
	{
		$result = 1;
	}
	else
	{
		$result = 0;
	}

	return $result;
}




=head2

    SUBROUTINE      login

    PURPOSE

        AUTHENTICATE USER USING ONE OF TWO WAYS:

		1. IF EXISTS 'LDAP_SERVER' ENTRY IN CONF FILE, USE THIS TO AUTHENTICATE

			THEN GENERATE A SESSION_ID IF SUCCESSFULLY AUTHENTICATED. STORE

			SESSION_ID IN sessions TABLE AND PRINT IT TO STDOUT

		2. OTHERWISE, CHECK INPUT PASSWORD AGAINST STORED PASSWORD IN users TABLE

			THEN GENERATE A SESSION_ID IF SUCCESSFULLY AUTHENTICATED. STORE

			SESSION_ID IN sessions TABLE AND PRINT IT TO STDOUT

	INPUTS

		1. JSON->USERNAME

		2. JSON->PASSWORD

	OUTPUTS

		1. SESSION ID

=cut

method login {

    my $dbobject = $self->dbobject();
	my $json = $self->json();

	my $username	=	$json->{username};
	my $password 	=	$json->{password};

    #### CHECK USERNAME AND PASSWORD DEFINED AND NOT EMPTY
    if ( not defined $username )    {   return; }
    if ( not defined $password )    {   return; }
    if ( not $username )    {   return; }
    if ( not $password )    {   return; }

	#### ADMIN VALIDATES AGAINST DATABASE, OTHER USERS VALIDATE
	#### AGAINST LDAP OR DATABASE IF LDAP NOT AVAILABLE
	my $is_admin = $self->isAdminuser($username);

	#### VALIDATE USING LDAP IF EXISTS 'LDAP_SERVER' ENTRY IN CONF FILE
	my $conf 		=	$self->conf();
	my $ldap_server = $conf->getKeyValue('agua', "LDAP_SERVER");
	my $match = 0;
	if ( not $is_admin and defined $ldap_server )
	{
		$match = $self->ldap();
	}

	#### OTHERWISE, GET STORED PASSWORD FROM users TABLE
	else
	{
		my $query = qq{SELECT password FROM aguausers
	WHERE username='$username'};
		my $stored_password = $dbobject->query($query);	

		#### CHECK FOR INPUT PASSWORD MATCHES STORED PASSWORD

		$match = $password =~ /^$stored_password$/; 
	}


	#### GENERATE SESSION ID
	my $session_id;

	#### IF PASSWORD MATCHES, STORE SESSION ID AND RETURN '1'
	my $exists;



#$match = 1;

	####
	my $now = "DATETIME('NOW')";
	my $dbtype = $self->conf()->getKeyValue('agua', 'DBTYPE');
	$now = "NOW()" if $dbtype eq "MySQL";

	if ( $match )
	{
		while ( not defined $session_id )
		{
			#### CREATE A RANDOM SESSION ID TO BE STORED IN dojo.cookie
			#### AND PASSED WITH EVERY REQUEST
			$session_id = time() . "." . $$ . "." . int(rand(1000));

			#### LATER: 
			# 1) create a file with that name
			# 2) keep data associated with the session in the sessionId file


			#### CHECK IF THIS SESSION ID ALREADY EXISTS
			my $exists_query = qq{
			SELECT username FROM sessions
			WHERE username = '$username'
			AND sessionid = '$session_id'};
			$exists = $dbobject->query($exists_query);

			if ( defined $exists )
			{
				$session_id = undef;
			}
			else
			{
			}
		}        

		#### IF IT DOES EXIST, UPDATE THE TIME
		if ( defined $exists )
		{
			my $update_query = qq{UPDATE sessions
			SET datetime = $now
			WHERE username = '$username'
			AND sessionid = '$session_id'};
			my $update_success = $dbobject->query($update_query);
		}

		#### IF IT DOESN'T EXIST, INSERT IT INTO THE TABLE
		else
		{
			my $query = qq{
INSERT INTO sessions
(username, sessionid, datetime)
VALUES
('$username', '$session_id', $now )};
			my $success = $dbobject->do($query);
			if ( $success )
			{
			}
		}		
	}

	#### LATER:: CLEAN OUT OLD SESSIONS
	# DELETE FROM sessions WHERE datetime < ADDDATE(NOW(), INTERVAL -48 HOUR)
	# DELETE FROM sessions WHERE datetime < DATE_SUB(NOW(), INTERVAL 1 DAY)
	my $delete_query = qq{
DELETE FROM sessions
WHERE datetime < DATETIME ('NOW', 'LOCALTIME', '-24 HOURS') };
	$dbobject->do($delete_query);

	if ( not defined $session_id and defined $ldap_server)
	{
		print "{ error: 'LDAP authentication failed for user: $username'  }";
		exit;
	}
	elsif ( not defined $session_id )
	{
		print "{ error: 'Authentication failed for user: $username'  }";
		exit;
	}

	print "{ sessionId : '$session_id' }";
    exit;
}



=head2

    SUBROUTINE      newuser

    PURPOSE

        CHECK 'admin' NAME AND PASSWORD AGAINST INPUT VALUES. IF

        VALIDATED, CREATE A NEW USER IN THE users TABLE

=cut

method newuser {
    my $dbobject = $self->dbobject();

    my $json;
    if ( not $self->validate() )
    {
        $json = "{ {validated: false} }";
        print $json;
        exit;
    }

	my $username	=	$self->cgiParam('username');
	my $password 	=	$self->cgiParam('password');
	my $newuser 	=	$self->cgiParam('newuser');
	my $newuserpassword 	=	$self->cgiParam('newuserpassword');


    ##### CREATE TABLE users IF NOT EXISTS
    #$self->create_table('users');

	#### CHECK IF USER ALREADY EXISTS
	my $query = qq{SELECT username FROM users WHERE username='$newuser' LIMIT 1};
	my $exists_already = $dbobject->query($query);
	if ( $exists_already )
	{
		print "User exists already: $username\n";
		exit;
	}

    $query = qq{INSERT INTO users VALUES ('$newuser', '$newuserpassword', NOW())};
    my $success = $dbobject->do($query);

    if ( $success )
    {
        print "New user created\n";
    }
    else
    {
        print "Failed to create new user\n";
    }
    exit;
}









}


