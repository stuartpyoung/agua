package Agua::Common::Group;
use Moose::Role;
use Moose::Util::TypeConstraints;

=head2

	PACKAGE		Agua::Common::Group

	PURPOSE

		GROUP METHODS FOR Agua::Common

=cut


use Data::Dumper;

=head2

    SUBROUTINE:     getGroups

    PURPOSE:

		RETURN AN ARRAY OF group HASHES

			E.G.:
			[
				{
				  'name' : 'NGS',
				  'desciption' : 'NGS analysis team',
				  'notes' : 'This group is for ...',
				},
				{
					...
			]

=cut


sub getGroups {
	my $self		=	shift;


    my $dbobject		=	$self->dbobject();
	my $json			=	$self->json();

	#### GET USERNAME AND SESSION ID
    my $username = $json->{'username'};

    #### VALIDATE    
    print "{ error: 'User session not validated' }" and exit unless $self->validate();

	#### GET ALL SOURCES
	my $query = qq{SELECT * FROM groups
WHERE username='$username'};
	my $groups = $dbobject->queryhasharray($query);

	$groups = [] if not defined $groups;

	return $groups;
}




=head2

    SUBROUTINE:     getGroupMembers

    PURPOSE:

		RETURN THESE GROUP-RELATED TABLES:

			groupmember 

	INPUTS

		1. JSON OBJECT CONTAINING username AND sessionId

	OUTPUTS	

			groups: JSON groupmember HASH

			E.G.:

			{
				'nextgen' : [],
				'bioinfo' : [
					{
					  'groupdesc' : '',
					  'owner' : 'syoung',
					  'location' : '/mihg/data/NGS/syoung/base/pipeline/human-chr-squashed',
					  'groupname' : 'bioinfo',
					  'name' : 'Eland squashed human chromosomes (36.1)',
					  'type' : 'source',
					  'description' : 'A directory containing human chromosome sequence files processed with squashGenome for input into Eland(Build 36.1)'
					},
					{
					  'groupdesc' : undef,
					  'owner' : 'syoung',
					  'location' : '/mihg/data/NGS/syoung/base/pipeline/human-chr/fa',
					  'groupname' : 'bioinfo',
					  'name' : 'Human chromosome FASTA (36.1)',
					  'type' : 'source',
					  'description' : 'Human chromosome FASTA files (Build 36.1)'
					},
				   ...
				]
			}

=cut


sub getGroupMembers {


	my $self		=	shift;

    my $dbobject		=	$self->dbobject();
	my $json			=	$self->json();

	#### GET USERNAME
    my $username = $json->{'username'};

	#### VALIDATE USER USING SESSION ID	
	print "{ error: 'User $username not validated' }" and exit unless $self->validate($username);

	#### GET ALL GROUPS AND USERS IN THIS USER'S groupmember TABLE
    my $query = qq{SELECT *
FROM groupmember 
WHERE owner='$username'
ORDER BY groupname, name};
    my $groupmember = $dbobject->queryhasharray($query);
	$groupmember = [] if not defined $groupmember;

	return $groupmember;
	#
	##### SORT SOURCES BY groupname INTO AN ARRAY OF HASHES
	#my $groups = {};
	#my $groupsarray = [];
	#my $current_groupname = $$groupmember[0]->{groupname};
	#foreach my $groupuser( @$groupmember)
	#{
	#	my $groupname = $groupuser->{groupname};
	#	my $name = $groupuser->{name};
	#
	#	if ( $groupname eq $current_groupname )
	#	{
	#
	#		
	#		push @$groupsarray, $groupuser;
	#	}
	#	else
	#	{
	#		
	#		#### ADD TO GROUP USERS HASH
	#		$groups->{$current_groupname} = $groupsarray;
	#
	#		#### REINITITALISE GROUP USERS ARRAY FOR NEW HASH KEY
	#		$groupsarray = [];
	#		#if ( $groupuser->{type} eq 'group' )
	#		#{
	#			push @$groupsarray, $groupuser;
	#		#}
	#		
	#		#### SET NEW HASH KEY
	#		$current_groupname = $groupname;
	#	}
	#}
	#
	##### ADD TO GROUP USERS HASH
	#$groups->{$current_groupname} = $groupsarray;
	#
	#return $groups;
}






=head2

	SUBROUTINE		deleteProject

	PURPOSE

		ADD A GROUP OBJECT TO THE projects TABLE

=cut

sub deleteProject {

	my $self		=	shift;

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

sub removeFromGroup {

	my $self		=	shift;

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

sub addToGroup {

	my $self		=	shift;

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

sub deleteSource {

	my $self		=	shift;

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

sub saveSource {

	my $self		=	shift;

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




1;