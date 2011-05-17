package Agua::Common::Privileges;
use Moose::Role;
use Moose::Util::TypeConstraints;

=head2

	PACKAGE		Agua::Common::Privileges

	PURPOSE

		AUTHENTICATION AND ACCESS PRIVILEGE METHODS FOR Agua::Common

=cut


use Data::Dumper;

##############################################################################
#				AUTHENTICATION METHODS
##############################################################################
=head2

    SUBROUTINE      isAdminuser

    PURPOSE

        CONFIRM USER HAS ADMIN PRIVILEGES

	INPUTS

		1. JSON->USERNAME

		2. CONF->ADMINS (COMMA-SEPARATED)

	OUTPUTS

		1. RETURN 1 ON SUCCESS, 0 ON FAILURE

=cut

sub isAdminuser {
	my $self		=	shift;

	my $json 		=	$self->json();
	my $conf 		=	$self->conf();
	my $username	=	$json->{username};
	my $isAdminusers 	= 	$conf->getKeyValue('agua', "ADMINS");

	if ( not defined $isAdminusers )
	{
		print "{  error: 'No admin users defined in conf file'  }";
		exit;
	}
	my @users = split ",", $isAdminusers;
	foreach my $isAdminuser ( @users )
	{
		return 1 if $isAdminuser eq $username;
	}

	return 0;
}

=head2

    SUBROUTINE      validate

    PURPOSE

        CHECK SESSION ID AGAINST STORED SESSION IN sessions TABLE

        TO VALIDATE USER

=cut

sub validate {
	my $self		=	shift;



	return 1 if  $self->validated();

	my $dbobject	=	$self->dbobject();
	my $json 		=	$self->json();



	my $username	=	$json->{username};
	my $requestor	=	$json->{requestor};
	my $session_id	=	$json->{sessionId};


	##### GET FROM CGI IF NOT DEFINED
	#$username =	$self->cgiParam('username') if not defined $username;
	#$requestor =	$self->cgiParam('requestor') if not defined $requestor;
	#$session_id =	$self->cgiParam('sessionId') if not defined $session_id;

	print "{ error: 'Agua::Common::Privileges::validate    username not defined' }" and exit if not defined $username;

	#### SET 'AS USER'
	my $as_user = $username;
	$as_user = $requestor if defined $requestor and $requestor;

	#### VALIDATE BASED ON PRESENCE OF SESSION ID IN sessions TABLE
	my $query = qq{
	SELECT username FROM sessions
	WHERE username = '$as_user'
	AND sessionid = '$session_id'};

    my $validated = $dbobject->query($query);

    #### IF IT DOES EXIST, UPDATE THE TIME
    if ( defined $validated )
    {
		my $conf = $self->conf();
		my $now = "DATETIME('NOW')";
		$now = "NOW()" if $conf->getKeyValue("database", "DBTYPE") =~ /^MYSQL$/i;

        my $update_query = qq{UPDATE sessions
        SET datetime = $now
        WHERE username = '$username'
        AND sessionid = '$session_id'};
        my $update_success = $dbobject->do($update_query);
    }

	if ( not defined $validated )	{	return	0;	}


	$self->validated(1);
	return 1;
}



##############################################################################
#				PRIVILEGE METHODS
##############################################################################
=head2

    SUBROUTINE:     _canAccess

    PURPOSE:

		INPUTS

			1. IDENTITY OF SHARER AND SHAREE

			2. LOCATION OF PROJECT OR SOURCE

			3. TYPE - WHETHER PROJECT OR SOURCE

		OUTPUTS

			1. RETURN 1 IF A USER CAN ACCESS ANOTHER USER'S PROJECT.

			2. RETURN 0 OTHERWISE

=cut

sub _canAccess {
    my $self        =   shift;    
    my $owner       =   shift;
    my $groupname     =   shift;
    my $requestor   =   shift;
	my $type		=	shift;



    my $dbobject = $self->dbobject();

    #### GET GROUP NAME FOR PROJECT
    my $query = qq{SELECT DISTINCT groupname from groupmember where owner = '$owner' and groupname = '$groupname' and type = '$type'};
    my $groupnames = $dbobject->queryarray($query);

	#### RETURN IF GROUP NAME NOT IN access TABLE
	if ( not defined $groupnames )
    {
        return;
    }

    #### CONFIRM THAT USER BELONGS TO THIS GROUP
	my $accessible = 0;
	for my $groupname ( @$groupnames )
	{
		$query = qq{SELECT owner FROM groupmember WHERE groupname = '$groupname' AND owner = '$owner' AND name = '$requestor' AND type = 'user'};
		my $access = $dbobject->query($query);
		if ( $access )
		{
			$accessible = 1;
			last;
		}
	}

	if ( not $accessible )
	{
		return;
	}

	return 1 if defined $accessible;

	return 0;
}


=head2

    SUBROUTINE:     projectPrivilege

    PURPOSE:

        1. Check the rights of a user to access another user's project

=cut

sub projectPrivilege {
    my $self        =   shift;    
    my $owner       =   shift;
    my $project     =   shift;
    my $requestor   =   shift;
    my $privilege   =   shift;



    my $dbobject = $self->dbobject();

    #### LATER: GET MAX PRIVILEGES ACROSS ALL GROUPS THE USER BELONGS TO
    my $max_privilege;

    #### GET GROUP NAME FOR PROJECT
    my $query = qq{SELECT groupname FROM  groupmember
WHERE owner = '$owner'
AND name = '$project'
AND type = 'project'};
    my $groupnames = $dbobject->queryarray($query);
	return if not defined $groupnames;

	#### RETURN THE MAXIMUM PRIVILEGE FOR THIS USER AMONG ALL OF THE
    #### GROUPS THEY BELONG TO (FOR THIS SHARER)
    my $max_write = 0;
    my $max_copy = 0;
    my $max_view = 0;
    foreach my $groupname ( @$groupnames )
    {
        $query = qq{SELECT owner
FROM groupmember
WHERE groupname = '$groupname'
AND owner = '$owner'
AND name = '$requestor'
AND type = 'user'};
        my $owner = $dbobject->query($query);
        next if not defined $owner;

        #### MAKE SURE THAT THE GROUP PRIVILEGES ALLOW ACCESS
        $query = qq{SELECT groupwrite, groupcopy, groupview
FROM access
WHERE owner = '$owner'
AND groupname='$groupname'};
        my $privileges = $dbobject->queryhash($query);
        next if not defined $privileges;
        $max_write = $privileges->{groupwrite} if $max_write < $privileges->{groupwrite};
        $max_copy = $privileges->{groupcopy} if $max_copy < $privileges->{groupcopy};
        $max_view = $privileges->{groupview} if $max_view < $privileges->{groupview};
    }


    return { groupwrite => $max_write, groupcopy => $max_copy, groupview => $max_view };
}

sub sumPrivilege () {
	my $self		=	shift;
	my $rights		=	shift;
	my $group		=	shift;

	my $privilege = 0;
	$privilege += $rights->{$group . "write"} if defined $rights->{$group . "write"};
	$privilege += $rights->{$group . "copy"}  if defined $rights->{$group . "copy"};
	$privilege += $rights->{$group . "view"}  if defined $rights->{$group . "view"};

	return $privilege;
}


sub canCopy() {
    my $self    =   shift;

    my $privileges = $self->projectPrivilege(@_);

    my $can_copy = $privileges->{groupcopy};

    return $can_copy;
}



1;