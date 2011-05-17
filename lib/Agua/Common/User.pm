package Agua::Common::User;
use Moose::Role;
use Moose::Util::TypeConstraints;

=head2

	PACKAGE		Agua::Common::User

	PURPOSE

		USER METHODS FOR Agua::Common

=cut


use Data::Dumper;

=head2

    SUBROUTINE     getUsers

    PURPOSE

		RETURN USER-RELATED INFO

=cut
sub getUsers {

	my $self		=	shift;

    my $dbobject		=	$self->dbobject();
	my $json			=	$self->json();

	#### GET USERNAME AND SESSION ID
    my $username = $json->{'username'};

	#### VALIDATE USER USING SESSION ID	
	print "{ error: 'User $username not validated' }" and exit unless $self->validate($username);

	#### GET ALL SOURCES
	my $query = qq{SELECT DISTINCT username, firstname, lastname, email, description FROM users ORDER BY username};
	my $users = $dbobject->querytwoDarray($query);
	$users = [] if not defined $users;

	return $users;
}


=head2

    SUBROUTINE     updateUser

	PURPOSE

		SAVE USER INFORMATION TO users TABLE IN TWO

		SCENARIOS:

		1. THE USER UPDATES THEIR OWN INFORMATION

			A PASSWORD CHECK IS REQUIRED.

		2. IF THE USER IS AN ADMIN USER. NO oldpassword

			IS REQUIRED. THIS ALLOWS ADMIN USERS TO CHANGE

			THE PASSWORD OF THE data->{username} USER.

=cut

sub updateUser {
	my $self		=	shift;


    my $json 			=	$self->json();
    my $dbobject        =	$self->dbobject();

 	print "{ error: 'Agua::Common::User::updateUser    json->{data}->{oldpassword} is not defined}' }"
        and exit unless defined $json->{data}->{oldpassword}
		or $self->isAdminuser();

    #### ONLY AN ADMIN USER CAN CHANGE THE PASSWORD WITHOUT
    #### PROVIDING THE OLD PASSWORD
    my $old_password = $json->{data}->{oldpassword};
    my $new_password = $json->{data}->{newpassword};

    if ( not defined $new_password and not $new_password )
    {
		print "{ error: 'Agua::Common::User::updateUser    newpassword not defined' }" and exit; 
    }

    my $is_admin = $self->isAdminuser();
    if ( not defined $old_password and not $old_password and not $is_admin )
    {
		print "{ error: 'Agua::Common::User::updateUser    oldpassword not defined' }" and exit; 
    }

	#### GET USERNAME TO BE CHANGED
	my $changer_username = $json->{username};
	my $changed_username = $json->{data}->{username};

	#### GET USER PASSWORD AND COMPARE TO INPUT oldpassword
	my $query = "SELECT password FROM users WHERE username='$changed_username'";

    my $stored_password = $dbobject->query($query);

	#### QUIT IF NO STORED PASSWORD AND NOT ADMIN USER
    if ( not defined $stored_password and not $stored_password and not $is_admin)
    {
		print "{ error: 'Agua::Common::User::updateUser    stored_password not defined' }" and exit; 
    }

	#### QUIT IF CHANGING SOMEONE ELSE'S PASSWORD AND NOT ADMIN USER
	if ( $changer_username ne $changed_username and not $is_admin )
	{
		print "{ error: 'Agua::Common::User::updateUser    User does not have sufficient privileges: $changer_username' }" and exit; 
	}

	#### QUIT IF THE PASSWORD DOES NOT MATCH AND NOT ADMIN USER
	if ( $stored_password ne $old_password and not $is_admin )
	{
		print "{ error: 'Agua::Common::User::updateUser    Password does not match: $old_password' }" and exit; 
	}

	#### QUIT IF THE PASSWORD DOES NOT MATCH AND ADMIN USER CHANGING OWN PASSWORD
	if ( $stored_password ne $old_password
		and $changer_username eq $changed_username and $is_admin )
	{
		print "{ error: 'Agua::Common::User::updateUser    Password does not match: $old_password' }" and exit; 
	}

	#### REMOVE FROM users TABLE IF EXISTS ALREADY
	my ($success, $user) = $self->_removeUser();

    #### ADD TO users TABLE
	$success = $self->_addUser();	
 	print "{ error: 'Agua::Common::User::updateUser    Could not save user $json->{data}->{username} to users table' }" and exit if not defined $success;

    #### CHANGE USER PASSWORD
    my $password = $json->{data}->{newpassword};

#	#### FEDORA:
#    (NB: CAN'T LOGIN ...)
#    my $passwd = "echo $password | passwd --stdin $changed_username";

	#### UBUNTU
	my $change_password = "echo '$changed_username:$password' | chpasswd";
	print `$change_password`;	

 	print "{ status: 'Agua::Common::User::updateUser    Updated user $changed_username in users table'}";
}


=head2

    SUBROUTINE     addUser

	PURPOSE

		ADD A REPORT TO THE report TABLE

=cut

sub addUser {
	my $self		=	shift;


	my $json 			=	$self->json();
	if ( not $self->isAdminuser() )
	{
		my $username = $json->{username};
		print "{ error: 'Insufficient privileges for user creation. Please supply old password for user: $username' }" and exit; 
	}

	#### GET USERNAME
	my $username = $json->{data}->{username};

	#### ADD USER TO DATABASE
	my $success = $self->_addUser();


    #### GET USER DIR
    my $userdir = $self->conf()->getKeyValue("agua", 'USERDIR');

	#### CREATE USER ACCOUNT AND HOME DIRECTORY

	#### REDHAT:
	#my $addgroup = "/usr/sbin/groupadd $username";
	#my $adduser = "useradd -g $username -s/bin/bash -p $username -d $userdir/$username -m  $username";

	#### UBUNTU
	my $set_home = "useradd -D -b $userdir";
	print `$set_home`;
	my $adduser = "useradd -m jgilbert";
	print `$adduser`;


	#### CREATE AGUA DIR
	my $aguadir = $self->conf()->getKeyValue("agua", 'AGUADIR');
	my $mkdir = "mkdir -p $userdir/$username/$aguadir";
	print `$mkdir`;

    #### SET CHOWN SO THAT OTHER USERS ARE E
	#### GET USERNAME, AND APACHE USER FOR PERMISSIONS LATER
    my $apache_user = $self->conf()->getKeyValue("agua", 'APACHE_USER');
    my $chown = "chown -R $username:$apache_user $userdir/$username $userdir/$username/$aguadir";
	print `$chown`;

    #### SET chmod SO THAT agua USER CAN ACCESS AND CREATE FILES
    my $chmod = "chmod 770 $userdir/$username $userdir/$username/$aguadir";
	print `$chmod`;

    #### SET USER UMASK (APACHE'S UMASK IS ALREADY SET AT 002)
    my $umask = qq{echo "umask 0002" >> $userdir/$username/.bashrc};
	print `$umask`;

    my $test = qq{echo "test file" >> $userdir/$username/$aguadir/testfile.txt};
	print `$test`;

 	print "{ error: 'Agua::Common::User::addUser    Could not add user $json->{data}->{username} to user table\n" and exit if not defined $success;
 	print "{ status: 'Agua::Common::User::addUser    Added user $json->{data}->{username} to user table' }";
	exit;
}




=head2

    SUBROUTINE     removeUser

	PURPOSE

		REMOVE A REPORT FROM THE report TABLE

=cut

sub removeUser {
	my $self		=	shift;


    my $json 			=	$self->json();

	#### DO THE REMOVE
	my $success = $self->_removeUser();

	#### CREATE USER ACCOUNT AND HOME DIRECTORY
	my $username = $json->{data}->{username};
	my $command = "userdel  $username";
	print "$command\n";
	print `$command`;

 	print "{ error: 'Agua::Common::User::removeUser    Could not remove user $json->{data}->{username} from user table\n" and exit if not defined $success;
 	print "{ status: 'Agua::Common::User::removeUser    Removed user $json->{data}->{username} from user table\n";
	exit;



}



=head2

    SUBROUTINE     _addUser

	PURPOSE

		ADD A REPORT TO THE report TABLE

=cut

sub _addUser {
	my $self		=	shift;


    my $dbobject        =	$self->dbobject();
    my $json 			=	$self->json();


	#### SET TABLE AND REQUIRED FIELDS	
	my $table = "users";
	my $required_fields = ["username", "email"];

	#### SET password FIELD
	$json->{data}->{password} = $json->{data}->{newpassword};

	#### CHECK REQUIRED FIELDS ARE DEFINED
	my $not_defined = $dbobject->notDefined($json->{data}, $required_fields);
    print "{ error: 'Agua::Common::User::_addUser    undefined values: @$not_defined' }" and exit if @$not_defined;

	#### LATER: FIX THIS SO THAT THE DATE IS ADDED PROPERLY 
	$json->{data}->{datetime} = $dbobject->now();

	#### DO THE ADD
	return $self->_addToTable($table, $json->{data}, $required_fields);	
}



=head2

    SUBROUTINE     _removeUser

	PURPOSE

		REMOVE A REPORT FROM THE report TABLE

=cut

sub _removeUser {
	my $self		=	shift;


    my $dbobject        =	$self->dbobject();
    my $json 			=	$self->json();


	#### SET TABLE AND REQUIRED FIELDS	
	my $table = "users";
	my $required_fields = ["username"];

	#### CHECK REQUIRED FIELDS ARE DEFINED
	my $not_defined = $dbobject->notDefined($json->{data}, $required_fields);
    print "{ error: 'Agua::Common::User::_removeUser    undefined values: @$not_defined' }" and exit if @$not_defined;

	#### DO THE REMOVE
	return $self->_removeFromTable($table, $json->{data}, $required_fields);
}



1;