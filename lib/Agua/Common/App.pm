package Agua::Common::App;
use Moose::Role;
use Moose::Util::TypeConstraints;

=head2

	PACKAGE		Agua::Common::App

	PURPOSE

		APPLICATION METHODS FOR Agua::Common

=cut


use Data::Dumper;

=head2

    SUBROUTINE:     getApps

    PURPOSE:

        RETURN A JSON LIST OF THE APPLICATIONS IN THE apps TABLE

=cut

sub getApps {

	my $self		=	shift;

    my $dbobject		=	$self->dbobject();
	my $json			=	$self->json();

	#### VALIDATE    
    my $username = $json->{username};
    my $session_id = $json->{sessionId};

    print "{ error: 'Agua::Common::App::getApps    User $username not validated' }"
		and exit unless $self->validate($username, $session_id);

	#### GET USER'S OWN APPS	
    my $query = qq{SELECT * FROM app
WHERE owner = '$username'
ORDER BY type,name};
    my $userapps = $dbobject->queryhasharray($query);
	$userapps = [] if not defined $userapps;

	#### GET admin USER'S APPS	
	my $admin = $self->conf()->getKeyValue("agua", 'ADMINUSER');
    $query = qq{SELECT * FROM app
WHERE owner = '$admin'
ORDER BY type};
    my $adminapps = $dbobject->queryhasharray($query);
	$adminapps = [] if not defined $adminapps;

	my $apps = {
		'common'	=>	$adminapps,
		'custom'	=>	$userapps
	};

	return $apps;
}

=head2

    SUBROUTINE:     deleteApp

    PURPOSE:

        VALIDATE THE admin USER THEN DELETE AN APPLICATION

=cut

sub deleteApp {
	my $self		=	shift;

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

sub saveApp {
	my $self		=	shift;

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

sub saveParameter {

	my $self		=	shift;

    my $dbobject		=	$self->dbobject();
	my $json			=	$self->json();

	#### VALIDATE    
    my $username = $json->{'username'};
	print "{ error: 'Admin::saveParameter    User $username not validated' }" and exit unless $self->validate($username);

	#### GET DATA FOR PRIMARY KEYS FOR parameters TABLE:
	####    name, type, location
	my $data = $json->{data};
	my $appname = $data->{appname};
	my $name = $data->{name};
	my $type = $data->{type};

	#### CHECK INPUTS
	print "{ error: 'Admin::saveParameter    Name $appname not defined or empty' }" and exit if not defined $appname or $appname =~ /^\s*$/;
	print "{ error: 'Admin::saveParameter    Name $name not defined or empty' }" and exit if not defined $name or $name =~ /^\s*$/;
	print "{ error: 'Admin::saveParameter    Name $type not defined or empty' }" and exit if not defined $type or $type =~ /^\s*$/;


	#### SET owner AS USERNAME IN data
	$data->{owner} = $username;

	#### EXIT IF ONE OR MORE PRIMARY KEYS IS MISSING	
	print "{   error: 'Admin::saveParameter    Either name, type or appname not defined'   }" and exit if not defined $name or not defined $type or not defined $appname;

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
		#### CONVERT ' TO \'
		$value =~ s/'/\\'/g;

		$value = '' if not defined $value;
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

sub deleteParameter {

	my $self		=	shift;

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




1;