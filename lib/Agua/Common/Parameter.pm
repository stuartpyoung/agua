package Agua::Common::Parameter;
use Moose::Role;
use Moose::Util::TypeConstraints;

=head2

	PACKAGE		Agua::Common::File

	PURPOSE

		APPLICATION PARAMETER AND STAGE PARAMETER METHODS FOR Agua::Common

=cut


use Data::Dumper;


##############################################################################
#				STAGEPARAMETER METHODS
##############################################################################
=head2

    SUBROUTINE:     getStageParameters

    PURPOSE:

        GET ARRAY OF HASHES:

			[
				appname1: [ parameters ], appname2 : [...], ...  ] 


=cut

sub getStageParameters {
	my $self		=	shift;
	my $owner		=	shift;


    my $dbobject		=	$self->dbobject();
	my $json			=	$self->json();

	#### SET OWNER IF NOT DEFINED
	$owner = $json->{username} if not defined $owner;

    #### VALIDATE    
    print "{ error: 'User session not validated' }" and exit unless $self->validate();

	my $query = qq{SELECT * FROM stageparameter
WHERE username='$owner'
ORDER BY appname, name};
	my $stageparameters = $dbobject->queryhasharray($query);
	$stageparameters = [] if not defined $stageparameters;


	return $stageparameters;
}










=head2

	SUBROUTINE		addStageParameter

	PURPOSE

		ADD A PARAMETER TO THE stageparameter TABLE

	INPUTS

		1. ENOUGH STAGE INFORMATION FOR UNIQUE ID

			AND TO SATISFY REQUIRED TABLE FIELDS

{"project":"Project1","workflow":"Workflow1","appname":"clusterMAQ","name":"cpus","appnumber":"1","description":"","discretion":"optional","type":"integer","paramtype":"integer","value":"optasdfasdional","username":"admin","sessionId":"9999999999.9999.999","mode":"addStageParameter"}

=cut

sub addStageParameter {
	my $self		=	shift;

    my $dbobject        =	$self->dbobject();
    my $json 			=	$self->json();


    #### VALIDATE
    print "{ error: 'Agua::Common::Parameter::addStageParameter    User session not validated' }" and exit unless $self->validate();

	#### SET TABLE AND REQUIRED FIELDS	
	my $table = "stageparameter";
	my $required_fields = ["username", "project", "workflow", "appname", "appnumber", "name"];

	#### CHECK REQUIRED FIELDS ARE DEFINED
	my $not_defined = $dbobject->notDefined($json, $required_fields);
    print "{ error: 'Agua::Common::Parameter::addStageParameter    undefined values: @$not_defined' }" and exit if @$not_defined;

	#### REMOVE IF EXISTS ALREADY
	$self->_removeFromTable($table, $json, $required_fields);

	#### ADD ALL FIELDS OF THE PASSED STAGE PARAMETER TO THE TABLE
	my $success = $self->_addToTable($table, $json, $required_fields);	
 	print "{ error: 'Agua::Common::Parameter::addStageParameter    Could not add parameter to stageparameter table'}" and exit if not defined $success;

 	print "{ status: 'Agua::Common::Parameter::addStageParameter    Added parameter $json->{name} to stageparameter table' }" and exit;
}



=head2

    SUBROUTINE     addParameters

    PURPOSE

		COPY parameter TABLE ENTRIES FOR THIS APPLICATION

		TO stageparameter TABLE

=cut

sub addParameters {
	my $self		=	shift;

    my $dbobject		=	$self->dbobject();
	my $json			=	$self->json();

    #### VALIDATE    
    print "{ error: 'User session not validated' }" and exit unless $self->validate();


	#### GET APPLICATION OWNER
	my $owner = $json->{owner};

	#### GET PRIMARY KEYS FOR parameters TABLE
    my $username = $json->{username};
	my $appname = $json->{name};
	my $project = $json->{project};
	my $workflow = $json->{workflow};
	my $number = $json->{number};


	#### CHECK INPUTS
	print "{ error: 'Agua::Common::Parameter::addParameters    appname $appname not defined or empty' }" and exit if not defined $appname or $appname =~ /^\s*$/;
	print "{ error: 'Agua::Common::Parameter::addParameters    username $username not defined or empty' }" and exit if not defined $username or $username =~ /^\s*$/;

	#### DELETE ENTRY IF EXISTS IN stageparameter TABLE
	my $success;
	my $query = qq{DELETE FROM stageparameter
WHERE username='$username'
AND project='$project'
AND workflow='$workflow'
AND appname='$appname'
AND appnumber='$number'};
	$success = $dbobject->do($query);


	#### GET ORIGINAL PARAMETERS FROM parameter TABLE
	$query = qq{SELECT * FROM parameter
WHERE owner='$owner'
AND appname='$appname'};
    my $parameters = $dbobject->queryhasharray($query);
	print "{ error: 'Agua::Common::Parameter::addParameters    no entries in parameter table'}" and exit if not defined $parameters;


	##### SET QUERY WITH PLACEHOLDERS
	my $table = "stageparameter";
	my $fields = $dbobject->fields($table);
	#my $fields_csv = $dbobject->fields_csv($table);

	$success = 1;
	foreach my $parameter ( @$parameters )
	{
		$parameter->{username} = $username;
		$parameter->{project} = $project;
		$parameter->{workflow} = $workflow;
		$parameter->{appnumber} = $number;

		#### INSERT %OPTIONAL% VARIABLES
		$parameter->{value} =~ s/%project%/$project/g;
		$parameter->{value} =~ s/%workflow%/$workflow/g;
		$parameter->{value} =~ s/%username%/$username/g;

		$parameter->{description} =~ s/'/\\'/g; 
		$parameter->{paramFunction} =~ s/'/\\'/g; 
		print Dumper $parameter;

		#### DO INSERT
		my $values_csv = $dbobject->fields_csvline($fields, $parameter);
		my $query = qq{INSERT INTO $table 
	VALUES ($values_csv) };
		my $do_result = $dbobject->do($query);

		$success = 0 if not $do_result;
	}

	print "{ error: 'Agua::Common::Parameter::addParameters    Could not copy parameters for application $appname to stageparameter table'}" and exit if not $success;

	print "{ status: 'Agua::Common::Parameter::addParameters    Copied parameters for application $appname to stageparameter table'}";

	exit;
}



##############################################################################
#				PARAMETER METHODS
##############################################################################
=head2

    SUBROUTINE:     getParameters

    PURPOSE:

        RETURN AN ARRAY OF PARAMETER KEY:VALUE PAIR HASHES

			[ { project: ..., workflow: ..., name: } , ... ]

=cut

sub getParameters {


	my $self		=	shift;
	my $owner		=	shift;

    my $dbobject		=	$self->dbobject();
	my $json			=	$self->json();

	#### SET DEFAULT OWNER IF NOT DEFINED
	$owner = $self->conf()->getKeyValue('agua', "ADMINUSER") if not defined $owner;

    #### VALIDATE    
    print "{ error: 'User session not validated' }" and exit unless $self->validate();

	#### GET USER'S OWN APPS	
    my $username = $json->{'username'};
    my $query = qq{SELECT * FROM parameter
WHERE owner = '$username'
ORDER BY apptype, appname};
    my $parameters = $dbobject->queryhasharray($query);
	$parameters = [] if not defined $parameters;

	return $parameters;
}





1;