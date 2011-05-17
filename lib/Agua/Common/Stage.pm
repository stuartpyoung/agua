package Agua::Common::Stage;
use Moose::Role;
use Moose::Util::TypeConstraints;

=head2

	PACKAGE		Agua::Common::Stage

	PURPOSE

		STAGE METHODS FOR Agua::Common

=cut


use Data::Dumper;

=head2

	SUBROUTINE		updateStage

	PURPOSE

		ADD A STAGE TO THE stage TABLE

	INPUTS

		1. STAGE IDENTIFICATION - PROJECT, WORKFLOW, STAGE AND STAGE NUMBER

		2. INFORMATION FOR RUNNING THE STAGE (name, location, etc.)

	OUTPUTS

		1. A NEW ENTRY IN THE stage TABLE

	NOTES

		THE PARAMETERS FOR THIS STAGE ARE ADDED TO THE stageparameter TABLE

		IN A SEPARATE CALL TO THE $self->addParameters() SUBROUTINE

=cut

sub updateStage {
	my $self		=	shift;


    my $dbobject        =	$self->dbobject();
    my $json 			=	$self->json();


	#### DO ADD STAGE WITH NEW NUMBER
	my $old_number = $json->{number};
	$json->{number} = $json->{newnumber};

	#### DO REMOVE STAGE WITH OLD NUMBER
	my $success = $self->_removeStage();
	print "{ status: 'Agua::Common::Stage::updateStage    Successfully removed stage $json->{name} from stage table' }\n" if $success;
	print "{ error: 'Agua::Common::Stage::updateStage    Could not remove stage $json->{name} from stage table' }\n" if not $success;

	$success = $self->_addStage();
	print "{ status: 'Agua::Common::Stage::updateStage    Successfully added stage $json->{name} into stage table' }\n" if $success;
	print "{ status: 'Agua::Common::Stage::updateStage    Could not add stage $json->{name} into stage table' }\n" if not $success;

	#### UPDATE THE APPNUMBER FIELD FOR ALL STAGE PARAMETERS
	#### BELONGING TO THIS PROJECT, WORKFLOW, APPNAME AND APPNUMBER
	#### NB: USE OLD NUMBER
	$json->{appnumber} = $old_number;
	$json->{appname} = $json->{name};
	my $unique_keys = ["username", "project", "workflow", "appname", "appnumber"];
	my $where = $dbobject->where($json, $unique_keys);
	my $query = qq{UPDATE stageparameter
SET appnumber='$json->{newnumber}'
$where};
	print "Agua::Common::Stage::updateStage    query: $query\n";
	$success = $dbobject->do($query);
	print "Agua::Common::Stage::updateStage    success: $success\n";

	print "{ status: 'Agua::Common::Stage::updateStage    Successful update of appnumber to $json->{newnumber} in stageparameter table' }\n" if $success;
	print "{ status: 'Agua::Common::Stage::updateStage    Could not update appnumber to $json->{newnumber} in stageparameter table' }\n" if not $success;
	exit;
}

=head2

    SUBROUTINE:     getStages

    PURPOSE:

        GET ARRAY OF HASHES:

			[
				appname1: [ s ], appname2 : [...], ...  ] 

=cut

sub getStages {
	my $self		=	shift;
	my $owner		=	shift;


    my $dbobject		=	$self->dbobject();
	my $json			=	$self->json();

	#### SET OWNER IF NOT DEFINED
	$owner = $json->{username} if not defined $owner;

    #### VALIDATE    
    print "{ error: 'User session not validated' }" and exit unless $self->validate();

	my $query = qq{SELECT * FROM stage
WHERE username='$owner'
ORDER BY project, workflow, number};
	my $stages = $dbobject->queryhasharray($query);
	$stages = [] if not defined $stages;


	return $stages;
}

=head2

	SUBROUTINE		addStage

	PURPOSE

		ADD A STAGE TO THE stage TABLE

	INPUTS

		1. STAGE IDENTIFICATION - PROJECT, WORKFLOW, STAGE AND STAGE NUMBER

		2. INFORMATION FOR RUNNING THE STAGE (name, location, etc.)

	OUTPUTS

		1. A NEW ENTRY IN THE stage TABLE

	NOTES

		NB: THE PARAMETERS FOR THIS STAGE ARE ADDED TO THE stageparameter TABLE

		IN A SEPARATE CALL TO THE $self->addParameters() SUBROUTINE

=cut

sub addStage {
	my $self		=	shift;


    my $dbobject        =	$self->dbobject();
    my $json 			=	$self->json();


	my $success = $self->_removeStage();
	$success = $self->_addStage();
	print "{ status: 'Agua::Common::Stage::addStage    Successful insert of stage $json->{name} into stage table' }" if $success;
	exit;
}



=head2

	SUBROUTINE		insertStage

	PURPOSE

		ADD A STAGE TO THE stage TABLE

	INPUTS

		1. STAGE OBJECT CONTAINING FIELDS: project, workflow, name, number

		2. INFORMATION FOR RUNNING THE STAGE (name, location, etc.)

	OUTPUTS

		1. A NEW ENTRY IN THE stage TABLE

	NOTES

		NB: THE PARAMETERS FOR THIS STAGE ARE ADDED TO THE stageparameter TABLE

		IN A SEPARATE CALL TO THE $self->addParameters() SUBROUTINE

=cut

sub insertStage {
	my $self		=	shift;


    my $dbobject        =	$self->dbobject();
    my $json 			=	$self->json();


	#### GET THE STAGES BELONGING TO THIS WORKFLOW
	my $where_fields = ["username", "project", "workflow"];
	my $where = $dbobject->where($json, $where_fields);
	my $query = qq{SELECT * FROM stage
$where};
	my $stages = $dbobject->queryhasharray($query);

	#### GET THE STAGE NUMBER 
	my $number = $json->{number};

	#### CHECK IF REQUIRED FIELDS ARE DEFINED
	my $required_fields = ["username", "owner", "project", "workflow", "name", "number"];
	my $not_defined = $dbobject->notDefined($json, $required_fields);
    print "{ error: 'Agua::Common::Stage::insertStage    undefined values: @$not_defined' }" and exit if @$not_defined;

	#### JUST ADD THE STAGE AND ITS PARAMETERS TO stage AND stageparameters
	#### IF THERE ARE NO EXISTING STAGES FOR THIS WORKFLOW
	if ( not defined $stages or scalar(@$stages) == 0 )
	{
		my $success = $self->_addStage();
		print "{ error: 'Agua::Common::Stage::insertStage    Could not insert stage $json->{name} into stage table' }" and exit if not $success;

		$success = $self->addParameters();
		print "{ error: 'Agua::Common::Stage::insertStage    Could not insert stage $json->{name} into stage table' }" and exit if not $success;

		print "{ status: 'Agua::Common::Stage::insertStage    Inserted stage $json->{name} into stage and stageparameter tables\n" if $success;

		exit;
	}

	#### INCREMENT THE number FOR DOWNSTREAM STAGES IN THE stage TABLE
	for ( my $i = @$stages - 1; $i > $number - 2; $i-- )
	{
		my $stage = $$stages[$i];
		my $new_number = $i + 2;
		my $where_fields = ["username", "project", "workflow", "number"];
		my $where = $dbobject->where($stage, $where_fields);
		my $query = qq{UPDATE stage SET
number='$new_number'
$where};
		my $success = $dbobject->do($query);
	}

	#### INCREMENT THE appnumber FOR DOWNSTREAM STAGES IN THE stageparameter TABLE
	for ( my $i = @$stages - 1; $i > $number - 2; $i-- )
	{
		my $stage = $$stages[$i];
		$stage->{appnumber} = $stage->{number};
		my $new_number = $i + 2;
		my $where_fields = ["username", "project", "workflow", "appnumber"];
		my $where = $dbobject->where($stage, $where_fields);
		my $query = qq{UPDATE stageparameter SET
appnumber='$new_number'
$where};
		my $success = $dbobject->do($query);
	}

	my $success = $self->_addStage();
	print "{ error: 'Agua::Common::Stage::insertStage    Could not insert stage $json->{name} into stage table' }" if not $success;

	$success = $self->addParameters();
	print "{ error: 'Agua::Common::Stage::insertStage    Could not insert stage parameters $json->{name} into stageparameters table' }" if not $success;


 	print "{ status: 'Agua::Common::Stage::insertStage    Inserted stage $json->{name} into stage and stageparameter tables\n" if $success;


	exit;
}


=head2

	SUBROUTINE		_addStage

	PURPOSE

		INTERNAL USE ONLY: ATOMIC ADDITION OF A STAGE TO THE stage TABLE

=cut

sub _addStage {
	my $self		=	shift;


    my $dbobject        =	$self->dbobject();
    my $json 			=	$self->json();


    #### VALIDATE
    print "{ error: 'Agua::Common::Stage::_addStage    User session not validated' }" and exit unless $self->validate();

	#### SET TABLE AND REQUIRED FIELDS	
	my $table = "stage";
	my $required_fields = ["username", "owner", "project", "workflow", "workflownumber", "name", "number", "type"];
	my $inserted_fields = ["username", "owner", "project", "workflow", "workflownumber", "name", "number", "location", "type", "executor", "cluster"];

	#### CHECK REQUIRED FIELDS ARE DEFINED
	my $not_defined = $dbobject->notDefined($json, $required_fields);
    print "{ error: 'Agua::Common::Stage::_addStage    undefined values: @$not_defined' }" and exit if @$not_defined;

	#### DO ADD
	my $success = $self->_addToTable($table, $json, $required_fields, $inserted_fields);	

	#### ADD IT TO THE report TABLE IF ITS A REPORT
	if ( $success and defined $json->{type} and $json->{type} eq "report" )
	{
		$json->{appname} = $json->{name};
		$json->{stagename} = $json->{name};
		$self->json($json);

		$success = $self->_addReport();
	}

	return $success;
}


=head2

	SUBROUTINE		_removeStage

	PURPOSE

		SAVE THE JSON (INPUTS, OUTPUTS AND ARGUMENTS) FOR EACH

        APPLICATION IN THE WORKFLOW, SENT AS A JSON STRING

=cut

sub _removeStage {
	my $self		=	shift;


    my $dbobject        =	$self->dbobject();
    my $json 			=	$self->json();



    #### VALIDATE
    print "{ error: 'Agua::Common::Stage::_removeStage    User session not validated' }" and exit unless $self->validate();

	#### CHECK UNIQUE FIELDS ARE DEFINED
	#### NB: ALSO CHECK name THOUGH NOT NECCESSARY FOR UNIQUE ID
	#### NNB: type IS NEEDED TO IDENTIFY IF ITS A REPORT
	my $required_fields = ["username", "project", "workflow", "number", "name", "type"];
	my $not_defined = $dbobject->notDefined($json, $required_fields);
    print "{ error: 'Agua::Common::Stage::_removeStage    undefined values: @$not_defined' }" and exit if @$not_defined;

	my $table = "stage";
	my $success = $self->_removeFromTable($table, $json, $required_fields);

	#### REMOVE IT FROM THE report TABLE TOO IF ITS A REPORT
	if ( $success and defined $json->{type} and $json->{type} eq "report" )
	{
		$json->{appname} = $json->{name};
		$json->{appnumber} = $json->{number};
		$self->json($json);

		$success = $self->_removeReport();
	}

	return $success;
}



=head2

	SUBROUTINE		removeStage

	PURPOSE

		REMOVE STAGE FROM this.stages

        REMOVE ASSOCIATED STAGE PARAMETER ENTRIES FROM this.stageparameter

=cut

sub removeStage {
	my $self		=	shift;


    my $dbobject        =	$self->dbobject();
    my $json 			=	$self->json();


    #### VALIDATE
    print "{ error: 'Agua::Common::Stage::removeStage    User session not validated' }" and exit unless $self->validate();

	#### GET THE STAGES BELONGING TO THIS WORKFLOW
	my $where_fields = ["username", "project", "workflow"];
	my $where = $dbobject->where($json, $where_fields);
	my $query = qq{SELECT * FROM stage
$where};
	my $stages = $dbobject->queryhasharray($query);
	if ( not defined $stages or scalar(@$stages) == 0 )
	{
	 	print "{ error: 'Agua::Common::Stage::removeStage    No stages in this workflow $json->{workflow} \n" and exit;
	}

#	#### GET THE STAGE PARAMETERS BELONGING TO THIS WORKFLOW
#	$query = qq{SELECT * FROM stageparameter
#$where};
#	my $stageparameters = $dbobject->queryhasharray($query);
#	if ( not defined $stageparameters or scalar(@$stageparameters) == 0 )
#	{
#	}

	#### REMOVE STAGE FROM stage TABLE
	my $success = $self->_removeStage();
 	print "{ error: 'Agua::Common::Stage::removeStage    Could not delete stage $json->{name} from stage table\n" and exit if not defined $success;

	#### REMOVE STAGE FROM stageparameter TABLE
	my $table2 = "stageparameter";
	$json->{appname} = $json->{name};
	$json->{appnumber} = $json->{number};
	my $required_fields2 = ["username", "project", "workflow", "appname", "appnumber"];
	$success = $self->_removeFromTable($table2, $json, $required_fields2);
 	print "{ error: 'Agua::Common::Stage::removeStage    Could not remove stage $json->{name} from $table2 table\n" and exit if not defined $success;

	#### QUIT IF THIS WAS THE LAST STAGE IN THE WORKFLOW
	my $number = $json->{number};
	if ( $number > scalar(@$stages) )
	{
	 	print "{ status: 'Agua::Common::Stage::addStage    Removed stage $json->{name} from stage and stageparameter tables\n";
		exit;
	}

	#### OTHERWISE, DECREMENT THE number FOR DOWNSTREAM STAGES IN THE stage TABLE
	for ( my $i = $number; $i < @$stages; $i++ )
	{
		my $stage = $$stages[$i];

		my $where_fields = ["username", "project", "workflow", "number"];
		my $where = $dbobject->where($stage, $where_fields);
		my $query = qq{UPDATE stage SET
number='$i'
$where};
		my $success = $dbobject->do($query);

		#### UPDATE report TABLE IF ITS A REPORT
		if ( $stage->{type} eq "report" )
		{
			$stage->{appname} = $stage->{name};
			$stage->{appnumber} = $stage->{number};
			my $where_fields = ["username", "project", "workflow", "appname", "appnumber"];
			my $where = $dbobject->where($stage, $where_fields);

			my $query = qq{UPDATE report SET
appnumber='$i'
$where};
			my $success = $dbobject->do($query);

		}
	}

	#### DECREMENT THE appnumber FOR DOWNSTREAM STAGES IN THE stageparameter TABLE
	for ( my $i = $number; $i < @$stages; $i++ )
	{
		my $stage = $$stages[$i];
		$stage->{appnumber} = $stage->{number};

		my $where_fields = ["username", "project", "workflow", "appnumber"];
		my $where = $dbobject->where($stage, $where_fields);
		my $query = qq{UPDATE stageparameter SET
appnumber='$i'
$where};
		my $success = $dbobject->do($query);
	}

 	print "{ status: 'Agua::Common::Stage::addStage    Removed stage $json->{name} from stage and stageparameter tables' }";

	exit;
}



1;