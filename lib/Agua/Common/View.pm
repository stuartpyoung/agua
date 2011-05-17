package Agua::Common::View;

=head2

	PACKAGE		Agua::Common::View

	PURPOSE

		VIEW METHODS FOR Agua::Common

=cut


use Moose::Role;
use Moose::Util::TypeConstraints;

#### EXTERNAL MODULES
use Data::Dumper;



=head2

    SUBROUTINE     getViews

    PURPOSE

        RETURN AN ARRAY OF VIEW HASHES FOR THIS USER

=cut

sub getViews {
    my $self        =   shift;


    my $dbobject = $self->dbobject();
    my $json = $self->json();    

	#### VALIDATE USER USING SESSION ID	
	print "{ error: 'Agua::Common::View::getView    User not validated' }" and exit unless $self->validate();

	#### GET VIEWS    
    my $username = $json->{'username'};
    my $query = qq{SELECT * FROM view
WHERE username = '$username'};
    my $views = $dbobject->queryhasharray($query);    

	$views = [] if not defined $views;

	return $views;
}

=head2

    SUBROUTINE     getViewFeatures

    PURPOSE

        RETURN AN ARRAY OF VIEW HASHES FOR THIS USER

=cut

sub getViewFeatures {
    my $self        =   shift;


    my $dbobject = $self->dbobject();
    my $json = $self->json();    

	#### VALIDATE USER USING SESSION ID	
	print "{ error: 'Agua::Common::View::getViewFeatures    User not validated' }" and exit unless $self->validate();

	#### GET VIEWS    
    my $username = $json->{'username'};
    my $query = qq{SELECT * FROM viewfeature
WHERE username = '$username'};
    my $viewfeatures = $dbobject->queryhasharray($query);    

	$viewfeatures = [] if not defined $viewfeatures;

	return $viewfeatures;
}

=head2

    SUBROUTINE     getFeatures

    PURPOSE

        RETURN AN ARRAY OF VIEW HASHES FOR THIS USER

=cut

sub getFeatures {
    my $self        =   shift;


    my $dbobject = $self->dbobject();
    my $json = $self->json();    

	#### VALIDATE USER USING SESSION ID	
	print "{ error: 'Agua::Common::View::getFeature    User not validated' }" and exit unless $self->validate();

	#### GET VIEWS    
    my $username = $json->{'username'};
    my $query = qq{SELECT project, workflow, feature, species, build
FROM feature
WHERE username = '$username'
AND workflow !=''};
    my $viewfeatures = $dbobject->queryhasharray($query);    

	$viewfeatures = [] if not defined $viewfeatures;

	return $viewfeatures;
}




=head2

    SUBROUTINE     updateView

	PURPOSE

		ADD A VIEW TO THE view TABLE

=cut

sub updateView {
	my $self		=	shift;


    my $json 			=	$self->json();


	#### REMOVE IF EXISTS ALREADY
	my $success = $self->_removeView();
 	print "{ error: 'Agua::Common::View::updateView    Could not remove view $json->{name} from view table\n" and exit if not defined $success and not $success;

	$success = $self->_addView();	
 	print "{ error: 'Agua::Common::View::updateView    Could not update view $json->{name} to view table\n" and exit if not defined $success;
 	print "{ status: 'Agua::Common::View::updateView    Updated view $json->{name} in view table\n";
	exit;
}

=head2

    SUBROUTINE     removeView

	PURPOSE

		REMOVE A VIEW FROM THE view TABLE

=cut

=head2

    SUBROUTINE     _addView

	PURPOSE

		ADD A VIEW TO THE view TABLE

=cut

sub _addView {
	my $self		=	shift;

	print "Agua::Common::View::_addView    Agua::Common::View::_addView()\n";

    my $dbobject        =	$self->dbobject();
    my $json 			=	$self->json();

	#### SET TABLE AND REQUIRED FIELDS	
	my $table = "view";
	my $required_fields = ["username", "project", "view"];

	#### CHECK REQUIRED FIELDS ARE DEFINED
	my $not_defined = $dbobject->notDefined($json, $required_fields);
    print "{ error: 'Agua::Common::View::_addView    undefined values: @$not_defined' }" and exit if @$not_defined;

	#### LATER: FIX THIS SO THAT THE DATE IS ADDED PROPERLY 
	$json->{datetime} = "NOW()";
	if ( $self->conf()->getKeyValue('database', 'DBTYPE') eq "SQLite" )
	{
		$json->{datetime} = "DATETIME('NOW');";
	}

	#### DO THE ADD
	return $self->_addToTable($table, $json, $required_fields);	
}



=head2

    SUBROUTINE     _removeView

	PURPOSE

		REMOVE A VIEW FROM THE view TABLE

=cut

sub _removeView {
	my $self		=	shift;


    my $dbobject        =	$self->dbobject();
    my $json 			=	$self->json();


	#### SET TABLE AND REQUIRED FIELDS	
	my $table = "view";
	my $required_fields = ["username", "project", "view"];

	#### CHECK REQUIRED FIELDS ARE DEFINED
	my $not_defined = $dbobject->notDefined($json, $required_fields);
    print "{ error: 'Agua::Common::View::_removeView    undefined values: @$not_defined' }" and exit if @$not_defined;

	#### DO THE REMOVE
	return $self->_removeFromTable($table, $json, $required_fields);
}








sub _addFeature {
	my $self		=	shift;
	my $object		=	shift;

	print "Agua::Common::View::_addFeature    Agua::Common::View::_addFeature(object)\n";
	print "Agua::Common::View::_addFeature    object:\n";
	print Dumper $object;

    my $dbobject        =	$self->dbobject();

	#### SET TABLE AND REQUIRED FIELDS	
	my $table = "feature";
	my $required_fields = ["username", "project", "workflow", "feature", "location"];

	#### CHECK REQUIRED FIELDS ARE DEFINED
	my $not_defined = $dbobject->notDefined($object, $required_fields);
    print "{ error: 'Agua::Common::View::_addFeature    undefined values: @$not_defined' }" and exit if @$not_defined;

	#### DO THE ADD
	return $self->_addToTable($table, $object, $required_fields);	
}

sub _removeFeature {
	my $self		=	shift;
	my $object		=	shift;

	print "Agua::Common::View::_removeFeature    Agua::Common::View::_removeFeature(object)\n";
	print "Agua::Common::View::_removeFeature    object:\n";
	print Dumper $object;

    my $dbobject        =	$self->dbobject();

	#### SET TABLE AND REQUIRED FIELDS	
	my $table = "feature";
	my $required_fields = ["username", "project", "workflow", "feature"];

	#### CHECK REQUIRED FIELDS ARE DEFINED
	my $not_defined = $dbobject->notDefined($object, $required_fields);
    print "{ error: 'Agua::Common::View::_removeFeature    undefined values: @$not_defined' }" and exit if @$not_defined;

	#### DO THE REMOVE
	return $self->_removeFromTable($table, $object, $required_fields);	
}


=head2

    SUBROUTINE     _addViewFeature

	PURPOSE

		ADD A VIEW TO THE view TABLE

=cut

sub _addViewFeature {
	my $self		=	shift;
	my $object		=	shift;

	print "Agua::Common::View::_addViewFeature    Agua::Common::View::_addViewFeature(featureObject)\n";

    my $dbobject        =	$self->dbobject();
    my $json 			=	$self->json();

	#### SET TABLE AND REQUIRED FIELDS	
	my $table = "viewfeature";
	my $required_fields = ["username", "project", "view", "feature", "location"];
	#my $required_fields = $dbobject->required_fields($table);

	#### CHECK REQUIRED FIELDS ARE DEFINED
	my $not_defined = $dbobject->notDefined($object, $required_fields);
    print "{ error: 'Agua::Common::View::_addViewFeature    undefined values: @$not_defined' }" and exit if @$not_defined;

	#### DO THE ADD
	return $self->_addToTable($table, $object, $required_fields);	
}




=head2

    SUBROUTINE     removeViewFeature

	PURPOSE

		REMOVE A VIEW FROM THE view TABLE

=cut

sub removeViewFeature {
	my $self		=	shift;


    my $json 			=	$self->json();

	#### DO THE REMOVE
	my $success = $self->_removeViewFeature();

 	print "{ error: 'Agua::Common::View::removeViewFeature    Could not remove feature $json->{feature} from viewtable table' }" and exit if not defined $success;
 	print "{ status: 'Agua::Common::View::removeViewFeature    Successfully removed feature $json->{feature} from viewfeature table' }";
	exit;
}

=head2

    SUBROUTINE     _removeViewFeature

	PURPOSE

		REMOVE A VIEW FROM THE view TABLE

=cut

sub _removeViewFeature {
	my $self		=	shift;


    my $dbobject        =	$self->dbobject();
    my $json 			=	$self->json();


	#### SET TABLE AND REQUIRED FIELDS	
	my $table = "viewfeature";
	my $required_fields = ["username", "project", "view", "feature"];

	#### CHECK REQUIRED FIELDS ARE DEFINED
	my $not_defined = $dbobject->notDefined($json, $required_fields);
    print "{ error: 'Agua::Common::View::_removeViewFeature    undefined values: @$not_defined' }" and exit if @$not_defined;

	#### DO THE REMOVE
	return $self->_removeFromTable($table, $json, $required_fields);
}



1;
