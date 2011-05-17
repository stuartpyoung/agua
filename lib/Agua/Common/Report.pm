package Agua::Common::Report;
use Moose::Role;
use Moose::Util::TypeConstraints;

=head2

	PACKAGE		Agua::Common::Report

	PURPOSE

		REPORT METHODS FOR Agua::Common

=cut


use Data::Dumper;

##############################################################################
#				REPORT METHODS
##############################################################################
=head2

    SUBROUTINE     getReports

    PURPOSE

        RETURN AN ARRAY OF REPORT HASHES FOR THIS USER

sub getReports
{
    my $self        =   shift;



    my $dbobject = $self->dbobject();
    my $json = $self->json();    

	#### VALIDATE USER USING SESSION ID	
	print "{ error: 'User not validated' }" and exit unless $self->validate();

	#### GET REPORTS    
    my $username = $json->{'username'};
    my $query = qq{SELECT * FROM report
WHERE username = '$username'};
    my $reports = $dbobject->queryhasharray($query);    

	$reports = [] if not defined $reports;

	return $reports;
}

=cut


=head2

    SUBROUTINE     updateReport

	PURPOSE

		ADD A REPORT TO THE report TABLE

=cut

sub updateReport {
	my $self		=	shift;


    my $json 			=	$self->json();


	#### REMOVE IF EXISTS ALREADY
	my $success = $self->_removeReport();
 	print "{ error: 'Agua::Common::Report::updateReport    Could not remove report $json->{name} from report table\n" and exit if not defined $success and not $success;

	$success = $self->_addReport();	
 	print "{ error: 'Agua::Common::Report::updateReport    Could not update report $json->{name} to report table\n" and exit if not defined $success;
 	print "{ status: 'Agua::Common::Report::updateReport    Updated report $json->{name} in report table\n";
	exit;
}


=head2

    SUBROUTINE     addReport

	PURPOSE

		ADD A REPORT TO THE report TABLE

=cut

sub addReport {
	my $self		=	shift;



    my $json 			=	$self->json();

	#### DO THE ADD
	my $success = $self->_addReport();

 	print "{ error: 'Agua::Common::Report::addReport    Could not add report $json->{name} to report table\n" and exit if not defined $success;
 	print "{ status: 'Agua::Common::Report::addReport    Added report $json->{name} to report table\n";
	exit;
}




=head2

    SUBROUTINE     removeReport

	PURPOSE

		REMOVE A REPORT FROM THE report TABLE

=cut

sub removeReport {
	my $self		=	shift;


    my $json 			=	$self->json();

	#### DO THE REMOVE
	my $success = $self->_removeReport();

 	print "{ error: 'Agua::Common::Report::removeReport    Could not add report $json->{name} to report table\n" and exit if not defined $success;
 	print "{ status: 'Agua::Common::Report::removeReport    Added report $json->{name} to report table\n";
	exit;
}



=head2

    SUBROUTINE     _addReport

	PURPOSE

		ADD A REPORT TO THE report TABLE

=cut

sub _addReport {
	my $self		=	shift;


    my $dbobject        =	$self->dbobject();
    my $cgi 			=	$self->cgi();
    my $json 			=	$self->json();


	#### SET TABLE AND REQUIRED FIELDS	
	my $table = "report";
	my $required_fields = ["username", "project", "workflow", "appname", "appnumber", "name"];

	#### CHECK REQUIRED FIELDS ARE DEFINED
	my $not_defined = $dbobject->notDefined($json, $required_fields);
    print "{ error: 'Agua::Common::Report::_addReport    undefined values: @$not_defined' }" and exit if @$not_defined;

	#### LATER: FIX THIS SO THAT THE DATE IS ADDED PROPERLY 
	$json->{datetime} = "NOW()";
	if ( $self->conf()->getKeyValue('agua', 'DBTYPE') eq "SQLite" )
	{
		$json->{datetime} = "DATETIME('NOW');";
	}

	#### DO THE ADD
	return $self->_addToTable($table, $json, $required_fields);	
}



=head2

    SUBROUTINE     _removeReport

	PURPOSE

		REMOVE A REPORT FROM THE report TABLE

=cut

sub _removeReport {
	my $self		=	shift;


    my $dbobject        =	$self->dbobject();
    my $cgi 			=	$self->cgi();
    my $json 			=	$self->json();


	#### SET TABLE AND REQUIRED FIELDS	
	my $table = "report";
	my $required_fields = ["username", "project", "workflow", "appname", "appnumber", "name"];

	#### CHECK REQUIRED FIELDS ARE DEFINED
	my $not_defined = $dbobject->notDefined($json, $required_fields);
    print "{ error: 'Agua::Common::Report::_removeReport    undefined values: @$not_defined' }" and exit if @$not_defined;

	#### DO THE REMOVE
	return $self->_removeFromTable($table, $json, $required_fields);
}



1;