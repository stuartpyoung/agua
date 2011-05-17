package View;



=head2

	PACKAGE		View

	PURPOSE

		THE View OBJECT PERFORMS THE FOLLOWING TASKS:

			1. ADD/REMOVE Views

			2. ALTER Views


perl View.cgi < View-ViewApplications.json


=cut

use strict;
use warnings;
use Carp;

#### USE LIB FOR INHERITANCE
use FindBin qw($Bin);
use lib "$Bin";
use lib "$Bin/external";


#### INHERIT FROM App CLASS
require Exporter;
our @ISA = 	qw(Common);

#### INTERNAL MODULES
use Admin;
use Common;
use Monitor::PBS;

#### EXTERNAL MODULES
use Data::Dumper;

#### FLUSH BUFFER
$| = 1;


#### SET SLOTS
our @DATA = qw(
    ROOT
	USERNAME
	SESSIONID
	PASSWORD
	DATABASE
	TYPE

	JSON
    MODE
    DBOBJECT
	CGI

	CLUSTER
	QSUB
	QSTAT
	SETUID

    CONF
    FIELDS
    View_APP
);
our $DATAHASH;
foreach my $key ( @DATA )	{	$DATAHASH->{lc($key)} = 1;	}

our $ROOT = 'admin';



#=head2
#
#	SUBROUTINE		views
#	
#	PURPOSE
#	
#		RETURN A REFERENCE TO THE HASHARRAY OF View NUMBER AND NAMES
#        
#        IN THE collectionView TABLE OF THE myEST DATABASE
#		
#=cut
#
#sub views
#{
#    my $self            =   shift;
#
#    my $dbh     = $self->{_dbh};
#	my $query = qq{SELECT DISTINCT project, Viewnumber, View, number, name
#    FROM stage ORDER BY Viewnumber, number};
#	my $views = Database::simple_queryhasharray($dbh, $query);
#	
#	return $views;
#}


################################################################################
##################			HOUSEKEEPING SUBROUTINES			################
################################################################################

=head2

	SUBROUTINE		initialise

	PURPOSE

		INITIALISE THE self OBJECT:

			1. LOAD THE DATABASE, USER AND PASSWORD FROM THE ARGUMENTS

			2. FILL OUT %VARIABLES% IN XML AND LOAD XML

			3. LOAD THE ARGUMENTS

=cut


sub initialise
{
    my $self		=	shift;
	my $arguments	=	shift;

    ######## SET DEFAULT 'ROOT' USER
    ####$self->value('root', $ROOT);

    #### VALIDATE USER-PROVIDED ARGUMENTS
	($arguments) = $self->validate_arguments($arguments, $DATAHASH);	

	#### LOAD THE USER-PROVIDED ARGUMENTS
	foreach my $key ( keys %$arguments )
	{		
		#### LOAD THE KEY-VALUE PAIR
		$self->value($key, $arguments->{$key});
	}

#use Data::Dumper;

    if ( not defined $self->{_dbobject} )
    {
        #### PRINT CONTENT TYPE
        print "View::Content-type: text/html\n\n";
        die "Database object not defined\n";
    }
}





1;


