package Agua::Common::Access;
use Moose::Role;
use Moose::Util::TypeConstraints;

=head2

	PACKAGE		Agua::Common::Access

	PURPOSE

		ACCESS METHODS FOR Agua::Common

=cut


use Data::Dumper;

##############################################################################
#				ACCESS METHODS
##############################################################################

=head2

    SUBROUTINE:     getAccess

    PURPOSE:

        RETURN A JSON STORE WITH THE FOLLOWING FORMAT    

			{
				//identifier: 'name',
				label: 'name',
				items: [
					{
						name: 'Project1',
						type: 'project',
						fullname: 'Next Generation Analysis Tools',
						groupname: 'bioinfo',
						groupdesc: 'Bioinformatics Team'
					},
					{   name: 'skhuri',
						type: 'user',
						fullname: 'Sawsan Khuri',
						groupname: 'bioinfo',
						groupdesc: 'Bioinformatics Team'
					},
					...
				]	
			}

=cut

sub getAccess {

	my $self		=	shift;

    my $dbobject		=	$self->dbobject();
	my $cgi 			=	$self->cgi();
	my $json			=	$self->json();

	my $jsonParser = JSON->new();

	#### VALIDATE    
    my $username = $json->{'username'};
    my $session_id = $json->{'sessionId'};
    if ( not $self->validate($username, $session_id) )
    {
        print "{ error: 'User $username not validated' }";
        exit;
    }

    my $query = qq{SELECT * FROM access WHERE owner='$username' ORDER BY groupname};
    my $access = $dbobject->queryhasharray($query);
    $access = [] if not defined $access;

    return $access;
}




1;