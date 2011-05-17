package Agua::Common::Admin;
use Moose::Role;
use Moose::Util::TypeConstraints;

=head2

	PACKAGE		Agua::Common::Admin

	PURPOSE

		ADMIN METHODS FOR Agua::Common

=cut


use Data::Dumper;

=head2

    SUBROUTINE		getHeadings

    PURPOSE

        VALIDATE THE USER, DETERIMINE IF THEY ARE THE

		ADMIN USER, THEN SEND THE APPROPRIATE LIST OF

		ADMIN PANES

=cut

sub getHeadings {
	my $self		=	shift;

    my $dbobject		=	$self->dbobject();
	my $json			=	$self->json();

	#### VALIDATE    
    my $username = $json->{username};
	print "{ error: 'Agua::Common::Admin::getHeadings    User $username not validated' }"
		and exit unless $self->validate($username);

	#### CHECK REQUESTOR
	print qq{ error: 'Agua::Common::Admin::getHeadings    Access denied to requestor: $json->{requestor}' } if defined $json->{requestor};

	my $headings = {

		#leftPane => [],
		#leftPane => [ "Users" ],
		##leftPane => [ "Access" ],
		#leftPane => ["Apps"],
		#leftPane => [ "Parameter" ],
		leftPane => [ "Clusters" ],
		#leftPane => [ "Apps", "Users" ],
		#leftPane => [ "Settings" ],
		#leftPane => [ "Settings", "Users" ],
		#leftPane => ["GroupProjects"],
		#leftPane => [ "Apps", "GroupProjects", "Settings" ],
		#leftPane => [ "Apps", "GroupProjects", "Settings", "Clusters" ],

		#middlePane => [],
		#middlePane => [ "Access" ],
		#middlePane => [ "Inputs" ],
		#middlePane => [ "Parameter" ],
		#middlePane => [ "Projects" ],
		middlePane => [ "Settings" ],
		#middlePane => [ "Groups" ],
		#middlePane => [ "Access" ],
		#middlePane => [ "Inputs" ],
		#middlePane => [ "Groups", "Access" ],
		#middlePane => [ "Groups", "Projects" ],
		##middlePane => [ "Groups", "Projects", "Access"],
		#middlePane => [ "Inputs", "Groups", "Projects", "Access"],
		#middlePane => [ "Inputs", "Groups", "Access" ],

		rightPane => []
		#rightPane => ["Users"]
		#rightPane => [ "Sources" ]
		#rightPane => ["GroupUsers"]
		#rightPane => [ "Sources", "GroupSources", "GroupUsers" ]
		#rightPane => [ "Sources", "GroupSources", "GroupUsers", "Users" ]
	};

	#### ADD ADMIN TABS: Users
	if ( $self->isAdminuser($username))
	{
		#push @{$headings->{rightPane}}, "Users";	
	}

    return $headings;
}




1;