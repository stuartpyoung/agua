package Agua::Common::Source;
use Moose::Role;
use Moose::Util::TypeConstraints;

=head2

	PACKAGE		Agua::Common::Source

	PURPOSE

		SOURCE METHODS FOR Agua::Common

=cut


use Data::Dumper;

=head2

    SUBROUTINE:     getSources

    PURPOSE:

		RETURN ALL ENTRIES IN THE source TABLE

			E.G.:
			[
				{
				  'location' : '/mihg/data/NGS/syoung/base/pipeline/human-genome-fa',
				  'name' : 'Human genome FASTA (36.1)',
				  'type' : 'fasta',
				  'description' : 'A single human genome file with one FASTA record pe
			r chromosome (Build 36.1)'
				},
				{
				  'location' : '/mihg/data/NGS/syoung/base/pipeline/human-chr/fa',
				  'name' : 'Human chromosome FASTA (36.1)',
				  'type' : 'fasta',
				  'description' : 'Human chromosome FASTA files (Build 36.1)'
				},
				{
					...
			]

=cut

sub getSources {

	my $self		=	shift;



    my $dbobject		=	$self->dbobject();
	my $json			=	$self->json();

	#### GET USERNAME AND SESSION ID
    my $username = $json->{'username'};

    #### VALIDATE    
    print "{ error: 'User session not validated' }" and exit unless $self->validate();

	#### GET ALL SOURCES
	my $query = qq{
SELECT * FROM source
WHERE username='$username'};
	my $sources = $dbobject->queryhasharray($query);
	$sources = [] if not defined $sources;


	return $sources;
}



=head2

    SUBROUTINE     sharedSources

    PURPOSE

		GET FILESYSTEM SOURCES SHARED BY OTHER USERS:

			1. FIND THE SHARED SOURCES

			2. CHECK THE USER BELONGS TO THE GROUP

=cut

sub sharedSources {	
    my $self        		=   shift;
    my $username			=   shift;



    #### GET DATABASE OBJECT
    my $dbobject = $self->dbobject();
    my $json = $self->json();

	$username = $json->{username} if not defined $username;

	my $sources;

	#### 1. FIND THE SHARED SOURCES
	my $admin = $self->conf()->getKeyValue('agua', 'ADMINUSER');
	my $query = qq{SELECT DISTINCT owner, groupname
FROM groupmember
WHERE type='source'
ORDER BY owner, groupname};
	my $sourcegroups = $dbobject->queryhasharray($query);

	#### 2. CHECK THE USER BELONGS TO THE GROUP
	$query = qq{SELECT owner, groupname
FROM groupmember
WHERE name='$username'
AND type='user'
};
	my $sharedgroups = $dbobject->queryhasharray($query);

	return [] if not defined $sourcegroups;
	return [] if not defined $sharedgroups;

	foreach my $sourcegroup ( @$sourcegroups )
	{

		my $found = 0; 
		foreach my $sharedgroup ( @$sharedgroups )
		{


			if ( $sharedgroup->{owner} eq $sourcegroup->{owner}
				and $sharedgroup->{groupname} eq $sourcegroup->{groupname} )
			{
				$found = 1;
			}

			if ( $found )
			{

				$query = qq{SELECT * FROM groupmember
WHERE owner='$sourcegroup->{owner}'
AND groupname='$sourcegroup->{groupname}'
AND type='source'
};
				my $groupsources = $dbobject->queryhasharray($query);
				foreach my $groupsource ( @$groupsources )
				{
					push @$sources, $groupsource;
				}
				last;
			}
		}
	}

	$sources = [] if not defined $sources;


	#### REMOVE DUPLICATE SOURCES
	my $already_seen;
	for ( my $i = 0; $i < @$sources; $i++ )
	{
		if ( not exists $already_seen->{$$sources[$i]->{name}} )
		{
			$already_seen->{$$sources[$i]->{name}} = 1;
		}
		else
		{
			splice (@$sources, $i, 1);
			$i--;
		}
	}

	return $sources;	
}




1;