package Report;

#### DEBUG

=head2

		PACKAGE		Report

		PURPOSE

			THE Reports OBJECT PERFORMS THE FOLLOWING TASKS:

				1. RETURNS FILE AND DIRECTORY LIST OF A GIVEN PATH AS A

                    dojox.data.FileStore JSON OBJECT TO BE RENDERED USING

                    FilePicker INSIDE dojox.dijit.Dialog


PROJECTS

cd C:\DATA\base\cgi-bin\Bioptic0.2.5

perl reports.cgi "mode=reports&sessionId=1228791394.7868.158&username=admin"


=cut

use strict;

require Exporter;
our @ISA = qw(Exporter Common);
our @EXPORT_OK = qw();
our $AUTOLOAD;

#### INTERNAL MODULES
use Util;
use Common;  # FOR AUTHENTICATION
use DBase::SQLite;

#### EXTERNAL MODULES
use File::Basename; # ALREADY IN PERL
use File::Copy;     # ALREADY IN PERL
use File::stat;     # ALREADY IN PERL
use Data::Dumper;   # ALREADY IN PERL
use Carp;           # ALREADY IN PERL
use JSON -support_by_pp;
use DBD::SQLite;

#### DOWNLOADED FROM CPAN
use File::Remove;
use File::Copy::Recursive;

#### SET SLOTS
our @DATA = qw(
	USERNAME
	SESSIONID
	JSON
    MODE
    DBOBJECT
	CGI    
    CONF
);
our $DATAHASH;
foreach my $key ( @DATA )	{	$DATAHASH->{lc($key)} = 1;	}
our $ROOT = 'admin';
our $DEFAULT_BYTES = 80;

=head2

	SUBROUTINE		new

	PURPOSE

		CREATE THE NEW self OBJECT AND INITIALISE IT, FIRST WITH DEFAULT 

		ARGUMENTS, THEN WITH PROVIDED ARGUMENTS

=cut


sub new {
 	my $class 		=	shift;
	my $arguments 	=	shift;

	my $self = {};
    bless $self, $class;

	#### CHECK CONF->FILEROOT IS PRESENT AND NOT EMPTY
	if ( not defined $arguments-> {CONF}
        or not $arguments-> {CONF}
        or not defined $arguments->{CONF}->{FILEROOT}
        or not $arguments->{CONF}->{FILEROOT}
    ) 
	{
		print "Conf->FILEROOT is not defined\n";
        exit;
    }
	else
	{
		$self->{_conf} = $arguments->{CONF};
	}

	#### INITIALISE THE OBJECT'S ELEMENTS
	$self->initialise($arguments);


    return $self;
}



=head2

    SUBROUTINE     reports

    PURPOSE

        1. Return a list of reports for the user

        2. Each report has permissions

=cut

#### Return a list of reports for the user
sub reports {
    my $self        =   shift;

    #### PRINT HEADER 
	print "Content-type: text/html\n\n";


    #### VALIDATE
    my $validated = $self->validate();
    if ( not $validated )
    {
        print "{ error: 'User session not validated'    }";
        exit;
    }

    #### GET DATABASE OBJECT
    my $dbobject = $self->{_dbobject};

    #### GET USER NAME, SESSION ID AND PATH FROM CGI OBJECT
    my $username = $self->cgiParam("username");

    #### GET THE PROJECTS OWNED BY THIS USER
    my $query = qq{SELECT * FROM reports WHERE owner = '$username'};
    my $user_reports = $dbobject->queryhasharray($query);    

    #### SET RIGHTS
    foreach my $report ( @$user_reports )
    {
        $report->{rights} = $report->{ownerrights};
    }


    ####
    #
    #   GENERATE LIST OF SOURCES (READABLE SHARED PROJECTS, SOME OF WHICH MAY BE WRITEABLE)
    #
    #    .schema groups
    #    CREATE TABLE groups
    #    (
    #            owner                   VARCHAR(20),
    #            groupname               VARCHAR(20),
    #            name                    VARCHAR(20),
    #            type                    VARCHAR(50),
    #    
    #            PRIMARY KEY (owner, groupname, name, type)
    #    );
    #    
    #    
    #    1. GET OWNER AND GROUP OF ANY GROUPS THE USER BELONGS TO
    #    select owner, group from groupmember where name='admin';
    #    syoung|bioinfo
    #    syoung|nextgen
    #
    #    2. GET THE NAMES OF ANY PROJECTS IN THE GROUPS THE USER BELONGS TO
    #    select name from groupmember where owner = 'syoung' and groupname = 'bioinfo' and type='report';
    #    Report1
    #    
    #    select name from groupmember where owner = 'syoung' and groupname = 'nextgen' and type='report';
    #    Report2
    #    
    #    PUSH TO ARRAY AND HASH FOR CHECKING IN WORLD-READABLE IF ALREADY USED
    #    
    #    
    #    3. GET THE PERMISSIONS AND LOCATIONS OF ANY SHARED PROJECTS IN THE USER'S GROUPS
    #
    #    .schema reports
    #    CREATE TABLE reports
    #    (
    #            report                 VARCHAR(20),
    #            owner                   VARCHAR(20),
    #            ownerrights             INT(1),
    #            grouprights             INT(1),
    #            worldrights             INT(1),
    #            location                TEXT,
    #    
    #            PRIMARY KEY (report, owner, ownerrights, grouprights, worldrights, loca
    #    tion)
    #    );
    #    
    #    select * from reports where owner = 'syoung' and report = 'Report1';
    #    Report1|syoung|7|5|1|syoung/Report1
    #    
    #    select * from reports where owner = 'syoung' and report = 'Report2';
    #    Report2|syoung|7|5|1|syoung/Report2
    #    
    #    

    #   1. GET OWNER AND GROUP OF ANY GROUPS THE USER BELONGS TO
    $query = qq{select owner, groupname from groupmember where name = '$username'};
    my $ownergroups = $dbobject->queryhasharray($query);

    my $owner_reports;         #### ARRAY
    my $used_owner_reports;    #### HASH 

    foreach my $hash( @$ownergroups )
    {
        my $owner = $hash->{owner};
        my $group = $hash->{groupname};


        #   2. GET THE NAMES OF ANY PROJECTS IN THE GROUPS THE USER BELONGS TO
        $query = qq{select name from groupmember where owner = '$owner' and groupname = '$group' and type='report'};
        my $queryarray = $dbobject->queryarray($query);
        foreach my $report ( @$queryarray )
        {
            $hash->{report} = $report;
            push @$owner_reports, $hash;

            #### PUSH TO 'USED OWNER PROJECTS' FOR CHECKING IN WORLD-READABLE LIST LATER
            $used_owner_reports->{"$owner-$report"} = 1;
        }
    }




    #   3. GET THE PERMISSIONS AND LOCATIONS OF ANY SHARED PROJECTS IN THE USER'S GROUPS
    my $sources;
    for my $owner_report ( @$owner_reports )
    {
        my $owner = $owner_report->{owner};
        my $report = $owner_report->{report};

        my $query = qq{select * from reports where owner = '$owner' and report = '$report'};
        my $hash = $dbobject->queryhash($query);
        $hash->{rights} = $hash->{grouprights};
        push @$sources, $hash;
    }


    ####
    #
    #   GENERATE LIST OF WORLD READABLE SOURCES (SOME OF WHICH MAY HAVE APPEARED
    #
    #   AMONG THE READABLE SHARED PROJECTS LIST ABOVE)
    #
    #    select * from reports where worldrights >= 1;
    #  
    #    Report1|syoung|7|5|1|syoung/Report1
    #    Report2|syoung|7|5|1|syoung/Report2
    #    Report0|admin|7|5|1|admin/Report0
    #    Report1|admin|7|5|1|admin/Report1
    #
    #

    $query = qq{select * from reports where worldrights >= 1 and owner != "$username"};
    my $world_sources = $dbobject->queryhasharray($query);
    for my $world_source ( @$world_sources )
    {
        my $key = "$world_source->{owner}-$world_source->{report}";
        if ( exists $used_owner_reports->{$key} )
        {
            next;
        }
        else
        {
            $world_source->{rights} = $world_source->{worldrights};

            push @$sources, $world_source;
        }
    }


    #### GENERATE PROJECTS JSON
    my $reports;
    $reports->{reports} = $user_reports;
    $reports->{sources} = $sources;
    my $jsonObject = JSON->new();
    #no warnings;
    my $json = $jsonObject->objToJson($reports, {pretty => 1, indent => 4});
    #use warnings;

    print $json;

    #### EXIT ON COMPLETION
    exit;    
}






=head2

    SUBROUTINE:     reportRights

    PURPOSE:

        1. Check the rights of a user to access another user's report

=cut
#### Check the rights of a user to access another user's report
sub reportRights {
    my $self        =   shift;    
    my $owner       =   shift;
    my $report     =   shift;
    my $username    =   shift;

    my $dbobject = $self->{_dbobject};

    #### GET GROUP NAME FOR PROJECT
    my $query = qq{SELECT groupname from groupmember where owner = '$owner' and name = '$report' and type = 'report'};
    print "$query\n";
    my $group = $dbobject->query($query);
    if ( not defined $group )
    {
        return;
    }

    #### CONFIRM THAT USER BELONGS TO THIS GROUP
    $query = qq{SELECT * FROM groupmember WHERE groupname = '$group' AND owner = '$owner' AND name = '$username' AND type = 'user'};
    my $access = $dbobject->query($query);
    if ( not defined $access )
    {
        return;
    }

    #### MAKE SURE THAT THE GROUP PRIVILEGES ALLOW ACCESS
    $query = qq{SELECT grouprights FROM reports WHERE owner = '$owner' AND report = '$report'};
    my $rights = $dbobject->query($query);
    if ( not defined $rights or not $rights )
    {
        return;
    }

    return $rights;
}


1;


