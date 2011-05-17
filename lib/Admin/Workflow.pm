package Admin::Workflow;


=head2

		PACKAGE		Admin::Workflow

		PURPOSE

			THE Admin::Workflow OBJECT PERFORMS THE FOLLOWING	TASKS:

				1. CREATE, MODIFY OR DELETE USERS

				2. AUTHENTICATE USER ACCESS BY PASSWORD AND SESSION ID

=cut

use strict;
use warnings;
use Carp;

#### USE LIB (FOR INHERITANCE)
use lib "../";

#### INHERIT FROM App CLASS
require Exporter;
our @ISA = 'Admin';
use Admin;

#### INTERNAL MODULES
#use Util;
#use App::Workflow;

#### EXTERNAL MODULES
use Data::Dumper;

#### FLUSH BUFFER
#$| = 1;

##### GET ROOTDIR AND BINDIR
#use Cwd qw(realpath);
#my $fullpath = realpath("../../../");
##my $libdir = Cwd::pwd('../../');
##print "Lib dir: $libdir\n";
#exit;
#


#### SET SLOTS
our @DATA = qw(
    ROOT
	USERNAME
	SESSIONID
	PASSWORD
	DATABASE
	TYPE
	MODE
    DBOBJECT
	DBH
	CGI
    CONF
    FIELDS
    WORKFLOW_APP
);
our $DATAHASH;
foreach my $key ( @DATA )	{	$DATAHASH->{lc($key)} = 1;	}

our $ROOT = 'admin';

=head2

	SUBROUTINE		new

	PURPOSE

		CREATE THE NEW self OBJECT AND INITIALISE IT, FIRST WITH DEFAULT 

		ARGUMENTS, THEN WITH PROVIDED ARGUMENTS

        THE PASSWORD FOR THE DEFAULT ROOT USER CAN BE SET ON THE MYSQL

        COMMAND LINE:

        use myEST;
        insert into users values('admin', 'mypassword', now());

=cut

sub new
{
 	my $class 		=	shift;
	my $arguments 	=	shift;

	my $self = {};
    bless $self, $class;

	#### SET THE MAIN .conf VALUES - ROOTDIR, BINDIR, PIPELINE_GENERATION, ETC.
	$self->conf();

	#### INITIALISE THE OBJECT'S ELEMENTS
	$self->initialise($arguments);

    return $self;
}



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

    #### SET DEFAULT 'ROOT' USER
    $self->value('root', $ROOT);

    #### VALIDATE USER-PROVIDED ARGUMENTS
	($arguments) = $self->validate_arguments($arguments, $DATAHASH);	

	#### LOAD THE USER-PROVIDED ARGUMENTS
	foreach my $key ( keys %$arguments )
	{		
		#### LOAD THE KEY-VALUE PAIR
		$self->value($key, $arguments->{$key});
	}

    if ( not defined $self->{_dbobject} )
    {
        #### PRINT CONTENT TYPE
        print "Content-type: text/xml\n\n";
        die "Database object not defined\n";
    }
}



=head2

	SUBROUTINE		workflow_applications

	PURPOSE

		RETURN THE JSON FOR ALL STAGES OF THE GIVEN WORKFLOW

=cut

sub workflow_applications
{
    my $self            =   shift;

    my $dbobject         =	$self->{_dbobject};

	#### PRINT CONTENT TYPE
	print "Content-type: text/xml\n\n";

    my $project = $self->cgiParam('project');	
    my $workflow = $self->cgiParam('workflow');	

    my $query = qq{SELECT * FROM stage WHERE projectname='$project' AND workflowname = '$workflow' ORDER BY stagenumber};

    my $stages = $dbobject->simple_queryhasharray($query);
    if ( not defined $stages )
    {
        print '';
        exit;
    }

	#### GET THE WORKFLOW NUMBER, NAME AND DESCRIPTION, AND DATABASE, USER & PASSWORD
	#### ADD THE WORKFLOW INPUT AND OUTPUT FILES, ARGUMENTS, STDERR AND STDOUT FILES
	#### AND TYPE (perl) AND EXECUTABLE (/usr/bin/perl) SO THAT IT CAN BE DOWNLOADED
	#### AND USED AS AN APPLICATION .xml FILE (AS A COMPONENT IN ANOTHER WORKFLOW)

    use JSON;
    my $jsonObject = JSON->new();
    my $json = $jsonObject->objToJson($stages, {pretty => 1, indent => 4});
    print $json;

    exit;
}




=head2

	SUBROUTINE		save_application

	PURPOSE

		RUN THE WORKFLOW

=cut

sub save_application
{
    my $self            =   shift;
    my $dbh         =	$self->{_dbh};
    my $cgi         =   $self->{_cgi};
    my $conf        =   $self->{_conf};

    my $workflow_number =   $cgi->param('workflow_number');
    my $stage_number    =   $cgi->param('stage_number');
    my $field           =   $cgi->param('field');
    my $value           =   $cgi->param('value');

    #### PRINT CONTENT TYPE
	print "Content-type: text/html\n\n";

    my $query = qq{UPDATE stage SET $field = '$value' WHERE workflownumber = '$workflow_number' AND stagenumber = '$stage_number'};
    print "$query\n";
    my $result = Database::do_query($dbh, $query);

    print $result;
}



=head2

	SUBROUTINE		run

	PURPOSE

		RUN THE WORKFLOW

=cut

sub run
{
    my $self            =   shift;
    my $dbh         =	$self->{_dbh};
    my $cgi         =   $self->{_cgi};
    my $conf        =   $self->{_conf};

	my $user		=	$conf->getKeyValue("database", "USER");
	my $password	=	$conf->getKeyValue("database", "PASSWORD");
    my $database	    =	$cgi->param('database');
    my $workflow_number =   $cgi->param('workflow_number');

	#### FLUSH BUFFER
	$| = 1;

	print "Running...\n";

	my $logfile = "/var/www/html/sandbox21/logs/workflow.log";
	#open(STDOUT, ">$logfile") or die "Can't redirect STDOUT to log file: $logfile\n";
	#open(STDERR, ">$logfile") or die "Can't redirect STDOUT to log file: $logfile\n";
	open(LOG, ">$logfile");

	my $date = `date`;
	print "Datetime: $date\n";
	print LOG "$date\n";

	my $query = qq{SELECT workflowxmlfile
	FROM stage
	WHERE workflownumber = '$workflow_number'};

	print LOG "$query\n";
	print "$query\n";

	my $xmlfile = Database::simple_query($dbh, $query);

	print "XMLFILE: $xmlfile\n";
	print LOG "XMLFILE: $xmlfile\n";


	my $workflow = App::Workflow->new(
		{
			'XMLFILE'	=>	$xmlfile,
			'DATABASE'	=>	$database,
			'USER'		=>	$user,
			'PASSWORD'	=>	$password
		}
	);

	print LOG `pwd`;
	print LOG `whoami`;

	my $env = sprintf Dumper %$ENV;
	print LOG $env;

	#my $dumper = sprintf Dumper $workflow;

	print LOG $workflow->run();

	#my $success = $workflow->run();
	#
}


=head2

	SUBROUTINE		stop

	PURPOSE

		GET THE PROCESS IDS OF ANY RUNNING STAGE OF THE WORKFLOW AND

		'kill -9' THEM (INCLUDES STAGE PID, App PARENT PID AND App CHILD PID)

=cut

sub stop
{
    my $self            =   shift;
    my $dbh         =	$self->{_dbh};
    my $cgi         =   $self->{_cgi};
    my $conf        =   $self->{_conf};

    my $database	    =	$cgi->param('database');
    my $workflow_number =   $cgi->param('workflow_number');

	my $message = '';

	if ( not defined $database or not $database )
	{
		return "Database not defined\n";
	}
	elsif ( not defined $workflow_number or not $workflow_number )
	{
		return "Workflow number not defined\n";
	}


	#### GET PROCESS IDS OF RUNNING STAGE(S)
	my $query = qq{SELECT stagenumber, workflowpid, stagepid, parentpid,childpid
	FROM stage
    WHERE workflownumber='$workflow_number'
	AND status = 'running'};
	my $running_stages = Database::simple_queryhasharray($dbh, $query);

	#### RETURN IF NO PIDS
	if ( not defined $running_stages )	{	return "No running stages."	}

	foreach my $stage ( @$running_stages )
	{

		#### OTHERWISE, KILL ALL PIDS
		$message .= $self->kill_pid($stage->{childpid});
		$message .= $self->kill_pid($stage->{parentpid});
		$message .= $self->kill_pid($stage->{stagepid});
		$message .= $self->kill_pid($stage->{workflowpid});

		#### SET THIS STAGE STATUS AS 'waiting'
		my $update_query = qq{UPDATE stage
		SET status = 'waiting'
		WHERE workflownumber='$workflow_number'
		AND stagenumber = '$stage->{stagenumber}'};
		print "$update_query\n";
		Database::do_query($dbh, $update_query);

	}

	return $message;
}

=head2

	SUBROUTINE		stop

	PURPOSE

		KILL THE PID AND RETURN THE kill COMMAND

=cut

sub kill_pid
{
    my $self            =   shift;
    my $pid		=	shift;

	my $kill_command = "kill -9 $pid";
	print $kill_command, "\n";

	`$kill_command &> /dev/null`;

	return $kill_command;
}



=head2

	SUBROUTINE		print_filepath

	PURPOSE

		PRINT THE CONTENTS OF THE FILE GIVEN BY THE filepath

		OR THE FIRST FILE IN THE ls OF THE filepath IF THERE

		IS A '*' IN THE filepath

=cut

sub print_filepath
{
    my $self            =   shift;
    my $filepath		=	shift;

	#### IF ITS A MULTI-FILE ADDRESS, PRINT THE FIRST FILE
	if ( $filepath =~ /\*/ )
	{
		my @files = split "\n", `ls $filepath`;
		if ( @files )
		{
			$filepath = $files[0];
		}
		else
		{
			print "No files found at: $filepath\n";
			return;
		}
	}	

	if ( -f $filepath )
	{
		if ( -z $filepath )
		{
			print "File is empty: $filepath\n";
		}
		else
		{
			open(FILE, $filepath) or die "Can't open file: $filepath\n";
			while ( <FILE> )	{	print $_;	}
			close(FILE);
		}
	}
	else
	{
		print "Can't find file: $filepath\n";
	}

}


=head2

	SUBROUTINE		workflow

	PURPOSE

		RETURN THE XML FOR ALL STAGES OF THE GIVEN WORKFLOW

=cut

sub workflowStatus
{
    my $self            =   shift;
	my $workflow_number	=	shift;

    my $dbh	    =	$self->{_dbh};

	my $query = qq{SELECT status FROM stage
	WHERE workflownumber='$workflow_number'
	LIMIT 1};
	my $workflowStatus = Database::simple_query($dbh, $query);

	return $workflowStatus;
}


=head2

	SUBROUTINE		workflow

	PURPOSE

		RETURN THE XML FOR ALL STAGES OF THE GIVEN WORKFLOW

=cut

sub workflow
{
    my $self            =   shift;

    return $self->stage_xml(@_);
}





#=head2
#
#	SUBROUTINE		stage_json
#	
#	PURPOSE
#	
#		RETURN THE JSON FOR ALL STAGES OF THE GIVEN WORKFLOW
#		
#=cut
#
#sub stage_json
#{
#    my $self            =   shift;
#    
#    my $dbh         =	$self->{_dbh};
#    my $cgi         =   $self->{_cgi};
#    my $conf        =   $self->{_conf};
#
#    #### GET fields AND rootdir
#    my $fields = $self->{_fields};
#    my $rootdir = $self->{_rootdir};
#    
#	#### PRINT CONTENT TYPE
#
#    my $workflow_number = $cgi->param('workflow_number');	
#    if ( not defined $workflow_number or not $workflow_number )
#    {
#        $workflow_number = 1;
#    }
#    
#    #print "workflow.cgi::workflow(dbh, workflow_number, fields_string)\n";
#    #print "Workflow number: $workflow_number\n";
#    #print "Fields: $fields\n";
#    
#	my $query = qq{SELECT * FROM stage
#    WHERE workflownumber='$workflow_number'
#    ORDER BY workflownumber,stagenumber};
#    #print "$query\n";	
#    my $stages = Database::simple_queryhasharray($dbh, $query);
#	#print Dumper $stages;
#    if ( not defined $stages )
#    {
#        exit;
#    }
#    
#	#### GET THE WORKFLOW NUMBER, NAME AND DESCRIPTION, AND DATABASE, USER & PASSWORD
#	#### ADD THE WORKFLOW INPUT AND OUTPUT FILES, ARGUMENTS, STDERR AND STDOUT FILES
#	#### AND TYPE (perl) AND EXECUTABLE (/usr/bin/perl) SO THAT IT CAN BE DOWNLOADED
#	#### AND USED AS AN APPLICATION .xml FILE (AS A COMPONENT IN ANOTHER WORKFLOW)
#
#    my @fields = split ",", $fields;
#    my $workflow = qq{  <stages>\n};
#	
#	for ( my $stage_counter = 0; $stage_counter < @$stages; $stage_counter++ )
#	{
#		my $stage = $$stages[$stage_counter];
#        $workflow .= "    <stage>\n";        
#        for ( my $field_counter = 0; $field_counter < $#fields + 1; $field_counter++ )
#    	{
#            my $field = $fields[$field_counter];
#            
#    		if ( defined $stage->{$field} )
#    		{
#				#### REMOVE THE ROOTDIR VALUE
#				my $cleaned_value = $stage->{$field};
#				$cleaned_value =~ s/$rootdir\///g;
#				$workflow .= qq{      <$field>$cleaned_value</$field>\n};
#            }
#    		else
#    		{
#				$workflow .= qq{      <$field> </$field>\n};
#            }
#		}
#        $workflow .= "    </stage>\n";        
#	}
#	
#	$workflow .= qq{  </stages>
#};
#	
#    exit;
#}
#
#
#=head2
#
#	SUBROUTINE		application_xml
#	
#	PURPOSE
#	
#		RETURN THE XML FOR ALL STAGES OF THE GIVEN WORKFLOW
#		
#=cut

sub application_xml
{
    my $self            =   shift;

    my $dbh         =	$self->{_dbh};
    my $cgi         =   $self->{_cgi};
    my $conf        =   $self->{_conf};

    #### GET fields AND rootdir
    my $fields = $self->{_fields};
    my $rootdir = $self->{_rootdir};

	#### PRINT CONTENT TYPE
	print "Content-type: text/xml\n\n";

    my $workflow_number = $cgi->param('workflow_number');	
    if ( not defined $workflow_number or not $workflow_number )
    {
        $workflow_number = 1;
    }
    my $stage_number = $cgi->param('stage_number');	
    if ( not defined $stage_number or not $stage_number )
    {
        $stage_number = 1;
    }

	my $query = qq{SELECT * FROM stage
    WHERE workflownumber='$workflow_number'
    AND stagenumber = '$stage_number'};
    my $application = Database::simple_queryhash($dbh, $query);
    if ( not defined $application )
    {
        print '';
        exit;
    }

	#### GET THE WORKFLOW NUMBER, NAME AND DESCRIPTION, AND DATABASE, USER & PASSWORD
	#### ADD THE WORKFLOW INPUT AND OUTPUT FILES, ARGUMENTS, STDERR AND STDOUT FILES
	#### AND TYPE (perl) AND EXECUTABLE (/usr/bin/perl) SO THAT IT CAN BE DOWNLOADED
	#### AND USED AS AN APPLICATION .xml FILE (AS A COMPONENT IN ANOTHER WORKFLOW)

    my @fields = split ",", $fields;
    my $application_xml = qq{  <applications>\n};
    $application_xml .= "    <application>\n";    
    for ( my $field_counter = 0; $field_counter < $#fields + 1; $field_counter++ )
    {
        my $field = $fields[$field_counter];

        if ( defined $application->{$field} )
        {
            #### REMOVE THE ROOTDIR VALUE
            my $cleaned_value = $application->{$field};
            $cleaned_value =~ s/$rootdir\///g;
            $application_xml .= qq{      <$field>$cleaned_value</$field>\n};
        }
        else
        {
            $application_xml .= qq{      <$field> </$field>\n};
        }
    }
    $application_xml .= "    </application>\n";        
    $application_xml .= qq{  </applications>\n};

	print $application_xml;
    exit;
}



=head2

	SUBROUTINE		workflow_menu

	PURPOSE

		RETURN THE HTML FOR THE DATABASE SELECT BOX CONTAINING ALL DATABASE

		NAMES IN THE collectionworkflow TABLE OF THE myEST DATABASE

=cut


sub workflow_menu
{
    my $self            =   shift;

    my $dbh     =	$self->{_dbh};

	#### PRINT CONTENT TYPE
	print "Content-type: text/html\n\n";
	my $workflows = $self->workflows();
    if ( not defined $workflows )
    {
        print qq{<ul><li><h2 id="workflow" class="workflow_menu">No workflows</h2></li></ul>\n};
        exit;
    }

	my $workflow_menu = qq{	<ul>
		<li>
			<h2 id="workflow" class="workflow_menu">Workflow</h2>
			<ul>\n};

	for ( my $i = 0; $i < @$workflows; $i++ )
	{
		my $workflow = $$workflows[$i];
		my $workflow_number = $workflow->{workflownumber};
		my $workflow_name = $workflow->{workflowname};
		if ( defined $workflow )
		{
			$workflow_menu .= qq{		<li>  <a>$workflow_number. $workflow_name</a> </li>\n};
		}
	}

	$workflow_menu .= qq{			
			 </ul>
		 </li>
	</ul>\n};

	print $workflow_menu;
    exit;
}




=head2

	SUBROUTINE		stage_menu

	PURPOSE

		RETURN THE HTML FOR THE DATABASE SELECT BOX CONTAINING ALL DATABASE

		NAMES IN THE collectionworkflow TABLE OF THE myEST DATABASE

=cut


sub stage_menu
{
    my $self            =   shift;

    my $dbh         =	$self->{_dbh};
    my $cgi         =   $self->{_cgi};
    my $conf        =   $self->{_conf};

	#### PRINT CONTENT TYPE
	print "Content-type: text/html\n\n";

    my $workflow_number = $cgi->param('workflow_number');	
    if ( not defined $workflow_number or not $workflow_number )
    {
        $workflow_number = 1;
    }


	my $query = qq{SELECT * FROM stage
    WHERE workflownumber='$workflow_number'
    ORDER BY workflownumber,stagenumber};
    my $stages = Database::simple_queryhasharray($dbh, $query);
    if ( not defined $stages )
    {
        print '';
        exit;
    }

	my $stage_menu = qq{	<ul>
		<li>
			<h2 id="stage" class="stage_menu">Stages</h2>
			<ul>\n};

	for ( my $i = 0; $i < @$stages; $i++ )
	{
		my $stage = $$stages[$i];
		my $stage_number = $stage->{stagenumber};
		my $stage_name = $stage->{stagename};
		if ( defined $stage )
		{
			$stage_menu .= qq{		<li>  <a>$stage_number. $stage_name</a> </li>\n};
		}
	}

	$stage_menu .= qq{			
			 </ul>
		 </li>
	</ul>\n};

	print $stage_menu;
    exit;
}



=head2

	SUBROUTINE		workflow_select

	PURPOSE

		RETURN THE HTML FOR THE DATABASE SELECT BOX CONTAINING ALL DATABASE

		NAMES IN THE collectionworkflow TABLE OF THE myEST DATABASE

=cut


sub workflow_select
{
    my $self            =   shift;


    my $dbh         =	$self->{_dbh};
    my $cgi         =   $self->{_cgi};
    my $conf        =   $self->{_conf};

	my $workflows = workflows();
	if ( not defined $workflows )
	{
		my $workflow_select = qq{<SELECT>
		<option> No workflows </options>
</SELECT>\n};

		return $workflow_select;	
	}

	my $workflow_select = qq{	<select>\n};

	for ( my $i = 0; $i < @$workflows; $i++ )
	{
		my $workflow_number = $$workflows[$i]->{workflownumber};
		my $workflow_name = $$workflows[$i]->{workflowname};
		if ( defined $workflow_number )
		{
			$workflow_select .= qq{		<option value='$workflow_number'> $workflow_number. $workflow_name </option>\n};
		}
	}

	$workflow_select .= qq{	</select>
};

	return $workflow_select;
}



=head2

	SUBROUTINE		workflow_xml

	PURPOSE

		RETURN THE XML FOR ALL WORKFLOWS IN THIS DATABASE FROM THE

        workflow TABLE

=cut


sub workflow_xml
{
    my $self            =   shift;

    my $dbh         =	$self->{_dbh};
    my $cgi         =   $self->{_cgi};
    my $conf        =   $self->{_conf};

	#### PRINT CONTENT TYPE
	print "Content-type: text/xml\n\n";

    #### GET rootdir
    my $rootdir = $self->{_rootdir};

	my $workflows = $self->workflows();
	if ( not defined $workflows )
	{
		my $workflow_xml = qq{<workflows></workflows>};	
		return $workflow_xml;	
	}

	my $workflow_xml = qq{<workflows>\n};
	for ( my $i = 0; $i < @$workflows; $i++ )
	{
        $workflow_xml .= qq{    <workflow>\n};
        my @keys = keys ( %{$$workflows[$i]} );
        foreach my $key ( @keys )
        {
            if ( defined $$workflows[$i]->{$key} )
            {
                #### REMOVE THE ROOTDIR VALUE
                my $cleaned_value = $$workflows[$i]->{$key};
                $cleaned_value =~ s/$rootdir\///g;
                $workflow_xml .= qq{        <$key> $cleaned_value </$key>\n};
            }
        }
        $workflow_xml .= qq{    </workflow>\n};
	}
    $workflow_xml .= qq{</workflows>\n};

	print $workflow_xml;
    exit;
}


=head2

	SUBROUTINE		workflows

	PURPOSE

		RETURN A REFERENCE TO THE HASHARRAY OF WORKFLOW NUMBER AND NAMES

        IN THE collectionworkflow TABLE OF THE myEST DATABASE

=cut

sub workflows
{
    my $self            =   shift;

    my $dbh     = $self->{_dbh};
	my $query = qq{SELECT DISTINCT workflownumber, workflowname, workflowdescription, workflowxmlfile
    FROM stage ORDER BY workflownumber};
	my $workflows = Database::simple_queryhasharray($dbh, $query);

	return $workflows;
}
















=head2

    SUBROUTINE      validate

    PURPOSE

        CHECK SESSION ID AGAINST STORED SESSION IN sessions TABLE

        TO VALIDATE USER

=cut

sub validate
{
	my $self		=	shift;

	my $dbh			=	$self->{_dbh};
	my $username	=	$self->{_username};
	my $session_id	=	$self->{_sessionid};
	my $cgi			=	$self->{_cgi};

	#### GET FROM CGI IF NOT DEFINED
	if ( not defined $username )
	{
		$username =	$cgi->param('username');
	}
	if ( not defined $session_id )
	{
		$session_id =	$cgi->param('session_id');
	}

	my $query = qq{
	SELECT username FROM sessions
	WHERE username = '$username'
	AND sessionid = '$session_id'};
	my $valid = Database::simple_query($dbh, $query);

	if ( not defined $valid )	{	return	0;	}

	return 1;
}




=head2

    SUBROUTINE      login

    PURPOSE

        CHECK INPUT PASSWORD AGAINST STORED PASSWORD TO VALIDATE SESSION

        AND STORE SESSION ID IN sessions TABLE

=cut

sub login
{
	my $self	=	shift;

	my $dbh = $self->{_dbh};
	my $cgi = $self->{_cgi};
	#my $dbh		=	shift;
	#my $cgi		=	shift;	

	my $username	=	$cgi->param('username');
	my $password 	=	$cgi->param('password');
	my $session_id 	=	$cgi->param('session_id');
	my $database 	=	$cgi->param('database');

	#### EXIT IF SESSION ID NOT DEFINED
	if ( not defined $session_id or not $session_id )
	{
		exit;
		return 0;
	}


	#### GET STORED PASSWORD
	my $query = qq{SELECT password FROM users
WHERE username='$username'};
	my $stored_password = Database::simple_query($dbh, $query);

	#### CHECK FOR INPUT PASSWORD MATCHES STORED PASSWORD
	my $match = $password =~ /^$stored_password$/; 

	#### IF PASSWORD MATCHES, STORE SESSION ID AND RETURN '1'
	if ( $match )
	{
		#### CHECK IF THIS SESSION ID ALREADY EXISTS
		my $exists_query = qq{
		SELECT username FROM sessions
		WHERE username = '$username'
		AND sessionid = '$session_id'};
		my $exists = Database::simple_query($dbh, $exists_query);
        if ( defined $exists )
        {
        }
        else
        {
        }


		#### IF IT DOES EXIST, UPDATE THE TIME
		if ( defined $exists )
		{
			my $now = "DATETIME('NOW')";
			$now = "NOW()" if $self->get_conf()->{DBTYPE} =~ /^MYSQL$/i;

			my $update_query = qq{UPDATE session
			SET datetime = $now
			WHERE username = '$username'
			AND sessionid = '$session_id'};
			my $update_success = Database::do_query($dbh, $update_query);
		}

		#### IF IT DOESN'T EXIST, INSERT IT INTO THE TABLE
		else
		{
			my $query = qq{INSERT INTO sessions
			(username, sessionid, datetime)
			VALUES
			('$username', '$session_id', NOW() )};
			my $success = Database::do_query($dbh, $query);
			if ( $success )
			{
			}
		}		
	}

	#### CLEAN OUT OLD SESSIONS
	# DELETE FROM sessions WHERE datetime < ADDDATE(NOW(), INTERVAL -48 HOUR)
	# DELETE FROM sessions WHERE datetime < DATE_SUB(NOW(), INTERVAL 1 DAY)
	my $delete_query = qq{
	DELETE FROM sessions WHERE datetime < ADDDATE(NOW(), INTERVAL -24 HOUR)};
	Database::do_query($dbh, $delete_query);

	return $match;
}


=head2

    SUBROUTINE      newuser

    PURPOSE

        CHECK 'ROOT' NAME AND PASSWORD AGAINST INPUT VALUES. IF

        VALIDATED, CREATE A NEW USER IN THE users TABLE

=cut

sub newuser
{
    my $self            =   shift;

	my $dbh = $self->{_dbh};
	my $cgi = $self->{_cgi};


	my $username    =	$cgi->param('username');
	my $password 	=	$cgi->param('password');
	my $newuser     =	$cgi->param('newuser');
	my $newuserpassword	=	$cgi->param('newuserpassword');


	#### CHECK IF USER ALREADY EXISTS
	my $query = qq{SELECT username FROM users WHERE username='$newuser' LIMIT 1};
	my $exists_already = Database::simple_query($dbh, $query);
	if ( $exists_already )
	{
		print "User exists already: $username\n";
		exit;
	}

    $query = qq{INSERT INTO users VALUES ('$newuser', '$newuserpassword', NOW())};
    my $success = Database::do_query($dbh, $query);

    if ( $success )
    {
        print "New user created\n";
    }
    else
    {
        print "Failed to create new user\n";
    }
    exit;
}


=head2

    SUBROUTINE      registerdatabase

    PURPOSE

        CREATE A NEW DATABASE FILE SYSTEM AND REGISTER IT IN THE collections

        TABLE OR RE-REGISTER AN EXISTING DATABASE

=cut

sub registerdatabase
{
    my $self            =   shift;

    my $dbh = $self->{_dbh};
	my $cgi = $self->{_cgi};

	my $username 		=	$cgi->param('username');
	my $password 	=	$cgi->param('password');
	my $session_id 	=	$cgi->param('session_id');
	my $database 	=	$cgi->param('database');
    my $type        =	$cgi->param('type');
	my $description	=	$cgi->param('description');

    if ( not defined $database  or not $database )
    {
        print "Database not defined\n";
        exit;
    }

    if ( not $self->validate() )
    {
        print "User $username not validated\n";
        exit;
    }
    else
    {
    }


	#### MAKE description MYSQL SAFE
	$description =~ s/'/\\'/g;

	#### CHECK IF DATABASE ALREADY EXISTS
	my $query = qq{SELECT databasename FROM collections WHERE databasename='$database' LIMIT 1};
	my $exists_already = Database::simple_query($dbh, $query);
	if ( $exists_already )
	{
		print "Database name exists already\n";
		exit;
	}

    $query = qq{INSERT INTO collections
    (username, databasename, type, description, created)
    VALUES ('$username', '$database', '$type', '$description', NOW() )};
    my $success = Database::do_query($dbh, $query);
    if ( $success )
    {
        print "Created new database '$database'\n";
    }
    else
    {
        print "Failed to create new database '$database'\n";
    }
    exit;
}


=head2

    SUBROUTINE      deletedatabase

    PURPOSE

        CREATE A NEW DATABASE FILE SYSTEM AND REGISTER IT IN THE collections

        TABLE OR RE-REGISTER AN EXISTING DATABASE

=cut

sub deletedatabase
{
    my $self            =   shift;

    my $dbh = $self->{_dbh};
	my $cgi = $self->{_cgi};

	my $username 		=	$cgi->param('username');
	my $password 	=	$cgi->param('password');
	my $database 	=	$cgi->param('database');

    if ( not defined $database  or not $database )
    {
        print "Database not defined\n";
        exit;
    }

    if ( not $self->validate() )
    {
        print "User $username not validated\n";
        exit;
    }
    else
    {
    }


	#### CHECK IF DATABASE EXISTS
	my $query = qq{SELECT databasename FROM collections
    WHERE databasename='$database'
    AND username='$username'
    LIMIT 1};
    print "$query\n";
	my $exists = Database::simple_query($dbh, $query);
	if ( not $exists )
	{
		print "Database does not exist or access denied exists\n";
		exit;
	}

    $query = qq{DELETE FROM collections
    WHERE username='$username'
    AND databasename = '$database'};
    my $success = Database::do_query($dbh, $query);
    if ( $success )
    {
        print "Deleted database '$database'\n";
    }
    else
    {
        print "Failed to deleted database '$database'\n";
    }
    exit;
}







=head2

    SUBROUTINE      registerworkflow

    PURPOSE

        CREATE A NEW WORKFLOW ENTRY IN THE workflows TABLE

        TABLE OR RE-REGISTER AN EXISTING DATABASE

=cut

sub registerworkflow
{
    my $self            =   shift;

    my $dbh = $self->{_dbh};
	my $cgi = $self->{_cgi};

	my $username 		=	$cgi->param('username');
	my $password 	=	$cgi->param('password');
	my $session_id 	=	$cgi->param('session_id');
	my $database 	=	$cgi->param('database');
    my $type        =	$cgi->param('type');
	my $description	=	$cgi->param('description');

    if ( not defined $database  or not $database )
    {
        print "Database not defined\n";
        exit;
    }

    if ( not $self->validate() )
    {
        print "User $username not validated\n";
        exit;
    }
    else
    {
    }


	#### MAKE description MYSQL SAFE
	$description =~ s/'/\\'/g;

	#### CHECK IF DATABASE ALREADY EXISTS
	my $query = qq{SELECT databasename FROM collections WHERE databasename='$database' LIMIT 1};
	my $exists_already = Database::simple_query($dbh, $query);
	if ( $exists_already )
	{
		print "Database name exists already\n";
		exit;
	}

    $query = qq{INSERT INTO collections
    (username, databasename, type, description, created)
    VALUES ('$username', '$database', '$type', '$description', NOW() )};
    my $success = Database::do_query($dbh, $query);
    if ( $success )
    {
        print "Created new database '$database'\n";
    }
    else
    {
        print "Failed to create new database '$database'\n";
    }
    exit;
}


=head2

    SUBROUTINE      deleteworkflow

    PURPOSE

        CREATE A NEW DATABASE FILE SYSTEM AND REGISTER IT IN THE collections

        TABLE OR RE-REGISTER AN EXISTING DATABASE

=cut

sub deleteworkflow
{
    my $self            =   shift;

    my $dbh = $self->{_dbh};
	my $cgi = $self->{_cgi};

	my $username 		=	$cgi->param('username');
	my $password 	=	$cgi->param('password');
	my $database 	=	$cgi->param('database');



    if ( not defined $database  or not $database )
    {
        print "Database not defined\n";
        exit;
    }

    if ( not $self->validate() )
    {
        print "User $username not validated\n";
        exit;
    }
    else
    {
    }


	#### CHECK IF DATABASE EXISTS
	my $query = qq{SELECT databasename FROM collections
    WHERE databasename='$database'
    AND username='$username'
    LIMIT 1};
    print "$query\n";
	my $exists = Database::simple_query($dbh, $query);
	if ( not $exists )
	{
		print "Database does not exist or access denied exists\n";
		exit;
	}

    $query = qq{DELETE FROM collections
    WHERE username='$username'
    AND databasename = '$database'};
    my $success = Database::do_query($dbh, $query);
    if ( $success )
    {
        print "Deleted database '$database'\n";
    }
    else
    {
        print "Failed to deleted database '$database'\n";
    }
    exit;
}






1;




##### CONVERT THE ARGS ARRAY (EITHER A HASH OF PARAMETER '-KEY->VALUE' PAIRS OR AN
##### ALREADY ORDERED ARRAY OF VALUES PASSED TO THE CONSTRUCTOR) INTO AN ORDERED
##### ARRAY OF DATA VALUES FOR ONLY THE KEYS LISTED IN THE @data ARRAY
#
#sub _rearrange
#{
#	my $self		=	shift;
#	my $order		=	shift;
#	my $keyvalues	=	shift;
#		
#	#### RETURN THE DATA IF IT'S ALREADY ORDERED (I.E., NO '-' AT START OF FIRST ELEMENT IN ARRAY )	
#	# return @_ if not index($$keyvalues[0],'-');
#	
#	#### MAKE IT AN EVEN-NUMBERED ARRAY
#	push @$keyvalues, undef unless (scalar(@$keyvalues) + 1) % 2;
#	
#	my ($key, $value, %param);
#	while ( @$keyvalues )
#	{
##			$key =~ s/\-//g; 		# REMOVE FIRST '-' 
#		$key = shift(@$keyvalues);
#		$value = shift(@$keyvalues);
#		$param{$key} = $value;
#	}
#		
#	return map { $_ => $param{$_} } @$order;
#}
#	
## NB:	
##	%hash = map { getkey($_) => $_ } @array;
##	
##	is just a funny way to write
##	
##	%hash = ();
##	foreach $_ (@array) {
##		$hash{getkey($_)} = $_;
##	}
##

=head2

    SUBROUTINE:     cgiParam

    PURPOSE:    Return the value of the input CGI parameter

=cut

sub cgiParam
{
    my $self    =   shift;
    my $param   =   shift;

    if ( not defined $self->{_cgi} )
    {
        return;
    }

    return $self->{_cgi}->{$param}[0];
}


