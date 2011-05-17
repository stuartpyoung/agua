use MooseX::Declare;

=head2

		PACKAGE		Project

		PURPOSE

			THE Project OBJECT PERFORMS THE FOLLOWING TASKS:

				1. RETURNS THE FILE AND DIRECTORY LISTINGS OF A GIVEN PATH AS A

                    dojox.data.FileStore JSON OBJECT (TO BE RENDERED USING

                    FilePicker INSIDE dojox.dijit.Dialog)

                2. MAINTAINS THE PERMISSIONS ON THE FILES AND FOLDERS TO

					PERMIT ACCESS BY THE USER AND THE APACHE USER

                3. PROVIDES THE FOLLOWING FUNCTIONALITY:

                    - PROJECT FOLDERS

                        - ADD
                        - RENAME
                        - DELETE 

                    - WORKFLOWS

                        - ADD
                        - RENAME
                        - DELETE
                        - MOVE TO PROJECT
                        - COPY TO PROJECT

                    - FILES

                        - ADD
                        - RENAME
                        - DELETE
                        - MOVE TO WORKFLOW OR LOWER DIRECTORY
                        - COPY TO WORKFLOW OR LOWER DIRECTORY



				NB: Agua::Project DOES NOT DO THE FOLLOWING:

					- KEEP TRACK OF FILESYSTEM SIZES IN EACH PROJECT

					- LIMIT FILE ADDITIONS BEYOND A PRESET USER QUOTA


	NOTES

		EXAMPLE OF dojox.data.FileStore JSON OBJECT

{
    "total":4,
    "items":
    [
        {
            "name":"dijit",
            "parentDir":".",
            "path":".\/dijit",
            "directory":true,
            "size":0,
            "modified" :1227075503,
            "children":
            [
                "ColorPalette.js",
                "Declaration.js",
                "Dialog.js",
                "dijit-all.js",
                "dijit.js","Editor.js","form","InlineEditBox.js","layout","LICENSE","Menu.081119-added-currentTarget.js","Menu.js","nls","ProgressBar.js","resources","robot.js","robotx.js","templates","tests","themes","TitlePane.js","Toolbar.js","Tooltip.js","Tree.js","_base","_base.js","_Calendar.js","_Container.js","_editor","_Templated.js","_TimePicker.js","_tree","_Widget.js"
            ]
        },
        {
            "name":"dojo",
            "parentDir":".",
            "path":".\/dojo",
            "directory":true,
            "size":0,
            "modified":1226987387,
            "children":
            [
                "AdapterRegistry.js","back.js","behavior.js","cldr","colors.js","cookie.js","currency.js","data","date","date.js","DeferredList.js","dnd","dojo.js","dojo.js.uncompressed.js","fx","fx.js","gears.js","html.js","i18n.js","io","jaxer.js","LICENSE","nls","NodeList-fx.js","NodeList-html.js","number.js","OpenAjax.js","parser.js","regexp.js","resources","robot.js","robotx.js","rpc","string.js","tests","tests.js","_base","_base.js","_firebug"
            ]
        },
        {
            "name":"dojox",
            "parentDir":".",
            "path":".\/dojox",
            "directory":true,
            "size":0,
            "modified":1228105583,
            "children":
            [
                "analytics","analytics.js","av","charting","collections","collections.js","color","color.js","cometd","cometd.js","data","date","dtl","dtl.js","editor","embed","encoding","flash","flash.js","form","fx","fx.js","gfx","gfx.js","gfx3d","gfx3d.js","grid","help","highlight","highlight.js","html","html.js","image","io","json","jsonPath","jsonPath.js","lang","layout","LICENSE","math","math.js","off","off.js","presentation","presentation.js","resources","robot","rpc","secure","sketch","sketch.js","sql","sql.js","storage","storage.js","string","testing","timing","timing.js","uuid","uuid.js","validate","validate.js","widget","wire","wire.js","xml","xmpp"
            ]
        }
        ,
        {
            "name":"util",
            "parentDir":".",
            "path":".\/util",
            "directory":true,
            "size":0,
            "modified":1226987022,
            "children":
            [
                "buildscripts", "docscripts","doh","jsdoc","LICENSE","maven","resources","shrinksafe"
            ]
        }
    ]
}

=cut


#### USE LIB FOR INHERITANCE
use FindBin qw($Bin);
use lib "$Bin";
use lib "$Bin/../";
use lib "$Bin/../..";
use lib "$Bin/../external";

use strict;
use warnings;

class Agua::Project with (Agua::Cluster::Checker,
	Agua::Common::Base,
	Agua::Common::Project,
	Agua::Common::Privileges,
	Agua::Common::Transport,
	Agua::Common::Util)
{


#### EXTERNAL MODULES
use File::Basename;
use File::Copy;    
use File::Path;		
use File::stat;     
use Data::Dumper;   
use Carp;           
use JSON -support_by_pp;

#### INTERNAL MODULES
use Agua::Common::Util;
use Agua::DBaseFactory;

# STRINGS
has 'json'		=> ( isa => 'HashRef|Undef', is => 'rw', default => undef );
has 'fileroot'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'username'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'project'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'workflow'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'validated'	=> ( isa => 'Bool|Undef', is => 'rw', default => undef );
has 'bytes'		=> ( isa => 'Int', is => 'rw', default => 200 );
# OBJECTS
has 'dbobject'	=> ( isa => 'Agua::DBase::MySQL', is => 'rw', required => 0 );
has 'conf' 	=> (
	is =>	'rw',
	'isa' => 'Conf::Agua',
	default	=>	sub { Conf::Agua->new(	backup	=>	1, separator => "\t"	);	}
);

####////}}

method BUILD ($hash) {

	#### SET DATABASE HANDLE
	$self->setDbh();		

    #### VALIDATE
    print "{ error: 'Agua::Project::BUILD    User session not validated' }" and exit unless $self->validate();
}

method newFolder {


    my $folderpath = $self->json()->{folderpath};

    #### VALIDATE
    print "{ error: 'User session not validated' }" and exit unless $self->validate();

    #### GET FULL PATH TO USER'S HOME DIRECTORY
    my $fileroot = $self->getFileroot();

   #### SET FULL PATH TO FILE TO BE COPIED
    $folderpath = "$fileroot/$folderpath";

    #### CONVERT folderpath AND destinationPath ON WINDOWS
    if ( $^O =~ /^MSWin32$/ )   {   $folderpath =~ s/\//\\/g;  }

    #### CHECK FILE IS PRESENT
    #### IF NOT, PRINT ERROR AND EXIT
	print "{ error: 'A folder already exists in path: $folderpath'    }" and exit if -d $folderpath;
	print "{ error: 'A file already exists in path: $folderpath'    }" and exit if -f $folderpath;

	#### CREATE NEW FOLDER
	File::Path::mkpath($folderpath);
	print "{ error: 'Could not create folder in path: $folderpath'    }" and exit if not -d $folderpath;

    #### EXIT ON COMPLETION
	print "{ status: 'Created folder in path: $folderpath' }" and exit;    
}

#### COPY FOLDER
method copyFolder {
	return $self->copyFile();
}

#### REMOVE A FOLDER
method deleteFolder {
	return $self->removeFile();
}

method moveFolder {
	return $self->renameFile();
}

#### RENAME A FOLDER
method renameFolder {
	return $self->renameFile();	
}

#### LATER: KEEP TRACK OF FILESYSTEM SIZES IN EACH PROJECT
method discUsage {

}

#### LATER: LIMIT FILE ADDITIONS BEYOND A PRESET USER QUOTA (DEFAULT = NONE)
method checkQuota {

}


=head2

    SUBROUTINE     renameWorkflow

    PURPOSE

        1. Return a list of getFolders for the user

        2. Each project has permissions

=cut

method renameWorkflow {


    #### VALIDATE
    my $validated = $self->validate();
    if ( not $validated )
    {
        print "{ error: 'User session not validated'    }";
        exit;
    }

    #### GET DATABASE OBJECT
    my $dbobject = $self->dbobject();

    #### GET USER NAME, SESSION ID AND PATH FROM CGI OBJECT
    my $username = $self->json()->{username};
    my $oldpath = $self->json()->{oldPath};
    my $newpath = $self->json()->{newPath};
	my ($new_workflow) = $newpath =~ /([^\/]+)$/;
	my ($old_workflow) = $oldpath =~ /([^\/]+)$/;
    my ($project) = $oldpath =~ /^([^\/]+)/;

    my $query = qq{UPDATE project
    SET workflow = '$new_workflow'
	WHERE username ='$username'
    AND project = '$project'
    AND workflow = '$old_workflow'};
    my $success = $dbobject->do($query);	
	if ( not $success )
	{
		print "{ error: 'Project::renameWorkflow. Could not update project workflow name from $old_workflow to $new_workflow' }";
		exit;
	}
	else
	{
		print "{ status: 'Project::renameWorkflow. Successfully updated project workflow name from $old_workflow to $new_workflow'}";
		exit;
	}
}


=head2

    SUBROUTINE     renameFile

    PURPOSE

        1. RENAME A FILE OR DIRECTORY IN THE USERS' HOME DIRECTORY

		2. A DIRECTORY MAY BE A WORKFLOW, IN WHICH CASE renameWorkflow

			WILL BE CALLED AFTER THIS CALL TO renameFile

cd c:\data\base\cgi-bin\agua

perl project.cgi "mode=renameFile&sessionId=1228791394.7868.158&username=syoung&oldpath=Project1/Workflow1&newpath=Project1/WorkflowA"

=cut

method renameFile {


    #### VALIDATE
    my $validated = $self->validate();
    if ( not $validated )
    {
        print "{ error: 'User session not validated'    }";
        exit;
    }

    #### GET DATABASE OBJECT
    my $dbobject = $self->dbobject();

    #### GET USER NAME, SESSION ID AND PATH FROM CGI OBJECT
    my $username = $self->json()->{username};
    my $oldpath = $self->json()->{oldPath};
    my $newpath = $self->json()->{newPath};

	#### GET FILEROOT 
	my $fileroot = $self->getFileroot($username);

	$oldpath = "$fileroot/$oldpath";
	$newpath = "$fileroot/$newpath";
	if ( $^O =~ /^MSWin32$/ )   {   $newpath =~ s/\//\\/g;  }
	if ( $^O =~ /^MSWin32$/ )   {   $oldpath =~ s/\//\\/g;  }

    #### BEGIN RENAME OLD TO NEW
	use File::Copy;
	my $success = File::Copy::move($oldpath,$newpath);
	if ( not $success )
	{
		print "{ error: 'Project::renameFile. Could not rename from $oldpath to $newpath' }";
		exit;
	}
	else
	{
		print "{ status: 'Project::renameFile. Successfully renamed from $oldpath to $newpath' }";
		exit;
	}
}


=head2

    SUBROUTINE     getFolders

    PURPOSE

        1. Return a list of getFolders for the user

        2. Each project has permissions

cd c:\data\base\cgi-bin\agua

perl project.cgi "mode=getFolders&sessionId=1228791394.7868.158&username=syoung"

=cut

method getFolders {


#### DEBUG -- TIMING
#use Time::HiRes qw[gettimeofday tv_interval];
#use Data::Dumper;
#my $time1=[gettimeofday()];


    #### VALIDATE
    print "{ error: 'User session not validated' }" and exit unless $self->validate();

    #### GET DATABASE OBJECT
    my $dbobject = $self->dbobject();

    #### GET USER NAME, SESSION ID AND PATH FROM CGI OBJECT
    my $username = $self->json()->{username};

	#### POPULATE THE FOLDERS HASH	
	####		folders = {
	####			projects => [],
	####			sharedprojects => [],
	####			worldprojects => [],
	####			sources => []
	####		}
	my $folders = {};

	#### USER'S PROJECTS
	my $query = qq{SELECT DISTINCT project, description FROM workflow WHERE username='$username'};
	my $my_projects = $dbobject->queryhasharray($query);

	#### IF PROJECTS NOT DEFINED, CREATE DEFAULT PROJECT 1, WORKFLOW 1
	if ( not defined $my_projects or not @$my_projects )
	{
		my $project = "Project1";
		my $workflow = "Workflow1";

		$my_projects = [
			{
				'project' => $project,
				'description' => ''
			}
		];

		my $fileroot = $self->getFileroot($username);

		my $projectdir = "$fileroot/$project";
		my $workflowdir = "$fileroot/$project/$workflow";
		$projectdir =~ s/\//\\/g if $^O =~ /^MSWin32$/;
		$workflowdir =~ s/\//\\/g if $^O =~ /^MSWin32$/;

		File::Path::mkpath($projectdir) or die "Can't create project dir: $projectdir\n" if not -d $projectdir;

		File::Path::mkpath($workflowdir) or die "Can't create workflow dir: $workflowdir\n" if not -d $workflowdir;

		#### ADD TO groupmember TABLE
		my $query = qq{INSERT INTO groupmember
(owner, groupname, groupdesc, name, type, description, location)
VALUES ('$username', 'First Group', '$project', 'project', '', '$projectdir')};
	}

    ####
    #
    #   GENERATE LIST OF USER-ACCESSIBLE SHARED PROJECTS
	#
    #    1. GET OWNER AND GROUP NAME OF ANY GROUPS THE USER BELONGS TO
    #
    #    2. GET NAMES OF PROJECTS IN THE GROUPS THE USER BELONGS TO
    #    
    #    3. GET THE PERMISSIONS AND LOCATIONS OF SHARED PROJECTS
    #
    #    .schema access
    #    CREATE TABLE access
    #    (
    #            project                 VARCHAR(20),
    #            owner                   VARCHAR(20),
    #            ownerrights             INT(1),
    #            grouprights             INT(1),
    #            worldrights             INT(1),
    #            location                TEXT,
    #    
    #            PRIMARY KEY (project, owner, ownerrights, grouprights, worldrights, loca
    #    tion)
    #    );
    #

	my $shared_projects = $self->sharedProjects($username);

    ####
    #
    #   GENERATE LIST OF WORLD READABLE PROJECTS
	#
	#### 

	my $world_projects = $self->worldProjects($username, $shared_projects);


	my $sources = $self->sources($username);

	#### DEFINE IF NOT DEFINED
	$my_projects = [] if not defined $my_projects;
	$shared_projects = [] if not defined $shared_projects;
	$world_projects = [] if not defined $world_projects;
	$sources = [] if not defined $sources;	

	$folders->{projects} = $my_projects;
	$folders->{sharedprojects} = $shared_projects;
	$folders->{worldprojects} = $world_projects;
	$folders->{sources} = $sources;

	#### PRINT JSON
    my $jsonObject = JSON->new();
    my $json = $jsonObject->objToJson($folders, {pretty => 1, indent => 4});
    print $json;

	#my $milliseconds = tv_interval($time1)*1000;

	exit;	

}


=head2

    SUBROUTINE     removeFile

    PURPOSE

        Remove a file system with the following functions/constraints:

        1. Return "removing" if remove is already underway

        2. Return "completed" if remove is done

        3. Return "file system busy" if File is being copied from

    NOTES

        Remove constraints and progress are implemented using the first of these

        two methodologies: 

            A) Presence/absence of a flag file ".removeFile.txt"

                in the top directory of the file system to be removed. File

                is removed on completion of remove action.

            B) Remove progress stored in SQL table and updated by perl executable

                that carried out the remove command once the remove is completed


            http://search.cpan.org/~adamk/File-Remove-1.42/lib/File/Remove.pm   
            use File::Remove 'remove';

            # removes (without recursion) several files
            remove( '*.c', '*.pl' );

            # removes (with recursion) several directories
            remove( \1, qw{directory1 directory2} ); 

            # removes (with recursion) several files and directories
            remove( \1, qw{file1 file2 directory1 *~} );

            # trashes (with support for undeleting later) several files
            trash( '*~' );

=cut

method removeFile {

    my $file = $self->json()->{file};
    my $modifier = $self->json()->{modifier};
    if ( not defined $modifier )    {   $modifier = ''; }

    #### VALIDATE
    my $validated = $self->validate();
    if ( not $validated )
    {
        print "{ error: 'User session not validated'    }";
        exit;
    }

	#### GET USER NAME
    my $username = $self->json()->{username};

    #### GET FULL PATH TO USER'S HOME DIRECTORY
    my $fileroot = $self->getFileroot($username);

    #### SET FULL PATH TO FILE TO BE COPIED
    my $filePath = "$fileroot/$file";

    #### CONVERT filePath AND destinationPath ON WINDOWS
    if ( $^O =~ /^MSWin32$/ )   {   $filePath =~ s/\//\\/g;  }

    #### CHECK FILE IS PRESENT
    #### IF NOT, PRINT ERROR AND EXIT
    if ( not $modifier =~ /^status$/ )
    {
        if ( not -f $filePath and not -d $filePath )
        {
            print "{ error: 'File not found'    }";
            exit;
        }
    }

    #### CHECK FOR COPY FLAG FILE 
    my $copyFlagFile = "$filePath~copyFile";
    if ( -f $copyFlagFile )
    {
        print "{    error: 'File being copied'  }\n";
        exit;
    }

    #### SET REMOVE FLAG FILE 
    my $flagFile = "$filePath~removeFile";

    #### RETURN 'removing' IF FLAG FILE FOUND
    if ( -f $flagFile )
    {
        #open(FILE, $flagFile) or die "Can't open flagFile: $flagFile\n";

        if ( $modifier =~ /^status$/ )
        {        
            print "{ status: 'initiated'  }\n";
            exit;
        }
        else
        {
            print "{ error: 'removing'    }\n";
        }
        exit;
    }
    else
    {
        if ( $modifier =~ /^status$/ )
        {
            print "{ status: 'completed'  }\n";
            exit;
        }
    }

    #### CREATE FLAG FILE
    open(FLAG, ">$flagFile") or die "Can't open copyFile flag file\n";
    print FLAG time();
    close(FLAG);

    #### PRINT COPY STATUS: INITIATED
    print "{ status: 'initiated'  }";

    #### FLUSH BUFFER TO MAKE SURE STATUS MESSAGE IS SENT
    $| = 1;


	use File::Remove;

    #### BEGIN REMOVING FILE
    if ( -f $filePath )
    {
        my $result = File::Remove::rm($filePath);
    }
    elsif ( -d $filePath )
    {
        my $result = File::Remove::rm(\1, $filePath);
    }

    #### REMOVE FLAG FILE
    File::Remove::rm($flagFile);

    #### EXIT ON COMPLETION
    exit;    
}




=head2

    SUBROUTINE     copyFile

    PURPOSE

        Copy a file system to different file system with

        the following functions/constraints:

        1. Return "copying" if copy is already underway and

            modifier = 'status' (i.e., its just a status check)

        2. Return "copying" if attempting to copy a

            file that has not been completely copied

    NOTES

        Copy constraints and progress are implemented using the first of these

        two methodologies: 

            A) Presence/absence of a flag file ".projectCopying.txt"

                in the Project directory

            B) Copy progress stored in SQL table and updated by perl executable

                that carried out the copy command once the copy is completed

=cut

method copyFile {



    #### GET USERNAME, FILE, DESTINATION AND MODIFIER
    my $username = $self->json()->{username};
    my $owner = $self->json()->{owner};
    my $file = $self->json()->{file};
    my $destination = $self->json()->{destination};
    my $modifier = $self->json()->{modifier};
    if ( not defined $modifier )    {   $modifier = ''; }

    #### VALIDATE USER SESSION
    my $validated = $self->validate();
    if ( not $validated )
    {
        print "{ error: 'User session not validated'    }";
        exit;
    }

    #### SET FULL PATH TO FILE TO BE COPIED
    my $filePath = '';

    #### SET FULL PATH TO DESTINATION DIRECTORY
    my $destinationPath = '';

    #### CHECK THAT USER HAS RIGHTS TO READ OTHER FILE SYSTEM
    if ( defined $owner and $owner ne '' )
    {

        $file =~ /^([^\/]+)\/(.*)$/;
        my $project = $1;
        my $rest = $2;

        my $can_copy = $self->projectPrivilege($owner, $project, $username, "groupcopy");
		print "{ error: 'User $username does not have sufficient privileges to copy this file' }" and exit if not $can_copy;

		#### GET FILE ROOT
		my $fileroot = $self->getFileroot($owner);

        #### SET FULL PATH TO FILE TO BE COPIED
        $filePath = "$fileroot/$project";
        if ( defined $rest )
        {
            $filePath .= "/$rest";
        }
    }
    else
    {
		#### GET FILE ROOT
		my $fileroot = $self->getFileroot($username);

        $filePath = "$fileroot/$file";
    }

	my $fileroot = $self->getFileroot($username);
	$destinationPath = "$fileroot/$destination";

    #### CHECK THAT USER HAS RIGHTS TO WRITE TO OTHER FILE SYSTEM
    if ( $destination =~ /^owner:/ )
    {

        $destination =~ /^owner:([^\/]+)\/([^\/]+)(.*)$/;
        my $owner = $1;
        my $project = $2;
        my $rest = $3;

        my $privilege = $self->projectPrivilege($owner, $project, $username, "groupwrite");
        if ( not defined $privilege or not $privilege)
        {
            print "{ error: 'User $username does not have sufficient privileges to write to this file system' }\n";
            exit;
        }

        #### SET FULL PATH TO FILE TO BE COPIED
        $destinationPath = "$fileroot/$project";
        if ( defined $rest )
        {
            $destinationPath .= "$rest";
        }
    }
    else
    {
        $destinationPath = "$fileroot/$destination";
    }

    #### CONVERT filePath AND destinationPath ON WINDOWS
    if ( $^O =~ /^MSWin32$/ )   {   $filePath =~ s/\//\\/g;  }
    if ( $^O =~ /^MSWin32$/ )   {   $destinationPath =~ s/\//\\/g;  }



    #### CHECK FILE IS PRESENT
    #### IF NOT, PRINT ERROR AND EXIT
    if ( not -f $filePath and not -d $filePath )
    {
        print "{ error: 'File path not found: $filePath'    }";
        exit;
    }

    #### CHECK IF DESTINATION DIRECTORY IS PRESENT
    #### IF NOT, PRINT ERROR AND EXIT
    if ( not -d $destinationPath )
    {
        print "{ error: 'Destination directory not found: $destination'    }";
        exit;
    }

    #### SET FLAG FILE 
    my $flagFile = "$filePath~copyFile";

    #### RETURN 'copying' IF FLAG FILE FOUND
    if ( -f $flagFile )
    {
        #open(FILE, $flagFile) or die "Can't open flagFile: $flagFile\n";

        if ( $modifier =~ /^status$/ )
        {        
            print "{ status: 'copying'  }";
            exit;
        }

        print "{ error: 'File copy under way'    }";
        exit;
    }

    if ( $modifier =~ /^status$/ )
    {
        print "{ status: 'completed'  }";
        exit;
    }

    #### SET DESTINATION FILE
    my ($filename) = $file =~ /([^\/]+)$/; 
    my $destinationFile = "$destinationPath/$filename";

    #### CONVERT filePath AND destinationPath ON WINDOWS
    #if ( $^O =~ /^MSWin32$/ )   {   $filePath =~ s/\//\\/g;  }
    if ( $^O =~ /^MSWin32$/ )   {   $destinationFile =~ s/\//\\/g;  }

    #### CHECK IF DESTINATION FILE EXISTS
    if ( -f $destinationFile or -d $destinationFile )
    {
        if ( not $modifier =~ /^overwrite$/ )
        {
            print "{ error: 'File exists'  }";
            exit;
        }
    }

    #### CREATE FLAG FILE
    open(FLAG, ">$flagFile") or die "Can't open copyFile flag file\n";
    print FLAG time();
    close(FLAG);

    #### PRINT COPY STATUS: INITIATED
    print "{ status: 'initiated'  }";

    #### FLUSH BUFFER TO MAKE SURE STATUS MESSAGE IS SENT
    $| = 1;

    #### CONVERT FILE PATH ON WINDOWS
    if ( $^O =~ /^MSWin32$/ )   {   $filePath =~ s/\//\\/g;  }
    if ( $^O =~ /^MSWin32$/ )   {   $destinationFile =~ s/\//\\/g;  }



    #### BEGIN COPY FILE TO DESTINATION

	use File::Copy::Recursive;


    #my $result = File::Copy::cp($filePath, $destinationFile);
    my $result = File::Copy::Recursive::rcopy($filePath, $destinationFile);


	use File::Remove;


    #### REMOVE FLAG FILE
    File::Remove::rm($flagFile);

    #### EXIT ON COMPLETION
    exit;    
}



=head2

    SUBROUTINE:     fileStats

    PURPOSE:        Return the following file statistics:
                    -   size (bytes)
                    -   directory ("true"|"false")
                    -   modified (seconds Unix time)

=cut

method fileStats ($filepath) {


    my $fileStats;

    my $filesize = -s $filepath;
    if ( not defined $filesize )
    {
        $filesize = '';
    }

    my $directory = -d $filepath;
    if ( not -f $filepath and not -d $filepath )
    {
        $directory = "''";
    }
    elsif ( not $directory )
    {
        $directory = "false";
    }
    else
    {
        $directory = "true";
    }

    my $modified = -M $filepath;
    if ( not defined $modified )
    {
        $modified = "''"
    }
    else
    {
        $modified = int(time() - $modified * 24 * 60);
    }

    $fileStats->{filesize} = "$filesize";
    $fileStats->{directory} = "$directory";
    $fileStats->{modified} = "$modified";

    return $fileStats;
}


=head2

    SUBROUTINE:     fileSystem

    PURPOSE:        Return the JSON details of all files and directories

                    in the file system.

                    JSON FORMAT:

                    "total":4,
                    "items":
                    [
                       {
                           "name":"dijit",
                           "parentDir":".",
                           "path":".\/dijit",
                           "directory":true,
                           "size":0,
                           "modified" :1227075503,
                           "":
                           [
                               "ColorPalette.js",


./project.cgi "mode=fileSystem&username=syoung&sessionId=9999999999.9999.999&location=/mihg/data/NGS/syoung/base/pipeline/human-chr/fa&query=%22owner%3Aadmin%2Fundefined%22"

perl project.cgi "mode=fileSystem&username=syoung&sessionId=1228791394.7868.158&query=%22Project2%22"


=cut

method fileSystem {


    #### GET USER NAME
    my $username = $self->json()->{username};

    #### FOR SECURITY, REMOVE ALL '..', '\', '/', ETC. FROM USERNAME
    $username =~ s/[><\\\/\.\.\(\)]//g;

    #### VALIDATE
    my $validated = $self->validate();
    if ( not $validated )
    {
        print "{ error: 'User session not validated'    }";
        exit;
    }

	#### GET QUERY, I.E., NAME OF PROJECT
    my $query = $self->json()->{query};
    $query =~ s/"//g if defined $query;    
    $query =~ s/\%22//g if defined $query;    
    if ( not defined $query or $query =~ /^{}$/ )
    {
        $query = '';
    }


	#### GET ADDITIONAL PATH IF DEFINED
    my $path = $self->json()->{path};


	#### CHECK IF WE ARE REQUESTING ACCESS TO SOMEONE ELSE'S FILES
	my $owner = $self->json()->{owner};

	#### EXIT IF NOT ALLOWED ACCESS TO FILES OWNED BY 'OWNER'
	my $fileroot;
	if ( defined $owner )
	{	
		$fileroot = $self->getFileroot($owner);

		#### GET THE PROJECT NAME, I.E., THE BASE DIRECTORY OF THE FILE NAME
		my ($project) = ($query) =~ /^([^\/^\\]+)/;
		my $can_access = $self->canAccessProject($owner, $project, $username);

		if ( not $can_access )
		{
			print "{ error: 'user $username is not permitted access to file path $path' }";
			exit;
		}
	}
	else
	{
		#### GET FILE ROOT FOR THIS USER
		$fileroot = $self->getFileroot();
	}


	#### SET USER NAME TO OWNER IF DEFINED
	$username = $owner if defined $owner and $owner;


    #### SET USERDIR AND CHANGE TO BACKSLASH FOR WINDOWS
    my $fullPath = '';


	#### IF ITS A SOURCE, LOCATION IS DEFINED
	my $location = $self->json()->{location};
	#$location = "/mihg/data/solexa";


    if ( defined $location and $location )
    {
	    $fullPath = "$location";

        #### ADD query IF DEFINED
        if ( defined $path and $path )
        {
            $fullPath .= "/$path";
        }
        elsif ( defined $query and $query and $query !~ /^\// )
        {
            $fullPath .= "/$query";
        }
        elsif ( defined $query and $query and $query =~ /^\// )
        {
            $fullPath = $query;
        }

    }
    else
    {
        $fullPath = "$fileroot";

        #### ADD query IF DEFINED
        if ( defined $path and $path )
        {
            $fullPath .= "/$path";
        }
        elsif ( defined $query and $query )
        {
            $fullPath .= "/$query";
        }
    }

    if ( $^O =~ /^MSWin32$/ )   {   $fullPath =~ s/\//\\/g;  }


    if ( not defined $path or not $path )
    {
        $path = $query;
    }

    if ( not defined $path and not defined $query )
    {
        #die "Neither query nor path are defined";
        #exit;
        $path = "";
    }


    #### GET THE FILE/DIR NAME FROM THE FILEPATH
    my $name = '';
    my $parentPath = '';
    if ( defined $path and $path )
    {
        ($name) = $path =~ /([^\/\\]+)$/;        
    }

    #### IF ITS THE BASE DIRECTORY, I.E., query IS '{}' (path IS NOT DEFINED),
    #### DO { 'total': 3, 'items': [ ... ]} WHERE '...' IS AN
    #### ARRAY OF JSONS FOR EACH FILE/DIR IN THE BASE DIRECTORY
    #### (THE SUB-DIRECTORIES JSON INCLUDES THEIR children)
    my $json = '';
	if ( -d $fullPath )
    {

        opendir(DIRHANDLE, $fullPath) or die "Can't open base directory: $fullPath: $!";
        my @filenames = sort readdir(DIRHANDLE);
        close(DIRHANDLE);

        #### REMOVE '.' AND '..'
        shift @filenames;
        shift @filenames;

		#### REMOVE .SVN DIRECTORY
		for ( my $i = 0; $i < $#filenames + 1; $i++ )
		{

			#### REMOVE ALL .DOT DIRECTORIES
			if ( $filenames[$i] =~ /^\./ )
			{
				splice @filenames, $i, 1;
				$i--;
			}
		}

        my $total = $#filenames + 1;

        #### START JSON
        $json = "/* {\n";

        $json .= "'name': '$name',\n";
        $json .= "'path': '$path',\n";

        #### PRINT THE TOTAL ITEMS
        $json .= "'total': '$total',\n";

        #### START THE items SQUARE BRACKETS
        $json .= "'items': [\n";

        foreach my $file (@filenames)
        {
            #### SET FILEPATH AND CHANGE TO BACKSLASH FOR WINDOWS
            my $filepath = "$fullPath/$file";
            if ( $^O =~ /^MSWin32$/ ) { $filepath =~ s/\//\\/g; }


            $json .= $self->fileJson($path, $filepath, $path, undef);
            $json .= ",\n";
        }
        $json =~ s/,$/\n/;

        $json .= "]\n";
        $json .= "} */\n";

		$json =~ s/\s+//g;

        print $json;
        exit;
    }
    elsif ( -f $fullPath )
    {
        $json = "/* " . $self->fileJson($path, $fullPath, $path, undef) . " */";
    }
	else
	{
        $json = "/* " . $self->fileJson($path, $fullPath, $path, undef) . " */";
	}

	$json =~ s/\s+//g;
	print $json;
	exit;
}


=head2

    SUBROUTINE:     fileJson

    PURPOSE:        Get the details for each file and its children if its a directory

=cut

method fileJson ($path, $filepath, $parentDir, $extra) {

#exit;

    #### START JSON    
    my $json = "{\n";

    #### GET FILE STATS (size, directory, modified)
    my $fileStats = $self->fileStats($filepath);

    #### GET THE FILE/DIR NAME FROM THE FILEPATH
    my ($name) = $filepath =~ /([^\/\\]+)$/;

    $json .= "    'name': '$name',\n";
    $json .= "    'path': '$name',\n";
    $json .= "    'parentPath': '$parentDir',\n";
    $json .= "    'parentDir': '$parentDir',\n";
    $json .= "    'directory': $fileStats->{directory},\n";
    $json .= "    'size': '$fileStats->{filesize}',\n";
    $json .= "    'modified': $fileStats->{modified},\n";

    #### DO CHILDREN 
    if ( $fileStats->{directory} =~ /^true$/ )
    {
        $json .= "    'children': [\n";

        opendir(DIRHANDLE, $filepath) or die "Can't open filepath: $filepath: $!";
        my @filenames = sort readdir(DIRHANDLE);
        close(DIRHANDLE);

        #### REMOVE '.' AND '..'
        shift @filenames;
        shift @filenames;

        foreach my $file ( @filenames )
        {
            $json .= "        '$file',\n";    
        }
        $json =~ s/,\n$/\n/;
        $json .= "    ]\n";
    }
    else
    {
		#$self->bytes(200);

        #### GET A SAMPLE FROM THE TOP OF THE FILE
        my $sample;
		my $bytes = $self->bytes();
        if ( -B $filepath )
        {
            $sample = "Binary file";
        }
        elsif ( -f $filepath and not -z $filepath )
        {
            open(FILEHANDLE, $filepath);
            seek(FILEHANDLE, 0, 0);
            read(FILEHANDLE, $sample, $bytes);
        }

        #### SET TO '' IF FILE IS EMPTY
        if ( not defined $sample or not $sample )
        {
            $sample = '';
        }
		else
		{
			$sample = $self->jsonSafe($sample, 'toJson');
		}
        $json .= qq{    'sample': '$sample', };
        $json .= qq{    'bytes' : '$bytes' };
    }

    #### END JSON
    $json .= "}";

    return $json;

}



}

