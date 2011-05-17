use MooseX::Declare;

=head2

	PACKAGE		View

	PURPOSE

		THIS MODULE ENABLES THE FOLLOWING USE SCENARIOS:

		1. USER GENERATES JBROWSE JSON FILES IN OWN FOLDER

			1.1 CREATE FEATURES FOLDER

				1.1.1 CREATE WORKFLOW-SPECIFIC FEATURES FOLDER

			   users/__username__/__project__/__workflow__/features

			   TO HOLD INDIVIDUAL FEATURE JSON FILE DIRECTORIES INPUT FILES, E.G.:

			   users/__username__/__project__/__workflow__/features/__feature_name__

			   ... GENERATED FROM INDIVIDUAL FEATURE FILES, E.G.:

			   users/__username__/__project__/__workflow__/somedir/__feature_name__.gff


			1.2 SET LOCATION OF refSeqs.js FILE IF NOT DEFINED (E.G., WOULD BE

				SPECIFIED BY THE USER IN CASE OF CUSTOM GENOME)


			1.3 GENERATE JBrowse FEATURES IN PARALLEL ON A CLUSTER

				FOR EACH FEATURE:

					1.2.1 GENERATE A UNIQUE NAME FOR THIS PROJECT

						BY ADDING AN INCREMENTED DIGIT TO THE END OF THE FEATURE IN

						THE CASE OF DUPLICATE FEATURE NAMES IN THE SAME PROJECT OR

						WORKFLOW

					1.2.2 REGISTER FEATURE IN features TABLE

					1.2.3 GENERATE FEATURE IN FEATURES SUBFOLDER (Agua::JBrowse), E.G.:

					   users/__username__/__project__/__workflow__/features/__feature_name__


		2. USER CREATES A NEW VIEW

			2.1 CREATE A NEW USER- AND PROJECT-SPECIFIC FOLDER INSIDE

				THE JBROWSE ROOT FOLDER:

			   plugins/view/jbrowse/users/__username__/__project__/__viewname__

				THE VIEW location WILL BE USED BY View.js TO LOCATE THE refSeqs.js AND

				trackInfo.json FILES FOR LOADING THE TRACKS.


			2.2 LINK (ln -s) ALL STATIC FEATURE TRACK SUBFOLDERS TO THE

				VIEW FOLDER (STATIC FEATURES ARE STORED IN THE feature TABLE)

				NB: LINK AT THE INDIVIDUAL STATIC FEATURE LEVEL, E.G.:

				jbrowse/species/__species__/__build__/data/tracks/chr*/FeatureDir 

				__NOT__ AT A HIGHER LEVEL BECAUSE WE WANT TO BE ABLE

				TO ADD/REMOVE DYNAMIC TRACKS IN THE USER'S VIEW FOLDER WITHOUT

				AFFECTING THE PARENT DIRECTORY OF THE STATIC TRACKS.

				THE STATIC TRACKS ARE MERELY INDIVIDUALLY 'BORROWED' AS IS.


		3. USER ADDS/REMOVES TRACKS TO/FROM VIEW

		    3.1 CREATE USER-SPECIFIC VIEW FOLDER:

		        plugins/view/jbrowse/users/__username__/__project__/__view__/data

			3.2 ADD/REMOVE FEATURE trackData.json INFORMATION TO/FROM data/trackInfo.js

			3.3 RERUN generate-names.pl AFTER EACH ADD/REMOVE

=cut


#### USE LIB FOR INHERITANCE
use FindBin qw($Bin);
use lib "$Bin";
use lib "$Bin/../";
use lib "$Bin/../..";
use lib "$Bin/../external";

use strict;
use warnings;

class Agua::View extends Agua::JBrowse with (Agua::Cluster::Checker,
	Agua::Cluster::Cleanup,
	Agua::Cluster::Jobs,
	Agua::Cluster::Loop,
	Agua::Cluster::Usage,
	Agua::Cluster::Util,
	Agua::Common::Util,
	Agua::Common::Base,
	Agua::Common::View,
	Agua::Common::Privileges)
{


#### EXTERNAL MODULES
use FindBin qw($Bin);
use Data::Dumper;
use File::Path;
use File::Copy::Recursive;
use File::Remove;

#### INTERNAL MODULES
use Agua::Common::Util;
use Agua::DBaseFactory;
use Agua::JBrowse;

# STRINGS
has 'json'		=> ( isa => 'HashRef|Undef', is => 'rw', default => undef );
has 'username'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'project'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'view'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'feature'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'workflow'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'viewdir'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'featuresdir'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'jbrowsedir'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'validated'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
# OBJECTS
has 'dbobject'	=> ( isa => 'Agua::DBase::MySQL', is => 'rw', required => 0 );
has 'conf' 	=> (
	is =>	'rw',
	'isa' => 'Conf::Agua',
	default	=>	sub { Conf::Agua->new(	backup	=>	1, separator => "\t"	);	}
);

	####/////}
method BUILD ($hash) {

	#### SET DATABASE HANDLE
	$self->setDbh();

    #### VALIDATE
    print "{ error: 'Agua::View::BUILD    User session not validated' }" and exit unless $self->validate();

	#### SET DIRECTORIES
	$self->setFeaturesDir();
}

method addView () {


	my $json 		=	$self->json();
	#exit;

	#### SET VARIABLES
	my $username 	= 	$self->username($json->{username});
	my $project 	= 	$self->project($json->{project});
	my $view 		= 	$self->view($json->{view});
	my $species 	= 	$self->species($json->{species});
	my $build 		= 	$self->build($json->{build});
	my $htmlroot	=	$self->htmlroot();

	#### CHECK INPUTS	
	print "Agua::View::addView    htmlroot not defined\n" and exit if not defined $htmlroot;	
	print "Agua::View::addView    project not defined\n" and exit if not defined $project;
	print "Agua::View::addView    view not defined\n" and exit if not defined $view;
	print "Agua::View::addView    username not defined\n" and exit if not defined $username;


	my $viewdir = $self->setViewDir();

	File::Path::mkpath($viewdir) if not -d $viewdir;
	print "{ error: 'Agua::View::addView    Could not create viewdir: $viewdir' }" and exit if not -d $viewdir;

	my $success = $self->_addView();
	print "{ error: 'Agua::View::addView    Could not add view $json->{view} to view table' }" and exit if not defined $success or not $success;
 	print "{ status: 'Agua::View::addView    Added view $json->{view} to view table' }";
	exit;
}

method removeView {
    my $json 			=	$self->json();

	#### SET VARIABLES
	my $username 	= 	$self->username($json->{username});
	my $project 	= 	$self->project($json->{project});
	my $view 		= 	$self->view($json->{view});
	my $species 	= 	$self->species($json->{species});
	my $build 		= 	$self->build($json->{build});
	my $htmlroot	=	$self->htmlroot();

	#### DELETE THE FILE SYSTEM
	my $viewdir = $self->setViewDir();
 	print "Agua::View::removeView    viewdir: $viewdir";

	my $success = File::Remove::rm(\1, $viewdir);
	print "{ error: 'Could not remove viewdir: $viewdir' }"
		and exit if not $success;

	#### DO THE REMOVE
	$success = $self->_removeView();
 	print "{ error: 'Agua::View::removeView    Could not remove view $json->{view} from view table' } " and exit if not defined $success or not $success;
 	print "{ status: 'Agua::View::removeView    Removed view $json->{view} from view table' } ";
	exit;
}

method addViewFeature () {

	my $json 		=	$self->json();

	#### CHECK SOURCE INPUTS
	print "Agua::View::addViewFeature    sourceproject not defined\n"
		and exit if not $json->{sourceproject};
	print "Agua::View::addViewFeature    sourceworkflow not defined\n"
		and exit if not $json->{sourceworkflow};
	print "Agua::View::addViewFeature    feature not defined\n"
		and exit if not $json->{feature};
	print "Agua::View::addViewFeature    species not defined\n"
		and exit if not $json->{species};
	print "Agua::View::addViewFeature    build not defined\n"
		and exit if not $json->{build};

	#### GET LOCATION 
	my $fields = ["username", "project", "workflow", "feature"]; 
	my $query = qq{SELECT location
FROM feature
WHERE project='$json->{sourceproject}'
AND workflow='$json->{sourceworkflow}'
AND feature='$json->{feature}'
AND species='$json->{species}'
AND build='$json->{build}'};
	my $location = $self->dbobject()->query($query);

	#### SET FEATURE OBJECT
	my $feature_object;
	$feature_object->{location} 	= $location;
	$feature_object->{username} 	= 	$self->username($json->{username});
	$feature_object->{project} 		= 	$self->project($json->{project});
	$feature_object->{view} 		= 	$self->view($json->{view});
	$feature_object->{feature} 		= 	$self->feature($json->{feature});
	$feature_object->{species} 		= 	$self->species($json->{species});
	$feature_object->{build} 		= 	$self->build($json->{build});

	#### CHECK HTMLROOT - NEED FOR DEFINING FILEPATHS
	my $htmlroot		=	$self->conf()->getKeyValue("agua", 'HTMLROOT');
	print "Agua::View::addViewFeature    htmlroot not defined\n"
		and exit if not defined $htmlroot;

	#### CHECK OTHER INPUTS
	print "Agua::View::addViewFeature    project not defined\n"
		and exit if not $self->project();
	print "Agua::View::addViewFeature    view not defined\n"
		and exit if not $self->view();
	print "Agua::View::addViewFeature    username not defined\n"
		and exit if not $self->username();
	print "Agua::View::addViewFeature    feature not defined\n"
		and exit if not $self->feature();

	#### LINK ALL FEATURE TRACKS TO TARGET VIEW DIR
	my $viewdir = $self->activateView($self->project(), $self->view(), $self->species(), $self->build());
	print "{ error: 'viewdir not defined: $viewdir' }" and exit if not defined $viewdir;
	print "{ error: 'Could not create viewdir: $viewdir' }" and exit if not -d $viewdir;
	print "Agua::View::addViewFeature    viewdir: $viewdir\n";

return;

	my $chromodirs;
	my $feature = $self->feature();
	my $target_tracksdir = $self->getTargetTracksdir();
	print "Agua::View::addViewFeature    target_tracksdir: $target_tracksdir\n";

	opendir(DIR, $target_tracksdir)
		or die "Can't open target_tracksdir: $target_tracksdir\n";
	@$chromodirs = readdir(DIR);
	closedir(DIR) or die "Can't close target_tracksdir: $target_tracksdir";
	foreach my $chromodir ( @$chromodirs )
	{
		my $chromopath = "$target_tracksdir/$chromodir";
		next if $chromodir =~ /^\./ or not -d $chromopath;
		#File::Path::mkpath($chromopath) if not -d $chromopath;
		print "Agua::View::addViewFeature    Linking feature in chromodir: $chromodir\n";
		$self->addLink("$location/data/tracks/$chromodir/$feature", "$target_tracksdir/$chromodir/$feature");
	}

	#### UPDATE trackInfo.js
	my $source_tracksdir = $self->getSourceTracksdir($self->species(), $self->build());
	print "Agua::View::addViewFeature    source_tracksdir: $source_tracksdir\n";


	my $view_infofile = $self->getViewInfofile($self->species(), $self->build());
	print "Agua::View::addViewFeature    view_infofile: $view_infofile\n";
	my $feature_infofile = $self->getFeatureInfofile($location);
	print "Agua::View::addViewFeature    feature_infofile: $feature_infofile\n";

	$self->addTrackinfo($feature_infofile, $view_infofile);

	my $success = $self->_addViewFeature($feature_object);
 	print "{ error: 'Agua::View::addViewFeature    Could not add feature $json->{feature} to viewfeature table' }" and exit if not defined $success or not $success;
 	print "{ status: 'Agua::View::addViewFeature    Added feature $json->{feature} to viewfeature table' }";

	exit;
}

method removeViewFeature () {

	my $json 		=	$self->json();

	#### CHECK INPUTS
	print "Agua::View::removeViewFeature    project not defined\n"
		and exit if not $json->{project};
	print "Agua::View::removeViewFeature    view not defined\n"
		and exit if not $json->{view};
	print "Agua::View::removeViewFeature    feature not defined\n"
		and exit if not $json->{feature};
	print "Agua::View::removeViewFeature    username not defined\n"
		and exit if not $json->{username};

	#### SET FEATURE OBJECT
	my $feature_object;
	$feature_object->{username} 	= 	$self->username($json->{username});
	$feature_object->{project} 		= 	$self->project($json->{project});
	$feature_object->{view} 		= 	$self->view($json->{view});
	$feature_object->{feature} 		= 	$self->feature($json->{feature});


	#### CHECK HTMLROOT - NEED FOR DEFINING FILEPATHS
	my $htmlroot		=	$self->conf()->getKeyValue("agua", 'HTMLROOT');
	print "Agua::View::removeViewFeature    htmlroot not defined\n"
		and exit if not defined $htmlroot;

	my $chromodirs;
	my $feature = $self->feature();
	my $target_tracksdir = $self->getTargetTracksdir();
	print "Agua::View::removeViewFeature    target_tracksdir: $target_tracksdir\n";

	opendir(DIR, $target_tracksdir)
		or die "Can't open target_tracksdir: $target_tracksdir\n";
	@$chromodirs = readdir(DIR);
	closedir(DIR) or die "Can't close target_tracksdir: $target_tracksdir";
	foreach my $chromodir ( @$chromodirs )
	{
		my $chromopath = "$target_tracksdir/$chromodir";
		next if $chromodir =~ /^\./ or not -d $chromopath;
		#File::Path::mkpath($chromopath) if not -d $chromopath;
		print "Agua::View::removeViewFeature    Removing link for feature in chromodir: $chromodir\n";
		$self->removeLink("$target_tracksdir/$chromodir/$feature");
	}

	#### UPDATE trackInfo.js
	my $view_infofile = $self->getViewInfofile($self->species(), $self->build());
	print "Agua::View::removeViewFeature    view_infofile: $view_infofile\n";

	my $featureinfo;
	$featureinfo->{label} =  $feature;
	$featureinfo->{key} =  $feature;

	$self->removeTrackinfo($featureinfo, $view_infofile);

	my $success = $self->_removeViewFeature($feature_object);
 	print "{ error: 'Agua::View::removeViewFeature    Could not add feature $json->{feature} to viewfeature table' }" and exit if not defined $success or not $success;
 	print "{ status: 'Agua::View::removeViewFeature    Removed feature $json->{feature} from viewfeature table' }";
	exit;
}


method getFeatureInfofile (Str $location) {
	return "$location/data/trackInfo.js";
}

method getViewInfofile (Str $species, Str $build) {
	my $target_tracksdir = $self->getTargetTracksdir();
	return	"$target_tracksdir/../trackInfo.js"; 
}

method getTargetTracksdir () {
	my $viewdir = $self->setViewDir();
	return "$viewdir/data/tracks";
}

method getSourceTracksdir (Str $species, Str $build) {
	my $htmlroot = $self->conf()->getKeyValue("agua", 'HTMLROOT');
	print "Agua::View::getSourceTracksDir    htmlroot not defined\n"
		and exit if not defined $htmlroot;
	my $source_tracksdir 	=	"$htmlroot/plugins/view/jbrowse/species/$species/$build/data/tracks";

	return $source_tracksdir;
}

method setJbrowseDir (Str $species, Str $build) {

	return $self->jbrowsedir() if $self->jbrowsedir();

	print "Agua::View::setJbrowseDir    species not defined or empty\n"
		and return if not defined $species or not $species;
	print "Agua::View::setJbrowseDir    build not defined or empty\n"
		and return if not defined $build or not $build;

	my $jbrowse_data = $self->conf()->getKeyValue(("data", 'JBROWSEDATA'));
	print "Agua::View::setJbrowseDir    jbrowse_data not defined or empty\n"
		and return if not defined $jbrowse_data or not $jbrowse_data;

	return $self->jbrowsedir("$jbrowse_data/$species/$build/jbrowse");
}

method addLink (Str $source, Str $target) {

	#### REMOVE LINK IF EXISTS
	$self->removeLink($target);

	#### ADD LINK
	my $command = "ln -s $source $target";
	`$command`;

	print "Agua::View::addLink    link exists\n" and return 1 if -l $target;
	print "Agua::View::addLink    link DOES NOT exist\n" and return 0;
}

method removeLink (Str $target) {
	my $command = "unlink $target";
	`$command`;

	print "Agua::View::removeLink    link DOES NOT exist\n" and return 1 if not -l $target;
	print "Agua::View::removeLink    link exists\n" and return 0;
}

method activateView ( Str $project, Str $view, Str $species, Str $build) {

	#### GET VIEWDIR (TO BE RETURNED)
	my $viewdir = $self->setViewDir();

	#### CREATE TARGET tracks DIR IF NOT EXISTS
	my $target_tracksdir = $self->getTargetTracksdir();
	if ( -d $target_tracksdir )
	{
		return $viewdir;
	}
	File::Path::mkpath($target_tracksdir) if not -d $target_tracksdir;
	print "{ error: 'Can't create target_tracksdir: $target_tracksdir' }" and exit if not -d $target_tracksdir;

	#### SET SOURCE tracks DIR
	my $source_tracksdir 	=	$self->getSourceTracksdir($species, $build);

	#### LINK STATIC TRACKS TO VIEW 
	$self->linkTracks($source_tracksdir, $target_tracksdir);	

	#### ADD STATIC refSeq.js AND trackInfo.js FILES TO VIEW DIR
	$self->addRefseqsFile($viewdir, $species, $build);
	$self->addTrackinfoFile($viewdir, $species, $build);

	#### ADD names.json FILE TO VIEW DIR
	$self->updateNames($project, $view, $species, $build);	

	#### ADD SPECIES AND BUILD INFO TO VIEW IF NOT PRESENT
	$self->addSpeciesToView($project, $view, $species, $build);	

	return $viewdir;
}

method updateNames (Str $project, Str $view, Str $species, Str $build) {
	my $jbrowse = $self->conf()->getKeyValue(("applications", 'JBROWSE'));
	print "Agua::View::updateNames    jbrowse: $jbrowse\n";
	my $executable = "$jbrowse/generate-names.pl";
	my $dir = $self->setViewDir();
	$dir.= "/data";

	my $command = "$executable --dir $dir";
	print "Agua::View::updateNames    command: $command\n";
	print `$command`;
}

method addRefseqsFile (Str $viewdir, Str $species, Str $build) {
	my $jbrowsedir = $self->setJbrowseDir($species, $build);
	my $refseqfile = "$jbrowsedir/chr1/data/refSeqs.js";
	my $trackinfofile = "$jbrowsedir/chr1/data/trackInfo.js";
	my $copy_refseq = File::Copy::Recursive::rcopy($refseqfile, "$viewdir/data");
	print "{ error: 'Failed to copy refSeqs.js file: $refseqfile' }\n" and exit if not $copy_refseq;	
}

method addTrackinfoFile (Str $viewdir, Str $species, Str $build) {
	my $jbrowsedir = $self->setJbrowseDir($species, $build);
	my $refseqfile = "$jbrowsedir/chr1/data/refSeqs.js";
	my $trackinfofile = "$jbrowsedir/chr1/data/trackInfo.js";
	my $copy_trackinfo = File::Copy::Recursive::rcopy($trackinfofile, "$viewdir/data");
	print "Agua::View::addTrackinfoFile    copy_trackinfo: $copy_trackinfo\n";
	print "{ error: 'Failed to copy trackInfo.js file: $trackinfofile' }" and exit if not $copy_trackinfo;
}

method linkTracks (Str $source_tracksdir, Str $target_tracksdir) {
	#### LINK ALL FEATURE DIRECTORIES IN SOURCE tracks DIR TO TARGET tracks DIR
	opendir(DIR, $source_tracksdir) or die "Can't open source_tracksdir: $source_tracksdir\n";
	my $subdirs;
	@$subdirs = readdir(DIR);
	closedir(DIR) or die "Can't close source_tracksdir: $source_tracksdir\n";
	print "Agua::View::linkTracks    subdirs: @$subdirs\n";

	my $link_tracks = 1;
	foreach my $subdir ( @$subdirs )
	{
		my $subdirpath = "$source_tracksdir/$subdir";
		next if not -d $subdirpath or $subdir =~ /^\./;
		$link_tracks = 0
		 if not $self->addLink("$source_tracksdir/$subdir", "$target_tracksdir/$subdir");
	}
	print "{ error: 'Agua::View::linkTracks    Error linking track dirs' }\n" and exit if not $link_tracks;
}

method addSpeciesToView ( Str $project, Str $view, Str $species, Str $build) {

	my $query = qq{SELECT 1 FROM view
WHERE project='$project'
AND view='$view'
AND species=''
OR build=''};
	my $missing = $self->dbobject()->query($query);
	return if not defined $missing or not $missing;

	$query = qq{UPDATE view
SET species = '$species',
build='$build'
WHERE project='$project'
AND view='$view'};
	my $success = $self->dbobject()->do($query);
	print "Agua::View::addSpeciesToView    success: $success\n";
	print "{ error: 'Failed to add species $species and build $build to $project:$view in view table' }" and return 0 if not $success;
	return 1;
}


=head2

	SUBROUTINE		createFeaturesDir

	PURPOSE

		CREATE USER-SPECIFIC FEATURES FOLDER

		users/__username__/__project__/__workflow__/features

		TO HOLD INDIVIDUAL FEATURE INPUT FILES, E.G.:

		users/__username__/__project__/__workflow__/features/_feature_name.gff

=cut
method createFeaturesDir () {


	$self->setFeaturesDir();

	my $featuresdir = $self->featuresdir();

	File::Path::mkpath($featuresdir);
	print "Agua::View::createFeaturesDir Can't create featuresdir: $featuresdir\n" and exit if not -d $featuresdir;
}	

method setFeaturesDir () {
	my $userdir = $self->conf()->getKeyValue("agua", 'USERDIR');
	print "Agua::View::setFeaturesDir userdir not defined\n" and exit if not defined $userdir;

	my $featuresdir = $userdir;
	$featuresdir .= "/" . $self->username() . "/agua";
	$featuresdir .= "/" . $self->project();
	$featuresdir .= "/" . $self->workflow();
	$featuresdir .= "/jbrowse";

	$self->featuresdir($featuresdir);
}


method setViewDir () {


	return $self->viewdir() if $self->viewdir();

	my $userdir = $self->conf()->getKeyValue("agua", 'USERDIR');
	print "Agua::View::setViewDir userdir not defined\n" and exit if not defined $userdir;

	my $htmlroot = $self->conf()->getKeyValue("agua", 'HTMLROOT');
	print "Agua::View::setViewDir    htmlroot not defined\n" and exit if not defined $htmlroot;

	my $jbrowseroot = "$htmlroot/plugins/view/jbrowse";
	my $viewdir = "$jbrowseroot";
	$viewdir .= "/users";
	$viewdir .= "/" . $self->username();
	$viewdir .= "/" . $self->project();	
	$viewdir .= "/" . $self->view();	

	$self->viewdir($viewdir);

	return $viewdir;
}

around generateFeatures () {
	print "Agua::View::generateFeatures    Agua.View.generateFeatures()\n";
	#### CHECK INPUTS
	print "Agua::View::generateFeatures username not defined: " and exit if not $self->username();
	print "Agua::View::generateFeatures project not defined\n" and exit if not $self->project();
	print "Agua::View::generateFeatures workflow not defined\n" and exit if not $self->workflow();


	##### CREATE DB OBJECT USING DBASE FACTORY
	my $dbtype = $self->conf()->getKeyValue(("database", 'DBTYPE'));
	my $dbobject = 	Agua::DBaseFactory->new( $dbtype,
		{
			'dbfile'	=>	$self->conf()->getKeyValue(("database", 'DBFILE')),
			'database'	=>	$self->conf()->getKeyValue(("database", 'DATABASE')),
			'user'      =>  $self->conf()->getKeyValue(("database", 'USER')),
			'password'  =>  $self->conf()->getKeyValue(("database", 'PASSWORD'))
		}
	) or print qq{ error: 'Cannot create database object $self->conf()->getKeyValue(("database", 'DATABASE')): $!' } and exit;
	print qq{ error: 'dbobject is not defined' } and exit if not defined $dbobject;
	$self->dbobject($dbobject);

	#### 1.1 CREATE USER-SPECIFIC FEATURES FOLDER
	$self->createFeaturesDir();

	#### SET LOCATION OF refSeqs.js FILE
	$self->setRefseqfile() if not $self->refseqfile();

	#### 1.2 GENERATE JBrowse FEATURES IN PARALLEL ON A CLUSTER (Agua::JBrowse)
	####
	$self->$orig();
}


method registerFeature (Str $feature) {




	my $location = $self->featuresdir();
	$location .= "/" . $feature;	

	my $feature_object;
	$feature_object->{username} = $self->username();
	$feature_object->{project} = $self->project();
	$feature_object->{workflow} = $self->workflow();
	$feature_object->{species} = $self->species();
	$feature_object->{build} = $self->build();
	$feature_object->{feature} = $feature;
	$feature_object->{location} = $location;
	$feature_object->{type} = "dynamic";

	my $success = $self->_addFeature($feature_object);
	print "Agua::View::registerFeature    success: $success\n";
};

#around addTrackdata ($trackdatafile, $trackinfofile) {
#
#	my $trackinfo = $self->readJson($trackinfofile);
#
#}	





}  #### END
