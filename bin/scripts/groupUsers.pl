#!/usr/bin/perl -w
use strict;

$DEBUG = 1;

=head2

APPLICATION     groupUsers

PURPOSE:

    PARSE USERS AND THEIR GROUPS FROM THE /etc/groups FILE AND

    OUTPUT THEM IN THE users TABLE FORMAT:




        1. LINK cgi-bin TO /var/www/cgi-bin
        -----------------------------------

        ln -s html /var/www/html/agua
        ln -s cgi-bin /var/www/cgi-bin/agua


        2. LINK lib AND conf TO INSIDE OF cgi-bin
        -----------------------------------------

        ln -s lib /var/www/cgi-bin/agua/lib
        ln -s conf /var/www/cgi-bin/agua/conf


        3. PATHS TO BASIC EXECUTABLES 
        -----------------------------

        CONF FILE:

            BIN                     /nethome/syoung/base/bin
            BINDIR                  /nethome/syoung/base/html/Bioptic0.2.5/bin
            FILEROOT                /nethome/syoung/base/html/Bioptic0.2.5/fileroot
            ROOTDIR		            /nethome/syoung/base/pipeline

            RUNMAPPING              /nethome/syoung/base/apps/454/2.0.00.20-64/bin/runMapping 

            GAPIPELINE_BIN          /mihg/analysis/GAPipeline-1.3.2/bin/
            ELAND                   /mihg/analysis/GAPipeline-1.3.2/bin/ELAND_standalone.pl
            BUSTARD                 /mihg/analysis/GAPipeline-1.3.2/bin/bustard.py
            GOAT_PIPELINE           /mihg/analysis/GAPipeline-1.3.2/bin/goat_pipeline.py
            SQUASH                  /mihg/analysis/GAPipeline-1.3.2/bin/squashGenome

            MIRA    	            /nethome/syoung/base/apps/mira/bin/mira
            VELVET                  /nethome/syoung/base/apps/velvet/velvet
            NUCMER                  /nethome/syoung/base/apps/mummer/nucmer
            DELTAFILTER             /store/nethome/syoung/base/apps/mummer/delta-filter
            SHOWCOORDS              /store/nethome/syoung/base/apps/mummer/show-coords
            VCAKE                   /nethome/syoung/base/apps/vcake/VCAKE_1.0.pl
            HTMLURL		            http://solexa01.med.miami.edu
            HTMLDIR		            /var/www/html
            PHRED_PARAMETER_FILE	/Users/young/FUNNYBASE/apps/phred/phredpar.dat
            MYSQLDATA		        /usr/local/mysql/data
            CLUSTER_MYSQLDATA       /private/var/mysql


INPUT:

    cgi-bin     CGI executables directory - link 'cgi-bin' to this 
    html        Web directory - link 'html' directory to this
    base        Base directory to install bin, lib, html, cgi-bin, etc. directories
    conf        Replace entries to additional executables in default.conf file with user-defined paths


OUTPUT:

    The following directory structure

     <base>/bin
            cgi-bin --> /var/www/cgi-bin/agua
            conf    --> /var/www/cgi-bin/agua/conf
            fileroot
            html    --> /var/www/html/agua
            lib     --> /var/www/cgi-bin/agua/lib
            t

NOTES:

    USES FindBin TO LOCATE THE lib DIRECTORY

    USES Util FOR USER PROMPTS

=cut

#### FLUSH BUFFER
$| = 1;

#### EXTERNAL MODULES
use FindBin qw($Bin);

#### USE LIB
use lib "$Bin/../lib";
use lib "$Bin/../lib/external";

#### INTERNAL MODULES
#use File::Copy;
use File::Copy::Recursive;
use File::Remove 'remove';
use Config::JSON;

#### INSTALL THESE DIRECTORIES
my @directories = qw(bin cgi-bin conf fileroot html lib t);
#my @directories = qw(bin cgi-bin conf db fileroot lib t);

#### URL OF DOJO RELEASE TO BE DOWNLOADED
my $dojo_release = "http://download.dojotoolkit.org/release-1.2.2/dojo-release-1.2.2.tar.gz";

#### DEFAULT DESTINATION INSTALLATION DIRECTORIES
my $basedir = "/usr/share/agua";
my $webdir = "/var/www/html/agua";
my $cgidir = "/var/www/cgi-bin/agua";
my $url = "http://localhost/agua";

#my $basedir = "/usr/share/agua-test";
#my $webdir = "/var/www/html/agua-test";
#my $cgidir = "/var/www/cgi-bin/agua-test";
#my $url = "http://localhost/agua-test";

#### START USER PROMPTS
print "Welcome the Agua groupUsersuration utility (groupUsers.pl)\n";
print "The standard install requires permissions to create\nthe following directories (e.g., root):\n";
print "    INSTALL DIR : $basedir\n";
print "    HTML DIR    : $webdir\n";
print "    CGI DIR     : $cgidir\n";
print "    URL         : $url\n";
print "Type 'Y' for the standard installation or type 'N' for custom installation\n";

($basedir, $webdir, $cgidir, $url) = installation_parameters($basedir, $webdir, $cgidir, $url);

#### REPORT INSTALLATION PARAMETERS
print "\n";
print "Using these parameters for installation:\n";
print "basedir: $basedir\n";
print "webdir: $webdir\n";
print "cgidir: $cgidir\n";
print "url: $url\n";
print "\n";

#### CHECK IF BASE DIRECTORY ALREADY EXISTS
print "Checking if base directory already exists...\n";
if ( -d $basedir )
{
    print "Base directory already exists: $basedir\n";
    print "Do you want to overwrite the existing installation (and archive all the existing data)?\n";
    my $overwrite = 0;
    if ( yes("Type 'Y' for default installation or 'N' for custom installation") )    {   $overwrite  = 1; }
    if ( $overwrite )
    {
        print "Archiving existing installation...\n";

        #### CREATE ARCHIVE DIRECTORY IF NOT EXISTS
        my $archivedir = "$basedir/archive";
        if ( not -d $archivedir  )
        {
            mkdir($archivedir) or die "Can't create archive directory: $archivedir\n";
        }

        #### GET PREVIOUS ARCHIVED INSTALLATIONS INSIDE ARCHIVE
        opendir(DIR, $archivedir) or die "Can't open archive directory: $archivedir\n";
        my @archives = readdir(DIR);
        @archives = sort @archives;
        shift @archives; shift @archives;   #### REMOVE '..' AND '.'
        close(DIR);
        my $number = $#archives + 1;
        my $archive = "$basedir/archive/agua$number";
        print "Moving files to archive: $archive\n";
        mkdir($archive) or die "Can't create archive: $archive\n";

        foreach my $directory ( @directories )
        {
            #### KEEP FILEROOT IN PLACE
            next if $directory =~ /^fileroot$/;

            #### KEEP DATABASES IN PLACE
            next if $directory =~ /^db$/;

            my $fullpath = "$basedir/$directory";

            #### PRINT '.' AS PROGRESS COUNTER
            print ".";

            if ( -d $fullpath )
            {
                File::Copy::Recursive::rmove("$basedir/$directory", "$archive/$directory") or die "Can't move directory $basedir/$directory to $archive/$directory\n";
            }
            else
            {
                print "Can't find directory to move to archive: $fullpath\n";
            }
        }
        print "\ndone\n";
    }
    else
    {
        die "Exiting installation.\n";
    }
}
else
{
    #### MAKE BASE DIRECTORY
    print "Base directory does not exist. Creating base directory...\n";
    mkdir($basedir) or die "Can't create base directory: $basedir $@\n";


    mkdir("$basedir/html") or die "Can't create base directory: $basedir/html, $@\n";

}


#### COPY ALL FILES TO BASE DIR
print "\nCopying folders to base directory: @directories\n";
foreach my $directory ( @directories )
{
    my $fullpath = "$Bin/../$directory";

    #### PRINT '.' AS PROGRESS COUNTER
    print ".";

    if ( -d $fullpath )
    {
        my $success = File::Copy::Recursive::rcopy("$fullpath", "$basedir/$directory");
        if ( not $success )
        {
            die "Could not copy directory '$fullpath' to '$basedir/$directory': $!\n";
        }
    }
    else
    {
        die "Directory is missing from Agua distribution: $fullpath\n";
    }
}
print "\ndone\n\n";

#### LINK WEB DIR AND CGI DIR
print "Linking directories: html cgi-bin cgi-bin/lib cgi-bin/conf\n";
linkdir("$basedir/html", $webdir);
linkdir("$basedir/cgi-bin", $cgidir);
linkdir("$basedir/lib", "$basedir/cgi-bin/lib");
linkdir("$basedir/conf", "$basedir/cgi-bin/conf");
print "\ndone\n\n";


#### SET APPLICATION PATHS IF NECESSARY
my $install_groupUsers = "$Bin/../conf/groupUsers.json";
my $default_groupUsers = "$basedir/conf/default.json";
print "Copying install groupUsers to default groupUsers: $default_groupUsers\n";
File::Copy::copy($install_groupUsers, $default_groupUsers);

use Data::Dumper;
my $groupUsers = Config::JSON->new("$default_groupUsers");
my $apps = $groupUsers->get("apps");
$apps = applications($apps);
$groupUsers->set("apps", $apps);
$groupUsers->set("base", $basedir);
$groupUsers->set("url", $url);

print "Downloading dojo...\n";
my $targetdir = "$basedir/html";
chdir($targetdir) or die "Could not change to download target directory: $targetdir\n";
`wget $dojo_release`;
print "Download completed.\n";
print "Unzipping dojo release...\n";
my ($release_file) = $dojo_release =~ /\/([^\/]+)$/;
`tar xvfz $release_file`;
print "Unzip completed\n";
my ($dojo_dir) = $release_file =~ /^(.+)\.tar\.gz$/;
`mv $dojo_dir dojo`;
`rm -fr $release_file`;
print "done.\n";


print "\n";
print "*******************************************************\n";
print "*******************************************************\n";
print "Agua has been installed to this directory:\n";
print "    $basedir\n";
print "Browse to Agua here:\n";
print "    $url\n";
print "*******************************************************\n";
print "*******************************************************\n";
print "\n";
exit;




#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### ####             SUBROUTINES                 #### #### #### ####  
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 

=head2

    SUBROUTINE      installation_parameters

    PURPOSE

        SET PATHS FOR INSTALL DIR, HTML DIR, CGI DIR AND URL

=cut

sub installation_parameters
{
    my $basedir     =   shift;
    my $webdir      =   shift;
    my $cgidir      =   shift;
    my $url         =   shift;

    my $mode;
    if ( yes("Type 'Y' for the standard installation or type 'N' for custom installation") )    {   $mode  = "default"; }
    else {  $mode = "custom";    }

    if ( $mode =~ /^custom$/ )
    {
        print "Custom install selected!\n";

        print "Base directory [$basedir]:\n";
        my $input;
        $input = <STDIN>;
        chop($input);
        if ( $input !~ /^\s*$/ )
        {
            $basedir = $input;
        }

        print "Web directory [$webdir]:\n";
        $input = <STDIN>;
        chop($input);
        if ( $input !~ /^\s*$/ )
        {
            $webdir = $input;
        }

        print "CGI directory [$cgidir]:\n";
        $input = <STDIN>;
        chop($input);
        if ( $input !~ /^\s*$/ )
        {
            $cgidir = $input;
        } 

        print "Web URL [$url]:\n";
        $input = <STDIN>;
        chop($input);
        if ( $input !~ /^\s*$/ )
        {
            $url = $input;
        } 
    }

    return $basedir, $webdir, $cgidir, $url;
}



=head2

    SUBROUTINE      applications

    PURPOSE

        SET THE PATHS TO APPLICATIONS OR ADD NEW APPLICATIONS
=cut

sub applications
{
    my $apps        =   shift;

    print "Type Y to use the default application paths or type N to set custom application paths\n";
    return if yes("Type Y to use the default application paths or type N to set custom application paths");

    my $continue;
    while ( not defined $continue )
    {
        my $counter = 1;
        foreach my $app ( @$apps )
        {
            my @keys = keys %$app;
            my $key = $keys[0];
            print "$counter $key\t:\t$app->{$key}\n";
            $counter++;
        }

        $continue = input("Type application number or 'Q' for quit");

        if ( $continue !~ /^\d+$/ and $continue !~ /^q$/i )
        {
            print "\n**** Invalid input: $continue\n\n";
            $continue = undef;
            next;
        }

        #### CONTINUE MUST BE A NUMBER OR 'Q'
        if ( $continue =~ /^q$/i )
        {
            last;
        }
        else
        {
            my $app = $$apps[$continue - 1];
            if ( not defined $app )
            {
                print "No application entry for index: $continue\n";
            }
            else
            {
                my @keys = keys %$app;
                my $key = $keys[0];

                my $value = input("Enter path to $key [$app->{$key}]");
                $$apps[$continue - 1] = { $key => $value };    
            }

            $continue = undef;
        }
    }    

    return $apps;
}


=head2

    SUBROUTINE      input

    PURPOSE

        PROMPT THE USER TO ENTER 'Y' OR 'N'

=cut

sub input
{
	my $message		=	shift;
	return if ( not defined $message );
    print "$message\n";	

	$/ = "\n";
	my $input = <STDIN>;
	while ( $input =~ /^\s*$/ )
	{
		print "$message\n";
		$input = <STDIN>;
	}	

    chop($input);
	return $input;
}



=head2

    SUBROUTINE      link

    PURPOSE

        REMOVE EXISTING LINK AND MAKE SYMBOLIC LINK IN ITS PLACE

=cut

sub linkdir
{
    my $source      =   shift;
    my $target      =   shift;

    #### PRINT '.' AS PROGRESS COUNTER
    print ".";


    my $success = remove($target);

    $success = symlink($source, $target);
    die "Could not link source $source to target $target\n" if not $success;

    return $success;
}



=head2

    SUBROUTINE      yes

    PURPOSE

        PROMPT THE USER TO ENTER 'Y' OR 'N'

=cut

sub yes
{
	my $message		=	shift;
	return if ( not defined $message );
    my $max_times = 10;

	$/ = "\n";
	my $input = <STDIN>;
	my $counter = 0;
	while ( $input !~ /^Y$/i and $input !~ /^N$/i )
	{
		if ( $counter > $max_times )	{	print "Exceeded 10 tries. Exiting...\n";	}

		print "$message\n";
		$input = <STDIN>;
		$counter++;
	}	

	if ( $input =~ /^N$/i )	{	return 0;	}
	else {	return 1;	}
}

