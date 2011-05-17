#!/usr/bin/perl -w
use strict;

#### DEBUG
$DEBUG = 1;

#### TIME
my $time = time();

=head2

    APPLICATION     svn.pl

    PURPOSE

        REMOVE .svn DIRS IN DIRECTORIES AND SUBDIRECTORIES

    INPUT

        1. INPUT DIRECTORY /full/path/to/file

    OUTPUT

        1. INPUT DIRECTORY WITH ALL .svn FILES REMOVED

    USAGE

    ./svn.pl <--workingdir String> [--help]

    --workingdir         :   /Full/path/to/input/directory
    --help              :   Print help info

    EXAMPLES


perl E:\agua\bin\utils\svn.pl \
--mode clean \
--workingdir C:\DATA\base\html\agua\html\plugins\project


/data/agua/0.3//bin/scripts/svn.pl \
--mode sync \
--workingdir /data/agua/0.3svn/html/plugins/project \
--repository file:///srv/svn/agua/trunk/html/plugins/project


/data/agua/0.3//bin/scripts/svn.pl \
--mode sync \
--workingdir /data/agua/0.3/html/plugins/project \
--repository file:///srv/svn/agua/trunk/html/plugins/workflow \
--message "Added modular CSS for OptionsTitlePane upload. Onworking poll for status"


#### SYNCHRONISE BETWEEN A WORKING DIRECTORY AND THE RESPOSITORY
#### INCLUDING AUTOMATIC ADD AND DELETE

/data/agua/0.3/bin/scripts/svn.pl \
--mode backup \
--base /data/agua/0.3svn \
--source html/plugins/project \
--repository file:///srv/svn/agua/trunk/html/plugins/project \
--message "Added modular CSS for OptionsTitlePane upload. Onworking poll for status"



=cut

#### FLUSH BUFFER
$| = 1;

#### EXTERNAL MODULES
use FindBin qw($Bin);
use Data::Dumper;
use Term::ANSIColor qw(:constants);
use Getopt::Long;

#### USE LIB
use lib "$Bin/../../lib";
use lib "E:/agua/lib/external";

#### INTERNAL MODULES  
use Timer;
use FileTools;

#### GET OPTIONS
my $base;
my $source;
my $workingdir;
my $repository;
my $message;
my $mode;
my $help;
GetOptions (
    'base=s' => \$base,
    'source=s' => \$source,
    'workingdir=s' => \$workingdir,
    'repository=s' => \$repository,
    'message=s' => \$message,
    'mode=s' => \$mode,
    'help' => \$help) or die "No options specified. Try '--help'\n";

#### PRINT HELP IF REQUESTED
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
print "Mode must be 'clean', 'sync' or 'backup'\n" and exit if $mode !~ /^(clean|sync|backup)$/;

print "workingdir not specified (option --workingdir)\n" and exit if not defined $workingdir and $mode eq 'clean';

print "Both workingdir and repository must be specified when running command: $mode\n"
	and exit
	if $mode eq 'sync'
	and ( not defined $workingdir or not defined $repository);

print "Base, source and repository must be specified when running command: $mode\n"
	and exit
	if $mode eq 'backup'
	and ( not defined $base or not defined $source or not defined $repository);


#### INSTANTIATE FileTools OBJECT
my $filetools = FileTools->new();

#### SET ARGS
my $args = {
	'base' 		=>	$base,
	'source' 		=>	$source,
	'workingdir' 		=>	$workingdir,
	'repository'		=>	$repository,
	'message'			=>	$message,
	'mode'				=>	$mode
};

#### ALTER MODE TO FIT FORMAT: svnModename
$mode =~ s/^(.)//;
$mode = "svn" . uc($1) . $mode;

#### RUN COMMAND
no strict;
my $result = $filetools->$mode($args);
use strict;

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "\nRun time: $runtime\n";
print "Completed $0\n";
print Timer::datetime(), "\n";
print "****************************************\n\n\n";
exit;

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### ####             SUBROUTINES                 #### #### #### ####  
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 


sub usage
{
	print GREEN;
    print `perldoc $0`;
	print RESET;

	exit;
}

