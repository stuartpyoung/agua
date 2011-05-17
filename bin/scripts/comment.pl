#!/usr/bin/perl -w
use strict;

#### DEBUG

#### TIME
my $time = time();

=head2

    APPLICATION     comment.pl

    PURPOSE

        UNCOMMENT, COMMENT OR REMOVE COMMENTS IN DIFFERENT FILETYPES

    INPUT

        1. INPUT FILE /full/path/to/file

        2. OUTPUT FILE (MUST BE DIFFERENT TO INPUT FILE)

        3. FILE TYPE (E.G., js, perl), WHICH IS USED TO SPECIFY 

            TARGETED LINES:

                TYPE        TARGET



        4. ACTION TYPE (uncomment|comment|clean)

            uncomment   UNCOMMENT ONCE ALL TARGETED LINES

        5. REGEX TO SELECT LINES FOR (UN)COMMENTING

        6. COMMENT TEXT

    OUTPUT

        1. OUTPUT FILE WITH CHANGED LINES WITH RESPECT TO THE

            INPUT FILE (I.E., UNCOMMENTED, COMMENTED OR REMOVED LINES)

    USAGE

./comment.pl <--inputdir String> <--outputdir String> [--regex String] [--regex String] [--help]

    --inputfile     :   /Full/path/to/inputfile
    --outputfile    :   /Full/path/to/outputfile
    --inputdir      :   /Full/path/to/input/directory
    --outputdir     :   /Full/path/to/output/directory
    --type          :   Type of file (js|perl)
    --action        :   Action to perform (add | uncomment | comment | clean)
                        uncomment   :   Remove one comment from start of line
                        comment     :   Add one comment to the start of line
                        clean       :   Remove all commented lines
                        add         :    Add debug output at the start of every method
    --regex         :   Regex for filtering comment files
    --namedepth     :   Depth of directories to include in namespace
    --nameseparator :   Separator symbol in namespace, e.g., "." in js, "::" in perl
    --comment       :   Comment text (DEFAULTS: "//" for javscript, "#" for perl)
    --help          :   Print help info


    EXAMPLES

cd E:\0.4\bin\scripts

perl comment.pl --inputdir E:\0.4\html\plugins\view\jbrowse\js  --outputdir E:\0.4\html\plugins\view\jbrowse\js-added --action add --type js


perl comment.pl --inputfile E:\0.4\html\plugins\core\ExpandoPane.js  --outputfile E:\0.4\html\plugins\core\ExpandoPane.js --action add --type js


perl E:\0.4\bin\scripts\comment.pl --inputfile E:\0.4\html\dojo.1.2.2\dijit\layout\BorderContainer.js  --outputfile E:\0.4\html\dojo.1.2.2\dijit\layout\BorderContainer.js --action add --type js




cd C:\DATA\base\html\agua\bin\scripts
perl comment.pl 
--inputfile C:\DATA\base\html\agua\html\plugins\admin\TEST.js 
--outputfile C:\DATA\base\html\agua\html\plugins\admin\TEST-commented.js 
--action comment 
--type js 

cd C:\DATA\base\html\agua\bin\scripts
perl comment.pl --inputfile C:\DATA\base\html\agua\html\plugins\admin\TEST-commented.js --outputfile C:\DATA\base\html\agua\html\plugins\admin\TEST-uncommented.js --action uncomment --type js

cd C:\DATA\base\html\agua\bin\scripts
perl comment.pl --inputfile C:\DATA\base\html\agua\html\plugins\admin\TEST-uncommented.js --outputfile C:\DATA\base\html\agua\html\plugins\admin\TEST-uncommented-twice.js --action uncomment --type js





# CHANGE DIR

cd C:\DATA\base\html\agua\bin\scripts

# RUN comment.pl

perl comment.pl --inputfile C:\DATA\base\html\agua\html\plugins\project\Project-uncommented.js --outputfile C:\DATA\base\html\agua\html\plugins\project\Project-commented.js --action comment --type js


cd C:\DATA\base\html\agua\bin\scripts

FILE ADD COMMENTS:

perl comment.pl 
--inputfile C:\DATA\base\svn\hp-laptop\trunk\agua\html\dojo.1.2.2\dojox\widget\RollingList.js
 --outputfile C:\DATA\base\svn\hp-laptop\trunk\agua\html\dojo.1.2.2\dojox\widget\RollingList-add.js
 --action add
 --type js


DIRECTORY ADD COMMENTS:

perl comment.pl 
--inputdir C:\DATA\base\svn\hp-laptop\trunk\agua\html\dojo.1.2.2\dojox\widget
 --outputdir C:\DATA\base\svn\hp-laptop\trunk\agua\html\dojo.1.2.2\dojox\widget-add 
 --action add
 --type js





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

#### INTERNAL MODULES  
use Timer;
use FileTools;

#### GET OPTIONS
my $inputfile;
my $outputfile;
my $inputdir;
my $outputdir;
my $type;
my $action;
my $recursive;
my $regex;
my $function;
my $namedepth;
my $nameseparator;
my $comment;
my $help;
GetOptions (
    'inputfile=s' => \$inputfile,
    'outputfile=s' => \$outputfile,
    'inputdir=s' => \$inputdir,
    'outputdir=s' => \$outputdir,
    'type=s' => \$type,
    'action=s' => \$action,
    'recursive' => \$recursive,
    'function=s' => \$function,
    'regex=s' => \$regex,
    'namedepth=i' => \$namedepth,
    'nameseparator=i' => \$nameseparator,
    'comment=s' => \$comment,
    'help' => \$help) or die "No options specified. Try '--help'\n";

#### PRINT HELP IF REQUESTED
if ( defined $help )	{	usage();	}

#### INSTANTIATE FileTools OBJECT
my $filetools = FileTools->new();
$filetools->comment(
    {
        'inputfile'     => $inputfile,
        'outputfile'    => $outputfile,
        'inputdir'      => $inputdir,
        'outputdir'     => $outputdir,
        'action'        => $action,
        'type'          => $type,
        'comment'       => $comment,
        'regex'         => $regex,
        'recursive'     => $recursive,
        'function'      => $function,
        'namedepth'     => $namedepth,
        'nameseparator' =>  $nameseparator
    }
);

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

