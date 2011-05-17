#!/usr/bin/perl -w

#### DEBUG
$DEBUG = 01;

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

APPLICATION

    workflow

PURPOSE

    Manage and run individual or nested workflow and application files


USAGE

    ./workflow.pl mode [switch] [args] [--help]


mode     :    Type of workflow object (work|app|param)
switch   :    Nested object (e.g., work app, app param)
args     :    Arguments for the selected mode
--help   :    print help info

EXAMPLES

# create workflow
perl workflow.pl work create --wkfile /workflows/workflowOne.wk --name workflowOne

# create application from command file
perl workflow.pl app loadCmd --cmdfile /workflows/applicationOne.cmd --appfile /workflows/applicationOne.app --name applicationOne

# add application to workflow
perl workflow.pl work addApp --wkfile /workflows/workflowOne.wk --appfile /workflows/applicationOne.app --name applicationOne

# run single application
perl workflow.pl work app run --wkfile /workflows/workflowOne.wk --name applicationOne

# run all applications in workflow
perl workflow.pl work run --wkfile /workflows/workflowOne.wk 

=cut

use strict;
#use diagnostics;

#### USE LIBRARY
use Scalar::Util qw(weaken);
use FindBin qw($Bin);
use lib "$Bin/../../../lib/external";
use lib "$Bin/../../../lib";
BEGIN {
    unshift @INC, "/nethome/syoung/base/pipeline/moose/tmp/lib64/perl5/site_perl/5.8.8/x86_64-linux-thread-multi";
    unshift @INC, "/nethome/syoung/0.5/lib/external/perl5-32/site_perl/5.8.8";
    unshift @INC, "/nethome/syoung/0.5/lib/external/perl5-64/site_perl/5.8.8/x86_64-linux-thread-multi";
    unshift @INC, "/nethome/syoung/0.5/lib/external/perl5-32/5.8.8";    
}

#### EXTERNAL MODULES
use Term::ANSIColor qw(:constants);
use Data::Dumper;
#use Getopt::Long qw(GetOptionsFromString);
#use Getopt::Simple;
use X::Parameter;
use X::App;
use X::Workflow;


#### INTERNAL MODULES
use Timer;

#### GET MODE AND ARGUMENTS
my @arguments = @ARGV;
print "workflow.pl    No arguments supplied\n" and exit if not @arguments;

#### GET FILE
my $file = shift @ARGV;

#### GET SWITCH
#my $switch = shift @ARGV;
##print "Argument '$switch' not supported. Please use 'param', 'app' or 'work'\n"
#    and exit if $switch !~ /^(param|app|work)$/;

#### GET MODE
my $mode = shift @ARGV;
usage() if $mode =~ /^-h$/ or $mode =~ /^--help$/;

#### MANAGE INDIVIDUAL OR NESTED WORKFLOW FILES
if ( $file =~ /\.param$/ ) {
    my $parameter = X::Parameter->new(
        paramfile => $file
    );
    $parameter->getopts();
    $parameter->$mode();    
}
elsif ( $file =~ /\.app$/ ) {
    my $app = X::App->new(
        appfile => $file
    );
    $app->getopts();
    $app->$mode();    
}
elsif ( $file =~ /\.wk$/ )
{
    my $workflow = X::Workflow->new(
        wkfile => $file
    );
    $workflow->getopts();
    my $success = $workflow->$mode();
    print "workflow.pl    mode '$mode' not recognised\n" and exit if $success != 1;
}
else
{
    print "workflow.pl    file type '$file' not recognised (must be .wk, .app or .param)\n";
}

##### PRINT RUN TIME
#my $runtime = Timer::runtime( $time, time() );
#exit;


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#                                    SUBROUTINES
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


sub usage
{
    print GREEN;
    print `perldoc $0`;
    print RESET;
    exit;
}



