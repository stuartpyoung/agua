#!/usr/bin/perl -w
use strict;

$DEBUG = 1;

#### TIME
my $time = time();

=head2

APPLICATION     userGroups

PURPOSE:

    PARSE USERS AND THEIR GROUPS FROM THE /etc/groups FILE AND

    OUTPUT THEM IN userGroups TABLE FORMAT:

        username    userid  groupname   groupid

INPUT:

    /etc/group FILE 


OUTPUT:

    1. TSV FILE IN userGroups TABLE FORMAT


=cut

#### FLUSH BUFFER
$| = 1;

#### EXTERNAL MODULES
use FindBin qw($Bin);

#### USE LIB
#use lib "$Bin/../lib";
#use lib "$Bin/../lib/external";
#
##### INTERNAL MODULES
##use File::Copy;
#use File::Copy::Recursive;
#use File::Remove 'remove';
#use Config::JSON;

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

