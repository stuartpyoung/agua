package Venn::Sam;


=head2

	PACKAGE		Venn::Sam

	PURPOSE

		PERFORM BASIC FILE FILTERING TASKS, SUCH AS:

			1. match: FIND MATCHING AND NON-MATCHING TRANSCRIPTS IN TWO

				.GTF FILES

			2. filterTranscripts: SELECT ONLY TRANSCRIPT FEATURES FROM A .GTF FILE

=cut

use strict;
use warnings;
use Carp;

#### USE LIB FOR INHERITANCE
use FindBin qw($Bin);
use lib "$Bin/..";

require Exporter;
our @ISA = qw(Venn);
use Venn;
our $AUTOLOAD;

#### INTERNAL MODULES
use Util;
use Filter;

#### EXTERNAL MODULES
use File::Remove;
use File::Copy;
use File::Path;

{
	use Data::Dumper;	
}

#### FLUSH BUFFER
#$| = 1;

#### SET SLOTS
our @DATA = qw(
	QUERYFILE
	TARGETFILE
	OUTPUTDIR
	QUERYLABEL
	TARGETLABEL
);
our $DATAHASH;
foreach my $key ( @DATA )	{	$DATAHASH->{lc($key)} = 1;	}


=head2

	SUBROUTINE		compare

	PURPOSE

		CALCULATE THE VENN INTERSECTION OF READ HITS IN TWO *.sam 

		FORMAT ALIGNMENT FILES

=cut

sub compare
{
	my $self		=	shift;

	my $queryfile	=	$self->get_queryfile();
	my $targetfile	=	$self->get_targetfile();
	my $outputdir	=	$self->get_outputdir();
	my $querylabel	=	$self->get_querylabel();
	my $targetlabel	=	$self->get_targetlabel();


	#### SET OUTPUT FILES
	my $queryonly;
	my $targetonly;
	my $both;
	($queryonly, $targetonly, $both) = $self->setOutputFiles($outputdir, $querylabel, $targetlabel, $queryonly, $targetonly, $both);

	#### DO QUERY AND BOTH
	my $querycount = 0;
	my $bothcount = 0;

	my ($query_only_lines, $target_only_lines, $query_both_lines, $target_both_lines);

	##### DO QUERY
	($query_both_lines, $query_only_lines) = $self->match($queryfile, $targetfile);

	#### DO TARGET
	($target_both_lines, $target_only_lines) = $self->match($targetfile, $queryfile);

	#### THESE TWO SHOULD BE IN AGREEMENT
	@$target_both_lines = sort @$target_both_lines;
	@$query_both_lines = sort @$query_both_lines;	
	#

	print "Venn::Sam::compare    target both != query both\n" if @$target_both_lines != @$query_both_lines;


	#### PRINT OUTPUT FILES
	open(QUERYONLY, ">$queryonly") or die "Can't open queryonly file: $queryonly\n";
	foreach my $line ( @$query_only_lines )	{	print QUERYONLY $line; }
	close(QUERYONLY);
	open(TARGETONLY, ">$targetonly") or die "Can't open targetonly file: $targetonly\n";
	foreach my $line ( @$target_only_lines )	{	print TARGETONLY $line; }
	close(TARGETONLY);
	open(BOTH, ">$both") or die "Can't open both file: $both\n";
	foreach my $line ( @$query_both_lines )	{	print BOTH $line; }
	close(BOTH);

	#### PRINT INFO
	print "Venn::Sam::compare    output files printed:\n";
	print "Venn::Sam::compare    $queryonly\n";
	print "Venn::Sam::compare    $targetonly\n";
	print "Venn::Sam::compare    $both\n";
	print "Venn::Sam::compare    \n";
}



#### SORT BY START POSITION
sub by_start
{
	my ($aa) = $a =~ /^\S+\s+\S+\s+transcript\s+(\d+)/;
	my ($bb) = $b =~ /^\S+\s+\S+\s+transcript\s+(\d+)/;

	$aa <=> $bb;	
}

#### SET DEFAULT OUTPUT FILES
sub setOutputFiles
{
	my $self		=	shift;
	my $outputdir	=	shift;
	my $querylabel	=	shift;
	my $targetlabel	=	shift;
	my $queryonly	=	shift;
	my $targetonly	=	shift;
	my $both		=	shift;

	#### SET OUTPUT FILES
	if ( not defined $queryonly )
	{
		$queryonly = "$outputdir/$querylabel-NOT-$targetlabel";
	}
	if ( not defined $targetonly )
	{
		$targetonly = "$outputdir/$targetlabel-NOT-$querylabel";
	}
	if ( not defined $both )
	{
		$both = "$outputdir/$querylabel-AND-$targetlabel";
	}


	return ($queryonly, $targetonly, $both);
}



=head2

	SUBROUTINE		match

	PURPOSE

		1. LOAD QUERY FILE INTO MEMORY AND SCROLL THROUGH

			TARGET FILE LOOKING FOR MATCHING AND NON-MATCHING

			TRANSCRIPTS 

	INPUTS

		1. QUERY AND TARGET .gtf FORMAT FILES WITH 'transcript' IN 3RD FIELD

			OF ALL LINES INDICATING THAT THE LINE IS A TRANSCRIPT

	OUTPUTS

		1. RETURN AN ARRAY OF MATCHING TRANSCRIPT LINES

		2. RETURN AN ARRAY OF THE NON-MATCHING TRANSCRIPT LINES

=cut



sub match
{
	my $self		=	shift;
	my $queryfile	=	shift;
	my $targetfile	=	shift;



	#### OPEN TARGET FILE
	my @targets = [];
	if ( $targetfile =~ /\.gz$/ )
	{
		my $pipe_command = "zcat $targetfile |";
		open(TARGET, $pipe_command);
		@targets = <TARGET>;
		close(TARGET);
	}
	else
	{
		open(TARGET, $targetfile) or die "Can't open target file: $targetfile\n";
		@targets = <TARGET>;
		close(TARGET);
	}

	my $matched = [];
	my $unmatched = [];

	#### OPEN QUERY FILE
	my $counter = 0;
	my $dot = 100000;
	open(QUERY, $queryfile) or die "Can't open target file: $targetfile\n";
	while ( <QUERY> )
	{
		next if $_ =~ /^\s*$/ or /^#/;

		$counter++;
		print "$counter\n" if $counter % $dot == 0;


		my ($queryid) = $_ =~ /^(\S+)/;
		$queryid =~ s/#.+$//;
		my $match = 0;
		foreach my $target ( @targets )
		{
			my ($targetid) = $target =~ /^(\S+)/;
			$targetid =~ s/#.+$//;


			if ( $targetid eq $queryid )
			{
				push @$matched, $target;
				$match = 1;
				last;
			}
		}

		push @$unmatched, $_ if not $match;
	}
	close(QUERY);

	return $matched, $unmatched;
}







################################################################################
##################			HOUSEKEEPING SUBROUTINES			################
################################################################################

=head2

	SUBROUTINE		initialise

	PURPOSE

		INITIALISE THE self OBJECT WITH USER-INPUT ARGUMENT VALUES

=cut

sub initialise
{
    my $self		=	shift;
	my $arguments	=	shift;

    #### VALIDATE USER-PROVIDED ARGUMENTS
	($arguments) = $self->validate_arguments($arguments, $DATAHASH);	

	#### LOAD THE USER-PROVIDED ARGUMENTS
	foreach my $key ( keys %$arguments )
	{		
		#### LOAD THE KEY-VALUE PAIR
		$self->value($key, $arguments->{$key});
	}
}


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

	#### INITIALISE THE OBJECT'S ELEMENTS
	$self->initialise($arguments);

    return $self;
}


=head2

	SUBROUTINE		value

	PURPOSE

		SET A PARAMETER OF THE self OBJECT TO A GIVEN value

    INPUT

        1. parameter TO BE SET

		2. value TO BE SET TO

    OUTPUT

        1. THE SET parameter INSIDE THE BioSVG OBJECT

=cut
sub value
{
    my $self		=	shift;
	my $parameter	=	shift;
	my $value		=	shift;

	$parameter = lc($parameter);

    if ( defined $value)
	{	
		$self->{"_$parameter"} = $value;
	}
}

=head2

	SUBROUTINE		validate_arguments

	PURPOSE

		VALIDATE USER-INPUT ARGUMENTS BASED ON

		THE HARD-CODED LIST OF VALID ARGUMENTS

		IN THE data ARRAY
=cut

sub validate_arguments
{
	my $self		=	shift;
	my $arguments	=	shift;
	my $DATAHASH	=	shift;

	my $hash;
	foreach my $argument ( keys %$arguments )
	{
		if ( $self->is_valid($argument, $DATAHASH) )
		{
			$hash->{$argument} = $arguments->{$argument};
		}
		else
		{
			warn "'$argument' is not a known parameter\n";
		}
	}

	return $hash;
}

=head2

	SUBROUTINE		is_valid

	PURPOSE

		VERIFY THAT AN ARGUMENT IS AMONGST THE LIST OF

		ELEMENTS IN THE GLOBAL '$DATAHASH' HASH REF

=cut

sub is_valid
{
	my $self		=	shift;
	my $argument	=	shift;
	my $DATAHASH	=	shift;

	#### REMOVE LEADING UNDERLINE, IF PRESENT
	$argument =~ s/^_//;

	#### CHECK IF ARGUMENT FOUND IN '$DATAHASH'
	if ( exists $DATAHASH->{lc($argument)} )
	{
		return 1;
	}

	return 0;
}

=head2

	SUBROUTINE		AUTOLOAD

	PURPOSE

		AUTOMATICALLY DO 'set_' OR 'get_' FUNCTIONS IF THE

		SUBROUTINES ARE NOT DEFINED.

=cut

sub AUTOLOAD {
    my ($self, $newvalue) = @_;


    my ($operation, $attribute) = ($AUTOLOAD =~ /(get|set)(_\w+)$/);

    # Is this a legal method name?
    unless($operation && $attribute) {
        croak "Method name $AUTOLOAD is not in the recognized form (get|set)_attribute\n";
    }
    unless( exists $self->{$attribute} or $self->is_valid($attribute) )
	{
        #croak "No such attribute '$attribute' exists in the class ", ref($self);
		return;
    }

    # Turn off strict references to enable "magic" AUTOLOAD speedup
    no strict 'refs';

    # AUTOLOAD accessors
    if($operation eq 'get') {
        # define subroutine
        *{$AUTOLOAD} = sub { shift->{$attribute} };

    # AUTOLOAD mutators
    }elsif($operation eq 'set') {
        # define subroutine4

        *{$AUTOLOAD} = sub { shift->{$attribute} = shift; };

        # set the new attribute value
        $self->{$attribute} = $newvalue;
    }

    # Turn strict references back on
    use strict 'refs';

    # return the attribute value
    return $self->{$attribute};
}


# When an object is no longer being used, this will be automatically called
# and will adjust the count of existing objects
sub DESTROY {
    my($self) = @_;

	#if ( defined $self->{_databasehandle} )
	#{
	#	my $dbh =  $self->{_databasehandle};
	#	$dbh->disconnect();
	#}

#    my($self) = @_;
}

1;
