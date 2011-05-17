package Converter;


=head2

	PACKAGE		Converter

	PURPOSE

		RUN FILE FORMAT CONVERSION TOOLS

=cut

use strict;
use warnings;
use Carp;

#### INTERNAL MODULES
use Sampler;

#### HAS A
use Cluster;

#### EXTERNAL MODULES
use Data::Dumper;
use File::Path;
#use MD5;


require Exporter;
our @ISA = qw(Exporter Cluster);
#our @EXPORT_OK = qw();
our $AUTOLOAD;


#### SET SLOTS
our @DATA = qw(

INPUTFILES
OUTPUTFILES
MATEFILES
REFERENCEDIR
OUTPUTDIR
TEMPDIR
MAXLINES
CLEAN
SPLITFILE
SOLEXA

PARAMS
LABEL
CLUSTER
SAMTOOLS
SAMTOOLSINDEX
SPECIES
VERBOSE
CONVERT

MAXJOBS
CPUS
CLUSTER
QUEUE
WALLTIME
QSTAT
QSUB
SLEEP
CLEANUP
VERBOSE
TEMPDIR
DOT

BATCHSTATS
COMMAND
);
our $DATAHASH;
foreach my $key ( @DATA )	{	$DATAHASH->{lc($key)} = 1;	}





################################################################################
##################			HOUSEKEEPING SUBROUTINES			################
################################################################################


=head2

	SUBROUTINE		new

	PURPOSE

		CREATE THE NEW self OBJECT AND INITIALISE IT, FIRST WITH DEFAULT 

		ARGUMENTS, THEN WITH PROVIDED ARGUMENTS

=cut

sub new
{
 	my $class 		=	shift;
	my $arguments 	=	shift;

	my $self = {};
    bless $self, $class;

	#### INITIALISE THE OBJECT'S ELEMENTS
	$self->initialise($arguments);

	#### SET REFERENCE
	if ( defined $self->{_referencefile} )
	{
		my ($reference) = $self->{_referencefile} =~ /([^\/]+)$/;
		$reference =~ s/\.[^\.]+$//;

		$self->{_reference} = $reference;
	}

	#### SET DEBUG IF verbose

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

    #### VALIDATE USER-PROVIDED ARGUMENTS
	($arguments) = $self->validate_arguments($arguments, $DATAHASH);	

	#### LOAD THE USER-PROVIDED ARGUMENTS
	foreach my $key ( keys %$arguments )
	{		
		#### LOAD THE KEY-VALUE PAIR
		$self->value($key, $arguments->{$key});
	}

	#### SET TIME
	$self->{_starttime} = time();

	#### REQUIRED: CALL THE PARENT FOR ADDITIONAL INITIALISATION
	$self->SUPER::initialise();		
}


1;


