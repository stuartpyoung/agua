package Hapmap;


=head2

		PACKAGE		Hapmap

		VERSION		0.02

		PURPOSE

	        WRAPPER SCRIPT FOR RUNNING Hapmap SNP PREDICTION

		HISTORY
					0.02 ADDED CHUNK-BY-CHROMOSOME IF referencedir SPECIFIED
						AND RUN CUFFLINKS COMMAND
					0.01 BASIC VERSION
=cut

use strict;
use warnings;
use Carp;

require Exporter;
our @ISA = qw(Exporter);
#our @EXPORT_OK = qw();
our $AUTOLOAD;

#### EXTERNAL MODULES
use Data::Dumper;
use File::Path;
use FindBin qw($Bin);

#### SET SLOTS
our @DATA = qw(

INPUTFILE
OUTPUTFILE
IDS

);
our $DATAHASH;
foreach my $key ( @DATA )	{	$DATAHASH->{lc($key)} = 1;	}


=head2

	SUBROUTINE		uniqueHeterozygotes

	PURPOSE

		SELECT HAPMAP SNPS FOR WHICH THE INDIVIDUAL IN QUESTION

		DIFFERS FROM ALL OTHER MEMBERS OF THE INPUT GROUP IDS

=cut

sub uniqueHeterozygotes
{
	my $self		=	shift;
	my $inputfile	=	shift;
	my $outputfile	=	shift;
	my $ids			=	shift;


	my $idhash = {};
	my @idarray = split ",", $ids;
	foreach my $entry ( @idarray )
	{
		$idhash->{$entry} = 1;
	}

	#### GET INPUT DATA
	open(FILE, $inputfile) or die "Can't open inputfile: $inputfile\n";
	my @lines = <FILE>;
	close(FILE);

	#### GET HAPMAPS THAT ARE HETEROZYGOTE FOR EACH INDIVIDUAL
	my $heteros = {};
	foreach my $id ( @idarray )
	{
		$heteros->{$id} = $self->filterHeteros(\@lines, $id, $idhash);
	}


	open(OUT, ">$outputfile") or die "Hapmap::uniqueHeterozygotes    Can't open outputfile: $outputfile\n";
	foreach my $id ( @idarray )
	{
		my $lines = $heteros->{$id};
		print OUT @$lines;
	}
	close(OUT) or die "Hapmap::uniqueHeterozygotes    Can't close outputfile: $outputfile\n";

	print "Hapmap::uniqueHeterozygotes    Printed outputfile:\n\n$outputfile\n\n";
}	#	uniqueHeterozygotes


sub filterHeteros
{
	my $self		=	shift;
	my $lines		=	shift;
	my $individual	=	shift;
	my $idhash		=	shift;

	#

	#### FILTER HETEROZYGOTES FOR OUR INDIVIDUAL OF INTEREST
	my @heteros = ();
	my $heteros = [];
	foreach my $line ( @$lines )
	{
		next if $line =~ /^#/;

		#### REM: FORMAT:
		#### marker id	chromosome	position	strand	list of sample IDs	genotype
		my @elements = split "\t", $line;
		my @idlist = split " ", $elements[4];
		my @genotypes = split " ", $elements[5];
		my @group_genotypes = [];

		#### HASHIFY ID-GENOTYPE PAIRS
		my $id_genotypes = {};
		for ( my $i = 0; $i < @idlist; $i++ ) 
		{
			my $id = $idlist[$i];
			my $genotype = $genotypes[$i];
			$id_genotypes->{$id} = $genotype;
		}

		#### CHECK IF OUR INDIVIDUAL IS HETEROZYGOUS
		my $individual_genotype = $id_genotypes->{$individual};

		#### WARN IF INDIVIDUAL DOES NOT HAVE A GENOTYPE.
		#### AND SKIP THIS LINE		
		next if not defined $individual_genotype;

		#### SKIP IF HOMOZYGOTE
		my $isHomozygote = $self->isHomozygote($individual_genotype);
		next if $isHomozygote;

		#### WE'RE WORKING WITH ONLY HETEROZYGOTES FROM NOW ON


		#### CHECK TO MAKE SURE OTHERS IN GROUP ARE HOMOZYGOTE
		my @idkeys = keys ( %$idhash );
		my $heterozygote_found = 0;
		foreach my $id ( @idkeys )
		{
			next if $id eq $individual;
			next if not defined $id_genotypes->{$id};

			$heterozygote_found = 1 if not $self->isHomozygote($id_genotypes->{$id});
		}
		next if $heterozygote_found == 1;

		my $heteroline = "$elements[0]\t$elements[1]\t$elements[2]\t$elements[3]\t$individual\t$individual_genotype\n";
		push @$heteros, $heteroline;
	}

	print @$heteros;

	return $heteros;	
}



=head 		SUBROUTINE		isHomozygote

=head2		PURPOSE

		RETURN 1 IF HOMOZYGOTE, 0 OTHERWISE

=cut

sub isHomozygote
{
	my $self		=	shift;
	my $genotype	=	shift;

	return if not defined $genotype;

	my ($a, $b) = split "", $genotype;


	return 1 if $a eq $b;
	return 0;
}


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

    return $self;
}

=head2

	SUBROUTINE		initialise

	PURPOSE

		INITIALISE THIS OBJECT:

			1. LOAD THE ARGUMENTS

			2. CALL PARENT initialise TO LOAD ADDITIONAL INFO

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
    unless( exists $self->{$attribute} )
	{
        #croak "No such attribute '$attribute' exists in the class ", ref($self);
		#return undef;
    }

    # Turn off strict references to enable "magic" AUTOLOAD speedup
    no strict 'refs';

    # AUTOLOAD accessors
    if($operation eq 'get') {

		return $self->{$attribute};

        # define subroutine
        *{$AUTOLOAD} = sub { shift->{$attribute} };

    # AUTOLOAD mutators
    }elsif($operation eq 'set') {


        ## define subroutine
        #*{$AUTOLOAD} = sub { shift->{$attribute} = shift; };

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

}





1;

