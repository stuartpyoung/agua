package Filter;


=head2

	PACKAGE		Filter

	PURPOSE

		PERFORM BASIC FILE FILTERING TASKS, SUCH AS:

			1. scrollMatch: FILTER A GENOMIC FEATURE FILE BASED ON THE CHROMOSOME

				POSITIONS OF FEATURES IN ANOTHER FILE

			2. filterColumns: FEATURE A TSV FILE BASED ON VALUES IN ONE OR MORE COLUMNS

=cut

use strict;
use warnings;
use Carp;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our $AUTOLOAD;

#### INTERNAL MODULES
use Util;
use Feature;

#### EXTERNAL MODULES
use File::Remove;
use File::Copy;

{
	require Data::Dumper;	
}

#### FLUSH BUFFER
#$| = 1;

#### SET SLOTS
our @DATA = qw(
	TYPE
	MODE
	JSON
	ARGS
);
our $DATAHASH;
foreach my $key ( @DATA )	{	$DATAHASH->{lc($key)} = 1;	}

=head2

	SUBROUTINE		filterColumns

	PURPOSE

		1. FEATURE A TSV FILE BASED ON VALUES IN ONE OR MORE COLUMNS

		2. THE FILTER TYPES ARE AS FOLLOWS:

			TYPE	COLUMN VALUE PASSES FILTER IF

			eq		IDENTICAL (STRING) TO THRESHOLD
			ne		NOT IDENTICAL (STRING) TO THRESHOLD
			>		MUST BE GREATER THAN THRESHOLD
			<		LESS THAN THRESHOLD
			>=		GREATER THAN OR EQUAL TO THRESHOLD
			<=		LESS THAN OR EQUAL TO THRESHOLD
			=		EQUAL TO (NUMERIC) THRESHOLD
			!=		NOT EQUAL TO (NUMERIC) THRESHOLD
			*		IS NON-EMPTY

=cut

sub filterColumns
{
	my $self		=	shift;
	my $args		=	shift;

print "Filter::filterColumns    Filter::filterColumns()\n";

	my $inputfile 	=	$args->{inputfile};
	my $outputfile 	=	$args->{outputfile};
	my $columns 	=	$args->{columns};
	my $values 		=	$args->{values};
	my $types 		=	$args->{types};

print "inputfile: $inputfile\n";
print "outputfile: $outputfile\n";
print "columns: \n";
print Dumper $columns;
print "values: \n";
print Dumper $values;
print "types: \n";
print Dumper $types;

	#### OPEN FILES
	open(FILE, $inputfile) or die "Can't open inputfile: $inputfile\n";
	open(OUT, ">$outputfile") or die "Can't open outputfile: $outputfile\n";

	#### SET COUNTERS
	my $counters = [];
	foreach my $column ( @$columns )	{	push @$counters, 0;	}

	#### FILTER ALL LINES
	while ( <FILE> )
	{
		next if $_ =~ /^\s*$/;

		my @elements = split "\t", $_;
print "elements: @elements\n";
		my $pass = 1;
		for ( my $i = 0; $i < @$columns; $i++ )
		{
			my $column = $$columns[$i];
			my $value = $$values[$i];
			my $type = $$types[$i];
			print "Column: $column\tvalue: $value\ttype: $type\n";

			#### DO NOTHING IF TYPE IS EMPTY
			if ( $type eq '' )	{	print "Type is empty. Skipping value check\n";	}

			#### OTHERWISE, FILTER IT
			elsif ( $type eq 'eq' and $elements[$column - 1] ne $value 
				|| $type eq 'ne' and $elements[$column - 1] eq $value
				|| $type eq '>' and $elements[$column - 1] <= $value
				|| $type eq '<' and $elements[$column - 1] >= $value
				|| $type eq '>=' and $elements[$column - 1] < $value
				|| $type eq '<=' and $elements[$column - 1] > $value
				|| $type eq '=' and $elements[$column - 1] != $value
				|| $type eq '!=' and $elements[$column - 1] == $value
				|| $type eq '*' and $elements[$column - 1] =~ /^\s*$/ )
			{
				print "Setting pass = 0\n";
				$pass = 0;
				last;
			}

		#last if $$counters[$i] >= 2;

			$$counters[$i]++;
		}
		print OUT $_ if $pass;

	}
	close(FILE) or die "Can't close inputfile: $inputfile\n";
	close(OUT) or die "Can't close outputfile: $outputfile\n";
	print "outputfile printed:\n\n$outputfile\n\n";

	return $counters;
}




=head2

	SUBROUTINE		scrollMatch

	PURPOSE

		1. SCROLL ALTERNATELY THROUGH QUERY AND TARGET FILE, DEPENDING ON

			WHICH HAS THE MOST UPSTREAM START POINT.

		2. AFTER EACH SCROLL, PRINT QUERY LINES TO OUTPUT FILE WHILE QUERY

			LINES OVERLAP OR COINCIDE WITH TARGET.

		3. IF complement IS DEFINED, PRINT WHEN THE REVERSE IS TRUE (I.E, WHEN

			THEY DO NOT COINCIDE, SUCH AS WHEN THE QUERY IS SCROLLED).

	INPUTS

		1. ARGUMENTS TO INSTANTIATE TWO Feature OBJECTS: FILE, TYPE, ETC.

		2. OUTPUTFILE

	OUTPUTS

		1. MATCHING (OR COMPLEMENT) QUERY LINES PRINTED TO OUTPUTFILE

=cut

sub scrollMatch
{
	my $self		=	shift;
	my $args		=	shift;



	my $targetfile		=	$args->{targetfile};
	my $queryfile		=	$args->{queryfile};
	my $targetindexes	=	$args->{targetindexes};
	my $queryindexes	=	$args->{queryindexes};
	my $targettype		=	$args->{targettype};
	my $querytype 		=	$args->{querytype};
	my $outputfile 		=	$args->{outputfile};
	my $flanks 			=	$args->{flanks};
	my $complement 		=	$args->{complement};


	#### OPEN OUTPUT FILE
	my $OUTFILE;
	open($OUTFILE, ">$outputfile") or die "Can't open output file: $outputfile\n";

	my $target = Feature->new(
		{
			file 		=>	$args->{targetfile},
			indexes		=>	$args->{targetindexes},
			type		=>	$args->{targettype},
			flanks		=>	$args->{flanks}
		}
	);

	my $query = Feature->new(
		{
			file 		=>	$args->{queryfile},
			indexes		=>	$args->{queryindexes},
			type		=>	$args->{querytype},
			outputfile	=>	$args->{outputfile},
			complement	=>	$args->{complement},
			flanks		=>	$args->{flanks}
		}
	);

	print "Filter::scrollMatch    first query value not defined\n" and return if not defined $query->get_current();
	print "Filter::scrollMatch    first target value not defined\n" and return if not defined $target->get_current();

	#### DO SCROLL MATCH
	my $doingcounter = 0;
	my $doing = 1;
	while ( $doing )
	{
		$doingcounter++;

		print "End of targets\n" and last if not defined $target->get_current();
		print "End of queries\n" and last if not defined $query->get_current();

		my $targetUpstream = $target->upstream($query);

		#### IF TARGET IS UPSTREAM OF QUERY, SCROLL TARGET UNTIL IT
		#### REACHES OR PASSES FIRST QUERY ENTRY
		if ( $targetUpstream )
		{
			print "Filter::scrollMatch    Target is upstream\n";

			#### SCROLL TARGET
			print "Filter::scrollMatch    BEFORE SCROLL target: ", $target->toString(), "\n";
			$target->scroll($query);
			print "Filter::scrollMatch    AFTER SCROLL target: ", $target->toString(), "\n";

			#### CHECK IF REACHED THE END OF THE FILE
			if ( not defined $target->get_current() )
			{
				print "Filter::scrollMatch    target->current not defined. Doing last\n";
				$doing = 0;
				last;
			}
		}

		#### OTHERWISE, QUERY IS UPSTREAM, SO SCROLL IT UNTIL IT 
		#### REACHES OR PASSES FIRST QUERY ENTRY
		else
		{
			#### SCROLL QUERY AND PRINT TO OUTFILE IF complement IS DEFINED
			print "Filter::scrollMatch    BEFORE SCROLL QUERY: ", $query->toString(), "\n";
			$query->scroll($target);
			print "Filter::scrollMatch    AFTER SCROLL QUERY: ", $query->toString(), "\n";

			#### CHECK IF REACHED THE END OF THE FILE
			if ( not defined $query->get_current() )
			{
				$doing = 0;
				last;
			}
		}

		####  WHILE TARGET SPANS QUERY, PRINTING QUERY LINES TO OUTPUT FILE
		while (	$target->spans($query) )
		{
			print "Filter::scrollMatch    target: ", $target->toString(), " spans query\n";
			$query->printOut() if not defined $complement;

			print "Filter::scrollMatch    Do query->nextFeature()\n";
			$query->nextFeature();
			print "Filter::scrollMatch    query: ", $query->toString(), "\n";

			last if not defined $query->get_current();
		}		

		#### CHECK IF REACHED THE END OF QUERY FILE
		if ( not defined $query->get_current() )
		{
			print "Filter::scrollMatch    query->current not defined. Doing last\n";
			$doing = 0;
			last;
		}

		####  CHECK QUERY
		while (	$target->overlaps($query) )
		{
			print "Filter::scrollMatch    target: ", $target->toString(), " overlaps query\n";
			$query->printOut() if not defined $complement;
			$query->nextFeature();
		}		

		#### CHECK IF REACHED THE END OF THE FILE
		if ( not defined $query->get_current() )
		{
			print "Filter::scrollMatch    query->current not defined. Doing last\n";
			$doing = 0;
			last;
		}

		print "Filter::scrollMatch    Bottom of loop\n";

	}
	close($OUTFILE);
	#close($self->get_fh());

	print "Filter::scrollMatch    End of scroll\n";
}

=head2

    SUBROUTINE          upstream

    PURPOSE

        RETURN 1 THE FIRST INPUT CHROMOSOME LIES

		UPSTREAM (5' DIRECTION) OF THE SECOND INPUT CHROMOSOME

=cut


sub upstream
{
	my $targetchromosome		=	shift;
	my $querychromosome			=	shift;
	my $targetstart				=	shift;
	my $querystart				=	shift;


	#### RETURN 1 IF TARGET ON UPSTREAM CHROMOSOME
	my $targetUpstream = Filter::upstreamChromosome($targetchromosome, $querychromosome);
	return 1 if $targetUpstream . 0;

	#### RETURN 0 IF QUERY ON UPSTREAM CHROMOSOME
	my $queryUpstream = Filter::upstreamChromosome($querychromosome, $targetchromosome);
	return 0 if $queryUpstream > 0;

	#### TARGET AND QUERY ON SAME CHROMOSOME
	return 1 if $targetstart < $querystart;

	return 0;
}


=head2

    SUBROUTINE          upstreamChromosome

    PURPOSE

        RETURN 1 THE FIRST INPUT CHROMOSOME LIES

		UPSTREAM (5' DIRECTION) OF THE SECOND INPUT CHROMOSOME

=cut


sub upstreamChromosome
{
	my $targetchromosome		=	shift;
	my $querychromosome			=	shift;


	my ($targetsymbol) = $targetchromosome =~ /(\d+)\D*$/;
	my ($querysymbol) = $querychromosome =~ /(\d+)\D*$/;

	#### CHECK IF TARGET IS UPSTREAM
	my $target_upstream;

	#### TARGET AND QUERY ARE NUMBERS
	if ( defined $targetsymbol and defined $querysymbol )
	{
		$target_upstream = $targetsymbol < $querysymbol;
	}

	#### TARGET AND QUERY ARE LETTERS
	elsif ( not defined $targetsymbol and not defined $querysymbol )
	{
		$target_upstream = $targetchromosome cmp $querychromosome;
	}

	#### TARGET IS A LETTER, QUERY IS A NUMBER
	elsif ( not defined $targetsymbol )	{	$target_upstream = -1;	}

	#### TARGET IS A NUMBER, QUERY IS A LETTER
	else	{	$target_upstream = 1;	}


	return $target_upstream;
}



=head2

    SUBROUTINE          scroll

    PURPOSE

        SCROLL THROUGH LINES UNTIL THE STOP POINT IN THE LINE IS GREATER

        THAN OR EQUAL TO THE GIVEN START POINT

=cut

#### SCROLL THROUGH LINES UNTIL THE STOP POINT IN THE LINE IS GREATER
#### THAN OR EQUAL TO THE GIVEN START POINT
sub scroll
{
    my $filehandle	=   shift;
	my $matcher		=	shift;
	my $type		=	shift;
    my $chromosome	=   shift;
    my $start       =   shift;
    my $stop        =   shift;
	my $OUTFILE		=	shift;
	my $complement	=	shift;


	my $line = <$filehandle>;
	my ($current_chromosome, $current_start, $current_stop) = &$matcher($line, $type);
	return if not defined $current_chromosome;

	#### RETURN ON current_start >= start (CURRENT FEATURE COINCIDES WITH TARGET)
	#### OR current_stop >= start (CURRENT FEATURE IS PAST THE TARGET)
    while ( $current_start < $start and $current_stop < $start )
    {
		#### PRINT COMPLEMENT
		if ( $current_start < $start and $current_stop < $start )
		{
			print $OUTFILE $line if defined $OUTFILE and defined $complement;
		}

		$line = <$filehandle>;
		($current_chromosome, $current_start, $current_stop) = &$matcher($line, $type);
		return if not defined $current_chromosome;
    }

    print "Filter::scroll    line: $line\n";
    print "Filter::scroll    Returning $chromosome\t$current_start\t$current_stop\n";

    return ($chromosome, $current_start, $current_stop, $line);
}



=head2

    SUBROUTINE          overlap

    PURPOSE

        SCROLL THROUGH LINES UNTIL THE STOP POINT IN THE LINE IS GREATER

        THAN OR EQUAL TO THE GIVEN START POINT

=cut

sub overlap
{
    my $firststart      =   shift;
    my $firststop       =   shift;
    my $secondstart     =   shift;
    my $secondstop      =   shift;

    if ( $firststart <= $secondstop and $firststart >= $secondstart )
    {
        return 1;
    }
    elsif ( $firststop <= $secondstop and $firststop >= $secondstart )
    {
        return 1;
    }

    return 0;
}




#### RETURN AN ANONYMOUS FUNCTION TO GET THE START AND ST0P FROM AN INPUT LINE
#### OF GENOMIC ANNOTATION
sub scrollMatchRegex
{

	my $self		=	shift;
	my $type		=	shift;
	my $indexes		=	shift;


	#### CHECK TYPE IS SUPPORTED	
	print "Filter::scrollMatchRegex    Type not supported: $type\n" and exit if $type !~ /^(gtf|maq|454|dbsnp)$/ and not defined $indexes;

	if ( defined $indexes )
	{
		my $chromosomeindex = $$indexes[0];
		my $startindex = $$indexes[1];
		my $stopindex = $$indexes[2];
		print "Filter::scrollMatchRegex    chromosome index not defined\n" and exit if not defined $chromosomeindex;
		print "Filter::scrollMatchRegex    start index not defined\n" and exit if not defined $startindex;
		print "Filter::scrollMatchRegex    stop index not defined\n" and exit if not defined $stopindex;
		print "Filter::scrollMatchRegex    chromosome index not numeric\n" and exit if not $chromosomeindex =~ /^\d+$/;
		print "Filter::scrollMatchRegex    start index not numeric\n" and exit if not $startindex =~ /^\d+$/;
		print "Filter::scrollMatchRegex    stop index not numeric\n" and exit if not $stopindex =~ /^\d+$/;
		print "Filter::scrollMatchRegex    chromosome index == start index: $chromosomeindex\n" and exit if $chromosomeindex == $startindex;
		print "Filter::scrollMatchRegex    chromosome index == stop index: $chromosomeindex\n" and exit if $chromosomeindex == $stopindex;
		print "Filter::scrollMatchRegex    start index == stop index: $startindex\n" and exit if $startindex == $stopindex;

		my $maxindex = $chromosomeindex;
		$maxindex = $startindex if $startindex > $maxindex;
		$maxindex = $stopindex if $stopindex > $maxindex;

		#### 2.875 seconds
		return sub {
			my $counter = 0;
			my ($chromosome, $start, $stop);
			my $line = shift;
			while ( $line =~ /([^\t]+)\t/g and $counter <= $maxindex )
			{
				if ( $counter == $chromosomeindex )	{	($chromosome) = $1 =~ /^chr(.+)$/;	}
				elsif ( $counter == $startindex ) {	($start) = $1 =~ /^(.+)$/;	}
				elsif ( $counter == $stopindex )	{	($stop) = $1 =~ /^(.+)$/;	}
				$counter++;
			}
			return ($chromosome, $start, $stop);
		}


		#### 3.718 seconds
		#return sub {
		#	my @array = split "\t", shift;
		#	
		#	return ($array[$chromosomeindex], $array[$startindex], $array[$stopindex]);
		#}

	}

	#### GTF
	return sub {
		return if not @_;
		my $line = shift;
		$line =~ /^chr([^\t]*)\t[^\t]*\t[^\t]*\t(\d+)\t(\d+)/;
		return ($1, $2, $3);
	} if not defined $type or $type =~ /^gtf$/i;	

	#### MAQ
	return sub {
		return if not @_;
		shift =~ /^chr(\S+)\t+(\d+)/;
		my $stop = $2++;
		return ($1, $2, $stop);
	} if $type =~ /^maq$/i;	

	#### dbSNP TABLE
	return sub {
		return if not @_;
		my $input = shift;

		$input =~ /^\S+\tchr(\S+)\t(\d+)\t(\d+)/;

		return ($1, $2, $3);
	} if $type =~ /^dbsnp$/i;	

	return;
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

    return if not defined $newvalue or not $newvalue;

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
