package Feature;


=head2

	PACKAGE		Feature

	PURPOSE

		A GENOMIC FEATURE OBJECT REPRESENTING THE FEATURES WITHIN THE

		USER-SPECIFIED INPUT FILE, WITH THE FOLLOWING BEHAVIOURS:

			-	COMPARISON WITH OTHER FEATURE OBJECTS TO DETERMINE WHICH IS

				UPSTREAM

			-	SCROLL THROUGH FEATURES WHILE UPSTREAM OF ANOTHER FEATURE

				OBJECT'S CURRENT FEATURE

			-	PRINT TO OUTPUT FILE IF FEATURE OVERLAPS TARGET FEATURE

				(OR THE CONVERSE WHEN complement ARGUMENT IS SPECIFIED)


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

#### EXTERNAL MODULES
use File::Remove;
use File::Copy;

{
	use Data::Dumper;	
}

#### FLUSH BUFFER
#$| = 1;

#### SET SLOTS
our @DATA = qw(
TYPE
FILE
INDEXES
OUTPUTFILE
COMPLEMENT
RIGHTFLANK
LEFTFLANK
FLANKS
CHROMOSOMEINDEX
STARTINDEX
STOPINDEX
MAXINDEX
FH
OUTFH
);
our $DATAHASH;
foreach my $key ( @DATA )	{	$DATAHASH->{lc($key)} = 1;	}


=head2

	SUBROUTINES		get_chromsomome, get_start, get_stop

	PURPOSE

		RETURN THE RESPECTIVE VALUES OF self->current

=cut

sub get_chromosome	{
	my $self		=	shift;
	return $self->get_current()->{chromosome};
}
sub get_start	{
	my $self		=	shift;
	return $self->get_current()->{start};
}
sub get_stop	{
	my $self		=	shift;
	return $self->get_current()->{stop};
}



=head2

    SUBROUTINE      fillAndPrint

    PURPOSE

        FILL OUT self->line COLUMNS TO USER-DEFINED NUMBER AND PRINT OUT

=cut

sub fillAndPrint
{
    my $self            =   shift;
    my $number            =   shift;


	my @elements = split "\t", $self->get_current()->{line};

	print "Feature::fillAndPrint    elements: \n";
	print join "**\n", @elements;
	print "\n";

	while ( $#elements < $number - 1 )
	{
		push @elements, "";
	}
	print "Feature::fillAndPrint    no. elements: " , $#elements + 1 , "\n";
	print "Feature::fillAndPrint    elements: \n";
	print join "**\n", @elements;
	print "Feature::fillAndPrint    \n";
	#exit;

	$self->get_current()->{line} = join "\t", @elements;

	$self->printOut();
}






=head2

	SUBROUTINE		printOut

	PURPOSE

		PRINT TO THE OUTPUT FILE:

			1. OUTPUT IF DEFINED, OR

			2. THE CURRENT LINE

=cut

sub printOut
{
	my $self		=	shift;
	my $output		=	shift;


	my $outfh = $self->get_outfh();
	print $outfh $output if defined $output;
	return if defined $output;

	my $line = $self->get_current()->{line};
	$line .= "\n" if not $line =~ /\n$/;

	print $outfh $line;
}


=head2

	SUBROUTINE		toString

	PURPOSE

		RETURN THE chromosome, start AND stop VALUES OF self->current IN A STRING

=cut

sub toString
{
	my $self		=	shift;


	my $current = $self->get_current();
	return "<EMPTY>" if not defined $current;
	return "$current->{chromosome}\t$current->{start}..$current->{stop}";
}


=head2

	SUBROUTINE		setFiles

	PURPOSE

		SET THE INPUT AND OUTPUT FILE HANDLES

=cut

sub setFiles
{
	my $self		=	shift;



	#### OPEN OUTPUT FILE
	my $outputfile = $self->get_outputfile();
	if ( defined $outputfile )
	{
		my $OUTFILE;
		open($OUTFILE, ">$outputfile") or die "Can't open output file: $outputfile\n";
		$self->{_outfh} = $OUTFILE;
	}

	#### CHECK INPUTS
	my $file = $self->get_file();
	print "Feature::setFiles    file not defined\n" and exit if not defined $file;
	print "Feature::setFiles    file not found: $file\n" and exit if not -f $file;
	print "Feature::setFiles    First file is empty: $file\n" and exit if -z $file;

	#### OPEN FILES
	my $handle;
	open($handle, $file) or die "Can't open file: $file\n";
	$self->{_fh} = $handle;	
}


=head2

	SUBROUTINE		overlaps

	PURPOSE

		DETERMINE WHETHER self->current FEATURE OVERLAPS OR SPANS THE CURRENT

		FEATURE OF ANOTHER Feature OBJECT

=cut

sub nextFeature
{
	my $self		=	shift;


	#### PARSE NEXT LINE OF INPUT FILE
	my ($chromosome, $start, $stop, $line);
	my $handle = $self->get_fh();
	my $matcher = $self->get_matcher();

	while ( not defined $chromosome )
	{
		if ( not defined $handle )
		{
			$self->{_current} = undef;
			return;
		}
		$line = <$handle>;

		if ( not defined $line )
		{
			$self->{_current} = undef;
			return;
		}

		#$self->{_current} = undef and return if not defined $line;
		next if $line =~ /^#/ or $line =~ /^\s*$/;

		($chromosome, $start, $stop) = &$matcher($line);
	}


#### **** ####
#### **** ####

	#### LATER: TEST RETURN UNDEF
	$self->{_current} = undef if not defined $chromosome;
	$self->{_current} = undef if not defined $start;
	$self->{_current} = undef if not defined $stop;

#### **** ####
#### **** ####

	$self->{_current} = (
		{
			chromosome 	=>	$chromosome,
			start		=>	$start,
			stop		=>	$stop,
			line		=>	$line
		}
	);

	return $self->{_current};
}


=head2

	SUBROUTINE		overlaps

	PURPOSE

		DETERMINE WHETHER self->current FEATURE OVERLAPS THE CURRENT

		FEATURE OF ANOTHER Feature OBJECT

=cut

sub overlaps
{
	my $self		=	shift;
	my $other		=	shift;

	return if not defined $self->get_current();
	return if not defined $other->get_current();

	my $selfstart = $self->get_current()->{start} || return;
	my $selfstop = $self->get_current()->{stop} || return;

	my $otherstart = $other->get_current()->{start} || return;
	my $otherstop = $other->get_current()->{stop} || return;

	my $leftflank = $self->get_leftflank();
	my $rightflank = $self->get_rightflank();

print "Feature::overlaps    Feature::overlaps(other)\n";
print "Feature::overlaps    self: $selfstart..$selfstop\n";
print "Feature::overlaps    leftflank: $leftflank, rightflank: $rightflank\n";
print "Feature::overlaps    other: $otherstart..$otherstop\n";



	return 1 if 
	(
		#### SELF STRADDLES OTHER START
		( ($selfstart - $leftflank) <= $otherstart and ($selfstop + $rightflank) >= $otherstart)
		or 
		#### SELF STRADDLES OTHER STOP
		( ($selfstart - $leftflank) <= $otherstop and ($selfstop + $rightflank) >= $otherstart )
	);

	return 0;
}



=head2

	SUBROUTINE		spans

	PURPOSE

		DETERMINE WHETHER self->current FEATURE SPANS THE CURRENT

		FEATURE OF ANOTHER Feature OBJECT

=cut

sub spans
{
	my $self		=	shift;
	my $other		=	shift;

	my $selfstart = $self->get_current()->{start} || return;
	my $selfstop = $self->get_current()->{stop} || return;

	my $otherstart = $other->get_current()->{start} || return;
	my $otherstop = $other->get_current()->{stop} || return;

	my $leftflank = $self->get_leftflank();
	my $rightflank = $self->get_rightflank();

	#### THE self->current FEATURE SPANS THE OTHER (OR VICE-VERSA)
	return 1 if
	(
		( ($otherstart - $leftflank) <= $selfstart and ($otherstop + $rightflank) >= $selfstop )
		or 
		( $otherstart >= ($selfstart - $leftflank) and $otherstop <= ($selfstop + $rightflank) )
	);

	return 0;
}



=head2

    SUBROUTINE          upstream

    PURPOSE

        RETURN 1 IF THE FIRST INPUT CHROMOSOME LIES UPSTREAM

		(5' DIRECTION) OF THE SECOND INPUT CHROMOSOME.

		IF THE CHROMOSOMES ARE THE SAME, RETURN 1 IF THE FIRST

		START IS UPSTREAM OF THE SECOND START.

		OTHERWISE, RETURN 0.

=cut

sub upstream
{
	my $self					=	shift;
	my $other					=	shift;

	my $selfchromosome			=	$self->get_current()->{chromosome};
	my $selfstart				=	$self->get_current()->{start};
	my $otherchromosome			=	$other->get_current()->{chromosome};
	my $otherstart				=	$other->get_current()->{start};


	#### RETURN 1 IF TARGET ON UPSTREAM CHROMOSOME
	my $selfUpstream = $self->upstreamChromosome($selfchromosome, $otherchromosome);
	return 1 if $selfUpstream > 0;

	#### RETURN 0 IF QUERY ON UPSTREAM CHROMOSOME
	my $otherUpstream = $self->upstreamChromosome($otherchromosome, $selfchromosome);
	return 0 if $otherUpstream > 0;

	#### TARGET AND QUERY ON SAME CHROMOSOME
	return 1 if $selfstart < $otherstart;

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
	my $self				=	shift;	
	my $selfchromosome		=	shift;
	my $otherchromosome		=	shift;


	my ($selfsymbol) = $selfchromosome =~ /(\d+)\D*$/;
	my ($othersymbol) = $otherchromosome =~ /(\d+)\D*$/;

	#### CHECK IF TARGET IS UPSTREAM
	my $self_upstream;

	#### TARGET AND QUERY ARE NUMBERS
	if ( defined $selfsymbol and defined $othersymbol )
	{
		$self_upstream = $selfsymbol < $othersymbol;
	}

	#### TARGET AND QUERY ARE LETTERS
	elsif ( not defined $selfsymbol and not defined $othersymbol )
	{
		$self_upstream = $selfchromosome cmp $otherchromosome;
	}

	#### TARGET IS A LETTER, QUERY IS A NUMBER
	elsif ( not defined $selfsymbol )	{	$self_upstream = -1;	}

	#### TARGET IS A NUMBER, QUERY IS A LETTER
	else	{	$self_upstream = 1;	}


	return $self_upstream;
}



=head2

    SUBROUTINE          scroll

    PURPOSE

        SCROLL THROUGH LINES UNTIL THE STOP POINT IN THE LINE IS GREATER

        THAN OR EQUAL TO THE GIVEN START POINT

=cut

sub scroll
{
	my $self		=	shift;
	my $other		=	shift;


	return if not defined $other;

	#### 	
    my $chromosome	=   $other->get_chromosome() || return;
    my $start       =   $other->get_start() || return;
    my $stop        =   $other->get_stop() || return;

    my $filehandle	=   $self->get_fh();
	my $matcher		=	$self->get_matcher();
	my $OUTFILE		=	$self->get_outfh();
	my $complement	=	$self->get_complement();
	my $flanks		=	$self->get_flanks();


	$self->printOut() if defined $complement;
	$self->nextFeature();
	my $current = $self->get_current();
	return if not defined $current;

	#### RETURN ON current_start >= start (CURRENT FEATURE COINCIDES WITH TARGET)
	#### OR current_stop >= start (CURRENT FEATURE IS PAST THE TARGET)
    while ( $current->{start} < $start and $current->{stop} < $start )
    {
		#### PRINT COMPLEMENT
		if ( $current->{start} < $start and $current->{stop} < $start )
		{
			$self->printOut() if defined $complement;
		}

		$current = $self->nextFeature();
		return if not defined $current;
    }

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





=head2

    SUBROUTINE          setMatcher

    PURPOSE

		RETURN AN ANONYMOUS FUNCTION TO GET THE CHROMOSOME, START AND ST0P

		FROM AN INPUT LINE OF GENOMIC FEATURE ANNOTATION

=cut

sub setMatcher
{

	my $self		=	shift;

	my $type		=	$self->get_type();
	my $indexes		=	$self->get_indexes();


	#### CHECK TYPE IS SUPPORTED	
	print "Feature::setMatcher    Type not supported: $type\n" and exit if $type !~ /^(gtf|maq|454|dbsnp|pileup)$/ and not defined $indexes;

	my $subroutine;

	#### INDEXES
	if ( defined $self->get_indexes() )
	{
		my $chromosomeindex = $self->get_chromosomeindex();
		my $startindex = $self->get_startindex();
		my $stopindex = $self->get_stopindex();
		my $maxindex = $self->get_maxindex();

		#### 2.875 seconds
		$subroutine = sub {
			return if not @_;

			my $counter = 0;
			my ($chromosome, $start, $stop);
			my $line = shift;
			while ( $line =~ /([^\t]+)\t/g and $counter <= $maxindex )
			{
				next if not defined $1;
				if ( $counter == $chromosomeindex )	{	($chromosome) = $1 =~ /^chr(.+)$/;	}
				elsif ( $counter == $startindex ) {	($start) = $1 =~ /^(.+)$/;	}
				elsif ( $counter == $stopindex )	{	($stop) = $1 =~ /^(.+)$/;	}
				$counter++;
			}
			return ($chromosome, $start, $stop);
		}

		#### 3.718 seconds
		#$subroutine = sub {
		#	my @array = split "\t", shift;
		#	
		#	return ($array[$chromosomeindex], $array[$startindex], $array[$stopindex]);
		#}
	}
	#### 454
	#### >chrY	10621724	10621724
	elsif ( $type =~ /^454$/ )
	{
		$subroutine = sub {
			return if not @_;

			my $line = shift;
			$line =~ /^.*?chr(\S+)\s+(\d+)\s+(\d+)\s+/;
			return ($1, $2, $3);
		} 
	}
	#### GTF
	elsif ( $type =~ /^gtf$/i )
	{
		$subroutine = sub {
			return if not @_;
			my $line = shift;
			$line =~ /^chr([^\t]*)\t[^\t]*\t[^\t]*\t(\d+)\t(\d+)/;
			return ($1, $2, $3);
		} 
	}
	#### MAQ
	elsif ( $type =~ /^maq$/i )
	{
		$subroutine = sub {
			return if not @_;
			my $line = shift;
			$line =~ /^chr(\S+)\s+(\d+)/;
			my $stop = $2;
			$stop++;			
			return ($1, $2, $stop);
		} 
	}
	#### dbSNP
	elsif ( $type =~ /^dbsnp$/i )
	{
		$subroutine = sub {
			return if not @_;
			my $input = shift;

			$input =~ /^\S+\tchr(\S+)\t(\d+)\t(\d+)/;

			return ($1, $2, $3);
		} 
	}
	#### pileup
	elsif ( $type =~ /^pileup$/i )
	{
		$subroutine = sub {
			return if not @_;
			my $line = shift;
			$line =~ /^chr(\S+)\s+(\d+)/;
			my $stop = $2;
			$stop++;			
			return ($1, $2, $stop);
		} 
	}
	else
	{
		print "Feature::setMatcher    type not supported: $type\n" and exit;
	}

	$self->{_matcher} = $subroutine;
}


=head2

	SUBROUTINE		setIndexes

	PURPOSE

		SET REGEX MATCH INDEXES FOR chromosome, start AND stop IF DEFINED

=cut

sub setIndexes
{
	my $self		=	shift;


	my $indexes = $self->get_indexes();
	return if not defined $indexes;


	my $chromosomeindex = $$indexes[0];
	my $startindex = $$indexes[1];
	my $stopindex = $$indexes[2];
	print "Feature::setMatcher    chromosome index not defined\n" and exit if not defined $chromosomeindex;
	print "Feature::setMatcher    start index not defined\n" and exit if not defined $startindex;
	print "Feature::setMatcher    stop index not defined\n" and exit if not defined $stopindex;
	print "Feature::setMatcher    chromosome index not numeric\n" and exit if not $chromosomeindex =~ /^\d+$/;
	print "Feature::setMatcher    start index not numeric\n" and exit if not $startindex =~ /^\d+$/;
	print "Feature::setMatcher    stop index not numeric\n" and exit if not $stopindex =~ /^\d+$/;
	print "Feature::setMatcher    chromosome index == start index: $chromosomeindex\n" and exit if $chromosomeindex == $startindex;
	print "Feature::setMatcher    chromosome index == stop index: $chromosomeindex\n" and exit if $chromosomeindex == $stopindex;


	#### SET INDEXES
	$self->{_chromosomeindex} = $chromosomeindex;
	$self->{_startindex} = $startindex;
	$self->{_stopindex} = $stopindex;


	#### SET MAX INDEX
	my $maxindex = $chromosomeindex;
	$maxindex = $startindex if $startindex > $maxindex;
	$maxindex = $stopindex if $stopindex > $maxindex;

	$self->{_maxindex} = $maxindex;
}

=head2

	SUBROUTINE		setFlanks

	PURPOSE

		LEFT AND RIGHT FLANKS FOR start AND stop, RESPECTIVELY 

=cut

sub setFlanks
{
	my $self		=	shift;


	my $flanks = $self->get_flanks();

	#### flanks OVERRIDES leftflank AND rightflank
	if ( defined $flanks )
	{
		$self->set_leftflank($flanks);
		$self->set_rightflank($flanks);
	}

	#### SET TO ZERO IF NOT DEFINED
	$self->set_leftflank(0) if not defined $self->get_leftflank();
	$self->set_rightflank(0) if not defined $self->get_rightflank();
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

	print "Feature::new    file not defined\n" and exit if not defined $self->get_file();
	print "Feature::new    Neither regex type nor regex indexes are defined\n" and exit if not defined $self->get_type() and not defined $self->get_indexes();

	#### SET FLANKS IF DEFINED
	$self->setFlanks();

	#### SET INDEXES IF DEFINED
	$self->setIndexes();

	#### SET MATCHER SUBROUTINE
	$self->setMatcher();

	#### SET INPUT AND OUTPUT FILES
	$self->setFiles();

	#### SET FIRST FEATURE
	$self->nextFeature();

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

    #return if not defined $newvalue or not $newvalue;

    my ($operation, $attribute) = ($AUTOLOAD =~ /(get|set)(_\w+)$/);

    # Is this a legal method name?
    unless($operation && $attribute) {
        croak "Method name $AUTOLOAD is not in the recognized form (get|set)_attribute\n";
    }
    unless( exists $self->{$attribute} or $self->is_valid($attribute, $DATAHASH) )
	{
        croak "No such attribute '$attribute' exists in the class ", ref($self);
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
