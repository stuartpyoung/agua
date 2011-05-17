package BinData;

#### DEBUG

=head2

	PACKAGE		BinData

    VERSION:    0.1

	PURPOSE

        A BinData OBJECT CAN BE USED TO BIN NUMERICAL DATA, FROM EITHER

        A SIMPLE ARRAY OF NUMBERS OR FROM NUMERIC FIELDS IN AN ARRAY

        OF HASHES. IN THE LATTER CASE, EACH bin WILL INCLUDE THE NUMERIC

        COUNT OF ITEMS, PLUS ALL OF THE HASHES THAT BELONG IN THAT BIN.



        YOU CAN USE THIS CLASS IN TWO WAYS:

            1)  INITIALISE IT WITH THE DATA AND NUMBER OF BINS TO GET THE BINNED DATA

                my $binner = BinData( {  'DATA' => $data, 'NUMBER_BINS' = 100;   });

                my $bins = $binner->get_bins();

            2)  INITIALISE IT EMPTY, ADD DATA AND THEN GET THE BINS

                #### data IS AN ARRAY OF NUMBERS
                my $binner = BinData( {  'DATA' => $data, 'NUMBER_BINS' = 100;   });
                $binner->add( $data );  

                #### data IS AN ARRAY OF HASHES WITH A NUMERIC field
                my $binner = BinData(
                    {
                        'DATA' => $data,
                        'NUMBER_BINS' => 100;
                        'FIELD' =>  'evalue'
                    }
                );
                #### ADD SOME DATA
                $binner->add( $data );  #### data IS AN ARRAY

                #### ADD SOME MORE DATA
                $binner->add( 1 );  #### data IS NUMBER

                my $bins = $binner->get_bins();
                my $data_bins = $binner->get_data_bins();

                #### IF YOU INPUT AN ARRAY OF OBJECTS, YOU CAN GET THE
                #### OBJECT BINS:
                my $binner = BinData( {  'DATA' => $objects, 'NUMBER_BINS' = 100;   });
                my $object_bins = $binner->get_object_bins();

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
use POSIX;
use Data::Dumper;

#### DEFAULT PARAMETERS
our @DATA = qw(
	DATA
    MIN
    MAX
    BINS
    NUMBER_BINS
    BIN_SIZE
    FIELD
    STEP
);
our $DATAHASH;
foreach my $key ( @DATA )	{	$DATAHASH->{lc($key)} = 1;	}


=head2

	SUBROUTINE		new

	PURPOSE

		CREATE A NEW self OBJECT

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

		INITIALISE THE self OBJECT

=cut

sub initialise
{
    my $self		=	shift;
	my $arguments	=	shift;


    #### VALIDATE USER-PROVIDED ARGUMENTS
	($arguments) = $self->validate_arguments($arguments);	

    #### PROCESS USER-PROVIDED ARGUMENTS
	foreach my $key ( keys %$arguments )
	{
		$self->value($key, $arguments->{$key});
	}

    #### CHECK REQUIRED ARGUMENTS
    if ( not defined $self->{_number_bins}
        and not defined $self->{_bin_size}
        and not defined $self->{_bins} )
    {
        print "Number of bins, bins and bin size (must use one of the three) not defined. Exiting...\n";
        exit;
    }

    #### DO min/max IF NOT PROVIDED BUT data IS PROVIDED
    if ( not defined $self->{_min} or not defined $self->{_max} )
    {
        #if ( not defined $self->{_data} )
        #{
        #    exit;
        #}

        $self->minmax();
    }


    #### CREATE BINS IF data ARGUMENT IS PROVIDED
    $self->_bins();

    #### BIN DATA IF data ARGUMENT IS PROVIDED
    $self->add($self->{_data}) if defined $self->{_data};
}



=head2

    SUBROUTINE      add

    PURPOSE

        ADD ONE DATA ITEM (A NUMBER OR A HASH CONTAINING A 'field' key

        WITH A NUMERIC value) OR AN ARRAY OF DATA ITEMS (NUMBERS OR HASHES)

=cut

sub add
{

    my $self            =   shift;
    my $data            =   shift;

    if (not defined $data )
    {
        print "BinData::add(self, data). Data not defined!\n";
        return;
    }

    my $field   =   $self->{_field};

    #### RECURSE IF INPUT DATA IS AN ARRAY
    if ( ref($data) eq "ARRAY" )
    {
        for ( my $i = 0; $i < @$data; $i++ )
        {
            $self->add($$data[$i]);
        }
    }

    #### OTHERWISE, ADD THE DATA
    else
    {        
        #### BIN USING THE data->{field} VALUE
        if ( defined $field )
        {        

            my $bin_number = $self->bin_number($data, $field);

            #### ADD IT TO A BIN IF ITS WITHIN THE RANGE OF THE bins
            if ( defined $self->{_bins}[$bin_number] )
            {
                my $objects = $self->{_object_bins}[$bin_number];
                my $object = Util::copy_hash($data);
                push @$objects, $object;
                $self->{_object_bins}[$bin_number] = $objects;

                $self->{_data_bins}[$bin_number] = $self->{_data_bins}[$bin_number] + 1;
            }
            else
            {
                return 0;
            }
        }

        #### BIN USING THE data VALUE 
        else
        {        
            my $bin_number = $self->bin_number($data, $field);
            if ( not defined $bin_number )
            {
                return;
            }

            #### ADD IT TO A BIN IF ITS WITHIN THE RANGE OF THE bins
            if ( defined $self->{_bins}[$bin_number] )
            {                                
                $self->{_data_bins}[$bin_number] = $self->{_data_bins}[$bin_number] + 1;
            }
            else
            {
                print "Data '$data' is not within the range of bins ($self->{_min}.. $self->{_max}, step: $self->{_step})\n";
                return 0;
            }
        }  
    }

#exit;



    return 1;

}




=head2

    SUBROUTINE      bin_number

    PURPOSE

        CLEAR THE bins, data_bins AND object_bins SLOTS

=cut

sub bin_number
{
    my $self            =   shift;
    my $data            =   shift;
    my $field           =   shift;



    my $bin_array = $self->{_bins};

    ### PRINT BINS FOR DEBUGGING
    for ( my $i = 0; $i < @$bin_array; $i++ )
    {
    }

    if ( not defined $bin_array )
    {
        print "BinData::bin_number. No bins in bin array. Exiting...\n";
        exit;
    }

    my $min= $self->{_min};
    my $max= $self->{_max};

    if ( defined $field )
    {

        if ( not defined $data->{$field} )
        {
            return;
        }
        elsif ( $data->{$field} < $min or $data->{$field} > $max )
        {
            print "Data value $data->{$field} less than minimum $min\n" if $data->{$field} < $min;    
            print "Data value $data->{$field} greater than maximum $max\n" if $data->{$field} < $max;
            return;
        }
    }
    else
    {
        if ( $data < $min or $data > $max )
        {
            return;
        }
    }

    #### GO THROUGH ALL BINS FROM LOW TO HIGH VALUES. IF THE DATA VALUE 
    #### IS LESS THAN A BIN, RETURN THE INDEX OF THE PREVIOUS BIN
    #### NB: WE KNOW THAT THE DATA VALUE IS GREATER THAN THE MINIMUM,
    #### SO AT ITS LEAST VALUE IT MUST BE LESS THAN THE VALUE OF THE
    #### SECOND BIN.

    for ( my $index = 0; $index < @$bin_array - 1; $index++ )
    {
        my $bin_value = $$bin_array[$index];

        #### MAKE SURE THE BOTTOM BIN IS EQUAL TO MIN
        #if ( $index == 0 )  {   $bin_value = $min;  }

        if ( defined $field )
        {

            if ( $data->{$field} <= $bin_value )
            {
                return $index;
            }
        }
        else
        {

            if ( $data <= $bin_value )
            {
                return $index;
            }
        }
    }

    #if ( defined $field )
    #{
    #}
    #else
    #{
    #}

    return scalar(@$bin_array) - 1;
}




=head2

    SUBROUTINE      _bins

    PURPOSE

        SET THE BINS: bins, data_bins, object_bins

=cut

sub _bins
{   
    my $self            =   shift;

    #### DEFINE min/max
    if ( not defined $self->{_min} or not defined $self->{_max} )
    {
        $self->minmax();
    }

    if ( defined $self->{_bins} )
    {
        my $bins = $self->{_bins};

        # SET DATA BIN COUNTS TO ZERO
        my $data_bins = ();
        for (my $j = 0; $j < @$bins; $j++) {   $$data_bins[$j] = 0; }

        #### INITIALISE SEQUENCE DATA
        my $object_bins;
        for (my $j = 0; $j < @$bins; $j++) {   $$object_bins[$j] = []; }



        $self->{_data_bins} = $data_bins;
        $self->{_object_bins} = $object_bins;

        return;
    }

    #### GET PARAMETERS
	my $field			=	$self->{_field};
    my $min             =   $self->{_min};
    my $max             =   $self->{_max};
    my $number_bins     =   $self->{_number_bins};
    my $bin_size        =   $self->{_bin_size};
    my $step            =   $self->{_step} ? $self->{_step} : 0;

    #### SET range AND bin_size
    my $range = $max - $min + $step;
    $self->{_range} = $range;

    #### SET bin_size IF number_bins DEFINED
    if ( defined $number_bins )
    {
        $bin_size = $range / ( $number_bins );
        $self->{_bin_size} = $bin_size;
    }

    #### OTHERWISE, SET number_bins
    else
    {
        $number_bins = int( $range/$bin_size );
        $self->{_number_bins} = $number_bins;
    }


	#### CREATE BIN LEVELS
    my $bins;
    my $level = $min;
    my $counter = 0;
    for (my $j = 0; $j < $number_bins; $j++) 
    {
        if ( int($level) == $level )
        {
            $$bins[$j] = $level;
        }
        else
        {
            $$bins[$j] = sprintf "%.2f", $level;
        }
        $level += $bin_size;
    }


    # SET DATA BIN COUNTS TO ZERO
    my $data_bins = ();
    for (my $j = 0; $j < $number_bins; $j++) {   $$data_bins[$j] = 0; }

	#### INITIALISE SEQUENCE DATA
	my $object_bins;
	for (my $j = 0; $j < $number_bins; $j++) {   $$object_bins[$j] = []; }



    $self->{_bins} = $bins;
    $self->{_data_bins} = $data_bins;
    $self->{_object_bins} = $object_bins;

}





=head2

	SUBROUTINE		minmax

	PURPOSE

		SET THE MINIMUM AND MAXIMUM VALUES FOR THE BINS


=cut

sub minmax
{
    my $self            =   shift;


    my $data    =   $self->{_data};
    my $field   =   $self->{_field};


    my $max;
    my $min;

    #### IF bins SUPPLIED IN ARGUMENTS, TAKE FIRST AND LAST AS MIN/MAX
    my $bins = $self->{_bins};

    if ( defined $bins )
    {
        $self->{_min}   =   $$bins[0] if not defined $self->{_min};
        $self->{_max}   =   $$bins[@$bins - 1]  if not defined $self->{_max};
        return;
    }

    #### SET INITIAL min/max
    if ( defined $field )
    {
        $min = $$data[0]->{$field};
        $max = $$data[0]->{$field};
    }
    else
    {
        $min = $$data[0];
        $max = $$data[0];
    }

    #### GET MINMAX
    for ( my $i = 0; $i < @$data; $i++ )
    {
        if ( defined $field )
        {
            if ( defined $$data[$i]->{$field} )
            {
                if ( $$data[$i]->{$field} > $max )    {   $max = $$data[$i]->{$field};  }
                if ( $$data[$i]->{$field} < $min )    {   $min = $$data[$i]->{$field};  }
            }
        }
        elsif ( defined $$data[$i] )
        {
            if ( $$data[$i] > $max )    {   $max = $$data[$i];  }
            if ( $$data[$i] < $min )    {   $min = $$data[$i];  }
        }
        else
        {
            if ( defined $field )
            {
                print Dumper $$data[$i];
            }
            else
            {
            }
        }
    }

    if ( not defined $self->{_min} )    {   $self->{_min} = $min;   }
    if ( not defined $self->{_max} )    {   $self->{_max} = $max;   }




}




=head2

    SUBROUTINE      clear

    PURPOSE

        CLEAR THE bins, data_bins AND object_bins SLOTS

=cut

sub clear
{
    my $self            =   shift;

    $self->{_bins} = undef;    
    $self->{_data_bins} = undef;    
    $self->{_object_bins} = undef;    
}




=head2

	SUBROUTINE		value

	PURPOSE

		SET A PARAMETER OF THE self OBJECT TO A GIVEN value

    INPUT

        1. parameter TO BE SET

		2. value TO BE SET TO

    OUTPUT

        1. THE SET parameter INSIDE THE self OBJECT

=cut

sub value
{
    my $self		=	shift;
	my $parameter	=	shift;
	my $value		=	shift;

	$parameter = lc($parameter);

    if ( not defined $value)	{	return;	}
	$self->{"_$parameter"} = $value;
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

	my $hash;
	foreach my $argument ( keys %$arguments )
	{
		if ( $self->is_valid($argument) )
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
}



1;


