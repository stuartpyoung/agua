package SolexaUtil;


=head2

		PACKAGE		SolexaUtil

		PURPOSE

			A COLLECTION OF ROUTINES FOR PARSING SOLEXA OUTPUT FILES FOR INPUT

			INTO DIFFERENT ASSEMBLY APPLICATIONS AND OTHER ANCILLARY TASKS

			RELATED TO SOLEXA DATA ANALYSIS

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
use Data::Dumper;
#use overload '==' => 'identical';
#use overload 'eq' => 'equal';


#### SET SLOTS
our @DATA = qw(
	ELAND
);
our $DATAHASH;
foreach my $key ( @DATA )	{	$DATAHASH->{lc($key)} = 1;	}

sub avqualfile
{
	my $self	=	shift;
	my $file	=	shift;
	my $min		=	shift;
	my $max		=	shift;
	my $skip	=	shift;
	my $length	=	shift;

	$file =~ s/\.gz$//;
	$file =~ s/\.zip$//;
	$file =~ s/\.fastq$//;
	$file .= ".min$min" if defined $min;
	$file .= ".max$max" if defined $max;
	$file .= ".skip$skip" if defined $skip;
	$file .= ".length$length" if defined $length;
	$file .= ".avqual";

	return $file;
}


=head2



=cut

sub build_directories
{
	my $self			=	shift;
	my $makefile_out	=	shift;	#### goat_pipeline.pl MAKEFILE GENERATION STDOUT

	#### DIRECTORIES TO BE FOUND
	my ($firecrest_directory, $bustard_directory, $gerald_directory);

	#### EXTRACT Firecrest AND Bustard DIRECTORIES FROM MAKEFILE PATH
	my $contents = Util::contents($makefile_out);

	#### PARSE THROUGH LINES TO FIND DIRECTORY INFORMATION LINE
	my @stdout_lines = split "\n", $contents;
	foreach my $line ( @stdout_lines )
	{
		#### MATCH THIS LINE:

		#### [/nethome/syoung/base/pipeline/solexa-reruns/run3/080814_HWI-EAS185_0001_SeqCapture_Barcoding_RNA_JH/Data/C1-42_Firecrest1.3.2_17-03-2009_syoung/Bustard1.3.2_17-03-2009_syoung/GERALD_17-03-2009_syoung]
		#### Analysis folder: /mihg/data/solexa_xfer/pipeline_out/090130_HWI-EAS185_0014_20DNUAAXX_mito1-4_DH_Jia_Uz_ENJH/Data/C1-72_Firecrest1.3.2_19-03-2009_syoung
		#### Sequence folder: /mihg/data/solexa_xfer/pipeline_out/090130_HWI-EAS185_0014_20DNUAAXX_mito1-4_DH_Jia_Uz_ENJH/Data/C1-72_Firecrest1.3.2_19-03-2009_syoung/Bustard1.3.2_19-03-2009_syoung


		#### GA Pipeline v1.4.0 FORMAT:
		#### 2 [/mihg/data/NGS/syoung/base/pipeline/run15/1.3.2/090528_HWI-EAS185_0007_42C1TAAXX_Flemington_Duan/Data/C1-108_Firecrest1.4.0_07-07-2009_syoung.3/Bustard1.4.0_07-07-2009_syoung/GERALD_07-07-2009_syoung]
		#### OR FROM THIS DEBUG OUTPUT:
		####
		#### GERALD::configureGerald     outdir: /nethome/syoung/base/pipeline/profile/image2eland2/090515_HWI-EAS185_0005_4299MAAXX_EFlemington_1-4/Data/C1-52_Firecrest1.4.0_13-08-2009_syoung/Bustard1.4.0_13-08-2009_syoung/GERALD_13-08-2009_syoung


		if ( $line =~ /(\/\S+Firecrest[^\/]+)\/(Bustard[^\/]+)\/(GERALD[^\/]+)/ )
		{
			$firecrest_directory = $1;
			$bustard_directory = $2;
			$gerald_directory = $3;
			last;
		}

		elsif ( $line =~ /^\[(\S+Firecrest[^\/]+)\/(Bustard[^\/]+)\/(GERALD[^\/]+)\]\s*$/ )
		{
			$firecrest_directory = $1;
			$bustard_directory = $2;
			$gerald_directory = $3;

			$firecrest_directory =~ s/^\]//;

			last;
		}
		elsif ( $line =~ /[\[]*(\S+Firecrest[^\/]+)\/(Bustard[^\/]+)\/(GERALD[^\/]+)\]\s*$/ )
		{
			$firecrest_directory = $1;
			$bustard_directory = $2;
			$gerald_directory = $3;

			$firecrest_directory =~ s/^\]//;

			last;
		}
	}

	return ($firecrest_directory, $bustard_directory, $gerald_directory);
}









sub fastq2fasta
{

	my $self		=	shift;
	my $record		=	shift;


	if ( not defined $record or not $record )
	{
		return;
	}

	my ($sequence_header, $sequence, $quality_header, $quality) = $record =~ /^([^\n]+)\n([^\n]+)\n([^\n]+)\n([^\n]+)$/ms; 
	if ( not defined $sequence_header or not defined $sequence )
	{
		print "FASTA not defined in record: $record\n";
		return;
	}

    $sequence_header =~ s/^@/>/;
	$sequence_header =~ s/^(\S+)[^\n]*$/$1/;
    $quality_header =~ s/^\+/>/;
	$quality_header =~ s/^(\S+)[^\n]*$/$1/;

	$quality_header = $sequence_header if $quality_header =~ /^>\s*$/;
    #$sequence =~ s/\s+$//ms;
    #$quality =~ s/\s+$//ms;


	return "$sequence_header\n$sequence", "$quality_header\n$quality";
}


sub sequence_read
{
	my $self		=	shift;
	my $record		=	shift;

	if ( not defined $record or not $record )
	{
		return;
	}

	my ($fasta, $quality) = $record =~ /^(.+?)\n\+(.+)$/ms; 
	if ( not defined $fasta )
	{
		print "FASTA not defined in record: $record\n";
		return;
	}

    $fasta =~ s/\s+$//ms;
    $quality =~ s/\s+$//ms;

	my $header;
    my $sequence;
    my $symbolic_quality;
    ($header, $sequence) = split "\n", $fasta;
    ($header, $symbolic_quality) = split "\n", $quality;
    if ( not defined $symbolic_quality )
    {
        print "Symbolic quality not defined for record: $record\n";
		print "fasta: $fasta\n";
		print "quality: $quality\n";
        return;
    }

    if ( not defined $sequence )
    {
        print "Sequence not defined for record: $record\n";
		print "fasta: $fasta\n";
		print "quality: $quality\n";
        return;
    }

	my $sequence_read;
	$sequence_read->{header} = $header;
	$sequence_read->{sequence} = $sequence;
	$sequence_read->{quality} = $symbolic_quality;

	return $sequence_read;
}



sub split_fasta
{
	my $self		=	shift;
	my $record		=	shift;
	print "Record: $record\n";
	if ( not defined $record or not $record )
	{
		return;
	}

	my ($fasta, $symbolic_quality) = $record =~ /^(.+)\+(.+)$/ms; 
    $fasta =~ s/\s+$//ms;
    $symbolic_quality =~ s/\s+$//ms;

	#### REMOVE FIRST 'G' RESIDUE AND ITS QUALITY
    #@SLXA-EAS1_89:1:1:672:654/1
    #GCTACGGAATAAAACCAGGAACAACAGACCCAGCA
    #+SLXA-EAS1_89:1:1:672:654/1
    #cccccccccccccccccccc]c``cVcZccbSYbY

	my $trash;
    ($trash, $fasta) = split "\n", $fasta;
    ($trash, $symbolic_quality) = split "\n", $symbolic_quality;
    $fasta =~ s/^.//;
    $symbolic_quality =~ s/^.//;

	return $fasta, $symbolic_quality;
}


sub symbolic2numeric
{
	my $self		=	shift;
	my $quality		=	shift;
	my $type 		=	shift;



	#### CHECK TYPE
	print "Type not supported: $type" if defined $type and not $type =~ /^(solexa|sanger)$/;

	#### SET DEFAULT TYPE
	$type = "solexa" if not defined $type;

	my $offset = 64;
	$offset = 33 if $type eq "sanger";

	my $ascii = $self->{_ascii};

    #### REMOVE 'G' RESIDUES AT EITHER END

    my @symbols = split "", $quality;

    my $numeric_quality = '';
    foreach my $symbol ( @symbols )
    {
		if ( not defined $ascii->{$symbol} )
		{
			print "ascii not defined for $symbol. Exiting.\n" and exit;
		}

		my $number = $ascii->{$symbol} - $offset;
		if ( $number < 0 )
		{
			$number = 0;
		}
        $numeric_quality .= "$number ";
    }
	$numeric_quality =~ s/\s+$//;


	return $numeric_quality;	
}


sub numeric2symbolic
{
	my $self		=	shift;
	my $quality		=	shift;

	my $ascii = $self->_ascii();
	%$ascii = reverse(%$ascii);
	#exit;

    #### REMOVE 'G' RESIDUES AT EITHER END

    my @numbers = split " ", $quality;
    my $symbolic_quality = '';
    foreach my $number ( @numbers )
    {
		if ( not defined $ascii->{$number} )
		{
			print "ascii not defined for $number\n";
		}

		$number += 64;

		my $symbol;
		if ( $number == 35 )	{	$symbol = "\#";	}
		elsif ( $number == 36 )	{	$symbol = "\$";	}
		elsif ( $number == 37 )	{	$symbol = "\%";	}
		elsif ( $number == 64 )	{	$symbol = "\@";	}
		else
		{
			$symbol = $ascii->{$number};
		}
        $symbolic_quality .= $symbol;
    }


	return $symbolic_quality;	
}

sub _reverse_ascii
{
	my $self		=	shift;
	my $key			=	shift;
	my $value		=	shift;

    my $ascii;    

	if ( not defined $key )
	{
	    $key = "dec";
	}
	if ( not defined $value )
	{
		$value = "symbol";
	}

	return $self->_ascii($key, $value);
}


sub _ascii
{
	my $self		=	shift;
	my $key			=	shift;
	my $value		=	shift;

    my $ascii;    

	if ( not defined $key )
	{
		$key = "symbol";
	}
	if ( not defined $value )
	{
	    $value = "dec";
	}

    my $data = $self->data();
    my @lines = split "\n", $data;

    #my @headings = split "\t", <DATA>;
    #chop $headings[$#headings];

    my @headings = split "\t", $lines[0];

	my %headings_hash;
    my $counter = 0;
    foreach my $heading ( @headings )
    {
        $headings_hash{$heading} = $counter;
        $counter++;
    }

    my $value_number = $headings_hash{$value};
    my $key_number = $headings_hash{$key};
    if ( not defined $value_number )
    {
        die "sub ascii(). Value $value number not defined in headings: @headings\n";
    }
    if ( not defined $key_number )
    {
        die "sub ascii(). Key $key number not defined in headings: @headings\n";
    }

#    while ( <DATA>  )
#    {
#	    $_ =~ s/\s+$//;
#		my @elements = split "\t", $_;

    for ( my $i = 1; $i < $#lines + 1; $i++ )
    {
        my $line = $lines[$i];
 		my @elements = split "\t", $line;

		if ( @elements )
		{
	        $ascii->{$elements[$key_number]} = $elements[$value_number];
		}
    }

    return $ascii;
}


#=head2
#
#	SUBROUTINE		identical
#	
#	PURPOSE
#	
#		OVERLOAD THE '==' IDENTITY OPERATOR
#
#	NOTES
#	
#		The module overload.pm is a standard part of the Perl distribution. It allows
#		your objects to define how they will react to a number of Perl's operators.
#		
#		For example, we can add code like this to Number::Fraction:
#	
#		use overload '+' => 'add';
#			
#		Whenever a Number::Fraction is used as one of the operands to the + operator,
#		the add method will be called instead. Code like:
#		
#		$three_quarters = $half + '3/4';
#		
#		is converted to:
#		
#		$three_quarters = $half->add('3/4');
#
#=cut
#
#
#sub identical
#{
#	my $self		=	shift;
#	my $other		=	shift;
#
#	my $identical = 1;
#	
#	foreach my $key (keys %$self)
#	{
#		if ( $key =~ /pid$/ )	{	next;	}
#		
#		if ( ref($self->{$key}) eq '' )
#		{
#			#print "Self->{$key} = $self->{$key} VS $other->{$key} = Other->{$key}\n";
#			if ( $self->{$key} ne $other->{$key} )
#			{
#				$identical = 0;
#			}
#		}
#
#		#print "Identical: $identical\n";
#	}
#	
#	return $identical;
#}
#
#
#=head2
#
#	SUBROUTINE		equal
#	
#	PURPOSE
#	
#		OVERLOAD THE 'eq' STRING EQUALITY OPERATOR
#
#	NOTES
#	
#		The module overload.pm is a standard part of the Perl distribution. It allows
#		your objects to define how they will react to a number of Perl's operators.
#		
#		For example, we can add code like this to Number::Fraction:
#	
#		use overload '+' => 'add';
#			
#		Whenever a Number::Fraction is used as one of the operands to the + operator,
#		the add method will be called instead. Code like:
#		
#		$three_quarters = $half + '3/4';
#		
#		is converted to:
#		
#		$three_quarters = $half->add('3/4');
#
#=cut
#
#
#sub equal
#{
#
#
#	my $self		=	shift;
#	my $other		=	shift;
#
#	my $equal = 1;
#	
#	if ( $self->{_application} ne $other->{_application} )
#	{
#		$equal = 0;
#	}
#
#	return $equal;
#}
#
#



################################################################################
################################################################################
##################			HOUSEKEEPING SUBROUTINES			################
################################################################################
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
    if ( defined $arguments )
    {
    	$self->initialise($arguments);
    }

	#### GET ASCII HASH VALUES
	my $ascii = $self->_ascii();
	$self->{_ascii} = $ascii;

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

	##### 1. LOAD DATABASE, USER AND PASSWORD ARGUMENTS INTO self IN
	##### ORDER TO FILL OUT THESE %VARIABLES% IN XML FILE ENTRIES
	#if (defined $arguments->{DATABASE})	{	$self->{_database} = $arguments->{DATABASE};	}
	#if (defined $arguments->{USER})	{	$self->{_user} = $arguments->{USER};	}
	#if (defined $arguments->{PASSWORD})	{	$self->{_password} = $arguments->{PASSWORD};	}
	#if (defined $arguments->{INSTANCE})	{	$self->{_instance} = $arguments->{INSTANCE};	}
	#
	##### 2. LOAD THE XML FILE: FILL OUT %VARIABLES% IN XML AND LOAD XML
	#if ( defined $arguments->{XMLFILE} )
	#{
	#	$self->{_xmlfile} = $arguments->{XMLFILE};
	#	$self->load();
	#}
	#else
	#{
	#	#die "No XML file in arguments\n";	
	#}

	#exit;	

	#### 3. LOAD THE USER-PROVIDED ARGUMENTS
	foreach my $key ( keys %$arguments )
	{
		#### FILL OUT ANY %VARIABLES% IN PATHS
		$arguments->{$key} =	$self->path($arguments->{$key});

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

    # my($self) = @_;

}



sub data
{
    my $data = qq{dec	octal	hex	binary	symbol
0	0	0	0	NUL  (Null char.)
1	1	1	1	SOH  (Start of Header)
2	2	2	10	STX  (Start of Text)
3	3	3	11	ETX  (End of Text)
4	4	4	100	EOT  (End of Transmiss	ion)
5	5	5	101	ENQ  (Enquiry)
6	6	6	110	ACK  (Acknowledgment)
7	7	7	111	BEL  (Bell)
8	10	8	1000	BS  (Backspace)
9	11	9	1001	HT  (Horizontal Tab)
10	12	00A	1010	LF  (Line Feed)
11	13	00B	1011	VT  (Vertical Tab)
12	14	00C	1100	FF  (Form Feed)
13	15	00D	1101	CR  (Carriage Return)
14	16	00E	1110	SO  (Shift Out)
15	17	00F	1111	SI  (Shift In)
16	20	10	10000	DLE  (Data Link Escape	)
17	21	11	10001	DC1  (XON) (Device Control 1)
18	22	12	10010	DC2      (Device Control 2)
19	23	13	10011	DC3  (XOFF)(Device Control 3)
20	24	14	10100	DC4      (Device Control 4)
21	25	15	10101	NAK  (Negative Acknowledgement)
22	26	16	10110	SYN  (Synchronous Idle)
23	27	17	10111	ETB  (End of Trans. Block)
24	30	18	11000	CAN  (Cancel)
25	31	19	11001	EM  (End of Medium)
26	32	01A	11010	SUB  (Substitute)
27	33	01B	11011	ESC  (Escape)
28	34	01C	11100	FS  (File Separator)
29	35	01D	11101	GS  (Group Separator)
30	36	01E	11110	RS  (Req to Send)(RecSep)
31	37	01F	11111	US  (Unit Separator)
32	40	20	100000	SP  (Space)
33	41	21	100001	!
34	42	22	00100010	"
35	43	23	100011	#
36	44	24	100100	\$
37	45	25	100101	\%
38	46	26	100110	&
39	47	27	100111	'
40	50	28	101000	(
41	51	29	101001	)
42	52	02A	101010	*
43	53	02B	101011	+
44	54	02C	00101100	,
45	55	02D	101101	-
46	56	02E	101110	.
47	57	02F	101111	/
48	60	30	110000	0
49	61	31	110001	1
50	62	32	110010	2
51	63	33	110011	3
52	64	34	110100	4
53	65	35	110101	5
54	66	36	110110	6
55	67	37	110111	7
56	70	38	111000	8
57	71	39	111001	9
58	72	03A	111010	:
59	73	03B	111011	;
60	74	03C	111100	<
61	75	03D	111101	=
62	76	03E	111110	>
63	77	03F	111111	?
64	100	40	1000000	\@
65	101	41	1000001	A
66	102	42	1000010	B
67	103	43	1000011	C
68	104	44	1000100	D
69	105	45	1000101	E
70	106	46	1000110	F
71	107	47	1000111	G
72	110	48	1001000	H
73	111	49	1001001	I
74	112	04A	1001010	J
75	113	04B	1001011	K
76	114	04C	1001100	L
77	115	04D	1001101	M
78	116	04E	1001110	N
79	117	04F	1001111	O
80	120	50	1010000	P
81	121	51	1010001	Q
82	122	52	1010010	R
83	123	53	1010011	S
84	124	54	1010100	T
85	125	55	1010101	U
86	126	56	1010110	V
87	127	57	1010111	W
88	130	58	1011000	X
89	131	59	1011001	Y
90	132	05A	1011010	Z
91	133	05B	1011011	[
92	134	05C	1011100	\\
93	135	05D	1011101	]
94	136	05E	1011110	^
95	137	05F	1011111	_
96	140	60	1100000	`
97	141	61	1100001	a
98	142	62	1100010	b
99	143	63	1100011	c
100	144	64	1100100	d
101	145	65	1100101	e
102	146	66	1100110	f
103	147	67	1100111	g
104	150	68	1101000	h
105	151	69	1101001	i
106	152	06A	1101010	j
107	153	06B	1101011	k
108	154	06C	1101100	l
109	155	06D	1101101	m
110	156	06E	1101110	n
111	157	06F	1101111	o
112	160	70	1110000	p
113	161	71	1110001	q
114	162	72	1110010	r
115	163	73	1110011	s
116	164	74	1110100	t
117	165	75	1110101	u
118	166	76	1110110	v
119	167	77	1110111	w
120	170	78	1111000	x
121	171	79	1111001	y
122	172	07A	1111010	z
123	173	07B	1111011	{
124	174	07C	1111100	|
125	175	07D	1111101	}
126	176	07E	1111110	~
127	177	07F	1111111	DEL};

    return $data;
}



1;

__DATA__
dec	octal	hex	binary	symbol
0	0	0	0	NUL  (Null char.)
1	1	1	1	SOH  (Start of Header)
2	2	2	10	STX  (Start of Text)
3	3	3	11	ETX  (End of Text)
4	4	4	100	EOT  (End of Transmiss	ion)
5	5	5	101	ENQ  (Enquiry)
6	6	6	110	ACK  (Acknowledgment)
7	7	7	111	BEL  (Bell)
8	10	8	1000	BS  (Backspace)
9	11	9	1001	HT  (Horizontal Tab)
10	12	00A	1010	LF  (Line Feed)
11	13	00B	1011	VT  (Vertical Tab)
12	14	00C	1100	FF  (Form Feed)
13	15	00D	1101	CR  (Carriage Return)
14	16	00E	1110	SO  (Shift Out)
15	17	00F	1111	SI  (Shift In)
16	20	10	10000	DLE  (Data Link Escape	)
17	21	11	10001	DC1  (XON) (Device Con	trol 1)
18	22	12	10010	DC2      (Device Contr	ol 2)
19	23	13	10011	DC3  (XOFF)(Device Con	trol 3)
20	24	14	10100	DC4      (Device Contr	ol 4)
21	25	15	10101	NAK  (Negative Acknowl	edgement)
22	26	16	10110	SYN  (Synchronous Idle	)
23	27	17	10111	ETB  (End of Trans. Bl	ock)
24	30	18	11000	CAN  (Cancel)
25	31	19	11001	EM  (End of Medium)
26	32	01A	11010	SUB  (Substitute)
27	33	01B	11011	ESC  (Escape)
28	34	01C	11100	FS  (File Separator)
29	35	01D	11101	GS  (Group Separator)
30	36	01E	11110	RS  (Req to Send)(Rec	Sep)
31	37	01F	11111	US  (Unit Separator)
32	40	20	100000	SP  (Space)
33	41	21	100001	!
34	42	22	00100010	"
35	43	23	100011	#
36	44	24	100100	$
37	45	25	100101	%
38	46	26	100110	&
39	47	27	100111	'
40	50	28	101000	(
41	51	29	101001	)
42	52	02A	101010	*
43	53	02B	101011	+
44	54	02C	00101100	,
45	55	02D	101101	-
46	56	02E	101110	.
47	57	02F	101111	/
48	60	30	110000	0
49	61	31	110001	1
50	62	32	110010	2
51	63	33	110011	3
52	64	34	110100	4
53	65	35	110101	5
54	66	36	110110	6
55	67	37	110111	7
56	70	38	111000	8
57	71	39	111001	9
58	72	03A	111010	:
59	73	03B	111011	;
60	74	03C	111100	<
61	75	03D	111101	=
62	76	03E	111110	>
63	77	03F	111111	?
64	100	40	1000000	@
65	101	41	1000001	A
66	102	42	1000010	B
67	103	43	1000011	C
68	104	44	1000100	D
69	105	45	1000101	E
70	106	46	1000110	F
71	107	47	1000111	G
72	110	48	1001000	H
73	111	49	1001001	I
74	112	04A	1001010	J
75	113	04B	1001011	K
76	114	04C	1001100	L
77	115	04D	1001101	M
78	116	04E	1001110	N
79	117	04F	1001111	O
80	120	50	1010000	P
81	121	51	1010001	Q
82	122	52	1010010	R
83	123	53	1010011	S
84	124	54	1010100	T
85	125	55	1010101	U
86	126	56	1010110	V
87	127	57	1010111	W
88	130	58	1011000	X
89	131	59	1011001	Y
90	132	05A	1011010	Z
91	133	05B	1011011	[
92	134	05C	1011100	\
93	135	05D	1011101	]
94	136	05E	1011110	^
95	137	05F	1011111	_
96	140	60	1100000	`
97	141	61	1100001	a
98	142	62	1100010	b
99	143	63	1100011	c
100	144	64	1100100	d
101	145	65	1100101	e
102	146	66	1100110	f
103	147	67	1100111	g
104	150	68	1101000	h
105	151	69	1101001	i
106	152	06A	1101010	j
107	153	06B	1101011	k
108	154	06C	1101100	l
109	155	06D	1101101	m
110	156	06E	1101110	n
111	157	06F	1101111	o
112	160	70	1110000	p
113	161	71	1110001	q
114	162	72	1110010	r
115	163	73	1110011	s
116	164	74	1110100	t
117	165	75	1110101	u
118	166	76	1110110	v
119	167	77	1110111	w
120	170	78	1111000	x
121	171	79	1111001	y
122	172	07A	1111010	z
123	173	07B	1111011	{
124	174	07C	1111100	|
125	175	07D	1111101	}
126	176	07E	1111110	~
127	177	07F	1111111	DEL	