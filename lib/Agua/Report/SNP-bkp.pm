package Report::SNP;

#### DEBUG
$DEBUG = 1;

=head2

		PACKAGE		Report

		PURPOSE

			THE Reports OBJECT PERFORMS THE FOLLOWING TASKS:

				1. RETURNS FILE AND DIRECTORY LIST OF A GIVEN PATH AS A

                    dojox.data.FileStore JSON OBJECT TO BE RENDERED USING

                    FilePicker INSIDE dojox.dijit.Dialog


saveReport

cd C:\DATA\base\cgi-bin\Bioptic0.2.5

perl report.cgi "class=Report::SNP&mode=saveReport&sessionId=1228791394.7868.158&username=admin"


filterReport

cd C:\DATA\base\cgi-bin\Bioptic0.2.5

perl report.cgi "class=Report::SNP&mode=filterReport&sessionId=1228791394.7868.158&username=admin"


loadReport



=cut

use strict;

#### USE LIB FOR INHERITANCE
use FindBin qw($Bin);
use lib "$Bin/..";

#### INHERIT FROM App CLASS
require Exporter;
our @ISA = 	qw(Report); 
use Report;
use Common;

our $AUTOLOAD;

#### INTERNAL MODULES
use Util;
use DBase::SQLite;


#### EXTERNAL MODULES
use File::Basename; # ALREADY IN PERL
use File::Copy;     # ALREADY IN PERL
use File::stat;     # ALREADY IN PERL
use Data::Dumper;   # ALREADY IN PERL
use Carp;           # ALREADY IN PERL


use JSON -support_by_pp;
use DBD::SQLite;

#### DOWNLOADED FROM CPAN
use File::Remove;
use File::Copy::Recursive;

#### SET SLOTS
our @DATA = qw(
	USERNAME
	SESSIONID
	JSON
    MODE
    DBOBJECT
	CGI    
    CONF
);
our $DATAHASH;
foreach my $key ( @DATA )	{	$DATAHASH->{lc($key)} = 1;	}

our $ROOT = 'admin';
our $DEFAULT_BYTES = 80;

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

	#### CHECK CONF->FILEROOT IS PRESENT AND NOT EMPTY
	if ( not defined $arguments->{CONF} or not $arguments->{CONF} )
	{
		print "CONF is not defined\n";
        exit;
    }
	else
	{
		$self->{_conf} = $arguments->{CONF};
	}

	#### INITIALISE THE OBJECT'S ELEMENTS
	$self->initialise($arguments);


    return $self;
}


=head2

    SUBROUTINE     filterReport

    PURPOSE

        1. FILTER A LIST OF SNPS BASED ON USER-DEFINED CRITERIA

        2. OUTPUT THE NUMBER OF SNPS PASSING EACH SUCCESSIVE FILTER

        3. GENERATE A SET OF DATABASE TABLES, ONE FOR EACH FILTER LEVEL

        4. FILTERS:

            1. chromosome
            2. variant
            3. depth
            4. sense
            5. exonic
            6. dbSNP

	NOTES

		PRINT JUST THE TABLE HEADERS AND THEN AN ARRAY OF
		LINES TO REDUCE TRANSPORT:
		   outputResult : {
			   headers: [ 'column1', 'column2', ... ],
			   data: [
						   [ 1,2,3,... ],
						   [ ... ],
						   ...
			   ]
		   }

		AT THE CLIENT END, CONVERT THIS INTO THE FOLLOWING
		FORMAT:
		PRINT OUTPUT DATA TO JSON FILE FOR THIS REPORT
		CONVERT output_rows INTO AN ARRAY CONTAINING ARRAYS OF HASHES

		   {
			   identifier :"id",
			   label : "id",
			   items : [ array of hashes ]
		   }


		454HCDiffs.txt FORMAT:

			name    chromosome      ccdsstart       ccdsstop        referencenucleotide     variantnucleotide       depth   variantfrequency        chromosomestart chromosomestop  sense   referencecodon    variantcodon    referenceaa     variantaa       strand  snp     score   strand
			CCDS3.1 chr1    770     770     C       T       3       100%    881093  881093  missense        GCC     GTC     Alanine Valine  -                       -
			CCDS3.1 chr1    780     780     G       A       3       100%    879243  879243  synonymous      CTG     CTA     Leucine Leucine -                       -
			CCDS3.1 chr1    790     790     C       G       3       100%    879233  879233  missense        CTG     GTG     Leucine Valine  -                       -

=cut

sub filterReport
{



    my $self        =   shift;

    use JSON;
    my $jsonObject = JSON->new();

    my $dbobject 	=	$self->{_dbobject};
    my $json 		=	$self->{_json};
    print "Report::SNP::filterReport     Json:\n";
	print Dumper $json;
    print "\n";

exit;


    #### GET ABSOLUTE PATH
    my $username = $json->{username};
    my $fileroot = $self->getFileroot($username);

    my $diffsfile = $json->{fileInput};
    $diffsfile = "$fileroot/$diffsfile";
    if ( $^O =~ /^MSWin32$/ )
    {
        $diffsfile =~ s/\//\\/g;
    }

    my $dbfile = $diffsfile;
    $dbfile =~ s/\.txt$//;
    $dbfile .= ".dbl";    

    if ( not -f $diffsfile )
    {
        print "{ error: 'Report::SNP::filterReport    Could not find diffs file: $diffsfile' }";
        exit 1;
    }

    #### GENERATE DIFFS DATABASE
    my ($random) = $self->diffs2db($diffsfile, $dbfile);

    my $filter_order = ["chromosome", "variant", "depth", "sense", "exonic", "dbsnp"];
    my $filters = $json->{filters};


	#### DELETE EXISTING TABLES
	foreach my $table ( @$filter_order )
	{
	    my $query = "DROP TABLE IF EXISTS $table$random";
		my $result = $dbobject->do($query);
	}

    #### SET UP QUERIES FOR EACH FILTER
    #### "filterOrder": [ "chromosome", "variant", "depth", "sense", "exonic", "dbsnp" ]
    my $queries;

	#### CHROMOSOME
    $queries->{chromosome} = "CREATE TABLE chromosome$random AS SELECT * FROM diffs$random";
    if ( $filters->{chromosomeCheckbox} eq "true" )  {   $queries->{chromosome} .= " WHERE chromosome = '$json->{chromosomeCombo}'"; }

	#### VARIANT FREQUENCY
    $queries->{variant} = "CREATE TABLE variant$random AS SELECT * FROM chromosome$random";
    if ( $filters->{variantCheckbox} eq "true" )
	{
		my $variantInput = $json->{variantInput};
		$variantInput =~ s/%$//;
		$queries->{variant} .= " WHERE variantfrequency >= $variantInput";
	}

	#### READ DEPTH
    $queries->{depth} = "CREATE TABLE depth$random AS SELECT * FROM variant$random";    
    if ( $filters->{depth} )  {   $queries->{depth} .= " WHERE depth >= $json->{depth}";    }

	#### SENSE
    my $sense = lc($json->{sense});
	$queries->{sense} = "CREATE TABLE sense$random AS SELECT * FROM depth$random";
    if ( $filters->{sense} and $filters->{sense} !~ /false/ )  {   $queries->{sense} .= " WHERE sense = '$sense'";    }

	#### EXONIC
    $queries->{exonic} = "CREATE TABLE exonic$random AS SELECT * FROM sense$random";

    #### SINCE THIS IS EXOME DATA, DO NOTHING RIGHT NOW
    #### LATER: ADD EXON/INTRON COLUMN INTO *headers-SNPs.txt FILE
    # if ( defined $json->{exonic} )  {   $queries->{exonic} .= " WHERE chr=$json->{exonic}";    }

	#### dbSNP
    $queries->{dbsnp} = "CREATE TABLE dbsnp$random AS SELECT * FROM exonic$random";

	#### DISABLED RIGHT NOW FOR DEBUGGING OF dbSNP CHECKER
    #if ( $filters->{dbsnp} eq "Only dbSNP" )  {   $queries->{dbsnp} .= " WHERE snp != ''";    }
    #elsif ( $filters->{dbsnp} eq "Only non-dbSNP" )  {   $queries->{dbsnp} .= " WHERE snp == ''";    }



    #### PUT RESULTS HERE, TO BE RETURNED TO CLIENT AS JSON STRING
    my $results;

    my $query = "SELECT COUNT(*) FROM diffs$random";
    my $result = $dbobject->query($query);
    $results->{totalResult} = $result;

    #### DO SUCCESSIVE FILTERS
    my $project = $json->{project};
    my $workflow = $json->{workflow};
    my $report = $json->{report};
    my $outputfile = "$fileroot/$username/$project/$workflow/$report.json";
    if ( $^O =~ /^MSWin32$/ )
    {
        $outputfile =~ s/\//\\/g;
    }

    my $output_rows;
    for ( my $counter = 0; $counter < @$filter_order; $counter++ )
    {
        my $filter = $$filter_order[$counter];

        my $query = $queries->{$filter};
        $dbobject->do($query);

        my $result_query = "SELECT COUNT(*) FROM $filter$random";
        my $result = $dbobject->query($result_query);

        #### GET OUTPUT RESULTS IN LAST FILTER
        if ( $counter == @$filter_order - 1)
        {
            my $query = qq{SELECT * FROM $filter$random ORDER BY chromosome, chromosomestart, ccdsstart};
			$output_rows = $dbobject->querytwoDarray($query);
        }        

		my $resultName = $filter . "Result";
        if ( not exists $filters->{$filter} )
        {
            $results->{$resultName} = '';
        }
        else
        {
            $results->{$resultName} = $result;
        }
    }
	if ( not defined $output_rows )
	{
		print "{}";
		exit(0);
	}

    my $dataJson = $jsonObject->encode($output_rows);

    $dataJson =~ s/"/'/g;
    $results->{outputResult} = $dataJson;




    #open(OUTFILE, ">$outputfile") or die "Can't open output file: $outputfile\n";
    #close(OUTFILE);    

    my $report_json = $jsonObject->encode($results);
    print "$report_json\n";
    exit(0);
}




=head2

    SUBROUTINE     diffs2db

    PURPOSE

		1. CREATE A NEW UNIQUE diffsxxx TABLE WHERE xxx IS A RANDOM NUMBER

		2. LOAD THE DATA INTO THE diffsxxx TABLE

		3. RETURN THE RANDOM NUMBER


	####### NB: COULD DO THIS BY CREATING RANDOMISED DATABASE NAMES
	####### BUT THIS WOULD REQUIRE THE AGUA USER TO HAVE SUFFICIENT
	####### PRIVILEGES TO CREATE AND DESTROY TABLES...

=cut

sub diffs2db
{
    my $self        =   shift;
    my $diffsfile   =   shift;
    my $dbfile      =   shift;



	my $dbtype = $self->{_conf}->{DBTYPE};

	#### SET UP VARIABLES
	my $table;
	my $sqlfile = "sql/diffs.sql";  
	my $ignore = 1;		#### IGNORE THIS NUMBER OF INITIAL LINES	
	my $dbobject;	

	$dbobject = $self->{_dbobject};

	my $random = sprintf "%06d", rand(1000000);



	#### SET RANDOM NUMBER TO BE ADDED AT THE END OF TABLE NAMES
	$table = "diffs" . $random;

	my $counter = 0;
	while ( $dbobject->is_table($table, $sqlfile) )
	{

		$counter++;
		$random = sprintf "%06d", rand(1000000);
		$table = "diffs" . $random;

		last if $counter == 3;
	}

$random = 244079;
$table = "diffs" . $random;



	#### CONVERT DIFFS INPUT FILE PATH TO LINUX FORMAT FOR MYSQL
	$diffsfile =~ s/\\/\//g;

	#### GET SQL FOR TABLE
	open(FILE, $sqlfile) or die "Can't open sqlfile $sqlfile: $!\n";
	my $temp = $/;
	undef $/;
	my $sql = <FILE>;
	close(FILE);
	$/ = $temp;
	$sql =~ s/\s*$/;/g;

	#### INSERT NEW RANDOMISED TABLE NAME
	if ( $sql =~ /^\s*CREATE TABLE IF NOT EXISTS/msi )
	{
		$sql =~ s/^\s*CREATE TABLE IF NOT EXISTS \S+/CREATE TABLE IF NOT EXISTS $table/;
	}
	else
	{
		$sql =~ s/^\s*CREATE TABLE \S+/CREATE TABLE $table/;
	}

	#### CREATE THE TABLE. THIS WILL RETURN 0 IF QUERY IS MALFORMED, 1 OTHERWISE
	print "Report::SNP::diffs2db    create table sql: $sql\n";
	my $success = $dbobject->do($sql);
	print "Report::SNP::diffs2db    create table success: $success\n";
	print "{ error: 'Report::SNP::diffs2db    Could not create table $table using sql query: $sql' }" and exit if not $success;


	my $query = qq{DELETE FROM $table};
	my $delete = $dbobject->do($query);

	#### GET PRELOAD COUNT FOR DEBUGGING
	$query = qq{SELECT COUNT(*) FROM $table};
	my $count = $dbobject->query($query);

	#### LOAD DATA INTO TABLE WITHOUT DUPLICATE LINES
	#### AND SKIPPING THE HEADER LINES
	$query = qq{LOAD DATA LOCAL INFILE '$diffsfile' INTO TABLE $table IGNORE $ignore LINES};

	$success = $dbobject->do($query);

	$count = $dbobject->query("SELECT COUNT(*) FROM $table");

    return $random;
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

    #### SET DEFAULT 'ROOT' USER
    $self->value('root', $ROOT);

    #### SET DEFAULT BYTES (TO BE seeked FROM FILE FOR SAMPLE)
    $self->value('bytes', $DEFAULT_BYTES);

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

        1. THE SET parameter INSIDE THE self OBJECT

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

	print "Report::SNP::AUTOLOAD(self, $newvalue)\n";
	print "New value: $newvalue\n";

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


