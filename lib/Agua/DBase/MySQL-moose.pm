use MooseX::Declare;


#### DEBUG

=head2

	PACKAGE		Agua::DBase::MySQL

    VERSION:        0.01

    PURPOSE

        1. UTILITY FUNCTIONS TO ACCESS A MYSQL DATABASE

=cut 

use strict;
use warnings;
use Carp;

#### INTERNAL MODULES
use DBase;

#### EXTERNAL MODULES
use DBI;
use DBD::mysql;
use POSIX;
use Data::Dumper;


class Agua::DBase::MySQL extends Agua::DBase {

sub BUILD { 
	my $self		=	shift;

	$self->initialise();
}

=head2

	SUBROUTINE		initialise

	PURPOSE

		INITIALISE THE self OBJECT

=cut

sub initialise {
    my $self		=	shift;
	my $arguments	=	shift;


	#### CONNECT TO DATABASE
	my $database = $self->database();
	$database = "mysql" if not defined $database;
	my $user = $self->user();
	my $password = $self->password();

	die "Agua::DBase::MySQL::initialise    user not defined\n" if not defined $user;
	die "Agua::DBase::MySQL::initialise    password not defined\n" if not defined $password;

	my $dsn = "DBI:mysql:$database";
    my $dbh = DBI->connect( $dsn, $user, $password, {	'PrintError' => 0, 'RaiseError' => 0 } );
	$self->dbh($dbh);
}






=head2

	SUBROUTINE		querytwoDarray

	PURPOSE

		RETURN AN ARRAY OF HASHES

=cut 

sub querytwoDarray {
    my $self        =   shift;
    my $query       =   shift;



    my $dbh		=	$self->dbh();

	my $twoDarray;

	#### GET FIELDS 
    my ($fieldstring) = $query =~ /^SELECT (.+) FROM/i;
    $fieldstring =~ s/\s+//g;

    my $fields;
    if ( $fieldstring =~ /^\*$/ )
    {
        my ($table) = $query =~ /FROM\s+(\S+)/i;
        $fields = $self->fields($table);
    }
    else
    {
        @$fields = split ",", $fieldstring;

    }

    my $response = $dbh->selectall_arrayref($query);
    foreach( @$response )
    {
		my $hash; 
        foreach my $i (0..$#$_)
        {
           $hash->{$fields->[$i]} = $_->[$i];
        }

        my $array = [];
        foreach my $field ( @$fields )
        {
    		push @$array, $hash->{$field};
        }

		push @$twoDarray, $array;
    }

#
#exit;
#
#
	return $twoDarray;
}


=head2

	SUBROUTINE		query

	PURPOSE

		RETURN A SINGLE QUERY RESULT

=cut

sub query {
    my $self        =   shift;
    my $query   	=   shift;

    my $dbh		=	$self->dbh();


    my $sth = $dbh->prepare($query);
	$sth->execute;

    return $sth->fetchrow;
}


=head2

	SUBROUTINE		has_field

	PURPOSE

		CHECK IF A TABLE HAS A PARTICULAR FIELD

=cut

sub has_field {
    my $self        =   shift;
    my $table	=	shift;
	my $field	=	shift;

    my $dbh		=	$self->dbh();	

	my $has_table = has_table($dbh, $table);
	if ( not $has_table )	{	return 0;	}

	my $fields = fields_arrayref($dbh, $table);
	if ( not defined $fields )	{	return 0;	}

	for ( my $i = 0; $i < @$fields; $i++ )
	{
		if ( $$fields[$i] =~ /^$field$/ )	{	return 1;	}
	}

	return 0;
}



=head2

	SUBROUTINE		fields_arrayref

	PURPOSE

		CHECK IF THERE ARE ANY ENTRIES FOR THE GIVEN FIELD IN

		THE TABLE (I.E., MUST BE NON-EMPTY AND NON-NULL). IF NONE

		EXIST, RETURN ZERO. OTHERWISE, RETURN 1.

=cut


sub fields_arrayref {
    my $self        =   shift;
    my $table	=	shift;

    my $dbh		=	$self->dbh();	

	my $query = qq{SHOW CREATE TABLE $table};
	my $result = $self->query($query);

	#### FORMAT:
	# PRIMARY KEY  (`collectionid`,`collectionaccession`,`collectionversion`),
	# FULLTEXT KEY `collectionannotation` (`targetannotation`,`targetsource`,`targetid`)
	#) TYPE=MyISAM |

	my @lines = split "\n", $result;
	my $fields;
	for ( my $i = 0; $i < $#lines + 1; $i++ )
	{
		$lines[$i] =~ s/^\s*\S+\s+KEY\s+.+$//;
		if ( defined $1 )	{	push @$fields, $1;	}
	}

	return $fields;	
}

=head2

	SUBROUTINE		has_entry

	PURPOSE

		CHECK IF THERE ARE ANY ENTRIES FOR THE GIVEN FIELD IN

		THE TABLE (I.E., MUST BE NON-EMPTY AND NON-NULL). IF NONE

		EXIST, RETURN ZERO. OTHERWISE, RETURN 1.

=cut

sub has_entry {
    my $self        =   shift;
	my $table	=	shift;
	my $field	=	shift;

    my $dbh		=	$self->dbh();

	if ( not defined $table or not $table )	{	print "Table not defined or empty\n";	return;	}
	if ( not defined $field or not $field )	{	print "Field not defined or empty\n";	return;	}

	my $has_table = has_table($dbh, $table);
	if ( not $has_table )	{	return 0;	}

	my $query = qq{SELECT COUNT(*) FROM $table WHERE $field!='' OR $field IS NOT NULL};
	my $count = simple_query($dbh, $query);
	if ( not defined $count )	{	$count = 0;	}

	if ( $count == 0 )	{	return 0;	}

	return 1;
}


=head2

	SUBROUTINE		has_table

	PURPOSE

		1. CHECK IF A TABLE EXISTS IN THE DATABASE

		2. RETURN 1 IF EXISTS, RETURN ZERO IF OTHERWISE

=cut

sub has_table {
    my $self        =   shift;
	my $table	=	shift;

    my $dbh		=	$self->dbh();

	if ( not defined $table or not $table )	{	print "Table not defined or empty\n";	return;	}

	my $query = qq{SHOW TABLES};
	my $tables = simple_queryarray($dbh, $query);

	if ( not defined $tables )	{	return 0;	}

	for ( my $i = 0; $i < @$tables; $i++ )
	{
		if ( $table =~ /^$$tables[$i]$/ )	{	return 1;	}
	}

	return 0;
}


=head2

	SUBROUTINE		update_fulltext

	PURPOSE

        UPDATE THE FULLTEXT INDEX OF A TABLE

=cut 

sub update_fulltext {
    my $self        =   shift;    
	my $table		=	shift;
	my $fulltext	=	shift;
	my $fields		=	shift;

    my $dbh		=	$self->dbh();

	if ( not defined $dbh or not defined $table or not defined $fulltext or not defined $fields )	{	return;	}

    $| = 1;
	if ( has_index($dbh, $table, $fulltext) )
	{
		print "Agua::DBase::MySQL    Dropping index '$fulltext'\n";
		drop_index($dbh, $table, $fulltext);
	}

	my $query = qq{ALTER TABLE $table ADD FULLTEXT $fulltext ($fields)};
	print "Agua::DBase::MySQL    $query\n";
	my $success = do_query($dbh, $query);

	return $success;
}


=head2

	SUBROUTINE		has_index

	PURPOSE

		RETURN 1 IF THE NAME IS AMONG THE LIST OF A TABLE'S INDICES 

=cut 

sub has_index {
    my $self        =   shift;
	my $table	=	shift;
	my $index	=	shift;

    my $dbh		=	$self->dbh();

	my $indices = indices($dbh, $table);
	if ( not defined $indices )	{	return 0;	}

	for ( my $i = 0; $i < @$indices; $i++ )
	{
		if ( $$indices[$i] =~ /^$index$/ )	{	return 1;	}
	}

	return 0;
}


=head2

	SUBROUTINE		indices

	PURPOSE

		RETURN A LIST OF A TABLE'S INDICES 

=cut 

sub indices {
    my $self        =   shift;
    my $table	=	shift;

    my $dbh		=	$self->dbh();


    my $query = qq{SHOW CREATE TABLE $table};
	my $result = simple_query($dbh, $query);

	#### FORMAT:
	# PRIMARY KEY  (`collectionid`,`collectionaccession`,`collectionversion`),
	# FULLTEXT KEY `collectionannotation` (`targetannotation`,`targetsource`,`targetid`)
	#) TYPE=MyISAM |

	my @lines = split "\n", $result;
	my $indices;
	for ( my $i = 0; $i < $#lines + 1; $i++ )
	{
		$lines[$i] =~ /^\s*\S+\s+KEY\s+`(\w+)`\s+\(/;
		if ( defined $1 )	{	push @$indices, $1;	}
	}

	return $indices;	
}


=head2

	SUBROUTINE		drop_index

	PURPOSE

		RETURN A HASH FROM A TSV LINE AND ORDERED LIST OF FIELDS

=cut 

sub drop_index {
    my $self        =   shift;
	my $table	=	shift;
	my $index	=	shift;

    my $dbh		=	$self->dbh();

    my $query = qq{DROP INDEX $index ON $table};
	my $success = do_query($dbh, $query);

	return $success;	
}



sub sqlite_dbh {
	my $sqlitefile		=	shift;

	my $dbh = DBI->connect("dbi:SQLite:$sqlitefile", '', '') || die "Can't open database file '$sqlitefile': $@\n";   	

	return $dbh;
}


=head2

	SUBROUTINE		tsvline2hash

	PURPOSE

		RETURN A HASH FROM A TSV LINE AND ORDERED LIST OF FIELDS

=cut 

sub tsvline2hash {
    my $self        =   shift;

	my $tsvline		=	shift;		
	my $fields		=	shift;

	my $hash;
	my @elements = split "\t", $tsvline;
	my @fields_array = split "," , $fields;
	for ( my $i = 0; $i < $#fields_array + 1; $i++ )
	{
		my $element = $elements[$i];
		$$hash{$fields_array[$i]} = $element;
	}

	return $hash;
}


=head2

	SUBROUTINE		query

	PURPOSE

		RETURN AN ARRAY OF SCALARS

=cut 

sub do {
    my $self        =   shift;
	my $query		=	shift;

    my $dbh		=	$self->dbh();

	$dbh->{RaiseError} = 0; # DISABLE TERMINATION ON ERROR
	my $success = $dbh->do($query); ### or warn "Query failed: $DBI::errstr ($DBI::err)\n";
	$dbh->{RaiseError} = 1; # ENABLE TERMINATION ON ERROR

	if ( $success )	{	return 1;	}

	return 0;
}



=head2

	SUBROUTINE		queryarray

	PURPOSE

		RETURN AN ARRAY OF SCALARS

=cut 

sub queryarray {
    my $self        =   shift;
    my $query       =   shift;

    my $dbh		=	$self->dbh();

	$dbh->{RaiseError} = 0; # DISABLE TERMINATION ON ERROR	
    my $sth = $dbh->prepare($query);
    $dbh->{RaiseError} = 1; # ENABLE TERMINATION ON ERROR

	$sth->execute;

	my $array;
    my $result = 0;
	while ( defined $result )
    {
        $result = $sth->fetchrow;
        if ( defined $result )	{	push @$array, $result;	}
    }

	return $array;
}



=head2

	SUBROUTINE		queryhash

	PURPOSE

		RETURN A HASH

=cut 

sub queryhash {
    my $self        =   shift;
    my $query       =   shift;

    my $dbh		=	$self->dbh();

	$dbh->{RaiseError} = 0; # DISABLE TERMINATION ON ERROR	
    my $sth = $dbh->prepare($query);
    $dbh->{RaiseError} = 1; # ENABLE TERMINATION ON ERROR
	$sth->execute;

    return $sth->fetchrow_hashref;
}



=head2

	SUBROUTINE		queryhasharray

	PURPOSE

		RETURN AN ARRAY OF HASHES

=cut 

sub queryhasharray {
    my $self        =   shift;
    my $query       =   shift;

    my $dbh		=	$self->dbh();

	$dbh->{RaiseError} = 0; # DISABLE TERMINATION ON ERROR
    my $sth = $dbh->prepare($query);
	$sth->execute;
	$dbh->{RaiseError} = 1; # ENABLE TERMINATION ON ERROR

	my $hasharray;
	while ( my $hash = $sth->fetchrow_hashref )	{	push @$hasharray, $hash;	}

	return $hasharray;
}


=head2

	SUBROUTINE		is_empty

	PURPOSE

		RETURN 1 IF TABLE IS EMPTY

=cut 

sub is_empty {
    my $self        =   shift;
    my $table       =   shift;

    my $dbh		=	$self->dbh();

    my $query = qq{SELECT * FROM $table LIMIT 1};
    my $sth = $dbh->prepare($query);
    $sth->execute;
    my $result = $sth->fetchrow;

    if ( $result )  {   return 0;   };
    return 1;
}


=head2

	SUBROUTINE		sql2fields

	PURPOSE

		RETURN LIST OF FIELDS FROM AN SQL FILE

=cut 

sub sql2fields {
    my $self        =   shift;
    my $sqlfile         =   shift;

    my $fields;

    open(FILE, $sqlfile) or die "Can't open sql file: $sqlfile\n";	
    $/ = "END OF FILE";
    my $contents = <FILE>;

	$contents =~ s/^.+?\(//ms;    
    $contents =~ s/\(.+?$//ms;
    $contents =~ s/PRIMARY KEY.+$//ims;
    $contents =~ s/KEY.+$//ims;
    $contents =~ s/\W+$//;
	my @lines = split "\n", $contents;
    for ( my $i = 0; $i < $#lines + 1; $i++ )
    {
        my $line = $lines[$i];
        if ( $line =~ /^\s*$/ or $line =~ /^\s+#/ )   {   next;   }
        $line =~ s/^\s*(\S+)\s+.+$/$1/;
        $fields .= $line . ",";
    }
    $fields =~ s/,$//;

	return $fields;
}



=head2

	SUBROUTINE		fields_insert

	PURPOSE

		CREATE AN INSERT QUERY .TSV LINE (FOR LOAD) BASED ON THE FIELDS

		AND THE HASH, MAKING DEFINED WHEN A FIELD key IS 

		NOT PRESENT IN THE HASH

=cut 

sub fields_insert {
    my $self        =   shift;
    my $fields          =   shift;
    my $hash        	=   shift;
	my $table			=	shift;

	$fields =~ s/\s//g;


	#### SET THE ARRAY OF TABLE FIELDS
    my @fieldarray = split ",", $fields;

	#### CHECK VALUES ARE THERE
	if ( not defined $table or not $table )
	{
		die "Database::fields_insert(fields, hash, table). Table not defined or empty. Exiting...\n";
	}
	if ( not @fieldarray)
	{
		die "Database::fields_insert(fields, hash, table). Fields is empty. Exiting...\n";
	}
	if ( not defined $hash )
	{
		die "Database::fields_insert(fields, hash, table). Hash is empty. Exiting...\n";
	}

	#### START insert
    my $insert = 'INSERT INTO $table (\n';

	#### DO THE FIELD NAMES
    for ( my $i = 0; $i < $#fieldarray + 1; $i++ )
    {
        $insert .= "\t" . $fieldarray[$i]. ",\n";        
    }
	$insert .= ")\n";

	#### DO THE FIELD VALUES
    for ( my $i = 0; $i < $#fieldarray + 1; $i++ )
    {
        no warnings;
        print "Agua::DBase::MySQL    Field $fieldarray[$i]: ", make_defined($$hash{$fieldarray[$i]}), "\n";
        $insert .= "\t" . make_defined($$hash{$fieldarray[$i]}) . ",\n";        
        use warnings;        
    }

    return $insert;
}



=head2

	SUBROUTINE		fields_tsvline

	PURPOSE

		CREATE A .TSV LINE (FOR LOAD) BASED ON THE FIELDS

		AND THE HASH, MAKING DEFINED WHEN A FIELD key IS 

		NOT PRESENT IN THE HASH

=cut 

sub fields_tsvline {
    my $self        =   shift;
    my $fields          =   shift;
    my $hash        =   shift;

    my $tsvline = '';

    for ( my $i = 0; $i < @$fields; $i++ )
    {
        $tsvline .= $$hash{$$fields[$i]} ? $$hash{$$fields[$i]} : '';

        $tsvline .= "\t";
    }
    $tsvline .= "\n";


    return $tsvline;
}


=head2

	SUBROUTINE		fields_csvline

	PURPOSE

		CREATE A .CSV LINE FOR TABLE INSERTS USING ''

		WHEN NO KEY VALUE IS PRESENT IN THE INPUT HASH

	INPUTS

		1. COMMA-SEPARATED fields

		2. HASH OF TABLE field KEY-VALUE PAIRS 

=cut 

sub fields_csvline {
    my $self        =   shift;
    my $fields          =   shift;
    my $hash        =   shift;




    my $csvline = '';
    for ( my $i = 0; $i < @$fields; $i++ )
    {
        no warnings;
        my $value = $$hash{$$fields[$i]} ? $$hash{$$fields[$i]} : '';

		$value = "'$value'," if $value !~ /^NOW\(\)$/i and $value !~ /^DATETIME\('NOW'\)$/i;
		$value = "$value," if $value =~ /^NOW\(\)$/i or $value =~ /^DATETIME\('NOW'\)$/i;
        use warnings;

        $csvline .= $value;
    }
	$csvline =~ s/,$//;


    return $csvline;
}




=head2

	SUBROUTINE		set

	PURPOSE

		CREATE A 'SET ...' LINE FOR ONE OR MORE FIELD VALUES

	INPUTS

		1. HASH OF TABLE field KEY-VALUE PAIRS 

		2. ARRAY OF FIELDS TO BE USED

=cut 

sub set {
    my $self        =   shift;
    my $hash        =   shift;
    my $fields          =   shift;

	return if not defined $hash or not defined $fields or not @$fields;

	my $set = '';
    for ( my $i = 0; $i < @$fields; $i++ )
    {
        my $value = $$hash{$$fields[$i]} ? $$hash{$$fields[$i]} : '';
		$set .= qq{$$fields[$i] = '$value'\n} if $i == (scalar(@$fields) - 1);
		$set .= qq{$$fields[$i] = '$value',\n} if $i < (scalar(@$fields) - 1);
    }
	return if not $set;

	$set = "SET " . $set;

    return $set;
}




=head2

	SUBROUTINE		where

	PURPOSE

		CREATE A 'WHERE ...' LINE FOR ONE OR MORE FIELD VALUES

	INPUTS

		1. HASH OF TABLE field KEY-VALUE PAIRS 

		2. ARRAY OF FIELDS TO BE USED

=cut 

sub where {
    my $self        =   shift;
    my $hash        =   shift;
    my $fields          =   shift;

	return if not defined $hash or not defined $fields or not @$fields;

	my $where = '';
    for ( my $i = 0; $i < @$fields; $i++ )
    {
        my $value = $$hash{$$fields[$i]} ? $$hash{$$fields[$i]} : '';
		$where .= qq{WHERE $$fields[$i] = '$value'\n} if $i == 0;
		$where .= qq{AND $$fields[$i] = '$value'\n} if $i != 0;
    }

    return $where;
}




=head2

	SUBROUTINE		notDefined

	PURPOSE

		RETURN A HASH OF UNDEFINED FIELD KEY-PAIRS IN A HASH

	INPUTS

		1. HASH OF TABLE field KEY-VALUE PAIRS 

		2. ARRAY OF FIELDS TO BE USED

=cut 

sub notDefined {

    my $self        =   shift;
    my $hash        =   shift;
    my $fields          =   shift;



	return if not defined $hash or not defined $fields or not @$fields;

	my $notDefined = [];
    for ( my $i = 0; $i < @$fields; $i++ )
    {
        push( @$notDefined, $$fields[$i]) if not defined $$hash{$$fields[$i]};
    }

    return $notDefined;
}






=head2

	SUBROUTINE		make_defined

	PURPOSE

		RETURN '' IF NOT DEFINED

=cut 

sub make_defined {
    my $self        =   shift;
    my $input   =   shift;

    if ( not defined $input )   {   return '';  }

    return $input;
}


=head2

	SUBROUTINE		load

	PURPOSE

		RETURN '' IF NOT DEFINED

=cut 

sub load {
    my $self        =   shift;
	my $table			=	shift;
	my $file			=	shift;
    my $extra           =   shift;

    my $dbh		=	$self->dbh();


	#### LOAD DATA INTO TABLE
	my $query = qq{LOAD DATA INFILE '$file' INTO TABLE $table};
    if ( defined $extra )   {   $query .= " $extra";    }

	$dbh->{RaiseError} = 0; # DISABLE TERMINATION ON ERROR
	my $count = $dbh->do($query) or warn "Query failed: $DBI::errstr ($DBI::err)\n";
	if ( not $count ) {   print "Agua::DBase::MySQL    DID NOT LOAD!\n";    }
	$dbh->{RaiseError} = 1; # ENABLE TERMINATION ON ERROR
}



=head2

	SUBROUTINE		load_nodups

	PURPOSE

		RETURN '' IF NOT DEFINED

=cut 

sub load_nodups {
    my $self        =   shift;
	my $table			=	shift;
	my $file			=	shift;

    my $dbh		=	$self->dbh();

	#### LOAD DATA INTO TABLE
	my $query = qq{LOAD DATA LOCAL INFILE '$file' INTO TABLE $table};
	$dbh->{RaiseError} = 0; # DISABLE TERMINATION ON ERROR
	my $count = $dbh->do($query) or warn "Query failed: $DBI::errstr ($DBI::err)\n";
	if ( not $count ) {   print "Agua::DBase::MySQL    DID NOT LOAD!\n";    }
	$dbh->{RaiseError} = 1; # ENABLE TERMINATION ON ERROR
}


=head2

	SUBROUTINE		field_present

	PURPOSE

		RETURN '' IF NOT DEFINED

=cut 

sub field_present {
    my $self        =   shift;
	my $table			=	shift;
	my $field				=	shift;

    my $dbh		=	$self->dbh();

	my $sth = $dbh->prepare("SELECT * FROM $table") or die "Can't prepare statement!\n";
	$sth->execute or die "Can't execute statement!\n";
	my $field_present = 0;
	my @fields = @{$sth->{NAME}};
	for ( my $i = 0; $i < $#fields + 1; $i++ )
	{
		if ( $fields[$i] =~ /^$field$/ )	{	$field_present = 1;	}	
	}

	return $field_present;
}


=head2

	SUBROUTINE		fields

	PURPOSE

		RETURN AN ARRAY OF THE TABLES FIELDS OR '' IF NOT DEFINED

=cut 

sub fields {
    my $self        =   shift;
	my $table		=	shift;	

    my $dbh	= $self->dbh();
    my $sth = $dbh->prepare("SELECT * FROM $table LIMIT 1") or die "Can't prepare field names query!\n";
    $sth->execute or die "Can't execute field names query!\n";
    my $fields;
    $fields = $sth->{NAME};


    return $fields;   
}


=head2

	SUBROUTINE		create_database

	PURPOSE

        CREATE DATABASE IF NOT EXISTS

=cut 

sub create_database {
    my $self                    =   shift;
    my $database				=	shift;

	my $query = qq{CREATE DATABASE $database};
	my $result = $self->do($query);

	return 0 if not defined $result;
	return 1;
}


=head2

	SUBROUTINE		use_database

	PURPOSE

        USE DATABASE IF EXISTS

=cut 

sub use_database {
    my $self                    =   shift;
    my $database				=	shift;

	my $query = qq{USE $database};
	my $result = $self->do($query);

	return 0 if not defined $result;
	return 1;
}


=head2

	SUBROUTINE		create_table

	PURPOSE

		CREATE A TABLE FROM AN .SQL FILE CONTAINING PROPER 

		SYNTAX, E.G.: "CREATE TABLE IF NOT EXISTS tablename..."

=cut 

sub create_table {
	my $self		=	shift;
	my $sqlfile		=	shift;


	#### GET CONTENTS OF SQL FILE
    open(FILE, $sqlfile) or die "Can't open sqlfile $sqlfile: $!\n";
    my $temp = $/;
    undef $/;
    my $sql = <FILE>;
    close(FILE);
    $/ = $temp;

	$sql =~ s/\n//g;

	#$sql =~ s/\s*;\s*$//g;

	my $dbh = $self->dbh();
    my $result = $dbh->do($sql);

	return $result;
}



=head2

	SUBROUTINE		is_database

	PURPOSE

		RETURN 1 IF A DATABASE EXISTS WITH THE GIVEN NAME

		RETURN 0 OTHERWISE.

=cut 

sub is_database {
	my $self			=	shift;
	my $database 		=	shift;

	my $query = qq{SHOW DATABASES};
	my $databases = $self->queryarray($query);
	print "databases: @$databases\n";
	for my $db ( @$databases )	{	return 1 if $db eq $database; }

	return 0;
}

=head2

	SUBROUTINE		is_table

	PURPOSE

		RETURN 1 IF A TABLE EXISTS WITH THE GIVEN NAME IN THIS DATABASE.

		RETURN 0 OTHERWISE.

=cut

sub is_table {
	my $self			=	shift;
	my $table 		=	shift;


	my $query = qq{SHOW TABLES};
	my $tables = $self->queryarray($query);
	for my $tablename ( @$tables )	{	return 1 if $tablename eq $table; }

	return 0;
}


sub drop_database {
	my $self			=	shift;
	my $database 		=	shift;

	my $query = qq{DROP DATABASE $database};
	my $result = $self->do($query);

	return 0 if not defined $result;	
	return 1;
}


sub is_reserved_word {
	my $self			=	shift;
	my $word 			=	shift;

	my @reserved_words = qw(
		ACCESSIBLE
		ALTER
		AS
		BEFORE
		BINARY
		BY
		CASE
		CHARACTER
		COLUMN
		CONTINUE
		CROSS
		CURRENT_TIMESTAMP
		DATABASE
		DAY_MICROSECOND
		DEC
		DEFAULT
		DESC
		DISTINCT
		DOUBLE
		EACH
		ENCLOSED
		EXIT
		FETCH
		FLOAT8
		FOREIGN
		GRANT
		HIGH_PRIORITY
		HOUR_SECOND
		IN
		INNER
		INSERT
		INT2
		INT8
		INTO
		JOIN
		KILL
		LEFT
		LINEAR
		LOCALTIME
		LONG
		LOOP
		MATCH
		MEDIUMTEXT
		MINUTE_SECOND
		NATURAL
		NULL
		OPTIMIZE
		OR
		OUTER
		PRIMARY
		RANGE
		READ_WRITE
		REGEXP
		REPEAT
		RESTRICT
		RIGHT
		SCHEMAS
		SENSITIVE
		SHOW
		SPECIFIC
		SQLSTATE
		SQL_CALC_FOUND_ROWS
		STARTING
		TERMINATED
		TINYINT
		TRAILING
		UNDO
		UNLOCK
		USAGE
		UTC_DATE
		VALUES
		VARCHARACTER
		WHERE
		WRITE
		ZEROFILL
		ALL
		AND
		ASENSITIVE
		BIGINT
		BOTH
		CASCADE
		CHAR
		COLLATE
		CONSTRAINT
		CREATE
		CURRENT_TIME
		CURSOR
		DAY_HOUR
		DAY_SECOND
		DECLARE
		DELETE
		DETERMINISTIC
		DIV
		DUAL
		ELSEIF
		EXISTS
		FALSE
		FLOAT4
		FORCE
		FULLTEXT
		HAVING
		HOUR_MINUTE
		IGNORE
		INFILE
		INSENSITIVE
		INT1
		INT4
		INTERVAL
		ITERATE
		KEYS
		LEAVE
		LIMIT
		LOAD
		LOCK
		LONGTEXT
		MASTER_SSL_VERIFY_SERVER_CERT
		MEDIUMINT
		MINUTE_MICROSECOND
		MODIFIES
		NO_WRITE_TO_BINLOG
		ON
		OPTIONALLY
		OUT
		PRECISION
		PURGE
		READS
		REFERENCES
		RENAME
		REQUIRE
		REVOKE
		SCHEMA
		SELECT
		SET
		SPATIAL
		SQLEXCEPTION
		SQL_BIG_RESULT
		SSL
		TABLE
		TINYBLOB
		TO
		TRUE
		UNIQUE
		UPDATE
		USING
		UTC_TIMESTAMP
		VARCHAR
		WHEN
		WITH
		YEAR_MONTH
		ADD
		ANALYZE
		ASC
		BETWEEN
		BLOB
		CALL
		CHANGE
		CHECK
		CONDITION
		CONVERT
		CURRENT_DATE
		CURRENT_USER
		DATABASES
		DAY_MINUTE
		DECIMAL
		DELAYED
		DESCRIBE
		DISTINCTROW
		DROP
		ELSE
		ESCAPED
		EXPLAIN
		FLOAT
		FOR
		FROM
		GROUP
		HOUR_MICROSECOND
		IF
		INDEX
		INOUT
		INT
		INT3
		INTEGER
		IS
		KEY
		LEADING
		LIKE
		LINES
		LOCALTIMESTAMP
		LONGBLOB
		LOW_PRIORITY
		MEDIUMBLOB
		MIDDLEINT
		MOD
		NOT
		NUMERIC
		OPTION
		ORDER
		OUTFILE
		PROCEDURE
		READ
		REAL
		RELEASE
		REPLACE
		RETURN
		RLIKE
		SECOND_MICROSECOND
		SEPARATOR
		SMALLINT
		SQL
		SQLWARNING
		SQL_SMALL_RESULT
		STRAIGHT_JOIN
		THEN
		TINYTEXT
		TRIGGER
		UNION
		UNSIGNED
		USE
		UTC_TIME
		VARBINARY
		VARYING
		WHILE
		XOR);

	for my $reserved_word ( @reserved_words )
	{
		return 1 if $reserved_word eq $word;
	}

	return 0;
}


	Agua::DBase::MySQL->meta->make_immutable(inline_constructor => 0)


}

1;




