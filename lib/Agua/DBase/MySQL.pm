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
use FindBin qw($Bin);
use lib "$Bin/../../";
use Agua::DBase;

class Agua::DBase::MySQL extends Agua::DBase {

has 'dbtype'	=> ( isa => 'Str|Undef', is => 'ro', default => 'mysql' );
has 'dbh'		=> ( isa => 'Any', is => 'rw', default => '' );
has 'database'	=> ( isa => 'Str', is => 'rw', default => '' );
#has 'user'		=> ( isa => 'Str', is => 'rw', default => '' );
#has 'password'	=> ( isa => 'Str', is => 'rw', default => '' );

#### EXTERNAL MODULES
use DBI;
use DBD::mysql;
use POSIX;
use Data::Dumper;

	#////}

method BUILD ($hash) {
	$self->initialise($hash);
}

=head2

	SUBROUTINE		initialise

	PURPOSE

		INITIALISE THE self OBJECT

=cut

method initialise ($hash) {


	#### CONNECT TO DATABASE
	my $database = $self->database();
	$database = "mysql" if not $database;
	my $user = $self->user();
	my $password = $self->password();

	die "Agua::DBase::MySQL::initialise    user not defined\n" if not defined $user;
	die "Agua::DBase::MySQL::initialise    password not defined\n" if not defined $password;

	my $dsn = "DBI:mysql:$database";

    my $dbh = DBI->connect($dsn, $user, $password, { 'PrintError' => 1, 'RaiseError' => 1 });

	$self->dbh($dbh);
}


=head2
GRANT ALL PRIVILEGES ON agua TO $user\@localhost IDENTIFIED BY '$password';	
	SUBROUTINE		querytwoDarray

	PURPOSE

		RETURN AN ARRAY OF HASHES

=cut 

method querytwoDarray ($query) {    


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

    my $response = $self->dbh()->selectall_arrayref($query);
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
	return $twoDarray;
}


=head2

	SUBROUTINE		has_entry

	PURPOSE

		CHECK IF THERE ARE ANY ENTRIES FOR THE GIVEN FIELD IN

		THE TABLE (I.E., MUST BE NON-EMPTY AND NON-NULL). IF NONE

		EXIST, RETURN ZERO. OTHERWISE, RETURN 1.

=cut

method has_entry ($table, $field) {

	if ( not defined $table or not $table )	{	print "Table not defined or empty\n";	return;	}
	if ( not defined $field or not $field )	{	print "Field not defined or empty\n";	return;	}

	my $has_table = has_table($self->dbh(), $table);
	if ( not $has_table )	{	return 0;	}

	my $query = qq{SELECT COUNT(*) FROM $table WHERE $field!='' OR $field IS NOT NULL};
	my $count = simple_query($self->dbh(), $query);
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

method has_table ($table) {

	if ( not defined $table or not $table )	{	print "Table not defined or empty\n";	return;	}

	my $query = qq{SHOW TABLES};
	my $tables = simple_queryarray($self->dbh(), $query);

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

method update_fulltext ($table, $fulltext, $fields) {

	if ( not defined $self->dbh() or not defined $table or not defined $fulltext or not defined $fields )	{	return;	}

    $| = 1;
	if ( has_index($self->dbh(), $table, $fulltext) )
	{
		print "Agua::DBase::MySQL    Dropping index '$fulltext'\n";
		drop_index($self->dbh(), $table, $fulltext);
	}

	my $query = qq{ALTER TABLE $table ADD FULLTEXT $fulltext ($fields)};
	print "Agua::DBase::MySQL    $query\n";
	my $success = do_query($self->dbh(), $query);

	return $success;
}


=head2

	SUBROUTINE		has_index

	PURPOSE

		RETURN 1 IF THE NAME IS AMONG THE LIST OF A TABLE'S INDICES 

=cut 

method has_index ($table, $index) {

	my $indices = indices($self->dbh(), $table);
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

method indices ($table) {

    my $query = qq{SHOW CREATE TABLE $table};
	my $result = simple_query($self->dbh(), $query);

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

method drop_index ($table, $index) {

    my $query = qq{DROP INDEX $index ON $table};
	my $success = do_query($self->dbh(), $query);

	return $success;	
}


=head2

	SUBROUTINE		tsvline2hash

	PURPOSE

		RETURN A HASH FROM A TSV LINE AND ORDERED LIST OF FIELDS

=cut 

method tsvline2hash ($tsvline, $fields) {
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

method do ($query) {

	$self->dbh()->{RaiseError} = 0; # DISABLE TERMINATION ON ERROR
	my $success = $self->dbh()->do($query); ### or warn "Query failed: $DBI::errstr ($DBI::err)\n";
	$self->dbh()->{RaiseError} = 1; # ENABLE TERMINATION ON ERROR

	if ( $success )	{	return 1;	}

	return 0;
}



=head2

	SUBROUTINE		query

	PURPOSE

		RETURN A SINGLE QUERY RESULT

=cut

method query ($query) {

	#use Data::Dumper;

    my $sth = $self->dbh()->prepare($query);
	$sth->execute;

    return $sth->fetchrow;
}


=head2

	SUBROUTINE		queryarray

	PURPOSE

		RETURN AN ARRAY OF SCALARS

=cut 

method queryarray ($query) {

	$self->dbh()->{RaiseError} = 0; # DISABLE TERMINATION ON ERROR	
    my $sth = $self->dbh()->prepare($query);
    $self->dbh()->{RaiseError} = 1; # ENABLE TERMINATION ON ERROR

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

method queryhash ($query) {

	$self->dbh()->{RaiseError} = 0; # DISABLE TERMINATION ON ERROR	
    my $sth = $self->dbh()->prepare($query);
    $self->dbh()->{RaiseError} = 1; # ENABLE TERMINATION ON ERROR
	$sth->execute;

    return $sth->fetchrow_hashref;
}



=head2

	SUBROUTINE		queryhasharray

	PURPOSE

		RETURN AN ARRAY OF HASHES

=cut 

method queryhasharray ($query) {


	$self->dbh()->{RaiseError} = 0; # DISABLE TERMINATION ON ERROR
    my $sth = $self->dbh()->prepare($query);
	$sth->execute;
	$self->dbh()->{RaiseError} = 1; # ENABLE TERMINATION ON ERROR

	my $hasharray;
	while ( my $hash = $sth->fetchrow_hashref )	{	push @$hasharray, $hash;	}

	return $hasharray;
}


=head2

	SUBROUTINE		is_empty

	PURPOSE

		RETURN 1 IF TABLE IS EMPTY

=cut 

method is_empty ($table) {

    my $query = qq{SELECT * FROM $table LIMIT 1};
    my $sth = $self->dbh()->prepare($query);
    $sth->execute;
    my $result = $sth->fetchrow;

    if ( $result )  {   return 0;   };
    return 1;
}

=head2

	SUBROUTINE		fields

	PURPOSE

		RETURN AN ARRAY OF THE TABLES FIELDS OR '' IF NOT DEFINED

=cut 

method fields ($table) {
    my $sth = $self->dbh()->prepare("SELECT * FROM $table LIMIT 1") or die "Can't prepare field names query!\n";
    $sth->execute or die "Can't execute field names query!\n";
    my $fields;
    $fields = $sth->{NAME};

    return $fields;   
}



=head2

	SUBROUTINE		has_field

	PURPOSE

		CHECK IF A TABLE HAS A PARTICULAR FIELD

=cut

method has_field ($table, $field) {

	my $has_table = has_table($self->dbh(), $table);
	if ( not $has_table )	{	return 0;	}

	my $fields = fields_arrayref($self->dbh(), $table);
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


method fields_arrayref ($table) {

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

	SUBROUTINE		fields_insert

	PURPOSE

		CREATE AN INSERT QUERY .TSV LINE (FOR LOAD) BASED ON THE FIELDS

		AND THE HASH, MAKING DEFINED WHEN A FIELD key IS 

		NOT PRESENT IN THE HASH

=cut 

method fields_insert ($fieldstring, $hash, $table) {
	$fieldstring =~ s/\s//g;


	#### SET THE ARRAY OF TABLE FIELDS
	my $fields;
    @$fields = split ",", $fieldstring;

	#### CHECK VALUES ARE THERE
	if ( not defined $table or not $table )
	{
		die "Database::fields_insert(fields, hash, table). Table not defined or empty. Exiting...\n";
	}
	if ( not @$fields)
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
    for ( my $i = 0; $i < @$fields; $i++ )
    {
        $insert .= "\t" . $$fields[$i]. ",\n";        
    }
	$insert .= ")\n";

	#### DO THE FIELD VALUES
    for ( my $i = 0; $i < @$fields; $i++ )
    {
        no warnings;
        print "Agua::DBase::MySQL    Field $$fields[$i]: ", make_defined($$hash{$$fields[$i]}), "\n";

        my $value = $$hash{$$fields[$i]} ? $$hash{$$fields[$i]} : '';
		$value =~ s/'/\\'/g;        
		$value =~ s/\\\\\\'/\\'/g;

        $insert .= "\t$value,\n";        
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

method fields_tsvline ($fields, $hash) {
    my $tsvline = '';

    for ( my $i = 0; $i < @$fields; $i++ )
    {
        my $value = $$hash{$$fields[$i]} ? $$hash{$$fields[$i]} : '';
		$value =~ s/'/\\'/g;        
		$value =~ s/\\\\\\'/\\'/g;
        $tsvline .= "$value\t";
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

method fields_csvline ($fields, $hash) {



    my $csvline = '';
    for ( my $i = 0; $i < @$fields; $i++ )
    {
        no warnings;
        my $value = $$hash{$$fields[$i]} ? $$hash{$$fields[$i]} : '';
		$value =~ s/'/\\'/g;
		$value =~ s/\\\\\\'/\\'/g;

		$value = "'$value'," if $value !~ /^NOW\(\)$/i and $value !~ /^DATETIME\('NOW'\)$/i;
		$value = "$value," if $value =~ /^NOW\(\)$/i or $value =~ /^DATETIME\('NOW'\)$/i;
        use warnings;

        $csvline .= $value;
    }
	$csvline =~ s/,$//;


    return $csvline;
}




=head2

	SUBROUTINE		field_present

	PURPOSE

		RETURN '' IF NOT DEFINED

=cut 

method field_present ($table, $field) {

	my $sth = $self->dbh()->prepare("SELECT * FROM $table") or die "Can't prepare statement!\n";
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

	SUBROUTINE		required_fields

	PURPOSE

		RETURN AN ARRAY OF THE TABLES FIELDS OR '' IF NOT DEFINED


=cut 

method required_fields ($table) {
	print "Agua::DBase::MySQL::required_fields    table not defined or empty.\n"
		and exit if not defined $table or not $table;

	my $structure = $self->get_structure($table);
	my $required_fields = [];
	foreach my $field ( @$structure )
	{
		push @$required_fields, $field->{f_name} if $field->{f_nullable} == 0 
	}

    return $required_fields;   
}

=head2

	SUBROUTINE		get_structure

	PURPOSE

		RETURN A HASH ARRAY OF THE TABLE'S FIELD DETAILS

		[
			{
			  'f_name' => 'username',
			  'f_autoinc' => 0,
			  'f_pri_key' => 1,
			  'f_nullable' => 0,
			  'f_length' => 30,
			  'f_type' => 'varchar'
			},
			...
		]


=cut

method get_structure ($table) {	
	my $query = "SELECT * FROM $table";
	my $sth 	=	$self->dbh()->prepare($query) or die "Can't prepare statement!\n";
	die "Error: " . $self->dbh()->errstr . "\n" unless ($sth);
	die "Error: " . $sth->errstr . "\n" unless ($sth->execute);
	#my $output = $self->query($query);

	my $names         = $sth->{'NAME'};
	my $type          = $sth->{'mysql_type_name'};
	my $length        = $sth->{'mysql_length'};
	my $is_nullable   = $sth->{'NULLABLE'};
	my $is_pri_key    = $sth->{'mysql_is_pri_key'};
	my $is_autoinc    = $sth->{'mysql_is_auto_increment'};

	my $structure;	
	@$structure = map {
	   { 
		  f_name      => $names->[$_],
		  f_type      => $type->[$_],
		  f_length    => $length->[$_], 
		  f_nullable  => $is_nullable->[$_] ? 1 : 0,
		  f_pri_key   => $is_pri_key->[$_]  ? 1 : 0,
		  f_autoinc   => $is_autoinc->[$_]  ? 1 : 0,
	   };
	} (0..$#{$names});

	#use Data::Dumper;

	return $structure;
}

=head2

	SUBROUTINE		load

	PURPOSE

		LOAD DATA FROM .TSV FILE INTO TABLE

=cut 

method load ($table, $file, $extra) {

	#### LOAD DATA INTO TABLE
	my $query = qq{LOAD DATA INFILE '$file' INTO TABLE $table};
    if ( defined $extra )   {   $query .= " $extra";    }

	$self->dbh()->{RaiseError} = 0; # DISABLE TERMINATION ON ERROR
	my $count = $self->dbh()->do($query) or warn "Query failed: $DBI::errstr ($DBI::err)\n";
	if ( not $count ) {   print "Agua::DBase::MySQL    DID NOT LOAD!\n";    }
	$self->dbh()->{RaiseError} = 1; # ENABLE TERMINATION ON ERROR
}



=head2

	SUBROUTINE		load_nodups

	PURPOSE

		LOAD DATA FROM FILE INTO TABLE AND IGNORE DUPLICATES

=cut 

method load_nodups ($table, $file) {

	#### LOAD DATA INTO TABLE
	my $query = qq{LOAD DATA LOCAL INFILE '$file' INTO TABLE $table};
	$self->dbh()->{RaiseError} = 0; # DISABLE TERMINATION ON ERROR
	my $count = $self->dbh()->do($query) or warn "Query failed: $DBI::errstr ($DBI::err)\n";
	if ( not $count ) {   print "Agua::DBase::MySQL    DID NOT LOAD!\n";    }
	$self->dbh()->{RaiseError} = 1; # ENABLE TERMINATION ON ERROR
}


=head2

	SUBROUTINE		create_database

	PURPOSE

        CREATE DATABASE IF NOT EXISTS

=cut 

method create_database ($database) {
	my $query = qq{CREATE DATABASE $database};
	my $result = $self->do($query);

	return 0 if not defined $result;
	return 1;
}
=head2 create database REMOTE CALL

	SUBROUTINE		create_database

	PURPOSE

        CREATE DATABASE IF NOT EXISTS

sub create_database {
    my $self                    =   shift;

    my $database				=	shift;
	my $user					=	shift;
	my $password				=	shift;

    my $drh = DBI->install_driver("mysql");
	my $remotecall = $drh->func('createdb', $database, 'localhost', $user, $password, 'admin');

    return $remotecall;
}

=cut

=head2

	SUBROUTINE		use_database

	PURPOSE

        USE DATABASE IF EXISTS

=cut 

method use_database ($database) {
	my $query = qq{USE $database};
	my $result = $self->do($query);

	return 0 if not defined $result;
	return 1;
}


=head2

	SUBROUTINE		is_database

	PURPOSE

		RETURN 1 IF A DATABASE EXISTS WITH THE GIVEN NAME

		RETURN 0 OTHERWISE.

=cut 

method is_database ($database) {
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

method is_table ($table) {

	my $query = qq{SHOW TABLES};
	my $tables = $self->queryarray($query);
	for my $tablename ( @$tables )	{	return 1 if $tablename eq $table; }

	return 0;
}




method drop_database ($database) {
	my $query = qq{DROP DATABASE $database};
	my $result = $self->do($query);

	return 0 if not defined $result;	
	return 1;
}


method is_reserved_word ($word) {
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

method now {
	return "NOW()";
}


} #### END


