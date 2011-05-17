#!/usr/bin/perl -w

package Database;

##############################################################################################################
#
#	NAME:	Database.pm
#
#
# 	PURPOSE:	1. UTILITY FUNCTIONS TO ACCESS THE DATABASE, OUTPUT RESULTS, ETC.
#
#
##############################################################################################################

# LOAD MODULES
use strict;
use DBI;
#use DBD::SQLite;
use DBD::mysql;
use Util;
use Data::Dumper;
use warnings;
use Carp;





=head2

	SUBROUTINE		reload_table

	PURPOSE

        1. DROP THE TABLE

        2. LOAD THE DUMPFILE

=cut

sub reload_table
{
    my $dbh         =   shift;
    my $user        =   shift;
    my $password    =   shift;
    my $database    =   shift;
    my $table       =   shift;
    my $dumpfile    =   shift;

    #### COPY DUMPFILE OVER TO LOCAL COPY
    if ( not -f $dumpfile )
    {
        die "Can't find dump file: $dumpfile\n";
    }

    #### LOAD THE orthologues TABLE
    print "Dropping '$table' table... ";
    Database::drop($dbh, 'table', $table);
    print "done.\n";
    print "Loading the '$table' table from dump file... ";
    Database::load_dumpfile($user, $password, $database, $dumpfile);
    print "done.\n";
}


=head2

	SUBROUTINE		dumpfile

	PURPOSE

		CREATE A DUMP FILE FROM A DATABASE TABLE

=cut


sub dumpfile
{
	my $user	    =	shift;
	my $password	=	shift;
    my $database    =   shift;
    my $table       =   shift;
    my $dumpfile    =   shift;

    if ( not defined $database )
    {
        print "++++ Database::dumpfile. Database not defined. Exiting...\n";
        exit;
    }
    if ( not defined $dumpfile )
    {
        print "++++ Database::dumpfile. Dump file not defined. Exiting...\n";
        exit;
    }
    if ( not defined $table )
    {
        print "++++ Database::dumpfile. Table not defined. Dumping whole of database '$database'...\n";
        $table = '';
    }

    my $command = "mysqldump -u $user -p$password $database $table > $dumpfile";
    print "$command\n";
    my $result = `$command`;

    return $result;
}



=head2

	SUBROUTINE		load_dumpfile

	PURPOSE

		LOAD A DUMP FILE INTO A DATABASE

=cut


sub load_dumpfile
{
	my $user	    =	shift;
	my $password	=	shift;
    my $database    =   shift;
    my $dumpfile    =   shift;

    my $command = "mysql -u $user -p$password $database < $dumpfile";
    print "$command\n";
    my $result = `$command`;

    return $result;
}




# GET DATABASE HANDLE
sub has_database
{
    my $dsn         = shift;    # data source name
    my $user		= shift;
    my $password	= shift;

    #### CHECK IF CAN CONNECT TO DATABASE
    my $dbh = DBI->connect( $dsn, $user, $password, {	'PrintError' => 0, 'RaiseError' => 0 } );
	###or warn "Can't connect to database: $dsn\n";

    # RETURN DATABASE HANDLE
    if ( defined $dbh )
    {
        return 1;
    }

    return 0;
}


=head2

	SUBROUTINE		has_fields

	PURPOSE

		CHECK IF A TABLE HAS A PARTICULAR FIELD

=cut

sub has_field
{
	my $dbh		=	shift;
	my $table	=	shift;
	my $field	=	shift;

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


sub fields_arrayref
{
	my $dbh		=	shift;
	my $table	=	shift;

	my $query = qq{SHOW CREATE TABLE $table};
	my $result = simple_query($dbh, $query);

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

sub has_entry
{
	my $dbh		=	shift;
	my $table	=	shift;
	my $field	=	shift;

	if ( not defined $table or not $table )	{	croak "Table not defined or empty\n";	return;	}
	if ( not defined $field or not $field )	{	croak "Field not defined or empty\n";	return;	}

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

sub has_table
{

	my $dbh		=	shift;
	my $table	=	shift;

	if ( not defined $table or not $table )	{	croak "Table not defined or empty\n";	return;	}

	my $query = qq{SHOW TABLES};
	my $tables = simple_queryarray($dbh, $query);

	if ( not defined $tables )	{	return 0;	}

	for ( my $i = 0; $i < @$tables; $i++ )
	{
		if ( $table =~ /^$$tables[$i]$/ )
        {
            return 1;
        }
	}

	return 0;
}



sub update_fulltext
{
	my $dbh			=	shift;
	my $table		=	shift;
	my $fulltext	=	shift;
	my $fields		=	shift;

	if ( not defined $dbh or not defined $table or not defined $fulltext or not defined $fields )	{	return;	}

    $| = 1;
	if ( has_index($dbh, $table, $fulltext) )
	{
		print "Dropping index '$fulltext'\n";
		drop_index($dbh, $table, $fulltext);
	}

	my $query = qq{ALTER TABLE $table ADD FULLTEXT $fulltext ($fields)};
	print "$query\n";
	my $success = do_query($dbh, $query);

	return $success;
}


sub has_index
{
	my $dbh		=	shift;
	my $table	=	shift;
	my $index	=	shift;

	my $indices = indices($dbh, $table);
	if ( not defined $indices )	{	return 0;	}

	for ( my $i = 0; $i < @$indices; $i++ )
	{
		if ( $$indices[$i] =~ /^$index$/ )	{	return 1;	}
	}

	return 0;
}


sub indices
{
	my $dbh		=	shift;
	my $table	=	shift;

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


sub drop_index
{
	my $dbh		=	shift;
	my $table	=	shift;
	my $index	=	shift;

	my $query = qq{DROP INDEX $index ON $table};
	my $success = do_query($dbh, $query);

	return $success;	
}



sub create_custom_table
{
    my $source_dbh             =   shift;
    my $table           =   shift;
    my $sqlfile         =   shift;

    my $query = Util::contents($sqlfile);
    if ( not defined $query or not $query ) {   print "No/empty sql file '$sqlfile'\n"; exit;   }

    $query =~ s/^(.+?)\s+\w+\s*\(/$1 $table (/;

    $source_dbh->{RaiseError} = 0;
    my $count = $source_dbh->do($query);
    if ( not $count )   {   print "Did not create table '$table'\n";    }
    $source_dbh->{RaiseError} = 1;
}


sub file2array
{
	my $file	=	shift;

	my $array;
	open(FILE, $file) or die "Can't open file to read array: $file\n";
	$/ = "\n";
	while ( my $entry = <FILE> )	{ push @$array, $entry;	}

	return $array;
}


sub array2file
{
	my $array	=	shift;
	my $file	=	shift;

	if ( not defined $array )	{	return 0;	}
	if ( not -f $file )	{	return 0;	}

	open(OUTFILE, ">$file") or die "Can't open file to write array: $file\n";
	for ( my $i = 0; $i < @$array; $i++ )
	{
		print OUTFILE $$array[$i], "\n";
	}
	close OUTFILE;

	return 1;
}

sub sqlite_dbh
{
	my $sqlitefile		=	shift;

	my $dbh = DBI->connect("dbi:SQLite:$sqlitefile", '', '') || die "Can't open database file '$sqlitefile': $@\n";   	

	return $dbh;
}


sub do_query
{
	my $dbh			=	shift;
	my $query		=	shift;

	$dbh->{RaiseError} = 0; # DISABLE TERMINATION ON ERROR
	my $success = $dbh->do($query); ### or warn "Query failed: $DBI::errstr ($DBI::err)\n";
	$dbh->{RaiseError} = 1; # ENABLE TERMINATION ON ERROR

	if ( $success )	{	return 1;	}

	return 0;
}

sub tsvline2hash
{
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


sub simple_insert
{
	my $dbh			=	shift;
	my $query		=	shift;

	$dbh->{RaiseError} = 0; # DISABLE TERMINATION ON ERROR
	my $success = $dbh->do($query) or warn "Query failed: $DBI::errstr ($DBI::err)\n";
	$dbh->{RaiseError} = 1; # ENABLE TERMINATION ON ERROR

	return $success;
}

sub twoDhash2file
{
    my $filename    =   shift;
    my $hash    	=	shift;

	open(OUTFILE, ">$filename") or die "Can't open hasharray output file: $!\n";	
	foreach my $key ( keys %$hash )
	{
		my $line;
		my $subhash = $$hash{$key};

        #exit;
        #
        $line .= "$key\t";

        foreach my $subkey ( keys %$subhash)	{	$line .= "$subkey => $$subhash{$subkey}\t";	}

		$line =~ s/\t$//;
		$line .= "\n";
		print OUTFILE $line;
	}
	close(OUTFILE);
}


sub file2twoDhash
{
    my $filename    =   shift;

	my $hash;
	$/ = "\n";

    open(FILE, $filename) or die "Can't open result file '$filename': $!\n";
	while ( my $line = <FILE> )
    {
        if ( $line =~ /^\s*$/ )    {   next;   }   

	    my @elements = split "\t", $line;
		my $key = shift @elements;

        #### GET SUBHASH
        my $subhash;
        for ( my $i = 0; $i < $#elements + 1; $i++ )
		{
			my ($subkey, $value) = $elements[$i] =~ /^(.+?)\s+=\>\s+(.+)$/;
			$$subhash{$subkey} = $value;
		}
		$$hash{$key} = $subhash;
    }

    return $hash;
}


sub hasharray2file
{
    my $filename    =   shift;
    my $hasharray	=	shift;

	open(OUTFILE, ">$filename") or die "Can't open hasharray output file: $!\n";	
	for ( my $i = 0; $i < @$hasharray; $i++ )
	{
		my $hash = $$hasharray[$i];
		my $line;
		foreach my $key ( keys %$hash)	{	$line .= "$key => $$hash{$key}\t";	}
		$line =~ s/\t$//;
		$line .= "\n";
		print OUTFILE $line;
	}
	close(OUTFILE);
}

sub file2hasharray
{
    my $filename    =   shift;

	my $array;
	$/ = "\n";
    open(FILE, $filename) or die "Can't open result file '$filename': $!\n";
	while ( my $line = <FILE> )
    {
        if ( $line =~ /^\s*$/ )    {   next;   }   
		my $hash;
	    my @elements = split "\t", $line;
		for ( my $i = 0; $i < $#elements + 1; $i++ )
		{
			my ($key, $value) = $elements[$i] =~ /^(.+?)\s+=\>\s+(.+)$/;
			if ( defined $key )
			{
				if ( defined $value )	{	$$hash{$key} = $value;	}
				else	{	$$hash{$key} = '';	}
			}
		}
		push @$array, $hash;
    }

    return $array;
}

sub simple_query
{
    my $dbh     =   shift;
    my $query   =   shift;

    my $sth = $dbh->prepare($query);
	$sth->execute;

    return $sth->fetchrow;
}


sub simple_queryarray
{
    my $dbh     =   shift;
    my $query   =   shift;


    my $sth = $dbh->prepare($query);
	$sth->execute;

	my $array;
    my $result = 0;
	while ( defined $result )
    {
        $result = $sth->fetchrow;
        if ( defined $result )	{	push @$array, $result;	}
    }

    if ( not defined $array )
    {
        my $result = Database::simple_queryhash($dbh, $query);
        push @$array, $result;
    }


	return $array;
}


sub simple_queryhash
{
    my $dbh     =   shift;
    my $query   =   shift;

    my $sth = $dbh->prepare($query);
	$sth->execute;

    return $sth->fetchrow_hashref;
}

sub simple_queryhasharray
{
    my $dbh     =   shift;
    my $query   =   shift;


	$dbh->{RaiseError} = 0; # DISABLE TERMINATION ON ERROR
    my $sth = $dbh->prepare($query);
	$sth->execute;
	$dbh->{RaiseError} = 1; # ENABLE TERMINATION ON ERROR

	my $hasharray;
	while ( my $hash = $sth->fetchrow_hashref )
    {
        push @$hasharray, $hash;
    }

    #### DO SIMPLE HASH QUERY
    if ( not defined $hasharray )
    {
        my $result = Database::simple_queryhash($dbh, $query);
        push @$hasharray, $result;
    }


	return $hasharray;
}

sub quickie
{
	my $dbh	=	shift;
	my $query	=	shift;

	return simple_query($dbh, $query);
}

sub empty
{
    my $dbh         =   shift;
    my $table       =   shift;
    my $query = qq{SELECT * FROM $table LIMIT 1};
    my $sth = $dbh->prepare($query);
    $sth->execute;
    my $result = $sth->fetchrow;

    if ( $result )  {   return 0;   };
    return 1;
}

sub fields
{
    my $sqlfile         =   shift;

    my $fields;

    my $contents = Util::contents($sqlfile);
	$fields = content_fields($contents);

    return $fields;
}

sub content_fields
{
	my $contents		=	shift;

	my $content_fields = '';

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
        $content_fields .= $line . ",";
    }
    $content_fields =~ s/,$//;

	return $content_fields;
}



=head2

	SUBROUTINE		fields_insert

	PURPOSE

		CREATE THE AN INSERT QUERY .TSV LINE (FOR LOAD) BASED ON THE FIELDS

		AND THE HASH, MAKING DEFINED WHEN A FIELD key IS 

		NOT PRESENT IN THE HASH

=cut 

sub fields_insert
{
    my $fields          =   shift;
    my $hash        	=   shift;
	my $table			=	shift;

	$fields =~ s/\s//g;

	print "Database::fields_insert(fields, hash, table).\n";
	print " Hash:\n";
	print Dumper $hash;

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
        print "Field $fieldarray[$i]: ", make_defined($$hash{$fieldarray[$i]}), "\n";
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

sub fields_tsvline
{
    my $fields          =   shift;
    my $hash        =   shift;


	$fields =~ s/\s//g;

    my $tsvline = '';

    my @fieldarray = split ",", $fields;
    for ( my $i = 0; $i < $#fieldarray + 1; $i++ )
    {
        no warnings;
        $tsvline .= make_defined($$hash{$fieldarray[$i]});

        use warnings;

        $tsvline .= "\t";
    }
    $tsvline .= "\n";


    return $tsvline;
}






sub make_defined
{
    my $input   =   shift;

    if ( not defined $input )   {   return '';  }

    return $input;
}


sub load
{
	my $dbh				=	shift;
	my $table			=	shift;
	my $file			=	shift;
    my $extra           =   shift;

	#### LOAD DATA INTO TABLE
	my $query = qq{LOAD DATA INFILE '$file' INTO TABLE $table};
    if ( defined $extra )   {   $query .= " $extra";    }
	$dbh->{RaiseError} = 0; # DISABLE TERMINATION ON ERROR
	my $count = $dbh->do($query) or warn "Query failed: $DBI::errstr ($DBI::err)\n";
	if ( not $count ) {   print "DID NOT LOAD!\n";    }
	$dbh->{RaiseError} = 1; # ENABLE TERMINATION ON ERROR
}


sub load_nodups
{
	my $dbh				=	shift;
	my $table			=	shift;
	my $file			=	shift;
    my $extra           =   shift;

	#### LOAD DATA INTO TABLE
	my $query = qq{LOAD DATA LOCAL INFILE '$file' INTO TABLE $table};
    if ( defined $extra )   {   $query .= " $extra";    }
	$dbh->{RaiseError} = 0; # DISABLE TERMINATION ON ERROR
	my $count = $dbh->do($query) or warn "Query failed: $DBI::errstr ($DBI::err)\n";
	if ( not $count ) {   print "DID NOT LOAD!\n";    }
	$dbh->{RaiseError} = 1; # ENABLE TERMINATION ON ERROR
}


sub drop
{
	my $dbh						=	shift;
	my $type					=	shift;
	my $name					= 	shift;

	my $query = qq{DROP $type IF EXISTS $name };
	$dbh->{RaiseError} = 0;						# DISABLE TERMINATION ON ERROR
	my $count = $dbh->do($query);
#	if ( not $count )   {   print "Did not drop $type $name (it may not exist yet)\n";    }
	$dbh->{RaiseError} = 1;						# ENABLE TERMINATION ON ERROR
}


sub rows2tsv
{
    my $rows        =   shift;
    if ( not defined $rows )    {   return; }

    my $tsv = '';
    for ( my $i = 0; $i < @$rows; $i++ )
    {
        my $string = join "\t", @{$$rows[$i]};
        $tsv .= $string;
        $tsv .= "\n";
    }

    return $tsv;
}


sub query_rows
{
    my $dbh         =   shift;
    my $query       =   shift;

    my $rows;
    my $sth = $dbh->prepare($query);
    $sth->execute;
    while ( my @row = $sth->fetchrow )
    {
        push @$rows, \@row;
    }

    return $rows;
}

sub selector
{
    my $dbh                 =   shift;
    my $table               =   shift;
    my $select              =   shift;
    my $extra               =   shift;

    if ( not defined $extra )   {   $extra = '';    }
    my $sequences;
    my $query = qq{
    SELECT $select
    FROM $table
    $extra
};    
    my $sth = $dbh->prepare($query);
    $sth->execute;    
    while ( my $result = $sth->fetchrow_hashref )   {   push @$sequences, $result;  }

    return $sequences;
}

sub select_where
{
    my $dbh                 =   shift;
    my $table               =   shift;
    my $select              =   shift;
    my $field               =   shift;
    my $value               =   shift;
    my $extra               =   shift;

    my $output;
    @$output = ();
    my $query = qq{
    SELECT $select
    FROM $table
    WHERE $field='$value'
    $extra
};


    my $sth = $dbh->prepare($query);
    $sth->execute;    
    while ( my $result = $sth->fetchrow_hashref )
    {
        push @$output, $result;
    }

    return $output;
}


sub field_present
{
	my $dbh					=	shift;
	my $table				=	shift;
	my $field				=	shift;

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


#### GET FIELD NAMES
sub field_names
{
    my $dbh         =   shift;
    my $table       =   shift;

    my $sth = $dbh->prepare("SELECT * FROM $table LIMIT 1") or die "Can't prepare field names query!\n";
    $sth->execute or die "Can't execute field names query!\n";
    my $fields;
    $fields = $sth->{NAME};

    return $fields;   
}


# CREATE DATABASE IF NOT EXISTS
sub create_database
{
	my $database				=	shift;
	my $user					=	shift;
	my $password				=	shift;

    my $drh = DBI->install_driver("mysql");
	#if ( defined $drh )	{	print "Driver handle defined: $drh\n";	}
	my $remotecall = $drh->func('createdb', $database, 'localhost', $user, $password, 'admin');


	# GET LIST OF DATABASES -- DOESN'T WORK
	#my @dbs = $drh->func("$host:$port", '_ListDBs');

	## DOESN'T WORK
	#my @databases = DBI->data_sources("mysql");
}


# GET DATABASE HANDLE
sub database_handle
{
    my $dsn			= shift;    # data source name
    my $user		= shift;
    my $password	= shift;

	my %error_attributes = (
	);
    my $dbh = DBI->connect( $dsn, $user, $password, {	'PrintError' => 0, 'RaiseError' => 0 } );
	###or warn "Can't connect to database: $dsn\n";

    # RETURN DATABASE HANDLE
    return $dbh;
}


sub create_table
{
	my $dbh						=	shift;
	my $table					=	shift;
	my $sqlfile					= 	shift;
	my $change_table            =   shift;

	my $contents = Util::contents($sqlfile);
    my $query = $contents; 
    if ( defined $change_table )
    {
        $contents =~ s/^.+?\(//ms;
        $query = "CREATE TABLE IF NOT EXISTS $table ($contents";
    }

	$dbh->{RaiseError} = 0;						# DISABLE TERMINATION ON ERROR
	my $count = $dbh->do($query);
	# if ( not $count )   {   print "Could not create table $table\n";    }
	$dbh->{RaiseError} = 1;						# ENABLE TERMINATION ON ERROR

    return $count;
}


sub delete_from
{
    my $dbh                     =   shift;
    my $table                   =   shift;
    my $field                   =   shift;
    my $value                   =   shift;

    my $query = qq{DELETE FROM $table };
	if ( defined $field and defined $value )	{	$query .= qq{WHERE $field='$value'};	}
    $dbh->{RaiseError} = 0; # DISABLE TERMINATION ON ERROR
    my $count = $dbh->do($query) or warn "Query failed: $DBI::errstr ($DBI::err)\n";
    $dbh->{RaiseError} = 1; # ENABLE TERMINATION ON ERROR    
}



sub quality_files
{
	my $dbh 					= shift;
	my $all_sequencefile_path 	= shift;

	my $count = 0;

#	my $query = qq{SELECT sequence_ID, phd_quality_file FROM master_seq2 WHERE Sequence_ID<86522};
	my $query = qq{SELECT distinct Sequence_name, phd_quality_file FROM master_seq2 WHERE phd_quality_file!='NULL'};
	my $sth = $dbh->prepare($query);
	$sth->execute();

	my $number_rows = $sth->rows;
	print "NO. .phd.1 FILES RETURNED FROM DATABASE: $number_rows\n"; 
	while ( my @row = $sth->fetchrow )
	{
		my $sequence_name    = $row[0] ;
		my $phd_quality_file	= $row[1];

		if ($phd_quality_file !~ /^\s*$/ and $phd_quality_file !~ /^[0\s]*$/ )
		{
			my $filename = "$sequence_name.phd.1";
			my ($experiment) = $sequence_name =~ /^(\d+)/;

			open(FASTAFILE, ">$all_sequencefile_path/$experiment/$filename");
			print FASTAFILE $phd_quality_file;
			close(FASTAFILE);
			$count++;
		}
	}

	$count++;
}

sub query
{
    my $dbh             = shift;
    my $query           = shift;

    my $sth = $dbh->prepare($query);
    $sth->execute();

    return $sth;
}

sub number_rows
{
    my $sth             = shift;
    my $number_rows = $sth->rows; 

    return $number_rows;
}

sub row
{
    my $sth             = shift;

    my @row = $sth->fetchrow();

	return wantarray ? @row : $row[0];	
}

sub rows
{
    my $sth             = shift;

    my @rows;

	no warnings;
    while ( my @this_row = $sth->fetchrow )
    {
		my $row = join "\t", @this_row;
		if ( defined $row )
		{
			push @rows, join "\t", @this_row;
		}
	}
	use warnings;

    return @rows;	
}

sub last_id
{
	my $dbh				=	shift;
	my $table			=	shift;
	my $field			=	shift;

	my $last_id_query = "SELECT MAX($field) FROM $table";
	my $sth = $dbh->prepare($last_id_query);
	$sth->execute();
	my $last_id = $sth->fetchrow();

    if ( not defined $last_id ) {   return 0;   }

	return $last_id;
}


sub funnybase_version
{
	my $dbh						=	shift;
	my $assembly_number			=	shift;

	my $funnybase_version;

	my $query = qq{
	SELECT max(funnybaseversion)
	FROM funnybasesequences
	WHERE assemblynumber!=$assembly_number
	};

	my $sth = $dbh->prepare($query);
	$sth->execute;
	$funnybase_version = $sth->fetchrow;

	return $funnybase_version;
}


1;
