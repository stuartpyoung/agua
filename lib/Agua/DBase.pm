use MooseX::Declare;

#### DEBUG

=head2

	PACKAGE		Agua::DBase

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

class Agua::DBase {

# STRINGS
has 'sqlite'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'user'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'password'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'database'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'dbfile'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'dbtype'	=> ( isa => 'Str|Undef', is => 'ro', default => '' );

#### EXTERNAL MODULES
use DBI;
use DBD::mysql;
use POSIX;
use Data::Dumper;

##///}


sub timestamp {
    my $self		=	shift;
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
    my $timestamp = sprintf
		"%4d-%02d-%02d %02d:%02d:%02d",
		$year+1900,$mon+1,$mday,$hour,$min,$sec;	

	return $timestamp;
}

=head2

	SUBROUTINE		set

	PURPOSE

		CREATE A 'SET ...' LINE FOR ONE OR MORE FIELD VALUES

	INPUTS

		1. HASH OF TABLE field KEY-VALUE PAIRS 

		2. ARRAY OF FIELDS TO BE USED

=cut 

method set ($hash, $fields) {
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

method where ($hash, $fields) {
	return if not defined $hash or not defined $fields or not @$fields;

	my $where = '';
    for ( my $i = 0; $i < @$fields; $i++ )
    {
        my $value = $$hash{$$fields[$i]} ? $$hash{$$fields[$i]} : '';
		$value =~ s/'/\\'/g;

		$where .= qq{WHERE $$fields[$i] = '$value'\n} if $i == 0;
		$where .= qq{AND $$fields[$i] = '$value'\n} if $i != 0;
    }

	#### MYSQL-SAFE
	#$where =~ s/'/\'/;

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

method notDefined ($hash, $fields) {


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

method make_defined ($input) {    
    return '' if not defined $input;    
    return $input;
}

=head2

	SUBROUTINE		query

	PURPOSE

		RETURN A SINGLE QUERY RESULT

=cut

sub query {}
=head2

	SUBROUTINE		has_field

	PURPOSE

		CHECK IF A TABLE HAS A PARTICULAR FIELD

=cut

sub has_field {}
=head2

	SUBROUTINE		fields_arrayref

	PURPOSE

		CHECK IF THERE ARE ANY ENTRIES FOR THE GIVEN FIELD IN

		THE TABLE (I.E., MUST BE NON-EMPTY AND NON-NULL). IF NONE

		EXIST, RETURN ZERO. OTHERWISE, RETURN 1.

=cut

sub has_entry {}
=head2

	SUBROUTINE		has_table

	PURPOSE

		1. CHECK IF A TABLE EXISTS IN THE DATABASE

		2. RETURN 1 IF EXISTS, RETURN ZERO IF OTHERWISE

=cut

sub has_table {}
=head2

	SUBROUTINE		update_fulltext

	PURPOSE

        UPDATE THE FULLTEXT INDEX OF A TABLE

=cut 

sub update_fulltext {}
=head2

	SUBROUTINE		has_index

	PURPOSE

		RETURN 1 IF THE NAME IS AMONG THE LIST OF A TABLE'S INDICES 

=cut 

sub has_index {}
=head2

	SUBROUTINE		indices

	PURPOSE

		RETURN A LIST OF A TABLE'S INDICES 

=cut 

sub indices {}
=head2

	SUBROUTINE		drop_index

	PURPOSE

		RETURN A HASH FROM A TSV LINE AND ORDERED LIST OF FIELDS

=cut 

sub drop_index {}
=head2

	SUBROUTINE		tsvline2hash

	PURPOSE

		RETURN A HASH FROM A TSV LINE AND ORDERED LIST OF FIELDS

=cut 

sub tsvline2hash {}
=head2

	SUBROUTINE		query

	PURPOSE

		RETURN AN ARRAY OF SCALARS

=cut 

sub do {}
=head2

	SUBROUTINE		queryarray

	PURPOSE

		RETURN AN ARRAY OF SCALARS

=cut 

sub queryarray {}
=head2

	SUBROUTINE		queryhash

	PURPOSE

		RETURN A HASH

=cut 

sub queryhash {}
=head2

	SUBROUTINE		queryhasharray

	PURPOSE

		RETURN AN ARRAY OF HASHES

=cut 

sub queryhasharray {}
=head2

	SUBROUTINE		is_empty

	PURPOSE

		RETURN 1 IF TABLE IS EMPTY

=cut 

sub is_empty {}
=head2

	SUBROUTINE		sql2fields

	PURPOSE

		RETURN LIST OF FIELDS FROM AN SQL FILE

=cut 

sub sql2fields {}
=head2

	SUBROUTINE		fields_insert

	PURPOSE

		CREATE AN INSERT QUERY .TSV LINE (FOR LOAD) BASED ON THE FIELDS

		AND THE HASH, MAKING DEFINED WHEN A FIELD key IS 

		NOT PRESENT IN THE HASH

=cut 

sub fields_insert {}
=head2

	SUBROUTINE		fields_tsvline

	PURPOSE

		CREATE A .TSV LINE (FOR LOAD) BASED ON THE FIELDS

		AND THE HASH, MAKING DEFINED WHEN A FIELD key IS 

		NOT PRESENT IN THE HASH

=cut 

sub fields_tsvline {}
=head2

	SUBROUTINE		load

	PURPOSE

		RETURN '' IF NOT DEFINED

=cut 

sub load {}
=head2

	SUBROUTINE		load_nodups

	PURPOSE

		LOAD DATA FROM FILE INTO TABLE AND IGNORE DUPLICATES

=cut 

sub load_nodups {}
=head2

	SUBROUTINE		field_present

	PURPOSE

		RETURN '' IF NOT DEFINED

=cut 

sub field_present {}
=head2

	SUBROUTINE		fields

	PURPOSE

		RETURN '' IF NOT DEFINED

=cut 

sub fields {}
=head2

	SUBROUTINE		create_database

	PURPOSE

        CREATE DATABASE IF NOT EXISTS

=cut 

sub create_database {}

} #### END

