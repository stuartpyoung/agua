package Agua::Common::Util;
use Moose::Role;
#use Moose::Util::TypeConstraints;

=head2

	PACKAGE		Agua::Common::Util

	PURPOSE

		UTILITY METHODS FOR Agua::Common

=cut

#use Agua::DBaseFactory;


#### USE LIB FOR INHERITANCE
use FindBin qw($Bin);
use lib "$Bin/../../";
use lib "$Bin/../../external";

has 'logfile'	=> ( isa => 'Str|Undef', is => 'rw', required => 0 );

sub setDbh {
	my $self		=	shift;


	my $conf = $self->conf();
	print "Agua::Util::setDbh    conf not defined\n" and exit if not $self->conf();
	use Data::Dumper;

	my $dbfile 		=	$conf->getKeyValue('database', 'DBFILE');
	my $dbtype 		=	$conf->getKeyValue('database', 'DBTYPE');
	my $user 		=	$conf->getKeyValue('database', 'USER');
	my $password 	=	$conf->getKeyValue('database', 'PASSWORD');
	my $database	=	$conf->getKeyValue('database', 'DATABASE');
	print "Agua::Util::setDbh    dbtype not defined\n" and exit if not $dbtype;
	print "Agua::Util::setDbh    user not defined\n" and exit if not $user;
	print "Agua::Util::setDbh    password not defined\n" and exit if not $password;
	print "Agua::Util::setDbh    database not defined\n" and exit if not $database;

	#### FOR DEBUGGING: SET DATABASE IF PROVIDED IN JSON
	if ( $self->can('json') )
	{
		my $json = $self->json();
		$database = $json->{database} if defined $json->{database} and $json->{database};
	}

	##### CREATE DB OBJECT USING DBASE FACTORY
	my $dbobject = 	Agua::DBaseFactory->new(
		$dbtype,
		{
			dbfile		=>	$conf->getKeyValue('database', 'DBFILE'),
			database	=>	$conf->getKeyValue('database', 'DATABASE'),
			user        =>  $conf->getKeyValue('database', 'USER'),
			password    =>  $conf->getKeyValue('database', 'PASSWORD')
		}
	) or print qq{ error: 'Agua::Util::setDbh    Cannot create database object $database: $!' } and exit;
	print "Agua::Util::setDbh    dbobject not defined\n" and exit if not defined $dbobject;

	$self->dbobject($dbobject);	
}

sub objectInArray {
	my $self	=	shift;

	return 1 if defined $self->_indexInArray(@_);

	return 0;
}

sub _indexInArray {
	my $self	=	shift;
	my $array	=	shift;
	my $object	=	shift;
	my $keys	=	shift;

	use Data::Dumper;

	return if not defined $array or not defined $object;

	for ( my $i = 0; $i < @$array; $i++ )
	{

		my $identified = 1;
		for ( my $j = 0; $j < @$keys; $j++ )
		{
			print "not identified\n" and $identified = 0 if $$array[$i]->{$$keys[$j]} ne $object->{$$keys[$j]};
			print "identified: $identified\n";
		}

		print "Agua::Common::Util::_indexInArray    Returning $i\n" if $identified;
		return $i if $identified;
	}

	return;
}


=head2

	SUBROUTINE		datetime

	PURPOSE

		RETURN THE CURRENT DATE AND TIME

=cut

sub datetime {
	my $self	=	shift;

	my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime();
	$sec = sprintf "%02d", $sec;
	$min = sprintf "%02d", $min;
	$hour = sprintf "%02d", $hour;
	$mday = sprintf "%02d", $mday;
	$mon = sprintf "%02d", $mon;
	$year -= 100;
	$year = sprintf "%02d", $year;
	my $datetime = "$year-$mon-$mday-$hour-$min-$sec";

	return $datetime;
}


=head2

	SUBROUTINE		parseHash

	PURPOSE

		HASH INTO ARRAY OF TAB-SEPARATED KEY PAIRS

=cut

sub parseHash {
	my $self	=	shift;
	my $json	= 	shift;


	#### INITIATE JSON PARSER
	use JSON -support_by_pp; 
	my $jsonParser = JSON->new();

	my $outputs = [];

	if ( defined $json )
	{
		my $hash = $jsonParser->decode($json);
		my @keys = keys %$hash;
		@keys = sort @keys;

		foreach my $key ( @keys )
		{
			my $value = $hash->{$key}->{value};
			if ( defined $value and $value )
			{
				push @$outputs, "$key\t$value\n";	
			}
		}
	}

	return $outputs;
}



sub hasharrayToHash {
	my $self		=	shift;
	my $hasharray	=	shift;
	my $key			=	shift;

	print "Agua::Commonn::Util::hasharrayToArray    hasharray not defined.\n"
		and exit if not defined $hasharray;
	print "Agua::Commonn::Util::keyToArray    key not defined.\n"
		and exit if not defined $key;

	my $hash = {};
	foreach my $entry ( @$hasharray )
	{
		#### GET PROJECT AND WORKFLOW
		my $key = $entry->{owner};
		if ( not exists $hash->{$key} )
		{
			$hash->{$key} = [ $entry ];
		}
		else
		{
			push @{$hash->{$key}}, $entry;
		}
	}

	return $hash;
}
=head2

	SUBROUTINE		json_parser

	PURPOSE

		RETURN A JSON PARSER OBJECT

=cut

sub json_parser {
	my $self		= 	shift;

	return $self->jsonParser() if $self->jsonParser();

	use JSON -support_by_pp; 
	my $jsonParser = JSON->new();

	$self->jsonParser($jsonParser);

	return $self->jsonParser();
}


sub cowCase {
    my $self    =   shift;
    my $string   =   shift;

	return uc(substr($string, 0, 1)) . substr($string, 1);
}


sub openLogfile {
	my $self		=   shift;
	my $logfile 	= 	shift;

	print "Agua::Common::Util::openLogfile    Agua::Common::Util::openLogfile(logfile)\n";
	print "Agua::Common::Util::openLogfile    logfile: $logfile\n";
	$logfile .= ".1";

	if ( -f $logfile )
	{
		my ($stub, $index) = $logfile =~ /^(.+?)\.(\d+)$/;
		$index++;
		$logfile = $stub . "." . $index;
	}

	open (STDOUT, "| tee -ai $logfile") or die "Can't split STDOUT to logfile: $logfile\n";
	print "Agua::Common::Util::openLogfile    Writing to logfile: $logfile\n";

	return $logfile;	
}

1;

