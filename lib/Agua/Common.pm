use MooseX::Declare;
#use Moose::Util::TypeConstraints;
#use Getopt::Simple;

##### SET THESE SLOTS IN INHERITING CLASS

#### USE LIB FOR INHERITANCE
use FindBin qw($Bin);
use lib "$Bin/../";
use lib "$Bin/../external";

use Agua::Common::Access;
use Agua::Common::Admin;
use Agua::Common::App;
use Agua::Common::Aws;
use Agua::Common::Base;
use Agua::Common::Cluster;
use Agua::Common::File;
use Agua::Common::Group;
use Agua::Common::Parameter;
use Agua::Common::Privileges;
use Agua::Common::Project;
use Agua::Common::Report;
use Agua::Common::SGE;
use Agua::Common::Shared;
use Agua::Common::Source;
use Agua::Common::Stage;
use Agua::Common::Transport;
use Agua::Common::User;
use Agua::Common::Util;
use Agua::Common::View;
use Agua::Common::Workflow;

=head2

	PACKAGE		Agua::Common

	PURPOSE 	
		O-O MODULE CONTAINING COMMONLY-USED Agua METHODS

=cut

class Agua::Common with (Agua::Common::Access,
	Agua::Common::Admin,
	Agua::Common::App,
	Agua::Common::Aws,
	Agua::Common::Base,
	Agua::Common::Cluster,
	Agua::Common::File,
	Agua::Common::Group,
	Agua::Common::Parameter,
	Agua::Common::Privileges,
	Agua::Common::Project,
	Agua::Common::Report,
	Agua::Common::SGE,
	Agua::Common::Shared,
	Agua::Common::Source,
	Agua::Common::Stage,
	Agua::Common::Transport,
	Agua::Common::User,
	Agua::Common::Util,
	Agua::Common::View,
	Agua::Common::Workflow)
{
	####///}

#### STRINGS/INTS	
has 'validated'	=> ( isa => 'Int', is => 'rw', default => 0 );
has 'cgi'		=> ( isa => 'Str|Undef', is => 'rw', default => 0 );
#### OBJECTS
has 'json'		=> ( isa => 'HashRef', is => 'rw', default => undef );
has 'conf'		=> ( isa => 'HashRef', is => 'rw', default => undef );
has 'dbobject'	=> ( isa => 'Agua::DBase::MySQL', is => 'rw', required => 0 );
has 'jsonParser'=> ( isa => 'JSON', is => 'rw');

use strict;
use warnings;
use Carp;

#### EXTERNAL MODULES
use FindBin qw($Bin);
use Data::Dumper;
use File::Path;
use File::Copy;
use File::Remove;
use File::stat;
use JSON;


}


