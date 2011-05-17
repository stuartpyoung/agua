#use Moose;
use MooseX::Declare;


=head2

	PACKAGE		Balancer

	PURPOSE

		MONITOR THE LOAD BALANCER RUNNING FOR EACH CLUSTER (I.E., SGE CELL)

		AND RESTART IT IF IT HAS BECOME INACTIVE (E.G., DUE TO A RECURRENT

		ERROR IN STARCLUSTER WHEN PARSING XML OUTPUT AND/OR XML ERROR MESSAGES

		CAUSED BY A BUG IN SGE)

	NOTES    

		1. CRON JOB RUNS THE SCRIPT checkBalancers.pl EVERY MINUTE

		* * * * * /agua/0.6/bin/scripts/checkBalancers.pl

=cut

use strict;
use warnings;
use Carp;

##### USE LIBS
#use FindBin qw($Bin);
#use lib "$Bin/../";
#use lib "$Bin/../external";

class Agua::Balancer with (Agua::Common::Aws, Agua::Common::Cluster, Agua::Common::Util) {

#### EXTERNAL MODULES
use Data::Dumper;

use FindBin qw($Bin);
use lib "$Bin/..";

#### INTERNAL MODULES	
use Agua::DBaseFactory;
use Conf::Agua;

# INTS
has 'keydir'	=> ( isa => 'Str|Undef', is => 'rw' );
has 'adminkey'	=> ( isa => 'Str|Undef', is => 'rw' );

# STRINGS
has 'fileroot'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'cluster'	=>  ( isa => 'Str', is => 'rw'	);
has 'username'  =>  ( isa => 'Str', is => 'rw'	);
# OBJECTS
has 'json'		=> ( isa => 'HashRef', is => 'rw', required => 0 );
has 'dbobject'	=> ( isa => 'Agua::DBase::MySQL', is => 'rw', required => 0 );

has 'conf' 	=> (
	is =>	'rw',
	'isa' => 'Conf::Agua',
	default	=>	sub { Conf::Agua->new(	backup	=>	1, separator => "\t"	);	}
);

####////}	

method BUILD ($hash) {
	print "Agua::Workflow::BUILD    Agua::Workflow::BUILD()\n";    

	#### SET DATABASE HANDLE
	$self->setDbh();
}

=head2

	SUBROUTINE		checkBalancers

	PURPOSE

		VERIFY THAT A LOAD BALANCER IS RUNNING FOR EACH ACTIVE CLUSTER

		AND RESTART THE LOAD BALANCER IF STOPPED

	NOTES    

		1. CHECKS clusterstatus DATABASE TABLE FOR RUNNING CLUSTERS

		2. CHECKS IF BALANCER PROCESS IS RUNNING USING ITS PID

		3. DOES THE FOLLOWING TO RESTART BALANCER:

			-   START NEW BALANCER PROCESS, CONCAT OUTPUT TO OUTPUT FILE

			- 	UPDATE clusterstatus TABLE WITH NEW PID

=cut

method checkBalancers {

	#### GET ALL status='running' CLUSTERS
	my $clusterobjects = $self->activeClusters();

	#### KEEP ONLY CLUSTERS WHERE THE LOAD BALANCER HAS FAILED	
	for ( my $i = 0; $i < @$clusterobjects; $i++ )
	{
		my $clusterobject = $$clusterobjects[$i];
		my $pid = $clusterobject->{pid};
		print "Agua::Workflow::checkBalancers    No pid for clusterobject $clusterobject->{clusterobject}\n" and next if not defined $pid or not $pid;
		my $running = kill 0, $pid;
		if ( $running )
		{
			splice @$clusterobjects, $i, 1;
			$i--;
		}
	}

	#### RESTART LOAD BALANCERS FOR REMAINING CLUSTERS
	foreach my $clusterobject (@$clusterobjects)
	{
		$self->restartBalancer($clusterobject);
	}
}

method restartBalancer ($clusterobject) {

	$clusterobject = $self->clusterInputs($clusterobject);

	#### INSTANTIATE STARCLUSTER OBJECT
	my $starcluster = Agua::StarCluster->new($clusterobject);
	$starcluster->launchBalancer();
}


}

1;


