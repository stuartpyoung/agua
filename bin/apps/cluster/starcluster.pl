#!/usr/bin/perl -w

#### DEBUG

=head2

=head3 B<APPLICATION>     starcluster.pl

	PURPOSE

		DRIVE TESTS OF StarCluster.pm, WHICH PERFORMS THE FOLLOWING TASKS:

			1. MOUNT BY DEFAULT /agua, /data AND /nethome ON STARCLUSTER NODES

			2. CLUSTER SHUTS DOWN WHEN ALL WORKFLOWS ARE COMPLETED

			3. ALLOW AQUARIUS USERS TO RUN JOBS ON small, medium OR large CLUSTERS
			    (ALL USERS USE admin USER'S CONFIG FILE)

			4. ALL WORKFLOWS USE ONE CLUSTER, SPECIFIED BY admin USER

	USAGE

		./starcluster.pl <mode> [additional_arguments]

/data/agua/0.5/bin/apps/cluster/starcluster.pl start \
--username admin \
--cluster smallcluster \
--privatekey /nethome/admin/.keypairs/private.pem \
--publiccert /nethome/admin/.keypairs/public.pem \
--keyname admin-key

=cut

use strict;

#### USE LIBS
use FindBin qw($Bin);
use lib "$Bin/../../../lib/";

my $configfile = "$Bin/../../../conf/default.conf";
use Conf::Agua;
my $conf = Conf::Agua->new(
	inputfile	=>	$configfile,
	backup		=>	1,
	separator	=>	"\t",
	spacer		=>	"\\s\+"
);

#### INTERNAL MODULES
use Agua::StarCluster;
my $starcluster = Agua::StarCluster->new(
	conf 		=>	$conf
);

#### GET MODE AND ARGUMENTS
my @arguments = @ARGV;
my $mode = shift @ARGV;

#### PRINT HELP
if ( $mode eq "-h" or $mode eq "--help" )	{	help();	}

#### FLUSH BUFFER
$| =1;

#### RUN QUERY
no strict;
eval { $starcluster->$mode() };
if ( $@ ){
	print "Error - mode '$mode' might not be supported\nDetailed error output:\n$@\n";
}

