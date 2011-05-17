package Agua::Common::Aws;
use Moose::Role;
#use Moose::Util::TypeConstraints;

=head2

	PACKAGE		Agua::Common::Aws

	PURPOSE

		AMAZON WEB SERVICES (AWS) METHODS FOR Agua::Common

=cut


use FindBin qw($Bin);
use lib "$Bin/../..";

use Data::Dumper;
use Agua::StarCluster;


=head2

	SUBROUTINE		getPublicfile

	PURPOSE

		SET LOCATION OF PUBLIC FILE

		NB: THIS SHADOWS AWS::get_publicfile

=cut

sub getPublicfile {
	my $self	=	shift;
	my $username=	shift;

	my $keydir = $self->getKeydir($username);
	my $publicfile = "$keydir/public.pem";

	return $publicfile;
}

=head2

	SUBROUTINE		getPrivatefile

	PURPOSE

		SET LOCATION OF PRIVATE FILE

		NB: THIS SHADOWS AWS::get_privatefile

=cut

sub getPrivatefile {
	my $self	=	shift;
	my $username=	shift;

	my $keydir = $self->getKeydir($username);
	my $privatefile = "$keydir/private.pem";

	return $privatefile;
}

=head2

	SUBROUTINE		getKeydir

	PURPOSE

		SET LOCATION OF PRIVATE FILE

		NB: THIS SHADOWS AWS::getKeydir BUT WITH

		AN IMPORTANT DIFFERENCE - THE KEYS ARE STORED

		IN USERDIR/USERNAME/.keypairs


=cut

sub getKeydir {
	my $self	=	shift;
	my $username=	shift;

	return $self->keydir() if $self->can('keydir') and defined $self->keydir();

	my $json	=	$self->json();
	$username = $json->{username} if not defined $username;
	$username = $self->username() if not defined $username;

	my $userdir = $self->conf()->getKeyValue('agua', 'USERDIR');

	print "Agua::Common::Aws::getKeydir    username not defined. Returning\n" and exit if not defined $username;

	my $keydir = "$userdir/$username/.keypairs";
	$self->keydir($keydir) if $self->can('keydir');

	return $keydir;
}


=head2

    SUBROUTINE:     getAws

    PURPOSE:

		RETURN THESE GROUP-RELATED TABLES:

			groupmember 

			etcgroup

	INPUTS

		1. JSON OBJECT CONTAINING username AND sessionId

	OUTPUTS	
			A HASH WITH TWO KEYS:

			groups: JSON groupmember HASH
			users : JSON etcusers ARRAY

=cut

sub getAws {
	my $self		=	shift;
	my $username	=	shift;

    my $dbobject		=	$self->dbobject();

	#### GET AWS INFO FOR USER
	my $adminkey = $self->getAdminKey($username);
	my $aws;
	if ( $adminkey )
	{
		my $query = qq{SELECT *
FROM aws
WHERE username='admin'};
		$aws = $dbobject->queryhash($query);
	}

	#### USE ADMIN KEYS IF USER HAS NO AWS CREDENTIALS
	if ( not defined $aws )
	{
		my $query = qq{SELECT *
	FROM aws
	WHERE username='$username'};
		print "$query\n";
		$aws = $dbobject->queryhash($query);
	}

	return $aws;
}


=head2

	SUBROUTINE		printKeyfiles

	PURPOSE

		PRINT THE PRIVATE KEY AND PUBLIC CERTIFICATE TO FILES

		NB: SHADOWS Agua::Common::Aws::printKeyfiles

=cut

sub printKeyfiles {
	my $self		=	shift;


	#### GET PRIVATE KEY AND PUBLIC CERTIFICATE
    my $json 		=	$self->json();
	print "Agua::Common::Aws::printKeyfiles    json:\n";
	print Dumper $json;

	my $publiccert	=	$json->{ec2publiccert};
	my $privatekey	=	$json->{ec2privatekey};

	#### PUBLIC CERT HAS 65-LETTER LINES
	my ($publicstring) = $publiccert =~ /^\s*\-{1,5}BEGIN CERTIFICATE\-{1,5}\s*(.+?)\s*\-{1,5}END CERTIFICATE\-{1,5}$/ms;
	$publicstring =~ s/\s+//g;
	my $public = $self->formatLines($publicstring, 64);
	$public = "-----BEGIN CERTIFICATE-----\n"
				. $public
				. "-----END CERTIFICATE-----";

	#### PUBLIC CERT HAS 65-LETTER LINES
	my ($privatestring) = $privatekey =~ /^\s*\-{1,5}BEGIN PRIVATE KEY\-{1,5}\s*(.+?)\s*\-{1,5}END PRIVATE KEY\-{1,5}$/ms;
	$privatestring =~ s/\s+//g;
	my $private = $self->formatLines($privatestring, 76);
	$private = "-----BEGIN PRIVATE KEY-----\n"
				. $private
				. "-----END PRIVATE KEY-----";

	#### GET KEYFILES
    my $publicfile = $self->getPublicfile();
	my $privatefile = $self->getPrivatefile();

	#### CREATE KEYDIR
	print "Agua::Common::Aws:printKeyfiles       Can't create keydir\n" and exit if not $self->createKeydir();

	open(PUBLICFILE, ">$publicfile") or die "AWS:printKeyfiles    Can't open publicfile: $publicfile\n";
	print PUBLICFILE $public;
	close(PUBLICFILE) or die "AWS:printKeyfiles    Can't close publicfile: $publicfile\n";
	#`chmod 600 $publicfile`;

	open(PRIVATEFILE, ">$privatefile") or die "AWS:printKeyfiles    Can't open privatefile: $privatefile\n";
	print PRIVATEFILE $private;
	close(PRIVATEFILE) or die "AWS:printKeyfiles    Can't close privatefile: $privatefile\n";
	`chmod 600 $privatefile`;

	print "Agua::Common::Aws::printKeyfiles    private.pem:\n";
	print `cat $privatefile`;
	print "Agua::Common::Aws::printKeyfiles    \n";
	print "Agua::Common::Aws::printKeyfiles    public.pem:\n";
	print `cat $publicfile`;
	print "Agua::Common::Aws::printKeyfiles    \n";

}


sub formatLines {
	my $self	=	shift;
	my $string	=	shift;
	my $length	=	shift;


	my $lines = '';
	my $offset = 0;
	while ( $offset < length($string) )
	{
		$lines .= substr($string, $offset, $length);
		$lines .= "\n";
		$offset += $length;
	}


	return $lines;	
}
sub createKeydir {
	my $self	=	shift;

	print "Agua::Common::Aws::createKeydir    AWS::createKeydir()\n";

	my $keydir = $self->getKeydir();
	print "Agua::Common::Aws::createKeydir    keydir: $keydir\n";
	File::Path::mkpath($keydir) if not -d $keydir;
	print "Agua::Common::Aws::createKeydir    Can't create keydir: $keydir\n" and return 0 if not -d $keydir;

	return 1;
}




sub generateKeypairfile {
	my $self		=	shift;

	#### PASSED VARIABLES
	my $conf 			= 	$self->conf();
	my $json 			= 	$self->json();

	#### SET KEYNAME
	my $username 		=	$json->{username};
	my $keyname = "$username-key";

	#### SET PRIVATE KEY AND PUBLIC CERT FILE LOCATIONS	
	my $privatekey = $self->getPrivatefile();
	my $publiccert = $self->getPublicfile();

	print "Agua::Common::Aws::generateKeypairfile    username: $username\n";
	print "Agua::Common::Aws::generateKeypairfile    keyname: $keyname\n";
	print "Agua::Common::Aws::generateKeypairfile    privatekey: $privatekey\n";
	print "Agua::Common::Aws::generateKeypairfile    publiccert: $publiccert\n";

	#### RUN starcluster.pl TO GENERATE KEYPAIR FILE IN .starcluster DIR
	my $starcluster = Agua::StarCluster->new(
		{
			privatekey	=>	$privatekey,
			publiccert	=>	$publiccert,
			username	=>	$username,
			keyname		=>	$keyname,
			conf		=>	$self->conf()
		}
	);

	print "Agua::Common::Aws::saveAws    Doing starcluster->generateKeypair()\n";
	$starcluster->generateKeypair();

#	my $starcluster = "$installdir/bin/apps/cluster/starcluster.pl";
#	my $command = qq{$starcluster generateKeypair \\
#--privatekey $privatekey \\
#--publiccert $publiccert \\
#--keyname $keyname \\
#--username $username
#};
}



=head2

    SUBROUTINE     addAws

	PURPOSE

		SAVE THE AWS AUTHENTICATION INFORMATION FOR THIS USER 

		TO THE aws TABLE

=cut

sub addAws {
	my $self		=	shift;

    my $json 		=	$self->json();
    my $dbobject 	=	$self->dbobject();


	#### CHECK IF THE USER ALREADY HAS STORED AWS INFO,
	#### IN WHICH CASE QUIT
	my $username = $json->{username};	
	my $query = qq{SELECT *
FROM aws
WHERE username='$username'};
	print "$query\n";
	my $aws = $dbobject->queryhash($query);
	print Dumper $aws;

	#### REMOVE IF EXISTS ALREADY
	if ( defined $aws )
	{
		my ($success, $user) = $self->_removeAws();
	}

	#### ADD TO TABLE
	my $success = $self->_addAws();
 	print "{ error: 'Agua::Common::Aws::addAws    Could not add entry $json->{name}' }" and exit if not defined $success;

	##### REMOVE WHITESPACE
	$json->{ec2publiccert} =~ s/\s+//g;
	$json->{ec2privatekey} =~ s/\s+//g;
	$aws->{ec2publiccert} =~ s/\s+//g;
	$aws->{ec2privatekey} =~ s/\s+//g;

	if ( $aws->{ec2privatekey} ne $json->{ec2privatekey}
		or $aws->{ec2publiccert} ne $json->{ec2publiccert} )
	{
		#### PRINT KEY FILES
		#$self->printKeyfiles();

		#### GENERATE KEYPAIR FILE FROM KEYS
		#$self->generateKeypairfile();
	}

 	print "{ status: 'Agua::Common::Aws::saveAws    Added AWS credentials for user $json->{username}' }";
	exit;
}


=head2

    SUBROUTINE     removeAws

	PURPOSE

		REMOVE A REPORT FROM THE entry TABLE

=cut

sub removeAws {
	my $self		=	shift;


    my $json 			=	$self->json();

	#### DO THE REMOVE
	my $success = $self->_removeAws();

 	print "{ error: 'Agua::Common::Aws::removeAws    Could not remove user $json->{username} from user table\n" and exit if not defined $success;
 	print "{ status: 'Agua::Common::Aws::removeAws    Removed user $json->{username} from user table\n";
	exit;
}



=head2

    SUBROUTINE     _addAws

	PURPOSE

		ADD A REPORT TO THE entry TABLE

=cut

sub _addAws {
	my $self		=	shift;


    my $dbobject        =	$self->dbobject();
    my $cgi 			=	$self->cgi();
    my $json 			=	$self->json();


	#### SET TABLE AND REQUIRED FIELDS	
	my $table = "aws";
	my $required_fields = ["username", "amazonuserid"];

	#### CHECK REQUIRED FIELDS ARE DEFINED
	my $not_defined = $dbobject->notDefined($json, $required_fields);
    print "{ error: 'Agua::Common::Aws::_addAws    undefined values: @$not_defined' }" and exit if @$not_defined;

	#### DO THE ADD
	return $self->_addToTable($table, $json, $required_fields);	
}

=head2

    SUBROUTINE     _removeAws

	PURPOSE

		REMOVE A REPORT FROM THE entry TABLE

=cut

sub _removeAws {
	my $self		=	shift;


    my $dbobject        =	$self->dbobject();
    my $cgi 			=	$self->cgi();
    my $json 			=	$self->json();


	#### SET TABLE AND REQUIRED FIELDS	
	my $table = "aws";
	my $required_fields = ["username"];

	#### CHECK REQUIRED FIELDS ARE DEFINED
	my $not_defined = $dbobject->notDefined($json, $required_fields);
    print "{ error: 'Agua::Common::Aws::_removeAws    undefined values: @$not_defined' }" and exit if @$not_defined;

	#### DO THE REMOVE
	return $self->_removeFromTable($table, $json, $required_fields);
}

=head2

	SUBROUTINE		mountVolume

	PURPOSE

		MOUNT AN EBS VOLUME FOR A GIVEN USER TO A WORKFLOW DIRECTORY

		IN THE USER'S HOME DIRECTORY

	INPUT

		1. THE USER'S AWS ACCESS KEY ID AND SECRET ACCESS KEY

		2. VOLUME (OR SNAPSHOT) ID, E.G., vol-e0e0e0e0, snap-e0e0e0e0

		3. DESTINATION MOUNT POINT, I.E., USER'S PROJECT AND WORKFLOW

=cut

sub mountUserVolume {
	my $self		=	shift;

	my $json		=	$self->get_json();
	print "Agua::Common::Aws::mountVolume    json:\n";
	print Dumper  $json;

	my $username	=	$json->{username};
	my $mountpoint	=	$json->{mountpoint};
	my $volume		=	$json->{volume};

	#### CHOOSE AN UNUSED DEVICE
	my $device		=	$json->{device};

	my $keys = $self->getKeys($username);
	print "Agua::Common::Aws::mountVolume    keys:\n";
	print Dumper $keys;

	my $accesskeyid = $keys->{accesskeyid};
	my $secretaccesskey = $keys->{secretaccesskey};



}


1;