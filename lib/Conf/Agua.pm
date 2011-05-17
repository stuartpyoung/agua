package Conf::Agua;

use Data::Dumper;

=head2

	PACKAGE		Conf::Agua

    PURPOSE

        1. READ AND WRITE ini-FORMAT CONFIGURATION FILES

			FOR STARCLUSTER, E.G.:

				[section name1]
				KEY1=VALUE1
				KEY2=VALUE2

				[section name2]
				KEY3=name2
				KEY4=VALUE3

=cut
use Conf;
use Moose;
with 'Conf';

#use MooseX::FollowPBP;
#use MooseX::Params::Validate;

has 'backup' 	=>	(	is	=>	'rw',	isa	=>	'Bool', default	=>	undef	);
has 'inputfile' =>	(	is	=>	'rw',	isa	=>	'Str'	);
has 'outputfile' => (	is	=>	'rw',	isa	=>	'Str'	);
has 'ignore' 	=>	(	is	=>	'rw',	isa	=>	'Str',	default	=>	"#;\\[?"	);
has 'spacer' 	=>	(	is	=>	'rw',	isa	=>	'Str',	default	=>	"\\t=\\s"	);
has 'separator'	=>	(	is	=>	'rw',	isa	=> 	'Str',	default	=>	"="		);
has 'comment'	=>	(	is	=>	'rw',	isa	=> 	'Str',	default	=>	"#"		);
has 'sections'	=>	(	is 	=>	'rw', 	isa =>	'ArrayRef' );

has 'separator'	=>	(	is	=>	'rw',	isa	=> 	'Str',	default	=>	"="		);
has 'comment'	=>	(	is	=>	'rw',	isa	=> 	'Str',	default	=>	"#"		);


#=head2
#
#	SUBROUTINE 		write
#	
#	PURPOSE
#
#		WRITE sections TO FILE WITH CONSTANT-WIDTH INDENT
#		
#=cut
#
#sub write {
#	my $self	=	shift;
#	my $file	=	shift;
#	my $sections=	shift;
#
#	
#	
#	
#	#### SANITY CHECK
#	
#	#### KEEP COPIES OF FILES IF backup DEFINED
#	$self->makeBackup($file) if defined $self->get_backup();
#	
#	open(OUT, ">$file") or die "Can't open file: $file\n";
#	foreach my $section ( @$sections )
#	{
#		my $key = $section->{key};
#		my $value = $section->{value};
#		
#		
#	}
#}
#
#
#
#sub addKeypair {
#	my $self	=	shift;
#	my $name 	=	shift;
#	my $location=	shift;
#	my $comment =	shift;
#	
#	my $header = "[key $name]";
#	
#
#	#### SANITY CHECK
#	
#
#}
#
#




##################			HOUSEKEEPING SUBROUTINES			################
#=head2
#
#	SUBROUTINE		new
#	
#	PURPOSE
#
#		INITIALISE OUR OBJECT AND VALIDATE INPUTS
#	
#=cut
#after 'new' => sub  {
#    my $self		=	shift;
#
#	validated_list(
#		\@_,
#		ignore		=>	{ isa	=>	'Str', 		optional	=>	1 },
#		backup		=>	{ isa	=>	'Str', 		optional	=>	1 },
#		spacer		=>	{ isa	=>	'Str', 		optional	=>	1 },
#		separator	=>	{ isa	=>	'Str', 		optional	=>	1 },
#		inputfile	=>	{ isa	=>	'Str', 		optional	=>	1 },
#		hash		=>	{ isa	=>	'HashRef',	optional	=>	1 },
#		outputfile	=>	{ isa	=>	'Str', 		optional	=>	1 },
#		command		=>	{ isa	=>	'Str', 		optional	=>	1 }
#	);
#};
#
##### DUMPER
#sub dump { 
#    my $self = shift;
#
#    require Data::Dumper;
#    $Data::Dumper::Maxdepth = shift if @_;
#}
#
#
=head1 LICENCE

This code is released under the GPL, a copy of which should
be provided with the code.

=end pod

=cut

1;

__END__

=head2

	SUBROUTINE 		_replaceKeyValue

	PURPOSE

		1. REPLACE AN EXISTING KEY-VALUE PAIR IN A SECTION

		2. RETAIN ANY EXISTING COMMENTS IF comment NOT DEFINED

		3. REPLACE COMMENTS IF comment DEFINED

=cut

sub _replaceKeyValue {
	my $self	=	shift;
	my $section	=	shift;
	my $key		=	shift;
	my $value	=	shift;
	my $comment =	shift;

	print "Conf::Agua::_replaceKeyValue    Conf::Agua::_replaceKeyValue(section, key, value)\n";	
	print "Conf::Agua::_replaceKeyValue    section: $section\n" if defined $section;
	print "Conf::Agua::_replaceKeyValue    key: $key\n";
	print "Conf::Agua::_replaceKeyValue    value: $value\n";

	#### SANTIY CHECK
	print "Conf::Agua::_replaceKeyValue    section not defined\n" and exit if not defined $section;
	print "Conf::Agua::_replaceKeyValue    key not defined\n" and exit if not defined $key;
	print "Conf::Agua::_replaceKeyValue    value not defined\n" and exit if not defined $value;

	$section->{keypairs}->{$key}->{value} = $value;
	$section->{keypairs}->{$key}->{comment} = $self->get_comment() . " " . $comment if defined $comment;
}



