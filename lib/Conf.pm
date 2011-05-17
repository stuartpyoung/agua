package Conf;

use Data::Dumper;


=head2

	PACKAGE		ConfigFile

    PURPOSE

        1. READ AND WRITE ini-FORMAT CONFIGURATION FILES:

			[section name]
			KEY<SPACER>VALUE

    NOTES

			PRINT USER CONFIG FILE EXAMPLES:

				admin       /home/admin/.starcluster/config
				jgilbert    /home/jgilbert/.starcluster/config

=cut

use Moose::Role;
#with 'MooseX::Getopt';
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

=head2

	SUBROUTINE 		read

	PURPOSE

		1. IF NO section IS PROVIDED, RETURN A HASH OF {section => {KEY-VALUE PAIRS}} PAIRS

		2. IF NO key IS PROVIDED, RETURN A HASH OF KEY-VALUE PAIRS FOR THE section

		3. RETURN THE KEY-VALUE PAIR FOR THE section AND key

		4. RETURN IF THE GIVEN section IS NOT IN THE ini FILE

		5. RETURN IF THE GIVEN key IS NOT IN THE ini FILE

=cut

sub read {
	my $self	=	shift;
	my $file	=	shift;

	$file = $self->inputfile() if not defined $file;	

	$file = $self->outputfile() if not defined $file;
	return if not defined $file;


	#### SKIP READ IF FILE ABSENT
	my $sections = [];
	my ($dir) = $file =~ /^(.+)\/([^\/]+)$/;
	if ( not -d $dir )
	{
		return $sections;	
	}
	return $sections if not -f $file;

	my $old = $/;
	$/ = undef;
	open(FILE, $file) or die "Conf::read    Can't open file: $file\n";
	my $contents = <FILE>;
	close(FILE) or die "Conf::read    Can't close file: $file\n";

	my @sections = split "\\n\\[", $contents;
	shift @sections if not $sections[0] =~ s/^\s*\[//;

	my $counter = 0;
	foreach my $section ( @sections )
	{
		$counter++;

		$section =~ /^(.+?)(\s+.+)?\s*\]/;
		my $section_key = $1;
		my $section_value = $2;
		$section_value =~ s/^\s+// if defined $section_value;
		my @lines = split "\n", $section;
		shift @lines;
		my $section_hash;
		$section_hash->{key} = $section_key;
		$section_hash->{value} = $section_value;

		my $line_comment = '';
		foreach my $line ( @lines )
		{
			my ($key, $value, $comment) = $self->parseline($line);

			#### STORE COMMENT IF PRESENT
			$line_comment .= "$comment\n" if defined $comment;

			#### ADD KEY-VALUE PAIR TO SECTION HASH
			if ( defined $key )
			{
				$section_hash->{keypairs}->{$key}->{value} = $value ;

				#### ADD COMMENT IF PRESENT
				$line_comment =~ s/\s+$//;
				$section_hash->{keypairs}->{$key}->{comment} = $line_comment if $line_comment;
				$line_comment = '';
			}
		}
		push @$sections, $section_hash;
	}

	$/ = $old;

	return $sections;
}

sub copy {
	my $self	=	shift;
	my $file	=	shift;


	$self->outputfile($file);
	$self->write();	
};

sub write {
	my $self	=	shift;
	my $sections=	shift;
	my $file	=	shift;


	$file = $self->outputfile() if not defined $file;	
	$file = $self->inputfile() if not defined $file;	

	$sections = $self->_getSections() if not defined $sections;

	#### SANITY CHECK
	print "Conf::write    file not defined\n" and exit if not defined $file;
	print "Conf::write    sections not defined\n" and exit if not defined $sections;

	#### KEEP COPIES OF FILES IF backup DEFINED
	$self->makeBackup($file) if defined $self->backup();

	#### CREATE FILE DIR IF NOT EXISTS	
	my ($dir) = $file =~ /^(.+)\/([^\/]+)$/;
	File::Path::mkpath($dir) if not -d $dir;
	print "Conf::write    Can't create dir: $dir\n" if not -d $dir;

	open(OUT, ">$file") or die "Conf::write    Can't open file: $file\n";
	foreach my $section ( @$sections )
	{
		my $section_key = $section->{key};
		my $section_value = $section->{value};

		my $header = "[" . $section_key;
		$header .= " $section_value" if defined $section_value;
		$header .= "]\n";
		print OUT $header;

		my $keypairs = $section->{keypairs};

		my @keys = keys %$keypairs;
		@keys = sort @keys;
		foreach my $key ( @keys )
		{
			my $comment = $keypairs->{$key}->{comment};
			my $value = $keypairs->{$key}->{value};
			my $entry = '';
			$entry .= "\n" . $comment . "\n" if defined $comment;
			$entry .= $key . $self->separator() . $value . "\n";
			print OUT $entry;
		}
		print OUT "\n";
	}
	close(OUT) or die "Can't close file: $file\n";

};

=head2

	SUBROUTINE 		parseline

	PURPOSE

		1. PARSE A KEY-VALUE PAIR FROM LINE

		2. RETURN A COMMENT AS A THIRD VALUE IF PRESENT

=cut

sub parseline {
	my $self	=	shift;
	my $line 	=	shift;


	my $ignore	=	$self->ignore();
	my $spacer	=	$self->spacer();

	return (undef, undef, $line) if $line =~ /^[$ignore]+/;

	my ($key, $value) = $line =~ /^(.+?)[$spacer\$]+(.+)\s*$/;

	#### DEFINE VALUE AS '' FOR BOOLEAN key (I.E., ITS A FLAG)
	$value = '' if not defined $value;

	return ($key, $value);	
};

=head2

=head3 C<SUBROUTINE 		conf>

	PURPOSE

		RETURN A SIMPLE HASH OF CONFIG KEY-VALUE PAIRS

=cut

sub conf {
	my $self	=	shift;
	my $file	=	shift;

	$file = $self->inputfile() if not defined $file;
	return if not defined $file;

	my $hash = {};
	open(FILE, $file) or die "Can't open file: $file\n";
	my @lines = <FILE>;
	close(FILE) or die "Can't close file: $file\n";
	foreach my $line ( @lines )
	{
		my ($key, $value) = $self->parseline($line);

		$hash->{$key} = $value if defined $key;
	}

	return $hash;
};


=head2

	SUBROUTINE 		makeBackup

	PURPOSE

		COPY FILE TO NEXT NUMERICALLY-INCREMENTED BACKUP FILE

=cut
sub makeBackup {
	my $self	=	shift;
	my $file	=	shift;

	print "Conf::makeBackup    file not defined\n" and exit if not defined $file;

	#### BACKUP FSTAB FILE
	my $counter = 1;
	my $backupfile = "$file.$counter";
	while ( -f $backupfile )
	{
		$counter++;
		$backupfile = "$file.$counter";
	}

	require File::Copy;
	File::Copy::copy($file, $backupfile);
};

=head2

	SUBROUTINE 		_getSections

	PURPOSE

		RETURN ALL SECTIONS IN THE INPUTFILE

=cut

sub _getSections {
	my $self	=	shift;


	my $sections = $self->read();
	$self->sections($sections);

	return $sections;
}

=head2

	SUBROUTINE 		_getSection

	PURPOSE

		1. RETURN A SECTION WITH THE GIVEN NAME (AND VALUE IF PROVIDED)

=cut

sub _getSection {
	my $self	=	shift;
	my $sections=	shift;
	my $name	=	shift;
	my $value	=	shift;


	#### SANITY CHECK
	print "Conf::_getSection    sections not defined\n" and exit if not defined $sections;
	print "Conf::_getSection    name not defined\n" and exit if not defined $name;

	#### SEARCH FOR SECTION AND ADD KEY-VALUE PAIR IF FOUND	
	foreach my $section ( @$sections )
	{

		next if not $section->{key} eq $name;
		next if defined $value and defined $section->{value}
			and not $section->{value} eq $value;


		return $section;
	}

	return {};
}

=head2

	SUBROUTINE 		_addSection

	PURPOSE

		1. ADD A SECTION

		2. ADD THE SECTION'S VALUE IF DEFINED

=cut

sub _addSection {
	my $self	=	shift;
	my $sections=	shift;
	my $name	=	shift;
	my $value	=	shift;


	#### SANTIY CHECK
	print "Conf::_addSection    sections not defined\n" and exit if not defined $sections;
	print "Conf::_addSection    name not defined\n" and exit if not defined $name;

	$value = '' if not defined $value;
	my $section;
	$section->{key} = $name;
	$section->{value} = $value;
	push @$sections, $section;

	return $section;
}


=head2

	SUBROUTINE 		_removeKeyValue

	PURPOSE

		1. REMOVE A KEY-VALUE PAIR FROM A SECTION

=cut

sub _removeKeyValue {
	my $self	=	shift;
	my $section	=	shift;
	my $key		=	shift;


	#### SANTIY CHECK
	print "Conf::_removeKeyValue    section not defined\n" and exit if not defined $section;
	print "Conf::_removeKeyValue    key not defined\n" and exit if not defined $key;

	delete $section->{keypairs}->{$key};
}





=head2

	SUBROUTINE 		_addKeyValue

	PURPOSE

		1. ADD A KEY-VALUE PAIR TO A SECTION

		2. RETAIN ANY EXISTING COMMENTS IF comment NOT DEFINED

		3. ADD COMMENTS IF comment DEFINED

=cut

sub _addKeyValue {
	my $self	=	shift;
	my $section	=	shift;
	my $key		=	shift;
	my $value	=	shift;
	my $comment =	shift;


	#### SANITY CHECK
	print "Conf::_addKeyValue    section not defined\n" and exit if not defined $section;
	print "Conf::_addKeyValue    key not defined\n" and exit if not defined $key;
	print "Conf::_addKeyValue    value not defined\n" and exit if not defined $value;

	$section->{keypairs}->{$key}->{value} = $value;
	$section->{keypairs}->{$key}->{comment} = $self->comment() . " " . $comment if defined $comment;
}

=head2

	SUBROUTINE 		_getKeyValue

	PURPOSE

		1. ADD A KEY-VALUE PAIR TO A SECTION

		2. RETAIN ANY EXISTING COMMENTS IF comment NOT DEFINED

		3. ADD COMMENTS IF comment DEFINED

=cut

sub _getKeyValue {
	my $self	=	shift;
	my $section	=	shift;
	my $key		=	shift;


	#### SANITY CHECK
	print "Conf::_getKeyValue    section not defined\n" and exit if not defined $section;
	print "Conf::_getKeyValue    key not defined\n" and exit if not defined $key;


	return $section->{keypairs}->{$key}->{value};
}


=head2

	SUBROUTINE 		_sectionIsEmpty

	PURPOSE

		1. RETURN 1 IF SECTION CONTAINS NO KEY-VALUE PAIRS

		2. OTHERWISE, RETURN 0

=cut

sub _sectionIsEmpty {
	my $self	=	shift;
	my $section =	shift;


	#### SANITY CHECK
	print "Conf::_sectionIsEmpty    section not defined\n" and exit if not defined $section;

	return 1 if not defined $section->{keypairs};
	my @keys = keys %{$section->{keypairs}};
	return 1 if not @keys;

	return 0;
}


=head2

	SUBROUTINE 		_removeSection

	PURPOSE

		1. REMOVE A SECTION FROM THE CONFIG OBJECT

=cut

sub _removeSection {
	my $self	=	shift;
	my $sections=	shift;
	my $section =	shift;


	#### SANITY CHECK
	print "Conf::_removeSection    sections not defined\n" and exit if not defined $sections;
	print "Conf::_removeSection    section not defined\n" and exit if not defined $section;

	#### SEARCH FOR SECTION AND REMOVE IF FOUND	
	for ( my $i = 0; $i < @$sections; $i++ )
	{
		my $current_section = $$sections[$i];
		next if not $current_section->{key} eq $section->{key};
		next if defined $section->{value} and not $current_section->{value} eq $section->{value};
		splice(@$sections, $i, 1);
		return 1;
	}

	return 0;
}




sub getKeyValue {
	my $self			=	shift;
	my $section_keypair	=	shift;
	my $key				=	shift;

	my ($name, $value) = $section_keypair =~ /^([^:]+):?(.*)$/;
	print "Conf::getKeyValue    name not defined\n" and exit if not defined $name;
	print "Conf::getKeyValue    key not defined\n" and exit if not defined $key;


	my $sections = $self->read();

	my $section = $self->_getSection($sections, $name, $value);

	return $self->_getKeyValue($section, $key);
}

=head2

	SUBROUTINE 		addKeyValue

	PURPOSE

		1. ADD A SECTION AND KEY-VALUE PAIR TO THE CONFIG

			FILE IF THE SECTION DOES NOT EXIST

		2. OTHERWISE, ADD A KEY-VALUE PAIR TO AN EXISTING SECTION

		3. REPLACE KEY ENTRY IF IT ALREADY EXISTS IN THE SECTION

=cut

sub addKeyValue {
	my $self	=	shift;
	my $section_name	=	shift;
	my $section_value	=	shift;
	my $key		=	shift;
	my $value	=	shift;
	my $comment	=	shift;


	#### SET SECTION VALUE TO UNDEF IF EMPTY
	$section_value = undef if defined $section_value and not $section_value;


	#### SANITY CHECK
	print "Conf::addKeyValue    section_name not defined\n" and exit if not defined $section_name;
	print "Conf::addKeyValue    key not defined\n" and exit if not defined $key;
	print "Conf::addKeyValue    value not defined for key: $key\n" and exit if not defined $value;

	#### SEARCH FOR SECTION AND ADD KEY-VALUE PAIR IF FOUND	
	my $sections = $self->_getSections();

	my $matched = 0;
	foreach my $current_section ( @$sections )
	{
		next if not $current_section->{key} eq $section_name;
		next if defined $section_value and not $current_section->{value} eq $section_value;
		$matched = 1;
		$self->_addKeyValue($current_section, $key, $value, $comment);
		last;
	}

	#### OTHERWISE, CREATE THE SECTION AND ADD THE 
	#### KEY-VALUE PAIR TO IT
	if ( not $matched )
	{
		my $section = $self->_addSection($sections, $section_name, $section_value);
		$self->_addKeyValue($section, $key, $value, $comment);
	}


	$self->write($sections);
}



=head2

	SUBROUTINE 		removeKeyValue

	PURPOSE

		1. ADD A SECTION AND KEY-VALUE PAIR IF SECTION DOES NOT EXIST

		2. ADD A KEY-VALUE PAIR TO AN EXISTING SECTION

		3. REPLACE KEY ENTRY IF IT ALREADY EXISTS IN THE SECTION

=cut

sub removeKeyValue {
	my $self		=	shift;
	my $section_name	=	shift;
	my $section_value	=	shift;
	my $key			=	shift;
	my $value		=	shift;

	#### SET SECTION VALUE TO UNDEF IF EMPTY
	$value = undef if defined $value and not $value;


	#### SANITY CHECK
	print "Conf::removeKeyValue    section_name not defined\n" and exit if not defined $section_name;
	print "Conf::removeKeyValue    key not defined\n" and exit if not defined $key;
	print "Conf::removeKeyValue    value not defined\n" and exit if not defined $value;

	#### SEARCH FOR SECTION AND ADD KEY-VALUE PAIR IF FOUND	
	my $sections = $self->read();
	my $matched = 0;
	foreach my $current_section ( @$sections )
	{
		next if not $current_section->{key} eq $section_name;
		next if defined $section_value and not $current_section->{value} eq $section_value;

		$matched = 1;
		$self->_removeKeyValue($current_section, $key, $value);

		#### REMOVE SECTION IF EMPTY
		if ( $self->_sectionIsEmpty($current_section) )
		{
			my $success = $self->_removeSection($sections, $current_section);
		}
		last;
	}


	$self->write($sections);
}

##################			HOUSEKEEPING SUBROUTINES			################
#=head2
#
#	SUBROUTINE		new
#	
#	PURPOSE
#
#		INITIALISE OUR OBJECT AND VALIDATE INPUTS
#=cut
#
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
#
#### DUMPER
sub dump { 
    my $self = shift;

    require Data::Dumper;
    $Data::Dumper::Maxdepth = shift if @_;
    print Data::Dumper::Dumper $self;
}


=head1 LICENCE

This code is released under the GPL, a copy of which should
be provided with the code.

=end pod

=cut

1;
