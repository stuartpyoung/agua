package Agua::Common::Transport;
use Moose::Role;
use Moose::Util::TypeConstraints;

=head2

	PACKAGE		Agua::Common::Transport

	PURPOSE

		TRANSPORT METHODS FOR Agua::Common

=cut


use Data::Dumper;

=head2

	SUBROUTINE		jsonSafe

	PURPOSE

		SUBSTITUTE OUT CHARACTERS THAT WOULD BREAK JSON PARSING

		WITH SAFE SYMBOL SETS, E.G., &quot; INSTEAD OF '"'

=cut

sub jsonSafe {
	my $self		=	shift;
	my $string		=	shift;
	my $mode		=	shift;



	#### CHECK MODE IS DEFINED
	die "Agua::Common::Transport::jsonSafe    mode not defined. Exiting.\n" if not defined $mode;


	#### SANITY CHECKS
	if ( not defined $string or not $string )
	{
		print "Agua::Common::Transport::jsonSafe    String not defined or empty. Returning ''\n";
		return ''; 
	}

	my $specialChars = [
		[ '&quot;',	"'" ],	
		[ '&quot;',	'"' ],	
		[ '&nbsp;', ' ' ],	
		#[ '&#35;', '#' ],	
		#[ '&#36;', '$' ],	
		#[ '&#37;', '%' ],	
		#[ '&amp;', '&' ],	
		#[ '&#39;', "'" ],	
		[ '&#40;', '\(' ],	
		[ '&#41;', '\)' ],	
		#[ '&frasl;', '\/' ],	
		[ '&#91;', '\[' ],	
		#[ '&#92;', '\\\\' ],	
		[ '&#93;', '\]' ],	
		#[ '&#96;', '`' ],	
		[ '&#123;', '\{' ],	
		#[ '&#124;', '|' ],	
		[ '&#125;', '\}' ]	
	];

	#### REMOVE LINE RETURNS AND '\s' 
	$string =~ s/\s/ /gms;
	$string =~ s/\n/\\n/gms;


	for ( my $i = 0; $i < @$specialChars; $i++)
	{

		if ( $mode eq 'toJson' )
		{
			no re 'eval';
			$string =~ s/$$specialChars[$i][1]/$$specialChars[$i][0]/msg;
			use re 'eval';
		}
		elsif ( $mode eq 'fromJson' )
		{
			no re 'eval';
			$string =~ s/$$specialChars[$i][0]/$$specialChars[$i][1]/msg;
			use re 'eval';
		}
	}


	return $string;	
}





=head2

	SUBROUTINE		downloadFile

	PURPOSE

		SEND A FILE DOWNLOAD

	INPUT

		1. USER VALIDATION: USERNAME, SESSION ID

		2. FILE PATH

	OUTPUT

		1. FILE DOWNLOAD Content-type: application/x-download

=cut

sub downloadFile {
	my $self		=	shift;

	my $json		=	$self->json();
	my $username	=	$json->{username};
	my $filepath	=	$json->{filepath};
	my $requestor	=	$json->{requestor};



	my $fileroot	=	$self->getFileroot($username);	

	#### ADD FILEROOT IF FILEPATH IS NOT ABSOLUTE
	if ( $filepath !~ /^\// )
	{
		#### GET THE FILE ROOT FOR THIS USER OR ANOTHER
		my $fileroot;

		if ( defined $requestor )
		{	
			#### GET THE PROJECT NAME, I.E., THE BASE DIRECTORY OF THE FILE NAME
			my ($project) = ($filepath) =~ /^([^\/^\\]+)/;

			#### EXIT IF PROJECT ACCESS NOT ALLOWED
			if ( not $self->_canAccess($username, $project, $requestor, "project") )
			{
				print "{ error: 'Agua::Common::Transport::downloadFiles    user $requestor is not permitted access to project $project owned by user $username' }";
				exit;
			}

			#### ELSE USE OTHER USER'S FILEROOT
			$fileroot = $self->getFileroot($username);
			print "{ error: 'Agua::Common::Transport::downloadFiles    Fileroot not found for username: $username' }" and exit unless defined $fileroot;

		}
		else
		{
			$fileroot = $self->getFileroot($username);
			print "{ error: 'Agua::Common::Transport::downloadFiles    Fileroot not found for user' }" and exit if not defined $fileroot;
		}

		$filepath = "$fileroot/$filepath";
		if ( $^O =~ /^MSWin32$/ )   {   $filepath =~ s/\\/\//g;  }
		if ( $^O =~ /^MSWin32$/ )   {   $filepath =~ s/\//\\/g;  }
		#if ( $^O =~ /^MSWin32$/ )   {   $filepath =~ s/\\/\//g;  }

	}

	#### GET FILE SIZE
	#my $filesize = (stat($filepath))[7] or die "Can't get size of file: $filepath\n";
	print "{ error: 'Cannot find file: $filepath' }\n"
		and exit if not -s $filepath;
	my $filesize = -s $filepath
		or print "{ error: 'Cannot get size of file: $filepath' }\n" and exit;

	#### GET FILE NAME
	my ($filename) = $filepath =~ /([^\/]+)$/;

	#### EXIT IF FILE IS EMPTY
	print "{ error : 'File not found: $filepath'}" and exit if not -f $filepath;
	print "{ error : 'File is empty: $filepath'}" and exit if -z $filepath;

	#### PRINT DOWNLOAD HEADER
	print qq{Content-type: application/x-download\n};
	print qq{Content-Disposition: attachment;filename="$filename\n};
	print qq{Content-Length: $filesize\n\n};

	open(FILE, $filepath) or die "Can't open file: $filepath\n";
	binmode FILE;
	$/ = undef;
	print <FILE>;
	close(FILE);
}




1;