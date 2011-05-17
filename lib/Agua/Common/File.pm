package Agua::Common::File;
use Moose::Role;
use Moose::Util::TypeConstraints;

=head2

	PACKAGE		Agua::Common::File

	PURPOSE

		FILE CHECKING AND INFO METHODS FOR Agua::Common

=cut


use Data::Dumper;

=head2

	SUBROUTINE		fileIntegrity

	PURPOSE

		THIS METHOD HAS TWO USES:

		1. DETERMINE WHETHER A FILE BELONGS TO A USER-DEFINED FILE TYPE

			(E.G., FASTA, FASTQ, BINARY FASTA, SQUASH, .ACE, .TSV, .CSV)

		2. IF IT DOES, DETERMINE WHETHER IT CONFORMS TO THE SPECIFICATIONS

			OF THAT FILE TYPE

=cut

sub fileIntegrity {
	#### PLACEHOLDER
}

=head2

	SUBROUTINE		fileType

	PURPOSE

		DETERMINE IF A FILE BELONGS TO ONE OF A LIST OF PREDETERMINED

		FILE TYPES (E.G., FASTA, FASTQ, BINARY FASTA, SQUASH, .ACE, .TSV, .CSV)

=cut

sub fileType {
    #### PLACEHOLDER	
}

=head2

	SUBROUTINE		filePeek

	PURPOSE

		PRINT LINES OF A FILE FROM A PARTICULAR OFFSET

=cut

sub filePeek {
	my $self		=	shift;

	my $filepath = $self->json()->{filepath};
	my $offset = $self->json()->{offset};
	my $lines = $self->json()->{lines};
	my $bytes = $self->json()->{bytes};

	return "{ error: 'Agua::Common::File::filePeek    offset not defined' }" and exit if not defined $offset;
	return "{ error: 'Agua::Common::File::filePeek    filepath not defined' }" and exit if not defined $filepath;
	return "{ error: 'Agua::Common::File::filePeek    lines and bytes not defined' }" and exit if not defined $lines and not defined $bytes;

	#### EXIT IF FILEPATH NOT FOUND
	if ( not -f $filepath and not -d $filepath )
	{
		return "{  error: 'Agua::Common::File::filepeek    file or directory not found: $filepath' }";
	}

	#### GET LINES 
	if ( defined $lines and defined $offset )
	{
		my $peek = `tail +$offset $filepath | head -n$lines`;
		if ( not defined $peek or not $peek )
		{
			return "{  error: 'Agua::Common::File::filepeek    peek at offset $offset and lines $lines not defined for filepath: $filepath' }";
		}		
		$peek =~ s/\n/\\n/g;
		$peek =~ s/'/\\'/;

		return $peek;
	}


	#### GET BYTES
	elsif ( defined $bytes and defined $offset )
	{
		my $peek = '';
		open(FILE, $filepath) or print "{ error: 'Agua::Common::File::filePeek    Could not find file: $filepath' }" and exit;
		seek FILE, $offset, 0;
		read(FILE, $peek, $bytes);
		close(FILE);

		if ( not defined $peek or not $peek )
		{
			return "{  error: 'Agua::Common::File::filepeek    peek at offset $offset and bytes $bytes not defined for filepath: $filepath' }";
		}

		print "{ peek: '$peek' }";
		exit;
	}
}


=head2

	SUBROUTINE		checkFiles

	PURPOSE

		PRINT FILE INFO FOR ONE OR MORE FILES OR DIRECTORIES

	INPUTS

		AN ARRAY OF FILE HASHES
		[
			{ filehash1 }, 
			{ filehash2 },
			...
		]

	OUTPUTS

		AN ARRAY (SAME ORDER AS INPUTS) OF FILE INFORMATION HASHES
		[
			{ fileinfo1 },
			{ fileinfo2 },
			...
		]		
=cut

sub checkFiles {
	my $self		=	shift;


    ### VALIDATE    
    print "{ error: 'Agua::Common::File::checkFiles    User session not validated' }" and exit unless $self->validate();

	#### GET ARRAY OF FILE HASHES
	my $files = $self->json()->{files};

	#### GET FILE INFO FOR EACH FILE/DIRECTORY
	my $fileinfos = [];
	foreach my $file ( @$files )
	{
		use Data::Dumper;
		my $filepath = $file->{value};

		push@$fileinfos, $self->_checkFile($filepath);	
	}

	#### PRINT FILE INFO	
	$self->printJSON($fileinfos);
	exit;
}


=head2

	SUBROUTINE		checkFile

	PURPOSE

		PRINT FILE INFO FOR A FILE OR DIRECTORY

=cut

sub checkFile {
	my $self		=	shift;


    #### GET JSON
    my $json  =	$self->json();

    #### VALIDATE    
    print "{ error: 'Agua::Common::File::checkFile    User session not validated' }" and exit unless $self->validate();

	my $fileinfo = $self->_checkFile();

	$self->printJSON($fileinfo);
	exit;
}

=head2

	SUBROUTINE		_checkFile

	PURPOSE

		RETURN A HASH CONTAINING FILE INFORMATION (EXISTS, SIZE,

		DATE MODIFIED, ETC.) FOR ONE OR MORE FILES AND DIRECTORIES

	INPUTS

		A FILE PATH STRING

	OUTPUTS

		A FILE INFORMATION HASH:
		{
			info1 : '...',
			info2 : '...',
			...
		}

		OR AN ERROR HASH:
		{
			error : '...'+
		}

	NOTES

		IF THE FILEPATH IS ABSOLUTE, THE USER SHOULD HAVE THE CORRECT PERMISSIONS
		TO ACCESS THE DATA BASED ON THE SCRIPT'S SET UID AND SET GID.

		IF THE requestor FIELD IS DEFINED, THIS MEANS ITS A REQUEST BY THE USER
		requestor TO VIEW FILE INFO IN PROJECTS OWNED BY THE USER username.

=cut

sub _checkFile {
	my $self		=	shift;
	my $filepath	=	shift;



    #### GET JSON
    my $json  =	$self->json();

    ### VALIDATE    
    print "{ error: 'Agua::Common::File::_checkFile    User session not validated' }" and exit unless $self->validate();

	#### GET FILE PATH AND USERNAME
	$filepath = $json->{filepath} if not defined $filepath;
	my $username = $json->{username};

	#### RETURN FILE NOT EXISTS IF FILEPATH IS EMPTY
	if ( $filepath =~ /^\s*$/)
	{
		$self->printJSON({ filepath => "$filepath", exists => "false"});
		exit;
	}

	#### IF THE FILEPATH IS ABSOLUTE, THE USER SHOULD HAVE THE CORRECT PERMISSIONS
	#### TO ACCESS THE DATA BASED ON THE SCRIPT'S SET UID AND SET GID
	if ( $filepath =~ /^\// )
	{
		if ( $^O =~ /^MSWin32$/ )   {   $filepath =~ s/\//\\/g;  }

		return $self->fileinfo($filepath);
	}

	#### ADD FILEROOT IF FILEPATH IS NOT ABSOLUTE.
	#### GET THE FILE ROOT FOR THIS USER IF requestor NOT DEFINED
	if ( not defined $json->{requestor} )
	{	

		#### DO SETUID OF WORKFLOW.CGI AS USER
		my $fileroot = $self->getFileroot($username);

		#### QUIT IF FILEROOT NOT DEFINED OR EMPTY
		return {error => 'Agua::Common::File::_checkFile    fileroot not defined for user: $username and filepath: $filepath' } if not defined $fileroot or not $fileroot;

		#### SET NEW FILEPATH
		$filepath = "$fileroot/$filepath";

		if ( $^O =~ /^MSWin32$/ )   {   $filepath =~ s/\//\\/g;  }

		return $self->fileinfo($filepath);
	}

	#### IF THE requestor FIELD IS DEFINED, THIS MEANS
	#### ITS A REQUEST BY THE USER requestor TO VIEW FILE INFO IN PROJECTS
	#### OWNED BY THE USER username.

	#### (THE CLIENT MAKES SEPARATE CALLS FOR SHARED FILES NOT OWNED
	#### BY THE LOGGED-IN USER WITH requestor=THIS USER AND username=OTHER USER, 
	#### SO THAT THE INSTANCE IS RUN AS SETUID THE OTHER USER AND THEREBY HAS
	#### ACCESS TO THE OTHER USER'S PROJECT DIRECTORIES

	my $requestor = $json->{requestor};
	my $fileroot = $self->getFileroot($username);
	return { error => 'Agua::Common::File::_checkFile    Fileroot not found for user' } if not defined $fileroot;

	#### GET THE PROJECT NAME, I.E., THE BASE DIRECTORY OF THE FILE NAME
	my ($project) = ($filepath) =~ /^([^\/^\\]+)/;
	my $can_access_project = $self->_canAccess($username, $project, $requestor, "project");

	#### IF PROJECT RIGHTS NOT DEFINED, USE DEFAULT PROJECT RIGHTS 
	if ( not $can_access_project )
	{
	}
	else
	{
		$fileroot = $self->getFileroot($username);
		return { error => 'Agua::Common::File::_checkFile    Fileroot not found for user' } if not defined $fileroot;

		#### QUIT IF FILEROOT NOT DEFINED OR EMPTY
		return  { error =>  "Agua::Common::File::_checkFile    fileroot not defined for user: $username and filepath: $filepath" } if not defined $fileroot or not $fileroot;

		$filepath = "$fileroot/$filepath";
		if ( $^O =~ /^MSWin32$/ )   {   $filepath =~ s/\\/\//g;  }
		if ( $^O =~ /^MSWin32$/ )   {   $filepath =~ s/\//\\\\/g;  }

		my ($file_manager) = $0;
		$file_manager =~ s/[^\/]+$/filemanager.pl/;
		my $jsonText = qq{{"filepath":"$filepath","mode":"checkFile"}};
		my $command = "echo $jsonText | C:\\perl\\bin\\perl $file_manager";

		my $fileinfo = `$command`;

		print $fileinfo;
		exit;
	}	
}


=head2

	SUBROUTINE		printJSON

	PURPOSE

		PRINT JSON FOR AN OBJECT USING THE FOLLOWING FLAGS:

			allow_nonref

=cut

sub printJSON {
	my $self		=	shift;
	my $data		=	shift;

	#### PRINT JSON AND EXIT
	my $jsonParser = JSON->new();
	#my $jsonText = $jsonParser->encode->allow_nonref->($data);
    my $jsonText = $jsonParser->objToJson($data, {pretty => 1, indent => 4});
	print $jsonText;
}



=head2

	SUBROUTINE		fileinfo

	PURPOSE

		RETURN FILE STATUS: EXISTS, SIZE, CREATED

=cut

sub fileinfo {
	my $self		=	shift;
	my $filepath		=	shift;


	if ( $^O =~ /^MSWin32$/ )   {   $filepath =~ s/\//\\/g;  }

	use File::stat;
	use Time::localtime;

	my $fileinfo;
	$fileinfo->{filepath} = $filepath;
	$fileinfo->{exists} = "false";
	return $fileinfo if $filepath =~ /^\s*$/;
	return $fileinfo if not -f $filepath and not -d $filepath;

    use File::stat;
    my $sb = stat($filepath);
	$fileinfo->{exists} = "true";
	$fileinfo->{size} = $sb->size;
	$fileinfo->{permissions} = $sb->mode & 0777;
	$fileinfo->{modified} = ctime($sb->mtime);
	#$fileinfo->{modified} = scalar localtime $sb->mtime;
	$fileinfo->{type} = "file";
	$fileinfo->{type} = "directory" if -d $filepath;

    return $fileinfo;
}


=head2

    SUBROUTINE:     fileStats

    PURPOSE:        Return the following file statistics:
                    -   size (bytes)
                    -   directory ("true"|"false")
                    -   modified (seconds Unix time)

=cut

sub fileStats {
    my $self        =   shift;
    my $filepath    =   shift;


    my $fileStats;

    my $filesize = -s $filepath;
    if ( not defined $filesize )
    {
        $filesize = '';
    }

    my $directory = -d $filepath;
    if ( not -f $filepath and not -d $filepath )
    {
        $directory = "''";
    }
    elsif ( not $directory )
    {
        $directory = "false";
    }
    else
    {
        $directory = "true";
    }

    my $modified = -M $filepath;
    if ( not defined $modified )
    {
        $modified = "''"
    }
    else
    {
        $modified = int(time() - $modified * 24 * 60);
    }

    $fileStats->{filesize} = "$filesize";
    $fileStats->{directory} = "$directory";
    $fileStats->{modified} = "$modified";

    return $fileStats;
}



1;