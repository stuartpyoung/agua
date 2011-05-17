package Sampler;

use strict;

#### DEBUG

=head2

	PACKAGE		Sampler

    VERSION:        0.01

    PURPOSE

        UTILITY FUNCTIONS TO SAMPLE READS FROM A COLLECTION

		OF FILES:

			1. fastaInfos.pl CREATE FIXED-WIDTH FASTA, QUAL AND INFO FILES

			2. listFiles.pl CREATE N listFiles PER ENTRY SEQUENCE CONTAINING

				BYTE-POSITIONS FOR RECORDS IN THE INPUT FASTQ FILE

			3. printFiles.pl PRINT FASTA, QUAL AND FASTQ FILES FOR EACH listFile

				OF EACH FASTQ FILE IN THE INPUT DIRECTORY

	EXAMPLES

/nethome/syoung/base/bin/comparison/fastaInfos.pl \
--inputdir /mihg/data/NGS/syoung/base/pipeline/SRA/test \
--outputdir /mihg/data/NGS/syoung/base/pipeline/SRA/test/fasta \
--fixed \
--id_length 21 \
--sequence_length 36 \
--dot 10 

/nethome/syoung/base/bin/comparison/listFiles.pl \
--inputdir /mihg/data/NGS/syoung/base/pipeline/SRA/test/fasta \
--outputdir /mihg/data/NGS/syoung/base/pipeline/SRA/test/fasta \
--size 10000


/nethome/syoung/base/bin/comparison/printFiles.pl \
--inputdir /mihg/data/NGS/syoung/base/pipeline/SRA/test/fasta \
--outputdir /mihg/data/NGS/syoung/base/pipeline/SRA/test/fasta/samples \
--mode fasta \
--cluster \
--queue "-q psmall" 

/nethome/syoung/base/bin/comparison/printFiles.pl \
--inputdir /mihg/data/NGS/syoung/base/pipeline/SRA/test/fasta \
--outputdir /mihg/data/NGS/syoung/base/pipeline/SRA/test/fasta/samples \
--mode qual \
--cluster \
--queue "-q psmall" 

/nethome/syoung/base/bin/comparison/printFiles.pl \
--inputdir /mihg/data/NGS/syoung/base/pipeline/SRA/test/fasta \
--outputdir /mihg/data/NGS/syoung/base/pipeline/SRA/test/fasta/samples \
--mode fastq \
--cluster \
--queue "-q psmall" 


/nethome/syoung/base/bin/comparison/mergeFiles.pl \
--inputdir /mihg/data/NGS/syoung/base/pipeline/SRA/test/fasta \
--outputdir /mihg/data/NGS/syoung/base/pipeline/SRA/test/samples \
--mode fastq



=cut 


################################################################################
######################		fastaInfos
################################################################################

=head2

	SUBROUTINE		compress

	PURPOSE

=cut

sub compress {
	my $filename	=	shift;
	my $compress	=	shift;

	my $compressfile = $filename;
	if ( $compress =~ /^gzip$/ )
	{
		$compressfile = $filename . ".gz";
		`rm -fr $compressfile`;
		my $command = "$compress $filename";
		print "$command\n";
		print `$command`;
	}
	elsif ( $compress =~ /^zip$/ )
	{
		$compressfile = $filename . ".zip";
		`rm -fr $compressfile`;
		my $command = "$compress $filename";
		print "$command\n";
		print `$command`;
	}	

	return $filename;
}


=head2

    SUBROUTINE      fixed_width

    PURPOSE

        MAKE FASTA RECORD OF FIXED LENGTH

=cut

sub fixed_width {
    my $args        =   shift;
	my $type		=	shift;


    $args->{$type} =~ s/>$//;
    my ($id, $sequence) = $args->{$type} =~ /^>*([^\n]+)\n(.+)$/;

    $id = substr($id, 0, $args->{id_length});
    $id = $id . " " x ($args->{id_length} - length($id));
    $sequence = substr($sequence, 0, $args->{sequence_length});
    $sequence = $sequence . " " x ($args->{sequence_length} - length($sequence));

    return ">$id\n$sequence";
}


################################################################################
######################		listFiles
################################################################################
=head2

	SUBROUTINE		listfiles

	PURPOSE

		EXTRACT FASTA RECORDS FROM FILES IN THE INPUT DIRECTORY

			1. GET READ COUNTS AND RECORD LENGTHS FROM .info FILES

			2. GENERATE RANDOM READ LISTS OF size READS PER FILE:

				1) PRINT total/size FILES PER .info FILE

				2) PRINT ONE READ LOCATION PER LINE

			3. FOR EACH ORIGINAL READ FILE, GENERATE SAMPLED .N.fasta

				FILES FROM .N.readlist FILES

			4. CONCAT FASTA FILES IN EACH FOLDER INTO SINGLE FASTA FILE

=cut

sub listfiles {
	my $args			=	shift;


use Data::Dumper;
print Dumper $args;

	#### GET INPUT DIRECTORIES	
	my $inputdir		=	$args->{inputdir};
	my $directories;
	@$directories = split ",", $inputdir;
	print "Sampler::listfiles    directories: \n";
	print join "\n", @$directories;
	print "\n";
	print "Sampler::listfiles    No. directories: ", scalar(@$directories), "\n";

	#### GET THE INFO FILES AND TOTAL READ COUNTS FOR EACH DIRECTORY
	my $samples = {};
	foreach my $directory ( @$directories )
	{
		my $temp_args = {};
		%$temp_args = %$args;
		$temp_args->{inputdir} = $directory;	
		my $infofiles = infofiles($temp_args);
		my $total = total($infofiles, $temp_args);
		$samples->{$directory}->{infofiles} = $infofiles;
		$samples->{$directory}->{total} = $total;
	}

	#### CALCULATE TOTAL READS IN ALL DIRECTORIES COMBINED	
	my $global_total = 0;
	foreach my $directory ( keys %$samples )
	{
		my $total = $samples->{$directory}->{total};
		die "Total not defined for directory: $directory\n" if not defined $total;
		$global_total += $total;
	}
	print "global_total: $global_total\n";
	print "\n";

	#### PRINT READ STATS TO STDOUT
	print "\ndirectory\ttotal reads\n";
	foreach my $directory ( keys %$samples )
	{
		print "$directory\t$samples->{$directory}->{total}\n";
	}
	print "\n";

	#### PRINT READ STATS TO EACH INPUT DIR
	foreach my $directory ( @$directories )
	{
		my $readstatsfile = "$directory/readstats.txt";
		open(OUT, ">$readstatsfile") or die "Sampler::listfiles     Can't open readstats file: $readstatsfile\n";

		foreach my $directory ( @$directories )
		{
			my $total = $samples->{$directory}->{total};
			print OUT "directory\treads\n";
			print OUT "$directory\t$total\n";
		}
		print OUT "Total\t$global_total\n";
	}

	#### GET FRACTION FOR EACH DIRECTORY
	print "directory\tfraction\n";
	foreach my $directory ( keys %$samples )
	{
		my $fraction = ($samples->{$directory}->{total})/ $global_total;
		print "$directory\t$fraction\n";
		$samples->{$directory}->{fraction} = $fraction;
		die "Fraction not defined for directory: $directory\n" if not defined $fraction;
	}

	return $samples;
}




=head2

	SUBROUTINE		printListFile

	PURPOSE

		GENERATE READ SUBFILES FOR A GIVEN INPUT FILE:

			1. GENERATE MULTIPLE *.list FILES CONTAINING

				RANDOM size NUMBER OF READ RECORD LOCATIONS

				IN INPUT FILE (RECORDS MUST BE FIXED WIDTH)

			2. PRINT total/size *.list FILES PER .info FILE

			3. PRINT ONE READ LOCATION PER LINE IN *.list FILE

	NOTES

		GO THROUGH FASTA FILE AS A STREAM:

			1. DEAL RECORDS LIKE CARDS, RANDOMLY TO LIST FILES IN GROUP

			2. DEAL UNTIL ALL LIST FILES HAVE RECEIVED ONE FASTA RECORD

			3. START A NEW ROUND, CONTINUE UNTIL END OF FILE

	INPUTS

		1. FRACTION (READS IN FILES IN INPUTFILE DIRECTORY / TOTAL READS)

		2. DESIRED SIZE OF OUTPUT FILE

	OUTPUTS

		1. MULTIPLE *.list FILES FOR EACH .info FILE, EACH CONTAINING 

			size NUMBER OF READ RECORD BYTE LOCATION ENTRIES

		2. REMAINDER *.list FILE CONTAINING REMAINING READ LOCATIONS

=cut

sub printListFile {
	my $args			=	shift;


	my $infofile 		=	$args->{inputfile};
	my $size			=	$args->{size};
	my $compress		=	$args->{compress};
	my $dot				=	$args->{dot};
	my $fraction 		=	$args->{fraction}; #### FRACTION THIS DIRECTORY MAKES UP OF THE TOTAL READS
	my $total 			=	$args->{total};

	print "Sampler::printListFile    infofile: $infofile\n";
	print "Sampler::printListFile    total (reads in directory): $total\n";
	print "Sampler::printListFile    fraction (this dir reads / all dir reads): $fraction\n";
	print "Sampler::printListFile    size (desired no. reads per file): $size\n";

	#### GET INFO FROM INFOFILE
	my $info = getInfo($infofile);

	my $record_length = $info->{totallength};
	print "Sampler::printListFile    record_length: $record_length\n";

	#### GET NUMBER OF READS IN INPUT FILE
	my $reads = $info->{reads};
	print "Sampler::printListFile    reads: $reads\n";

	#### ELUDE INPUTFILE FROM INFOFILE
	my $inputfile = inputfile($infofile);
	print "Sampler::printListFile    inputfile: $inputfile\n";

	#### TOTAL READS PER FILE ADJUSTED FOR FRACTION TOTAL READS
	my $readsnumber = int($fraction * $size);

	#### FRACTION OF TOTAL DIRECTORY READS IN THIS FASTA FILE
	my $directory_fraction = $reads/$total;		

	#### NUMBER OF READS IN THIS FILE FROM THIS DIRECTORY
	my $file_readsnumber = int($readsnumber * $directory_fraction);
	print "Sampler::printListFile    portion of reads from this directory = int(fraction * size)= $readsnumber\n";

	print "Sampler::printListFile    directory fraction (file reads/directory reads): $directory_fraction ($reads/$total)\n";


	#### CALCULATE THE NUMBER OF READ LIST FILES TO GENERATE FOR EACH FASTA FILE
	#### AND THE remainder NUMBER OF READS LEFT OVER
	my $numberfiles = int($reads / $file_readsnumber);
	my $remainder = $reads % $file_readsnumber;
	#$samples->{$inputdir}->{remainder}->{$infofilename} = $remainder;		
	print "Sampler::printListFile    number files: $numberfiles\n";
	print "Sampler::printListFile    remainder: $remainder\n";

	#### REMOVE LIST FILES
	my ($file, $subdir, $directory) = file_subdir($inputfile);
	$subdir = "$directory/$subdir";
	my $remove_command = "rm -fr $subdir/$file.*.list";
	print "Sampler::printListFile    remove_command: $remove_command\n";
	`$remove_command`;

	#### GENERATE FULL PATHS OF LIST FILES
	my $listfiles = set_listfiles($inputfile, $numberfiles);
	#my $remainderfile;
	#if ( defined $remainder and $remainder )
	#{
	#	#$remainderfile = pop(@$listfiles);
	#	$remainderfile->{reads} = $remainder;
	#}

	#### SKIP IF INFO FILE IS THE COMPLEMENTARY *_2.fastq* MATE FILE 
	#### (WILL USE read_1 LIST FILES TO GENERATE read_2 SUBFILES LATER)
	if ( $infofile =~ /_2\./ )
	{
		print "Sampler::printListFile    Sampler::printListFile    Skipping mate pair infofile: $infofile\n";
		return;
	}

	#### ADD TOTAL NUMBER OF READS IN THIS INPUT FILE
	#### TO EACH LIST FILE ENTRY
	foreach my $listfile ( @$listfiles )
	{
		$listfile->{reads} = $file_readsnumber;
	}

	#### GET FILEHANDLES FOR ALL LIST FILES
	$listfiles = open_listfiles($listfiles);

	print "Sampler::printListFile    Number listfiles: ", scalar(@$listfiles), "\n";
	foreach my $listfile ( @$listfiles )
	{
		print "Sampler::printListFile    $listfile->{filename}\n";
	}

	#### USE COPY OF LISTFILES TO PRESERVE ORDER OF FILES IN
	#### LISTFILES FOR USE WHEN PRINTING FASTA FILES LATER
	my $listfiles_copy;
	@$listfiles_copy = @$listfiles;

	#### DEAL RECORDS TO ALL LIST FILES
	my $record = 0;
	while ( $record < ($reads - $remainder) )
	{
		#### DO ONE ROUND, DEALING TO ALL FILES IN THE GROUP
		my $done = [];
		while ( @$listfiles_copy and $record < ($reads - $remainder) )
		{
			my $rand = int(rand(@$listfiles_copy));
			my $filehandle = $$listfiles_copy[$rand]->{filehandle};
			my $location = $record_length * $record;
			print $filehandle "$location\n";

			$record++;
			print "Sampler::printListFile    $record\n" if $record % $dot == 0;

			push @$done, splice(@$listfiles_copy, $rand, 1);
		}

		(@$listfiles_copy) = (@$listfiles_copy, @$done);
	}
	print "Sampler::printListFile    record index: $record\n";

	##### CLOSE FILE HANDLES FOR ALL LIST FILES
	$listfiles = close_listfiles($listfiles);

	#### PRINT REMAINING FASTA RECORD BYTE POSITIONS TO REMAINDER FILE
	#my $filenumber = $numberfiles + 1;
	#my $remainderfile = "$inputfile.$filenumber.list";

	#open(OUT, ">$remainderfile->{filename}") or die "Can't open remainder file: $remainderfile->{filename}\n";
	#while 	( $record < $reads )
	#{
	#	#print OUT $record_length * $record, "\n";
	#	my $location = $record_length * $record;
	#	$record++;
	#}
	#close(OUT);

	#### ADD REMAINDER FILE TO LIST FILES
	#push @$listfiles, $remainderfile;
}	

=head2

	SUBROUTINE		readlists

	PURPOSE

		GENERATE RANDOM READ LISTS OF size READS PER FILE:

			1. PRINT total/size FILES PER .info FILE

			2. PRINT ONE READ LOCATION PER LINE

	NOTES

		GO THROUGH FASTA FILE AS A STREAM:

			1. DEAL RECORDS LIKE CARDS, RANDOMLY TO LIST FILES IN GROUP

			2. DEAL UNTIL ALL LIST FILES HAVE RECEIVED ONE FASTA RECORD

			3. START A NEW ROUND, CONTINUE UNTIL END OF FILE

=cut

sub readlists {

	my $args			=	shift;
	my $samples			=	shift;
	#my $global_total	=	shift;

	my $inputdir		=	$args->{inputdir};
	my $size			=	$args->{size};
	my $compress		=	$args->{compress};
	my $dot				=	$args->{dot};
	my $paired 			=	$args->{paired};


	my $infofiles 		=	$samples->{$inputdir}->{infofiles};

	#### FRACTION THIS DIRECTORY MAKES UP OF THE TOTAL READS
	print "Doing readlists for directory: $inputdir\n";
	my $fraction 		=	$samples->{$inputdir}->{fraction};
	print "this directory contains this fraction of the total reads: $fraction\n";
	print "size (desired no. reads per file): $size\n";

	#### DIRECTORY READS = TOTAL READS PER FILE ADJUSTED FOR FRACTION TOTAL READS
	my $readsnumber = int($fraction * $size);
	print "portion of reads from this directory = int(fraction * size): $readsnumber\n";

	my $total 		=	$samples->{$inputdir}->{total};
	print "total reads in directory: $total\n";

	my $readlists;

	foreach my $infofile ( @$infofiles )
	{
		#### GET INPUT FILE NAME FROM INFO FILE NAME
		#use Data::Dumper;
		my @keys = keys %$infofile;
		my $infofilename = $keys[0];
		my $inputfile = inputfile($infofilename, $compress);

		#### NEXT IF THIS IS A MATE PAIR
		#### WE'LL COPY THE LIST FILES FROM read_1 TO THOSE FOR
		#### read_2 AT THE END OF THIS FOR LOOP
		if ( defined $paired )
		{
			if ( $infofilename =~ /_2[^\/]+$/ )
			{
				#exit;
				next;
			}
		}


		#### GET NUMBER OF READS IN INPUT FILE
		my $reads = $infofile->{$infofilename}->{reads};
		die "Reads not defined for info file: $infofilename\n" if not defined $reads;

		#### GET FRACTION OF TOTAL DIRECTORY READS IN THIS FASTA FILE
		my $directory_fraction = $reads/$total;
		my $file_readsnumber = int($readsnumber * $directory_fraction);

		print "directory fraction (file reads/directory reads): $directory_fraction ($reads/$total)\n";


		#### CALCULATE THE NUMBER OF READ LIST FILES TO GENERATE FOR EACH FASTA FILE
		#### AND THE remainder NUMBER OF READS LEFT OVER
		my $numberfiles = int($reads / $file_readsnumber);
		my $remainder = $reads % $file_readsnumber;
		#$samples->{$inputdir}->{remainder}->{$infofilename} = $remainder;		
		print "number files: $numberfiles\n";
		print "remainder: $remainder\n";

		#### GENERATE FULL PATHS OF LIST FILES
		my $listfiles = set_listfiles($inputfile, $numberfiles, $remainder);
		my $remainderfile;
		if ( defined $remainder and $remainder )
		{
			$remainderfile = pop(@$listfiles);
			$remainderfile->{reads} = $remainder;
		}


		#### ADD TOTAL NUMBER OF READS IN THIS INPUT FILE
		#### TO EACH LIST FILE ENTRY
		foreach my $listfile ( @$listfiles )
		{
			$listfile->{reads} = $file_readsnumber;
		}

		#### GET FILEHANDLES FOR ALL LIST FILES
		$listfiles = open_listfiles($listfiles);

		print "Number listfiles: ", scalar(@$listfiles), "\n";
		print "Listfiles:\n";
		foreach my $listfile ( @$listfiles )
		{
			print "$listfile->{filename}\n";
		}
		print "Remainderfile:\n";
		print "$remainderfile->{filename}\n";

		my $record_length = $infofile->{$infofilename}->{totallength};
		print "Record length: $record_length\n";

		#### USE COPY OF LISTFILES TO PRESERVE ORDER OF FILES IN
		#### LISTFILES FOR USE WHEN PRINTING FASTA FILES LATER
		my $listfiles_copy;
		@$listfiles_copy = @$listfiles;

		#### DEAL RECORDS TO ALL LIST FILES
		my $record = 0;
		while ( $record < ($reads - $remainder) )
		{
			#### DO ONE ROUND, DEALING TO ALL FILES IN THE GROUP
			my $done = [];
			while ( @$listfiles_copy and $record < ($reads - $remainder) )
			{
				my $rand = int(rand(@$listfiles_copy));
				my $filehandle = $$listfiles_copy[$rand]->{filehandle};
				my $location = $record_length * $record;
				print $filehandle "$location\n";


				$record++;
				print "$record\n" if $record % $dot == 0;

				push @$done, splice(@$listfiles_copy, $rand, 1);
			}

			(@$listfiles_copy) = (@$listfiles_copy, @$done);
		}
		print "record index: $record\n";

		##### CLOSE FILE HANDLES FOR ALL LIST FILES
		$listfiles = close_listfiles($listfiles);

		#### PRINT REMAINING FASTA RECORD BYTE POSITIONS TO REMAINDER FILE
		#my $filenumber = $numberfiles + 1;
		#my $remainderfile = "$inputfile.$filenumber.list";
		open(OUT, ">$remainderfile->{filename}") or die "Can't open remainder file: $remainderfile->{filename}\n";
		while 	( $record < $reads )
		{
			my $location = $record_length * $record;
			print OUT $location, "\n";
			$record++;
		}
		close(OUT);


		#### ADD REMAINDER FILE TO LIST FILES
		push @$listfiles, $remainderfile;


		#### IF PAIRED, COPY read_1 LIST FILES TO read_2 LIST FILES
		if ( defined $paired )
		{
			foreach my $listfile ( @$listfiles )
			{
				my $matefile = $listfile->{filename};
				$matefile =~ s/_1([^\/]+)$/_2$1/;
				print "matefile: $matefile\n";
				my $copy_command = "cp $listfile->{filename} $matefile";
				print "$copy_command\n";
				`$copy_command`;
			}			
		}


		#### ADD INFO FILE VS. LISTFILES HASH TO READ LIST
		my $readlist = {
			$infofilename => $listfiles
		};
		push @$readlists, $readlist;

	}	####	@$infofiles

	return $readlists;	
}	

=head2

    SUBROUTINE      listfiles

    PURPOSE

		PRINT 1..N FASTA FILES FOR EACH ORIGINAL SEQUENCE

		FILE IN THE DIRECTORY WHERE N IS THE NUMBER OF LIST FILES

=cut

sub get_listfiles {
	print "Sampler::get_listfiles    Sampler::get_listfiles(inputfile)\n";
	my $inputfile	=	shift;
	print "Sampler::get_listfiles    inputfile: $inputfile\n";

	my ($file, $subdir, $directory) = file_subdir($inputfile);

	my $listfiles;	
	@$listfiles = split "\n", `ls $directory/$subdir/*list`;

	my ($string) = $inputfile =~ /\/([^\/]+)$/;

	#### USE FIRST MATE LIST FILE FOR SECOND MATE
	$string =~ s/_2\./_1./;

	for ( my $i = 0; $i < @$listfiles; $i++ )
	{
		my ($list_string) = $$listfiles[$i] =~ /\/([^\/]+)$/;

		if ( $list_string!~ /^$string/ )
		{
			splice @$listfiles, $i, 1;
			$i--;
		}
	}

	return $listfiles;
}


=head2

	SUBROUTINE		regex_safe

	PURPOSE

		CONVERT THE INPUT STRING INTO A FORM THAT CAN BE USED AS A REGEX

=cut


sub regex_safe {
	my $string		=	shift;

	$string =~ s/\//\\\//g;
	$string =~ s/\./\\./g;
	$string =~ s/\-/\\-/g;
	$string =~ s/\[/\\[/g;
	$string =~ s/\]/\\]/g;
	$string =~ s/\(/\\(/g;
	$string =~ s/\)/\\)/g;
	$string =~ s/\*/\\*/g;
	$string =~ s/\+/\\+/g;

	return $string;
}


=head2

	SUBROUTINE		set_listfiles

	PURPOSE

		SET THE LIST OF INPUT FILES BASED ON:

			1. INFO FILE

			2. COMPRESS OPTION

=cut

sub set_listfiles {
	my $inputfile		=	shift;
	my $numberfiles		=	shift;
	my $remainder		=	shift;


	#$numberfiles++ if defined $remainder and $remainder;

	my ($file, $subdir, $directory) = file_subdir($inputfile);
	$subdir = "$directory/$subdir";

	mkdir($subdir) if not -d $subdir;

	my $listfiles = [];
	for ( my $i = 1; $i < $numberfiles + 1; $i++ )
	{
		my $filename = "$subdir/$file.$i.list";
		my $listfile;
		$listfile->{filename} = $filename;
		push @$listfiles, $listfile;		
	}

	return $listfiles;
}


=head2

	SUBROUTINE		file_subdir

	PURPOSE

		GENERATE A SUBDIR COMPOSED OF THE FIRST N CHARACTERS OF

		A FILENAME

=cut

sub file_subdir {
	my $inputfile		=	shift;

	my ($directory, $file) = $inputfile =~ /^(.+?)\/([^\/]+)$/;

	#### USE THIS NUMBER OF CHARACTERS FROM THE FILENAME TO CREATE THE SUBDIR NAME
	my $CHARS			=	8;	
	my ($subdir) = $file =~ /^(.{$CHARS})/;


	return $file, $subdir, $directory;
}


=head2

	SUBROUTINE		open_listfiles

	PURPOSE

		OPEN A FILE HANDLE FOR EACH LIST FILE AND ATTACH

		IT TO THE LIST FILE ENTRY

=cut

sub open_listfiles {
	my $listfiles	=	shift;

	foreach my $listfile ( @$listfiles )
	{
		my $filename = $listfile->{filename};
		my $filehandle;
		open($filehandle, ">$filename") or die "Can't open for writing list file: $filename\n";

		$listfile->{filehandle} = $filehandle;
	}

	return $listfiles;	
}


=head2

	SUBROUTINE		close_listfiles

	PURPOSE

		close THE FILE HANDLE FOR EACH LIST FILE

=cut

sub close_listfiles {
	my $listfiles	=	shift;

	foreach my $listfile ( @$listfiles )
	{
		my $filename = $listfile->{filename};
		my $filehandle = $listfile->{filehandle};
		close($filehandle) or die "Can't close for writing list file: $filename\n";
	}

	return $listfiles;	
}

=head2

	SUBROUTINE		inputfile

	PURPOSE

		RETURN INPUT FILE GIVEN INFO FILE AND COMPRESS OPTION

=cut

sub inputfile {
	my $infofile		=	shift;
	my $compress		=	shift;

	my $inputfile = $infofile;
	$inputfile =~ s/\.info//;
	if ( defined $compress and $compress =~ /^gzip$/ )
	{
		$inputfile .= ".gz";
	}
	if ( defined $compress and $compress =~ /^zip$/ )
	{
		$inputfile .= ".zip";
	}

	return $inputfile;
}



=head2

	SUBROUTINE		infofiles

	PURPOSE

		GET THE infofiles NUMBER OF READS IN ALL FILES IN THE INPUT DIRECTORY

=cut

sub infofiles {
	my $args			=	shift;

	my $inputdir		=	$args->{inputdir};
	my $compress		=	$args->{compress};
	my $dot				=	$args->{dot};

	#### GET ALL REQUIRED FILES IN THE INPUT DIRECTORY
	my $files = files($args);

	#### COUNT infofiles READS IN ALL FILES
	my $infofiles = [];
	foreach my $file ( @$files )
	{
		my $infofile = $file;
		$infofile =~ s/\.(gz|zip)$//;
		$infofile .= ".info";		
		$infofile = "$inputdir/$infofile";
		die "Sampler::infofiles    Can't find infofile: $infofile\n" if not -f $infofile;

		open(FILE, $infofile) or die "Can't open infofile: $infofile\n";
		$/ = undef;
		my ($reads, $headerlength, $sequencelength, $totallength) = <FILE> =~ /^(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/;
		close(FILE);
		next if not defined $reads or not $reads;

		push @$infofiles, {
			$infofile => {
				'reads' => $reads,
				'headerlength' => $headerlength,
				'sequencelength' => $sequencelength,
				'totallength' => $totallength
			}
		};
	}

	return $infofiles;
}


=head2

	SUBROUTINE		getInfo

	PURPOSE

		GET INFO FROM infofile

=cut

sub getInfo {
	my $infofile	=	shift;

	print "Sampler::getInfo    infofile: $infofile\n";
	die "Sampler::getInfo    Can't find info file: $infofile\n" if not -f $infofile;

	open(FILE, $infofile) or die "Can't open info file: $infofile\n";
	$/ = undef;
	my ($reads, $headerlength, $sequencelength, $totallength) = <FILE> =~ /^(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/;
	close(FILE);
	die "Sampler::getInfo    no reads info in infofile: $infofile\n" and exit if not defined $reads or not $reads;

	my $info = {
		'reads' => $reads,
		'headerlength' => $headerlength,
		'sequencelength' => $sequencelength,
		'totallength' => $totallength
	};

	return $info;
}


=head2

	SUBROUTINE		total

	PURPOSE

		GET THE TOTAL NUMBER OF READS IN ALL FILES IN THE INPUT DIRECTORY

=cut


sub total {

	my $infofiles		=	shift;
	my $args			=	shift;

	#### CALCULATE .info FILES IF NOT DEFINED
	$infofiles = infofiles($args) if not defined $infofiles;

	#### GET ALL REQUIRED FILES IN THE INPUT DIRECTORY
	my $files = files($args);

	#### COUNT TOTAL READS IN ALL FILES
	my $total = 0;
	foreach my $infofile ( @$infofiles )
	{
		my @keys = keys %$infofile;
		my $key = $keys[0];

		$total += $infofile->{$key}->{reads};
	}

	return $total;
}




=head2

	SUBROUTINE		files

	PURPOSE

		GET fastq OR fastq.gz FILES IN THE INPUT DIRECTORY

=cut

sub files {
	my $args			=	shift;


	use Data::Dumper;

	my $inputdir		=	$args->{inputdir};
	my $compress		=	$args->{compress};

	#### GET FILES
	my $files = [];
	my @dirs = split ",", $inputdir;

	foreach my $dir ( @dirs )
	{
		my $dirfiles = Util::files($dir);

		next if not defined $dirfiles;

		@$files = ( @$files, @$dirfiles );
	}
	die "No files in input directory: $inputdir\n" if not defined $files or scalar(@$files) == 0;

	#### FASTQ.GZ FILES ONLY
	if ( defined $compress and $compress =~ /^gzip$/ )
	{
		$files = Util::by_suffix($files, "\\fastq.\\.gz");
	}
	else
	{
		#### FASTQ FILES ONLY
		$files = Util::by_suffix($files, "\\.fastq");
	}
	die "No files in input directory: $inputdir\n" if not @$files or scalar(@$files) == 0;

	return $files;
}


=head2

	SUBROUTINE		filepaths

	PURPOSE

		GET REQUIRED FILES IN THE INPUT DIRECTORY

=cut

sub filepaths {
	my $args			=	shift;


	my $inputdir		=	$args->{inputdir};
	my $compress		=	$args->{compress};
	my $mode			=	$args->{mode};

	#### CHECK INPUTS
	die "No input directory in args: ", print Dumper $args, if not defined $inputdir;

	#### GET FILES
	my $files = [];
	my @dirs = split ",", $inputdir;
	foreach my $dir ( @dirs )
	{
		my $dirfiles = Util::files($dir);
		foreach my $dirfile ( @$dirfiles )
		{
			$dirfile = "$dir/$dirfile";
		}
		next if not defined $dirfiles;

		@$files = ( @$files, @$dirfiles );
	}
	die "No files in input directory: $inputdir\n" if not defined $files or scalar(@$files) == 0;

	if ( defined $compress and $compress =~ /^gzip$/ and $mode =~ /^fasta$/ )
	{
		$files = Util::by_suffix($files, "\\.fastq\\.gz");
	}
	else
	{
		$files = Util::by_suffix($files, "\\.fastq");
	}
	die "No files in input directory: $inputdir\n" if not @$files or scalar(@$files) == 0;

	#@$files = sort by_mostly_numeric ( @$files );
	@$files = sort { $a cmp $b } ( @$files );

	return $files;
}

sub by_mostly_numeric {
    ($a <=> $b) || ($a cmp $b);
}


################################################################################
######################		printFiles
################################################################################

=head2

    SUBROUTINE      subfiles

    PURPOSE

		PRINT 1..N FASTQ FILES FOR EACH ORIGINAL SEQUENCE

		FILE IN THE DIRECTORY WHERE N IS THE NUMBER OF LIST FILES

=cut

sub subfiles {
	my $args		=	shift;


	print "Sampler::subfiles    Sampler::subfiles(args)\n";	
	{
		print "Sampler::subfiles    hostname: ";
		print `hostname`;
		print "Sampler::subfiles    args: \n";
		use Data::Dumper;
		print Dumper $args;
	}

	my $compress 	=	$args->{compress};
	my $inputfile 	=	$args->{inputfile};
	my $width 		=	$args->{width};
	my $dot 		=	$args->{dot};
	my $outputdir	=	$args->{outputdir};
	my $paired		=	$args->{paired};

	my $listfiles	=	get_listfiles($inputfile);	
	my $number_listfiles = scalar(@$listfiles);

#exit;

	#### OPEN INPUT FILE AND SET RECORD SEPARATOR
	if( defined $compress and $compress =~ /^gzip$/ )
	{
		my $zcat_command = "zcat $inputfile.gz > $inputfile";
		print "$zcat_command\n";
		`$zcat_command`;
	}

	#### OPEN INPUT FILE
	open(FILE, $inputfile) or die "Sampler::subfiles    Can't open inputfile: $inputfile\n";

	#### LIST OF FASTQ FILES TO BE RETURNED
	my $outputfiles;
	my $counter = 0;
	foreach my $listfile ( @$listfiles )
	{		

		#### SET QUAL FILE AND ADD TO QUAL FILES ARRAY
		my ($filenumber) = $listfile =~ /(\d+)\.list$/;

		#### GET FILE AND SUBDIR
		my ($file, $subdir, $directory) = file_subdir($inputfile);

		#### SET FASTQ FILE
		my $outputfile = "$directory/$subdir/$file.$filenumber-$number_listfiles.reads-fastq";
		print "outputfile: $outputfile\n";
		push @$outputfiles, $outputfile;

		#### OPEN OUTPUT FASTA FILE
		open(OUTFILE, ">$outputfile") or die "Sampler::subfiles    Can't open outputfile: $outputfile\n";


		#### PRINT RECORDS IN THIS LIST FILE TO FILE FILE:
		#### 1. OPEN LIST FILE AND READ BYTE POSITION ON EACH LINE
		#### 2. SEEK INSIDE INPUT FILE FOR BYTE POSITION
		#### 3. EXTRACT FILE RECORD AND PRINT TO FILE FILE
		open(LISTFILE, $listfile) or die "Can't open list file: $listfile\n";
		$/ = "\n";
		while ( <LISTFILE> )
		{
			next if $_ =~ /^\s*$/;

			print "$counter\n" if $counter % $dot == 0;

			my $record;
			seek(FILE, $_, 0);
			read(FILE, $record, $width);

			print OUTFILE "$record";

			$counter++;


		}
		close(LISTFILE);
		close(OUTFILE);

#last;


		#### REMOVE UNCOMPRESSED FASTA FILE
		if( defined $compress and $compress =~ /^gzip$/ )
		{
			my $rm_command = "rm -fr $inputfile";
			print "$rm_command\n";
			`$rm_command`;
		}

	}	#### OUTFILE FILES

	close(FILE);

	return $outputfiles;	
}		



=head2

    SUBROUTINE      qual

    PURPOSE

		PRINT 1..N QUAL FILES FOR EACH ORIGINAL SEQUENCE

		FILE IN THE DIRECTORY WHERE N IS THE NUMBER OF LIST FILES

=cut

sub qual {
	print "Sampler::qual(args, filename)\n";

	my $args		=	shift;

	my $compress 	=	$args->{compress};
	my $inputfile 	=	$args->{inputfile};
	my $tempdir 	=	$args->{tempdir};
	my $width 		=	$args->{width};
	my $outputdir	=	$args->{outputdir};

	my $listfiles	=	get_listfiles($inputfile);	
	my $number_listfiles = scalar(@$listfiles);

	my $ascii		=	ascii();

	#### SET INPUTFILE AS .QUAL OR .QUAL.GZ FILE
	$inputfile =~ s/qual$/qual.qual/;
	$inputfile =~ s/qual\.gz$/qual.qual.gz/;

	my $qualfiles;
	foreach my $listfile ( @$listfiles )
	{		

		#### SET QUAL FILE AND ADD TO QUAL FILES ARRAY
		my ($filenumber) = $listfile =~ /(\d+)\.list$/;

		#### GET FILE AND SUBDIR
		my ($file, $subdir, $directory) = file_subdir($inputfile);

		#### SET QUAL FILE
		my $qualfile = "$directory/$subdir/$file.$filenumber-$number_listfiles.reads-fasta-qual";
		#if ( -f $qualfile )
		#{
		#	next;
		#}

		#### SET TEMPFILE
		my ($filename) = $inputfile =~ /([^\/]+)$/;
		my $tempfile = "$tempdir/$filename.$filenumber-$number_listfiles.reads-qual";

		push @$qualfiles, $qualfile;

		#### OPEN OUTPUT QUAL FILE
		open(QUAL, ">$tempfile") or die "Can't open for writing QUAL file: $tempfile\n";

		#### OPEN INPUT FILE AND SET RECORD SEPARATOR
		if( $inputfile =~ /\.gz$/ )
		{
			my $pipe_command = "zcat $inputfile |";
			open(FILE, $pipe_command) or die "Can't open input file: $inputfile\n";
		}
		else
		{
			open(FILE, $inputfile) or die "Can't open input file: $inputfile\n";
		}

		#### PRINT RECORDS IN THIS LIST FILE TO QUAL FILE:
		#### 1. OPEN LIST FILE AND READ BYTE POSITION ON EACH LINE
		#### 2. SEEK INSIDE INPUT FILE FOR BYTE POSITION
		#### 3. EXTRACT QUAL RECORD AND PRINT TO QUAL FILE
		open(LISTFILE, $listfile) or die "Can't open list file: $listfile\n";
		$/ = "\n";
		while ( <LISTFILE> )
		{
			next if $_ =~ /^\s*$/;

			my $output;
			seek(FILE, $_, 0);
			read(FILE, $output, $width);
			my ($id, $symbolic_quality) = $output =~ /^([^\n]+)\n(.+)$/;
			my $numeric_quality = symbolic2numeric($symbolic_quality, $ascii);
			$output = "$id\n$numeric_quality\n";
			print QUAL $output;

		}
		close(FILE);
		close(LISTFILE);
		close(QUAL);

		#### MOVE TEMP FILE TO QUAL FILE
		print "hostname: ";
		print `hostname`;
		my $move = "mv $tempfile $qualfile";
		print "move: $move\n";
		`$move`;

		print "QUAL file printed:\n$qualfile\n\n";

	}	#### QUAL FILE

	return $qualfiles;
}




=head2

    SUBROUTINE      fasta

    PURPOSE

		PRINT 1..N FASTA FILES FOR EACH ORIGINAL SEQUENCE

		FILE IN THE DIRECTORY WHERE N IS THE NUMBER OF LIST FILES

=cut

sub fasta {
	print "Sampler::fasta(args, filename)\n";

$DEBUG = 1;

	my $args		=	shift;

	my $compress 	=	$args->{compress};
	my $inputfile 	=	$args->{inputfile};
	my $tempdir 	=	$args->{tempdir};
	my $width 		=	$args->{width};
	my $outputdir	=	$args->{outputdir};

	my $listfiles	=	get_listfiles($inputfile);	
	my $number_listfiles = scalar(@$listfiles);

	my $fastafiles;
	my $filecounter = 0;
	foreach my $listfile ( @$listfiles )
	{
		$filecounter++;


		#### SET FASTA FILE AND ADD TO FASTA FILES ARRAY
		my ($filenumber) = $listfile =~ /(\d+)\.list$/;

		#### GET FILE AND SUBDIR
		my ($file, $subdir, $directory) = file_subdir($inputfile);

		#### SET FASTA FILE		
		my $fastafile = "$directory/$subdir/$file.$filenumber-$number_listfiles.reads-fasta";
		push @$fastafiles, $fastafile;

		#### SET TEMPFILE
		my ($filename) = $inputfile =~ /([^\/]+)$/;

		my $tempfile = $fastafile;
		$tempfile = "$tempdir/$filename.$filenumber-$number_listfiles.reads-fasta" if defined $tempdir and $tempdir;

		#### OPEN OUTPUT FASTA FILE
		open(FASTA, ">$tempfile") or die "Can't open for writing FASTA file: $tempfile\n";

		#### OPEN INPUT FILE AND SET RECORD SEPARATOR
		if( $compress =~ /^gzip$/ )
		{
			#my $unzip_command = "gunzip -f $inputfile.gz";
			#`$unzip_command`;
			my $zcat_command = "zcat $inputfile.gz > $inputfile";
			print "$zcat_command\n";
			`$zcat_command`;
		}

		#### OPEN INPUT FILE
		open(FILE, $inputfile) or die "Sampler::fasta    Can't open input file: $inputfile\n";

		#### PRINT RECORDS IN THIS LIST FILE TO FASTA FILE:
		#### 1. OPEN LIST FILE AND READ BYTE POSITION ON EACH LINE
		#### 2. SEEK INSIDE INPUT FILE FOR BYTE POSITION
		#### 3. EXTRACT FASTA RECORD AND PRINT TO FASTA FILE
		open(LISTFILE, $listfile) or die "Can't open list file: $listfile\n";
		$/ = "\n";

		my $counter = 0;
		while ( <LISTFILE> )
		{
			$counter++;

			next if $_ =~ /^\s*$/;

			my $seek = $_;
			$seek =~ s/\s+//g;


			my $output;
			seek FILE, $seek, 0;

			my $offset = tell(FILE);
			print "Sampler::fasta    offset: $offset\n";
			read(FILE, $output, $width);

			print "Sampler::fasta    counter $counter width: $width, seek $_, output: $output\n";			
			last if $counter >= 10;

			#exit;
		}
		close(FILE);
		close(LISTFILE);
		close(FASTA);

		#### REMOVE UNCOMPRESSED FASTA FILE
		if( $compress =~ /^gzip$/ )
		{
			my $rm_command = "rm -fr $inputfile";
			print "$rm_command\n";
			`$rm_command`;
		}

		#### MOVE TEMP FILE TO FASTQ FILE
		if ( defined $tempdir and $tempdir )
		{
			my $move = "mv $tempfile $fastafile";
			print "Sampler::fasta    move: $move\n";
			`$move`;
		}

		print "FASTA file printed:\n$fastafile\n\n";

		last if $filecounter >= 3;

	}	#### FASTA FILE

	return $fastafiles;
}


=head2

	SUBROUTINE		symbolic2numeric

	PURPOSE

		CONVERT SYMBOLIC QUALITY VALUES TO NUMERIC QUALITY VALUES

=cut

sub symbolic2numeric {
	my $quality		=	shift;
	my $ascii 		= 	shift;

    #### REMOVE 'G' RESIDUES AT EITHER END

    my @symbols = split "", $quality;
    my $numeric_quality = '';
    foreach my $symbol ( @symbols )
    {
		if ( not defined $ascii->{$symbol} )
		{
			print "ascii not defined for $symbol\n";
		}

		#### NB: MINUS 33 FOR STANDARD SANGER FASTQ (NOT 64 AS IN SOLEXA FASTQ) 
		#### http://seqanswers.com/forums/showthread.php?t=1110
		#### there are now three types of FASTQ files floating about; standard Sanger FASTQ with quality scores expressed as ASCII(Qphred+33), Solexa FASTQ with ASCII(Qsolexa+64) and Solexa FASTQ with ASCII(Qphred+64)		
		my $number = $ascii->{$symbol} - 33;
		#die "no number for symbol: $symbol\n" if not defined $ascii->{$symbol};


		if ( $number < 0 )
		{
			$number = 0;
		}
        $numeric_quality .= "$number ";
    }
	$numeric_quality =~ s/\s+$//;


	return $numeric_quality;	
}



=head2

	SUBROUTINE		ascii

	PURPOSE

		PROVIDE A HASH OF MAPPINGS BETWEEN ASCII SYMBOLS AND THEIR DECIMAL NUMBERS

=cut

sub ascii {
	my $self		=	shift;
	my $key			=	shift;
	my $value		=	shift;

    my $ascii;    

	if ( not defined $key )
	{
		$key = "symbol";
	}
	if ( not defined $value )
	{
	    $value = "dec";
	}

    my $data = data();
    my @lines = split "\n", $data;
    my @headings = split "\t", $lines[0];

	my %headings_hash;
    my $counter = 0;
    foreach my $heading ( @headings )
    {
        $headings_hash{$heading} = $counter;
        $counter++;
    }

    my $value_number = $headings_hash{$value};
    my $key_number = $headings_hash{$key};
    if ( not defined $value_number )
    {
        die "sub ascii(). Value $value number not defined in headings: @headings\n";
    }
    if ( not defined $key_number )
    {
        die "sub ascii(). Key $key number not defined in headings: @headings\n";
    }

#    while ( <DATA>  )
#    {
#	    $_ =~ s/\s+$//;
#		my @elements = split "\t", $_;

    for ( my $i = 1; $i < $#lines + 1; $i++ )
    {
        my $line = $lines[$i];
 		my @elements = split "\t", $line;

		if ( @elements )
		{
	        $ascii->{$elements[$key_number]} = $elements[$value_number];
		}
    }


	$ascii->{' '} = 32;


    return $ascii;
}



=head2

	SUBROUTINE		data

	PURPOSE

		PROVIDE A TABLE OF MAPPINGS BETWEEN ASCII SYMBOLS AND THEIR DECIMAL NUMBERS

		NB: MINUS 33 FOR STANDARD SANGER FASTQ (NOT 64 AS IN SOLEXA FASTQ) 

		http://seqanswers.com/forums/showthread.php?t=1110

		There are now three types of FASTQ files floating about:
			1. Standard Sanger FASTQ with quality scores expressed as ASCII(Qphred+33)
			2. Solexa FASTQ with ASCII(Qsolexa+64)
			3. Solexa FASTQ with ASCII(Qphred+64)		

=cut

sub data {
    my $data = qq{dec	octal	hex	binary	symbol
0	0	0	0	NUL  (Null char.)
1	1	1	1	SOH  (Start of Header)
2	2	2	10	STX  (Start of Text)
3	3	3	11	ETX  (End of Text)
4	4	4	100	EOT  (End of Transmiss	ion)
5	5	5	101	ENQ  (Enquiry)
6	6	6	110	ACK  (Acknowledgment)
7	7	7	111	BEL  (Bell)
8	10	8	1000	BS  (Backspace)
9	11	9	1001	HT  (Horizontal Tab)
10	12	00A	1010	LF  (Line Feed)
11	13	00B	1011	VT  (Vertical Tab)
12	14	00C	1100	FF  (Form Feed)
13	15	00D	1101	CR  (Carriage Return)
14	16	00E	1110	SO  (Shift Out)
15	17	00F	1111	SI  (Shift In)
16	20	10	10000	DLE  (Data Link Escape	)
17	21	11	10001	DC1  (XON) (Device Control 1)
18	22	12	10010	DC2      (Device Control 2)
19	23	13	10011	DC3  (XOFF)(Device Control 3)
20	24	14	10100	DC4      (Device Control 4)
21	25	15	10101	NAK  (Negative Acknowledgement)
22	26	16	10110	SYN  (Synchronous Idle)
23	27	17	10111	ETB  (End of Trans. Block)
24	30	18	11000	CAN  (Cancel)
25	31	19	11001	EM  (End of Medium)
26	32	01A	11010	SUB  (Substitute)
27	33	01B	11011	ESC  (Escape)
28	34	01C	11100	FS  (File Separator)
29	35	01D	11101	GS  (Group Separator)
30	36	01E	11110	RS  (Req to Send)(RecSep)
31	37	01F	11111	US  (Unit Separator)
32	40	20	100000	SP  (Space)
33	41	21	100001	!
34	42	22	00100010	"
35	43	23	100011	#
36	44	24	100100	\$
37	45	25	100101	\%
38	46	26	100110	&
39	47	27	100111	'
40	50	28	101000	(
41	51	29	101001	)
42	52	02A	101010	*
43	53	02B	101011	+
44	54	02C	00101100	,
45	55	02D	101101	-
46	56	02E	101110	.
47	57	02F	101111	/
48	60	30	110000	0
49	61	31	110001	1
50	62	32	110010	2
51	63	33	110011	3
52	64	34	110100	4
53	65	35	110101	5
54	66	36	110110	6
55	67	37	110111	7
56	70	38	111000	8
57	71	39	111001	9
58	72	03A	111010	:
59	73	03B	111011	;
60	74	03C	111100	<
61	75	03D	111101	=
62	76	03E	111110	>
63	77	03F	111111	?
64	100	40	1000000	\@
65	101	41	1000001	A
66	102	42	1000010	B
67	103	43	1000011	C
68	104	44	1000100	D
69	105	45	1000101	E
70	106	46	1000110	F
71	107	47	1000111	G
72	110	48	1001000	H
73	111	49	1001001	I
74	112	04A	1001010	J
75	113	04B	1001011	K
76	114	04C	1001100	L
77	115	04D	1001101	M
78	116	04E	1001110	N
79	117	04F	1001111	O
80	120	50	1010000	P
81	121	51	1010001	Q
82	122	52	1010010	R
83	123	53	1010011	S
84	124	54	1010100	T
85	125	55	1010101	U
86	126	56	1010110	V
87	127	57	1010111	W
88	130	58	1011000	X
89	131	59	1011001	Y
90	132	05A	1011010	Z
91	133	05B	1011011	[
92	134	05C	1011100	\\
93	135	05D	1011101	]
94	136	05E	1011110	^
95	137	05F	1011111	_
96	140	60	1100000	`
97	141	61	1100001	a
98	142	62	1100010	b
99	143	63	1100011	c
100	144	64	1100100	d
101	145	65	1100101	e
102	146	66	1100110	f
103	147	67	1100111	g
104	150	68	1101000	h
105	151	69	1101001	i
106	152	06A	1101010	j
107	153	06B	1101011	k
108	154	06C	1101100	l
109	155	06D	1101101	m
110	156	06E	1101110	n
111	157	06F	1101111	o
112	160	70	1110000	p
113	161	71	1110001	q
114	162	72	1110010	r
115	163	73	1110011	s
116	164	74	1110100	t
117	165	75	1110101	u
118	166	76	1110110	v
119	167	77	1110111	w
120	170	78	1111000	x
121	171	79	1111001	y
122	172	07A	1111010	z
123	173	07B	1111011	{
124	174	07C	1111100	|
125	175	07D	1111101	}
126	176	07E	1111110	~
127	177	07F	1111111	DEL};

    return $data;
}





################################################################################
######################		mergeFiles
################################################################################

=head2

	SUBROUTINE		mergeFiles

	PURPOSE

		MERGE A LIST OF FILES INTO A SINGLE FILE

=cut

sub mergeFiles {
	my $args			=	shift;


	my $inputdir		=	$args->{inputdir};	
	my $outputdir		=	$args->{outputdir};	
	my $label			=	$args->{label};
	my $paired			=	$args->{paired};
	my $indices			=	$args->{indices};
	my $maxindex		=	$args->{maxindex};
	my $tempdir 		=	$args->{tempdir};

    print "Sampler::mergeFiles    inputdir not defined\n" and exit if not defined $inputdir;
    print "Sampler::mergeFiles    outputdir not defined\n" and exit if not defined $outputdir;

	my $inputfiles = filepaths($args);
	die "No input files in directories: $inputdir\n" if not defined $inputfiles;
	print "Sampler::mergeFiles    Number inputfiles: ", scalar(@$inputfiles), "\n";

	#### LOOK UP MAX INDEX IN FILENAMES IF INDICES OR MAXINDEX NOT SUPPLIED
	if ( not defined $indices and not defined $maxindex )
	{
		my $inputfile = $$inputfiles[0];
		my ($file, $subdir, $directory) = file_subdir($inputfile);
		my $subfiles = Util::files("$directory/$subdir");
		$subfiles = Util::by_regex($subfiles, $file);
		$subfiles = Util::by_regex($subfiles, "\.reads-fastq\$");
		($maxindex) = $$subfiles[0] =~ /(\d+)\.reads-fastq$/;
	}

	#### GENERATE LISTFILES FROM INDICES
	$indices = "1-$maxindex" if not defined $indices;

	my $listfiles = [];
	my @array = split ",", $indices;
	foreach my $symbols ( @array )
	{
		die "Incorrect entry in indices option: $symbols\n" if $symbols !~ /^\d+$/ and $symbols !~ /^\d+\-\d+$/;
		if ( $symbols =~ /^(\d+)$/ )
		{
			push @$listfiles, $symbols;
			$maxindex = $symbols;
		}
		elsif ( $symbols =~ /^(\d+)\-(\d+)$/ )
		{
			$maxindex = $2;

			die "Second index is less than first index in indices option: $symbols\n" if $1 > $2;
			foreach my $number ( $1..$2 )
			{
				push @$listfiles, $number;
			}
		}
	}



	#### STORE IN HASH WITH KEY-VALUE PAIRS: targetfile => array of small files
	my $mergefiles = {};
	foreach my $listfile_index ( @$listfiles ) 
	{
		my $smallfiles = ();
		foreach my $inputfile ( @$inputfiles )
		{
			my ($file, $subdir, $directory) = file_subdir("$inputfile");

			my $smallfile = "$directory/$subdir/$file.$listfile_index-$maxindex.";
			$smallfile .= "reads-fastq";
			push @$smallfiles, $smallfile;
		}

		#### DEAL WITH PAIRED END
		if ( defined $paired )
		{
			my $files1 = Util::by_regex($smallfiles, "_1[^\/]+\$");
			my $files2 = Util::by_regex($smallfiles, "_2[^\/]+\$");

			#### SET READ 1 OUTPUT FILE
			my $singlefile1 = "$outputdir/$label-$listfile_index.reads_1.";
			$singlefile1 .= "fastq";

			#### SET READ 2 OUTPUT FILE
			my $singlefile2 = "$outputdir/$label-$listfile_index.reads_2.";
			$singlefile2 .= "fastq";

			#### ADD TO MERGE FILES
			$mergefiles->{$singlefile1} = $files1;
			$mergefiles->{$singlefile2} = $files2;
		}
		else
		{

			#### SET TEMPFILE
			my $tempfile = "$tempdir/reads.$listfile_index.";
			$tempfile .= "fastq";

			#### SET SINGLE FILE
			my $singlefile = "$outputdir/reads.$listfile_index.";
			$singlefile .= "fastq";

			$mergefiles->{$singlefile} = $smallfiles;
		}
	}	#	listfiles

	return $mergefiles;
}




################################################################################
######################		splitFiles
################################################################################

=head2

	SUBROUTINE		get_splitfiles

	PURPOSE

		PARSE SPLIT FILE NAMES FROM LIST FILE

=cut

sub get_splitfiles {    
    my $filename        =   shift;


    my $splitfiles;
	open(FILE, $filename) or die "Can't open split filename file: $filename\n";
	while ( <FILE> )
	{
		my ($filenumber, $index, $outputfile) = $_ =~ /^(\S+)\s+(\S+)\s+(\S+)$/;
		$splitfiles->[$filenumber][$index] = $outputfile;
	}
	close(FILE) or die "Can't close split filename: $filename\n";

    return $splitfiles;
}



=head2

	SUBROUTINE		set_splitfiles (KEEP AS LEGACY FOR run-maq FOR TIME BEING)

	PURPOSE

		SPLIT ONE OR MORE INPUT FILES INTO SMALLER FILES:

			1. SUBFILES CONTAIN A FIXED NUMBER OF LINES

			2. SUBFILES ARE NUMBERED

			3. SUBFILES ARE PRINTED TO NUMBERED SUBDIRS

				INSIDE OUTPUT DIRECTORY

=cut

sub set_splitfiles {
	return splitfiles(@_);
}


=head2

	***** DEPRECATED: USE splitFiles BELOW
	***** WILL REMOVE THIS LATER AFTER FIXING Comparison SCRIPTS TO
	***** USE splitFiles

	SUBROUTINE		splitfiles

	PURPOSE

		SPLIT ONE OR MORE INPUT FILES INTO SMALLER FILES:

			1. SUBFILES CONTAIN A FIXED NUMBER OF LINES

			2. SUBFILES ARE NUMBERED

			3. SUBFILES ARE PRINTED TO NUMBERED SUBDIRS

				INSIDE OUTPUT DIRECTORY

		INPUT

			A HASH CONTAINING THE FOLLOWING KEY PAIRS

				1. inputfiles:	COMMA-SEPARATED LIST OF INPUT FILES
				2. lines:		MAXIMUM NUMBER OF LINES IN EACH SPLITFILE
				3. splitfile:	PATH TO 'SPLITFILE' TO WHICH SPLIT FILE NAMES WILL BE WRITTEN
				4. outputdir:	OUTPUT DIRECTORY
				5. suffix:		(OPTIONAL) SUFFIX FOR SPLIT FILES

		OUTPUT

			1. 'SPLITFILE' CONTAINING NAMES OF SPLIT FILES

			2. SPLIT FILES IN OUTPUT DIRECTORY

		NOTES

			head /p/NGS/syoung/base/pipeline/SRA/eland/1/reads.1.fastq.split
				0       0       /p/NGS/syoung/base/pipeline/SRA/eland/1/1/reads_SANGER_1.1.fastq
				1       0       /p/NGS/syoung/base/pipeline/SRA/eland/1/2/reads_SANGER_1.2.fastq
				2       0       /p/NGS/syoung/base/pipeline/SRA/eland/1/3/reads_SANGER_1.3.fastq
				3       0       /p/NGS/syoung/base/pipeline/SRA/eland/1/4/reads_SANGER_1.4.fastq
				4       0       /p/NGS/syoung/base/pipeline/SRA/eland/1/5/reads_SANGER_1.5.fastq
				5       0       /p/NGS/syoung/base/pipeline/SRA/eland/1/6/reads_SANGER_1.6.fastq
				6       0       /p/NGS/syoung/base/pipeline/SRA/eland/1/7/reads_SANGER_1.7.fastq
				7       0       /p/NGS/syoung/base/pipeline/SRA/eland/1/8/reads_SANGER_1.8.fastq
				8       0       /p/NGS/syoung/base/pipeline/SRA/eland/1/9/reads_SANGER_1.9.fastq
				9       0       /p/NGS/syoung/base/pipeline/SRA/eland/1/10/reads_SANGER_1.10.fastq

			tail /p/NGS/syoung/base/pipeline/SRA/eland/1/reads.1.fastq.split
				90      1       /p/NGS/syoung/base/pipeline/SRA/eland/1/91/reads_SANGER_2.91.fastq
				91      1       /p/NGS/syoung/base/pipeline/SRA/eland/1/92/reads_SANGER_2.92.fastq
				92      1       /p/NGS/syoung/base/pipeline/SRA/eland/1/93/reads_SANGER_2.93.fastq
				93      1       /p/NGS/syoung/base/pipeline/SRA/eland/1/94/reads_SANGER_2.94.fastq
				94      1       /p/NGS/syoung/base/pipeline/SRA/eland/1/95/reads_SANGER_2.95.fastq
				95      1       /p/NGS/syoung/base/pipeline/SRA/eland/1/96/reads_SANGER_2.96.fastq
				96      1       /p/NGS/syoung/base/pipeline/SRA/eland/1/97/reads_SANGER_2.97.fastq
				97      1       /p/NGS/syoung/base/pipeline/SRA/eland/1/98/reads_SANGER_2.98.fastq
				98      1       /p/NGS/syoung/base/pipeline/SRA/eland/1/99/reads_SANGER_2.99.fastq
				99      1       /p/NGS/syoung/base/pipeline/SRA/eland/1/100/reads_SANGER_2.100.fastq

=cut

sub splitfiles {
	my $args			=	shift;


	my $inputfiles		=	$args->{inputfiles};
	my $lines			=	$args->{lines};
	my $splitfile		=	$args->{splitfile};
	my $outputdir		=	$args->{outputdir};
	my $suffix			=	$args->{suffix};


	{
		print "args:\n";
		use Data::Dumper;
		print Dumper $args;
	}

	die "Argument 'inputfiles' not defined\n" if not defined $inputfiles;
	die "Argument 'lines' not defined\n" if not defined $lines;
	die "Argument 'splitfile' not defined\n" if not defined $splitfile;
	die "Argument 'outputdir' not defined\n" if not defined $outputdir;

	#### CONVERT INPUTFILES INTO ARRAY
	my $files;
	@$files = split ",", $inputfiles;	

	#### CREATE OUTPUT DIRECTORY
	die "Can't create output directory: $outputdir\n" if not -d $outputdir and not mkdir($outputdir);

	#### OPEN SPLIT FILE FOR WRITING LIST OF SUBFILES
	open(SPLITFILE, ">$splitfile") or die "Can't open splitfile: $splitfile\n";

	my $splitfiles;
	for ( my $index = 0; $index < @$files; $index++ )
	{
		my $file = $$files[$index];

		my $filenumber = 1;

		my ($basedir, $filename) = $file =~ /^(.+?)\/([^\/]+)$/;
		$filename =~ s/\.([^\.]+?)$//;

		#### CREATE SUB DIRECTORY FOR SPLIT FILE
		my $subdir = "$outputdir/$filenumber";
		mkdir($subdir) or die "Can't create output directory: $subdir\n" if not -d $subdir;

		#### OPEN OUTPUT FILE 
		my $outputfile;
		if ( defined $suffix )
		{
			$filename =~ s/\.[^\.]+$//;
			$outputfile = "$subdir/$filename.$filenumber.$suffix";					
		}
		else
		{
			$outputfile = "$subdir/$filename.$filenumber";
		}
		$splitfiles->[$filenumber - 1][$index] = $outputfile;
		print SPLITFILE $filenumber - 1, "\t$index\t$outputfile\n";

		#### OPEN OUTPUT FILE
		open(OUTFILE, ">$outputfile") or die "Can't open output file: $outputfile\n";

		my $counter = 0;
		open(FILE, $file) or die "Can't open input file: $file\n";
		while ( <FILE> )
		{
			if ( $counter >= $lines )
			{
				#### CLOSE OLD OUTPUT FILE
				close(OUTFILE);
				$filenumber++;

				#### CREATE SUB DIRECTORY FOR SPLIT FILE
				$subdir = "$outputdir/$filenumber";
				mkdir($subdir) or die "Can't create output directory: $subdir\n" if not -d $subdir;

				#### SET NEW OUTPUT FILE NAME
				if ( defined $suffix )
				{
					$filename =~ s/\.[^\.]+$//;
					$outputfile = "$subdir/$filename.$filenumber.$suffix";					
				}
				else
				{
					$outputfile = "$subdir/$filename.$filenumber";
				}
				$splitfiles->[$filenumber - 1][$index] = $outputfile;

                #### PRINT FILENUMBER, INDEX AND OUTPUTFILE TO SPLIT FILE LIST
				print SPLITFILE $filenumber - 1, "\t$index\t$outputfile\n";

				#### OPEN NEW OUTPUT FILE
				open(OUTFILE, ">$outputfile") or die "Can't open output file: $outputfile\n";

				$counter = 0;
			}

			$counter++;
			print OUTFILE $_;
		}
		close(FILE);
	}
	close(SPLITFILE);

	return $splitfiles;
}



=head2

	SUBROUTINE		splitFiles

	PURPOSE

			1. SUBFILES CONTAIN A FIXED NUMBER OF LINES

			2. SUBFILES ARE NUMBERED

			3. SUBFILES ARE PRINTED TO NUMBERED SUBDIRS

				INSIDE OUTPUT DIRECTORY

		INPUT

			A HASH CONTAINING THE FOLLOWING KEY PAIRS

				1. inputfiles:	COMMA-SEPARATED LIST OF INPUT FILES (REQUIRED)
				1. matefiles:	COMMA-SEPARATED LIST OF MATE FILES (OPTIONAL)
				2. lines:		MAXIMUM NUMBER OF LINES IN EACH SPLITFILE
				3. splitfile:	PATH TO 'SPLITFILE' TO WHICH SPLIT FILE NAMES WILL BE WRITTEN
				4. outputdir:	OUTPUT DIRECTORY
				5. suffix:		(OPTIONAL) SUFFIX FOR SPLIT FILES

		OUTPUT

			1. 'SPLITFILE' CONTAINING NAMES OF SPLIT FILES

			2. SPLIT FILES IN OUTPUT DIRECTORY

		NOTES

			head /p/NGS/syoung/base/pipeline/SRA/eland/1/reads.1.fastq.split
				0       0       /p/NGS/syoung/base/pipeline/SRA/eland/1/1/reads_SANGER_1.1.fastq
				1       0       /p/NGS/syoung/base/pipeline/SRA/eland/1/2/reads_SANGER_1.2.fastq
				2       0       /p/NGS/syoung/base/pipeline/SRA/eland/1/3/reads_SANGER_1.3.fastq
				3       0       /p/NGS/syoung/base/pipeline/SRA/eland/1/4/reads_SANGER_1.4.fastq
				4       0       /p/NGS/syoung/base/pipeline/SRA/eland/1/5/reads_SANGER_1.5.fastq
				5       0       /p/NGS/syoung/base/pipeline/SRA/eland/1/6/reads_SANGER_1.6.fastq
				6       0       /p/NGS/syoung/base/pipeline/SRA/eland/1/7/reads_SANGER_1.7.fastq
				7       0       /p/NGS/syoung/base/pipeline/SRA/eland/1/8/reads_SANGER_1.8.fastq
				8       0       /p/NGS/syoung/base/pipeline/SRA/eland/1/9/reads_SANGER_1.9.fastq
				9       0       /p/NGS/syoung/base/pipeline/SRA/eland/1/10/reads_SANGER_1.10.fastq

			tail /p/NGS/syoung/base/pipeline/SRA/eland/1/reads.1.fastq.split
				90      1       /p/NGS/syoung/base/pipeline/SRA/eland/1/91/reads_SANGER_2.91.fastq
				91      1       /p/NGS/syoung/base/pipeline/SRA/eland/1/92/reads_SANGER_2.92.fastq
				92      1       /p/NGS/syoung/base/pipeline/SRA/eland/1/93/reads_SANGER_2.93.fastq
				93      1       /p/NGS/syoung/base/pipeline/SRA/eland/1/94/reads_SANGER_2.94.fastq
				94      1       /p/NGS/syoung/base/pipeline/SRA/eland/1/95/reads_SANGER_2.95.fastq
				95      1       /p/NGS/syoung/base/pipeline/SRA/eland/1/96/reads_SANGER_2.96.fastq
				96      1       /p/NGS/syoung/base/pipeline/SRA/eland/1/97/reads_SANGER_2.97.fastq
				97      1       /p/NGS/syoung/base/pipeline/SRA/eland/1/98/reads_SANGER_2.98.fastq
				98      1       /p/NGS/syoung/base/pipeline/SRA/eland/1/99/reads_SANGER_2.99.fastq
				99      1       /p/NGS/syoung/base/pipeline/SRA/eland/1/100/reads_SANGER_2.100.fastq

=cut

sub splitFiles {
	my $args			=	shift;
	my $inputfiles		=	$args->{inputfiles};
	my $matefiles		=	$args->{matefiles};
	my $lines			=	$args->{lines};
	my $splitfile		=	$args->{splitfile};
	my $outputdir		=	$args->{outputdir};
	my $suffix			=	$args->{suffix};
	my $label			=	$args->{label};



	{
		print "args:\n";
		use Data::Dumper;
		print Dumper $args;
	}

	die "Sampler::splitFiles    label not defined\n" if not defined $label;
	die "Sampler::splitFiles    inputfiles not defined\n" if not defined $inputfiles;
	die "Sampler::splitFiles    lines not defined\n" if not defined $lines;
	die "Sampler::splitFiles    splitfile not defined\n" if not defined $splitfile;
	die "Sampler::splitFiles    outputdir not defined\n" if not defined $outputdir;

	#### CREATE OUTPUT DIRECTORY
	die "Sampler::splitFiles    Can't create output directory: $outputdir\n" if not -d $outputdir and not mkdir($outputdir);

	#### OPEN SPLIT FILE FOR WRITING LIST OF SUBFILES
	open(SPLITFILE, ">$splitfile") or die "Sampler::splitFiles    Can't open splitfile: $splitfile\n";	

	my $index = 0;
	my @infiles = split ",", $inputfiles;

	my @mates = split ",", $matefiles if defined $matefiles;

	#### PRINT FIRST MATE FILES
	my $firstmate_label = $label . "_1";

	my $entries = printSplitfiles($outputdir, \@infiles, $lines, $index, $suffix, $label);
	foreach my $entry ( @$entries )
	{
		print SPLITFILE "$entry->{filenumber}\t$entry->{index}\t$entry->{outputfile}\n";
	}

	#### PRINT SECOND MATES FILES IF AVAILABLE	
	if ( defined $matefiles )
	{
		my $firstmate_label = $label . "_2";
		my $mate_index = 1;
		my @mates = split ",", $matefiles;
		my $mate_entries = printSplitfiles($outputdir, \@mates, $lines, $mate_index, $suffix, $label);
		foreach my $entry ( @$mate_entries )
		{
			print SPLITFILE "$entry->{filenumber}\t$entry->{index}\t$entry->{outputfile}\n";
		}

		@$entries = (@$entries, @$mate_entries);
	}
	close(SPLITFILE);


	my $splitfiles;
	foreach my $entry ( @$entries )
	{
		$$splitfiles[$entry->{filenumber}][$entry->{index}] = $entry->{outputfile};
	}



	return $splitfiles;
}





=head2

	SUBROUTINE		printSplitfiles

	PURPOSE

		FOR EACH INPUT FILE IN THE INPUT LIST OF FILES, WRITE CHUNKS TO SUBFILES UNTIL 

		YOU REACH THE END OF THE FILE.

		NB: THE LAST SUBFILE MAY CONTAIN LESS THAN THE REQUIRED READS PER CHUNK.

		FOR MATE PAIRS, THE INPUT FILE INDEX IS 0 AND THE MATE FILE INDEX IS 1.

	INPUTS

		outputdir	:	PRINT SPLIT FILES TO THIS DIRECTORY
		files		:	ARRAY OF INPUT FILES
		lines		:	HOW MANY LINES TO PRINT TO EACH SPLIT FILE
		index		:	USE THIS INDEX AS THE STARTING POINT FOR THE INDEXES
						OF THE SPLIT FILES (I.E., IF THE SPLIT FILES FOR THE
						FIRST MATE HAVE ALREADY BEEN PRINTED, START AT THE
						INDEX FOLLOWING THE LAST OF THEM)
		suffix		:	SUFFIX OF SPLIT FILES
		label		:	USE THIS AS THE STEM OF THE SPLIT FILE NAMES

=cut

sub printSplitfiles {


	my $outputdir	=	shift;
	my $files		=	shift;
	my $lines		=	shift;
	my $index		=	shift;
	my $suffix		=	shift;
	my $label		=	shift;



	#### RETURN THIS ARRAY OF SPLITFILE ENTRIES
	my $entries;

	#### FILE INDEX IS 1 OR 2, INDICATING FIRST OR SECOND MATE FILE
	my $fileindex = $index + 1;
	print "Sampler::printSplitfiles    fileindex: $fileindex\n";
	my $filename = $label . "_" . $fileindex;

	my $filenumber = 1;	
	foreach my $file ( @$files )
	{	

		#### CREATE SUB DIRECTORY FOR SPLIT FILE
		my $subdir = "$outputdir/$filenumber";
		mkdir($subdir) or die "Sampler::printSplitfiles    Can't create output directory: $subdir\n" if not -d $subdir;

		#### SET OUTPUT FILE 
		my $outputfile = splitfileName($subdir, $filename, $filenumber, $suffix);
		print "Sampler::printSplitfiles    $outputfile\n";

		#### PUSH ENTRY FOR FIRST SUBFILE FOR THIS FILE
		push @$entries, {
			filenumber 	=> 	$filenumber - 1,
			index		=> 	$index,
			outputfile	=>	$outputfile
		};

		#### OPEN OUTPUT FILE
		open(OUTFILE, ">$outputfile") or die "Sampler::printSplitfiles    Can't open output file: $outputfile\n";

		my $counter = 0;
		open(FILE, $file) or die "Sampler::printSplitfiles    Can't open input file: $file\n";
		while ( <FILE> )
		{
			if ( $counter >= $lines )
			{
				#### CLOSE OLD OUTPUT FILE
				close(OUTFILE);
				$filenumber++;

				#### CREATE SUB DIRECTORY FOR SPLIT FILE
				$subdir = "$outputdir/$filenumber";
				mkdir($subdir) or die "Sampler::printSplitfiles    Can't create output directory: $subdir\n" if not -d $subdir;

				#### SET NEW OUTPUT FILE NAME
				$outputfile = splitfileName($subdir, $filename, $filenumber, $suffix);
				print "Sampler::printSplitfiles    $outputfile\n";

				#### SAVE THIS ENTRY
				push @$entries, {
					filenumber 	=> 	$filenumber - 1,
					index		=> 	$index,
					outputfile	=>	$outputfile
				};

				#### OPEN NEW OUTPUT FILE
				open(OUTFILE, ">$outputfile") or die "Sampler::printSplitfiles    Can't open output file: $outputfile\n";

				$counter = 0;
			}
			$counter++;

			#### PRINT LINE TO OUTPUT SUBFILE
			print OUTFILE $_;

		}
		close(FILE);
		close(OUTFILE);

		$filenumber++;
	}

	return $entries;
}



=head2

	SUBROUTINE		splitfileName

	PURPOSE

		SET THE NAME OF A SPLITFILE

=cut

sub splitfileName {
	my $subdir		=	shift;
	my $filename		=	shift;
	my $filenumber		=	shift;
	my $suffix		=	shift;

#

	return "$subdir/$filename.$filenumber$suffix" if defined $suffix;

	return "$subdir/$filename.$filenumber" if not defined $suffix;
}



=head2

	SUBROUTINE		lines

	PURPOSE

		COUNT THE NUMBER OF LINES IN A FILE

=cut

sub lines {
	my $filename		=	shift;

	print "Sampler::lines(filename)\n";
	print "filename: $filename\n";

	open(FILE, $filename);
	my $lines = 0;
	while(<FILE>)
	{
		$lines++;
	}
	close(FILE);

	return $lines;
}


=head2

	SUBROUTINE		replaceAtQual

	PURPOSE

		REPLACE THE '@' QUALITY VALUE WITH 'A' (ONE VALUE HIGHER)

		BECAUSE IT INTERFERES WITH MAQ CONVERSION FROM solexa TO fastq

=cut

sub replaceAtQual {
	my $args			=	shift;

	my $inputfile		=	$args->{inputfile};
	my $outputfile		=	$args->{outputfile};
	my $dot				=	$args->{dot};



	{
		print "Sampler::replaceAtQual    Sampler::replaceAtQual    args:\n";
		use Data::Dumper;
		print Dumper $args;
	}

	#### SET OUTPUT FILE IF NOT DEFINED
	my $copyback = 0;
	$copyback = 1 if not defined $outputfile or $outputfile eq $inputfile;
	$outputfile = $inputfile . "-temp" if not defined $outputfile or $outputfile eq $inputfile;

print "Sampler::replaceAtQual    inputfile: $inputfile\n";
print "Sampler::replaceAtQual    outputfile: $outputfile\n";
#exit;

	#### OPEN FILES
	open(FILE, $inputfile) or die "Can't open input file: $inputfile\n";
	open(OUTFILE, ">$outputfile") or die "Can't open output file: $outputfile\n";

	#### PRINT ALTERED QUALITY RECORDS TO OUTPUT FILE
	print "Sampler::replaceAtQual    Printing to output file...\n";
	$/ = "\n";
	my $counter = 0;
	while ( <FILE> )
	{
		next if $_ =~ /^\s*$/;

		print "Sampler::replaceAtQual    $counter\n" if $counter % $dot == 0;

		my $sequence_header = $_;
		my $sequence 		= <FILE>;
		my $quality_header	= <FILE>;
		my $quality			= <FILE>;

		$quality =~ s/@/A/g;

		print OUTFILE $sequence_header, $sequence, $quality_header, $quality;

		$counter++;
	}
	close(FILE);
	close(OUTFILE);

	####
	if ( $copyback )
	{
		my $copy_command = "mv $outputfile $inputfile";
		`$copy_command`;
	}


	return $outputfile;
}


1;

