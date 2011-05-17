#!/usr/bin/perl -w

#### DEBUG
$DEBUG = 1;

#### TIME
my $time = time();

=head2

	APPLICATION     sampleReads

    PURPOSE

        1. SAMPLE READS FROM A LIBRARY 

    INPUT

		1. FASTA OR FASTA.GZ FILE

    OUTPUT

		OPTION: total

			1. TOTAL READS IN ALL READ FILES IN INPUT DIRECTORY

		OPTION: sample

			1. size NUMBER OF FASTA FILES, RECORDS RANDOMLY SELECTED FROM

				FILES IN THE INPUT DIRECTORY

    USAGE

    ./sampleReads.pl <--inputdir String> <--mode String> [-h]

    --inputdir             	:   /full/path/to/inputdir
    --mode                 	:   Functional mode (total, samples)
    --compress				:   Input files are compressed (gzip|zip)
    --help					:   print help info

    EXAMPLES

0600

/nethome/syoung/base/bin/comparison/sampleReads.pl \
--inputdir /mihg/data/NGS/syoung/base/pipeline/SRA/NA18507/SRX000600/fasta \
--compress gzip \
--mode total


1539

/nethome/syoung/base/bin/comparison/sampleReads.pl \
--inputdir /mihg/data/NGS/syoung/base/pipeline/SRA/NA18507/SRX001539/fasta \
--compress gzip \
--mode total



/nethome/syoung/base/bin/comparison/sampleReads.pl \
--inputdir /mihg/data/NGS/syoung/base/pipeline/SRA/NA18507/SRX000600/fasta \
--compress gzip \
--mode samples \
--size 10


=cut

use strict;

#### FLUSH BUFFER
$| = 1;

#### EXTERNAL MODULES
use Data::Dumper;
use File::Path;
use FindBin qw($Bin);
use Getopt::Long;
use Term::ANSIColor qw(:constants);

#### USE LIBRARY
use lib "$Bin/../../../lib";	

#### INTERNAL MODULES
use Timer;
use Util;
use Conf::Agua;

#### GET modeS
use Getopt::Long;
my $inputdir;
my $outputdir;
my $mode;
my $reads;
my $compress;
my $size;
my $cpus = 8;
my $dot = 1000000;
my $help;
GetOptions (
    'inputdir=s' => \$inputdir,
    'outputdir=s' => \$outputdir,
    'mode=s' => \$mode,
    'reads=i' => \$reads,
    'compress=s' => \$compress,
    'size=i' => \$size,
    'cpus=i' => \$cpus,
    'dot=i' => \$dot,
    'help' => \$help             
) or die "No modes specified. Use --help for usage\n";

usage() if defined $help;

#### CHECK INPUTS
die "Input directory not defined (use --help for usage)\n" if not defined $inputdir;
die "Output directory not defined (use --help for usage)\n" if not defined $outputdir;
die "Mode not defined (use --help for usage)\n" if not defined $mode;
die "Compress type must be 'gzip' or 'zip'\n" if defined $compress and $compress !~ /^(gzip|zip)$/;

#### CREATE OUTPUT DIR IF NOT PRESENT
File::Path::mkpath($outputdir) if not -d $outputdir;

no strict;
&$mode(
	{
		'inputdir'	=>	$inputdir,
		'outputdir'	=>	$outputdir,
		'reads'		=>	$reads,
		'compress'	=>	$compress,
		'size'		=>	$size,
		'cpus'		=>	$cpus,
		'dot'		=>	$dot
	}
);
use strict;


#### REPORT OUTPUT FILES PRINTED

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "\nRun time: $runtime\n";
print "Completed $0\n";
print Timer::datetime(), "\n";
print "****************************************\n\n\n";
exit;

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#									SUBROUTINES
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


=head2

	SUBROUTINE		samples

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

sub samples
{

	my $args			=	shift;

	my $inputdir		=	$args->{inputdir};
	my $outputdir		=	$args->{outputdir};
	my $compress		=	$args->{compress};
	my $size			=	$args->{size};

	die "size option (--size Integer) missing\n" if not defined $size;

	my $directories;
	@$directories = split ",", $inputdir;

	#### GET THE INFO FILES AND TOTAL READ COUNTS FOR EACH DIRECTORY
	my $samples = {};
	foreach my $directory ( @$directories )
	{
		$args->{inputdir} = $directory;	
		my $infofiles = infofiles($args);
		my $total = total($args, $infofiles);
		$samples->{$directory}->{infofiles} = $infofiles;
		$samples->{$directory}->{total} = $total;
	}

	#### GET GLOBAL TOTAL
	print "\ndirectory\ttotal reads\n";
	my $global_total = 0;
	foreach my $directory ( keys %$samples )
	{
		print "$directory\t$samples->{$directory}->{total}\n";
		my $total = $samples->{$directory}->{total};
		die "Total not defined for directory: $directory\n" if not defined $total;
		$global_total += $total;
	}
	print "global reads: $global_total\n";

	#### GET FRACTION FOR EACH DIRECTORY
	print "directory\tfraction\n";
	foreach my $directory ( keys %$samples )
	{
		my $fraction = ($samples->{$directory}->{total})/ $global_total;
		print "$directory\t$fraction\n";
		$samples->{$directory}->{fraction} = $fraction;
		die "Fraction not defined for directory: $directory\n" if not defined $fraction;
	}

	#### GENERATE READ LIST FILES	
	print "Generating read list files...\n";
	foreach my $directory ( keys %$samples )
	{
		print "directory: $directory\n";
		$args->{inputdir} = $directory;

		my $readlists = readlists($args, $samples, $global_total);
		$samples->{$directory}->{readlists} = $readlists;
	}
	print "Finished generating read list files\n";



	##### GENERATE FASTA FILES
	#my $numberfastafiles = 0;
	#foreach my $directory ( keys %$samples )
	#{
	#	$args->{inputdir} = $directory;
	#	
	#	my $fastafiles = print_files($args, $samples, 'fasta');
	#	$samples->{$directory}->{fastafiles} = $fastafiles;
	#	
	#	$numberfastafiles = scalar(@$fastafiles);
	#}
	#
	#
	##### GENERATE QUALITY FILES
	#my $ascii = ascii();
	#my $numberqualfiles = 0;
	#foreach my $directory ( keys %$samples )
	#{
	#	$args->{inputdir} = $directory;
	#	
	#	my $qualfiles = print_files($args, $samples, 'fasta.qual', $ascii);
	#	$samples->{$directory}->{qualfiles} = $qualfiles;
	#	
	#	$numberqualfiles = scalar(@$qualfiles);
	#}
	#




	##### GENERATE FASTQ FILES
	#my $numberfastqfiles = 0;
	#foreach my $directory ( keys %$samples )
	#{
	#	$args->{inputdir} = $directory;
	#	
	#	my $fastqfiles = fastqfiles($args, $samples);
	#	$samples->{$directory}->{fastqfiles} = $fastqfiles;
	#	
	#	$numberfastqfiles = scalar(@$fastqfiles);
	#}
	#
	##### MERGE FASTQ FILES
	#mergefiles($samples, $numberfastqfiles, 'fastqfiles');
	#


	##### MERGE FASTA FILES
	#mergefiles($samples, $numberfastafiles, 'fastafiles');
	#
	##### MERGE QUAL FILES
	#mergefiles($samples, $numberqualfiles, 'qualfiles');

	#### REPORT FASTA FILES AND INDIVIDUAL AND CUMULATIVE NUMBERS OF READS IN FILES

	print "Completed samples()\n";	
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

sub readlists
{

	my $args			=	shift;
	my $samples			=	shift;
	my $global_total	=	shift;

	my $inputdir		=	$args->{inputdir};
	my $size			=	$args->{size};
	my $compress		=	$args->{compress};
	my $dot				=	$args->{dot};
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
		#### GET INFO FILE NAME AND THEN INPUT FILE NAME
		my @keys = keys %$infofile;
		my $infofilename = $keys[0];
		my $inputfile = inputfile($args, $infofilename);

		#### GET NUMBER OF READS IN INPUT FILE
		my $reads = $infofile->{$infofilename}->{reads};
		die "Reads not defined for info file: $infofilename\n" if not defined $reads;

		#### GET FRACTION OF TOTAL DIRECTORY READS IN THIS FASTA FILE
		my $directory_fraction = $reads/$total;
		my $file_readsnumber = int($readsnumber * $directory_fraction);

		print "\n\ninfofilename: $infofilename\n";
		print "inputfile: $inputfile\n";
		print "reads in file: $reads (", $reads/$total, " total)\n";
		print "reads (all N listfiles) per directory: $readsnumber\n";
		print "reads (per N listfile) for this file: $file_readsnumber\n";
		print "directory fraction (file reads/directory reads): $directory_fraction ($reads/$total)\n";
		print "file readsnumber: $file_readsnumber\n";

		#### CALCULATE THE NUMBER OF READ LIST FILES TO GENERATE FOR EACH FASTA FILE
		#### AND THE remainder NUMBER OF READS LEFT OVER
		my $numberfiles = int($reads / $file_readsnumber);
		my $remainder = $reads % $file_readsnumber;
		$samples->{$inputdir}->{remainder}->{$infofilename} = $remainder;		
		print "number files: $numberfiles\n";
		print "remainder: $remainder\n";

		#### GENERATE FULL PATHS OF LIST FILES
		my $listfiles = listfiles($inputfile, $numberfiles);

		#### ADD TOTAL NUMBER OF READS IN THIS INPUT FILE
		#### TO EACH LIST FILE ENTRY
		foreach my $listfile ( @$listfiles )
		{
			$listfile->{reads} = $file_readsnumber;
		}

		#### GET FILEHANDLES FOR ALL LIST FILES
		$listfiles = open_listfiles($listfiles);

		print "Number listfiles: ", scalar(@$listfiles), "\n";
		#foreach my $listfile ( @$listfiles )
		#{
		#}

		my $record_length = $infofile->{$infofilename}->{totallength};
		print "Record length: $record_length\n";

		#### USE COPY OF LISTFILES TO PRESERVE ORDER OF FILES IN
		#### LISTFILES FOR USE WHEN PRINTING FASTA FILES LATER
		my $listfiles_copy;
		@$listfiles_copy = @$listfiles;

		#### DEAL RECORDS TO ALL LIST FILES
		my $record_index = 0;
		while ( $record_index < ($reads - $remainder) )
		{
			print "$record_index\n" if $record_index % $dot == 0;

			#### DO ONE ROUND, DEALING TO ALL FILES IN THE GROUP
			my $done = [];
			while ( @$listfiles_copy and $record_index < ($reads - $remainder) )
			{
				my $rand = int(rand(@$listfiles_copy));
				my $filehandle = $$listfiles_copy[$rand]->{filehandle};
				my $location = $record_length * $record_index;
				print $filehandle "$location\n";


				$record_index++;

				push @$done, splice(@$listfiles_copy, $rand, 1);
			}

			(@$listfiles_copy) = (@$listfiles_copy, @$done);
		}
		print "record index: $record_index\n";

		##### CLOSE FILE HANDLES FOR ALL LIST FILES
		$listfiles = close_listfiles($listfiles);

		#### PRINT REMAINING FASTA RECORD BYTE POSITIONS TO REMAINDER FILE
		my $filenumber = $numberfiles + 1;
		my $remainderfile = "$inputfile.$filenumber.list";
		open(OUT, ">$remainderfile") or die "Can't open remainder file: $remainderfile\n";
		while 	( $record_index < $reads )
		{
			my $location = $record_length * $record_index;
			print OUT $location, "\n";
			$record_index++;
		}
		close(OUT);

		#### ADD REMAINDER FILE TO LIST FILES
		my $listfile;
		$listfile->{filename} = $remainderfile;
		$listfile->{reads} = $remainder;
		push @$listfiles, $listfile;

		#### ADD INFO FILE VS. LISTFILES HASH TO READ LIST
		my $readlist = {
			$infofilename => $listfiles
		};
		push @$readlists, $readlist;
	}	#### 	foreach my $infofile ( @$infofiles )

	return $readlists;	
}	


=head2

	SUBROUTINE		fastafiles

	PURPOSE

		PRINT 1..N FASTA FILES FOR EACH ORIGINAL SEQUENCE

		FILE IN THE DIRECTORY

=cut

sub print_files
{

	my $args			=	shift;
	my $samples			=	shift;
	my $type			=	shift;
	my $ascii			=	shift;

	print "type: $type\n";

	my $inputdir		=	$args->{inputdir};
	my $compress		=	$args->{compress};

	my $infofiles		=	$samples->{$inputdir}->{infofiles};

	my $readlists 		=	$samples->{$inputdir}->{readlists};

	#### GET NUMBER OF LIST FILES PER INFO FILE FROM FIRST INFO FILE
	#### (SHOULD BE THE SAME FOR ALL INFO FILES IN THIS DIRECTORY)
	my $readlist = $$readlists[0];
	my @keys = keys %$readlist;
	my $infofilename = $keys[0];
	my $listfiles = $readlist->{$infofilename};
	my $number_listfiles = scalar(@$listfiles);

	my $printfiles;
	for( my $i = 0; $i < $number_listfiles; $i++ )
	{		
		#### SET FASTA FILE AND ADD TO FASTA FILES ARRAY
		my $filenumber = $i + 1;
		my $printfile = "$inputdir/reads.$filenumber-$number_listfiles.$type";
		print "printfile: $printfile\n";
		push @$printfiles, $printfile;

		#### OPEN OUTPUT FASTA FILE
		open(FASTA, ">$printfile") or die "Can't open for writing FASTA file: $printfile\n";

		for ( my $fileindex = 0; $fileindex < @$infofiles; $fileindex++ )
		{
			my $infofile = $$infofiles[$fileindex];
			my $readlist = $$readlists[$fileindex];

			#### GET INFO FILE NAME AND THEN INPUT FILE NAME
			my @keys = keys %$infofile;
			my $infofilename = $keys[0];
			print "infofilename: $infofilename\n";
			my $width = $infofile->{$infofilename}->{totallength};

			my $listfiles = $readlist->{$infofilename};
			my $listfile = $$listfiles[$i]->{filename};
			my $inputfile = inputfile($args, $infofilename);
			if ( $type =~ /qual$/ )
			{
				$inputfile =~ s/fasta$/fasta.qual/;
				$inputfile =~ s/fasta\.gz$/fasta.qual.gz/;
			}
			print "listfile: $listfile\n";
			print "inputfile: $inputfile\n";

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

			#### PRINT RECORDS IN THIS LIST FILE TO FASTA FILE:
			#### 1. OPEN LIST FILE AND READ BYTE POSITION ON EACH LINE
			#### 2. SEEK INSIDE INPUT FILE FOR BYTE POSITION
			#### 3. EXTRACT FASTA RECORD AND PRINT TO FASTA FILE
			open(LISTFILE, $listfile) or die "Can't open list file: $listfile\n";
			$/ = "\n";
			while ( <LISTFILE> )
			{
				next if $_ =~ /^\s*$/;

				my $output;
				seek(FILE, $_, 0);
				read(FILE, $output, $width + 3);

				if ( $type =~ /qual$/ )
				{
					my ($id, $symbolic_quality) = $output =~ /^([^\n]+)\n(.+)$/;
					my $numeric_quality = symbolic2numeric($symbolic_quality, $ascii);
					$output = "$id\n$numeric_quality\n";
					#exit;
				}

				print FASTA $output;

			}

		}	#### foreach my $infofile ( @$infofiles )

		close(FASTA);
		print "Fasta file printed:\n\n$printfile\n";

	}	#### FASTA FILE

	return $printfiles;	
}		





=head2

	SUBROUTINE		fastqfiles

	PURPOSE

		PRINT 1..N FASTA FILES FOR EACH ORIGINAL SEQUENCE

		FILE IN THE DIRECTORY

=cut

sub fastqfiles
{

	my $args			=	shift;
	my $samples			=	shift;

	my $inputdir		=	$args->{inputdir};
	my $compress		=	$args->{compress};
	my $cpus			=	$args->{cpus};

	my $infofiles		=	$samples->{$inputdir}->{infofiles};

	my $readlists 		=	$samples->{$inputdir}->{readlists};

	#### GET NUMBER OF LIST FILES PER INFO FILE FROM FIRST INFO FILE
	#### (SHOULD BE THE SAME FOR ALL INFO FILES IN THIS DIRECTORY)
	my $readlist = $$readlists[0];
	my @keys = keys %$readlist;
	my $infofilename = $keys[0];
	my $listfiles = $readlist->{$infofilename};
	my $number_listfiles = scalar(@$listfiles);

	#### LIST OF FASTQ FILES TO BE RETURNED
	my $fastqfiles;

	for( my $listfile_index = 0; $listfile_index < $number_listfiles; $listfile_index++ )
	{
		my $args;

		#### SET FASTA FILE AND ADD TO FASTA FILES ARRAY
		my $filenumber = $listfile_index + 1;
		my $fastqfile = "$inputdir/reads.$filenumber-$number_listfiles.fastq";
		print "fastqfile: $fastqfile\n";

		#### PUSH ONTO FASTQ FILES ARRAY TO BE RETURNED AT END OF SUBROUTINE
		push @$fastqfiles, $fastqfile;

		#### OPEN OUTPUT FASTA FILE
		open(FASTQ, ">$fastqfile") or die "Can't open for writing FASTA file: $fastqfile\n";

		for ( my $fileindex = 0; $fileindex < @$infofiles; $fileindex++ )
		{
			my $infofile = $$infofiles[$fileindex];
			my $readlist = $$readlists[$fileindex];

			#### GET INFO FILE NAME AND THEN INPUT FILE NAME
			my @keys = keys %$infofile;
			my $infofilename = $keys[0];
			my $width = $infofile->{$infofilename}->{totallength};

			my $listfiles = $readlist->{$infofilename};
			my $listfile = $$listfiles[$listfile_index]->{filename};

			#### SET FASTA FILE NAME
			my $fastafile = inputfile($args, $infofilename);
			$fastafile .= ".gz" if defined $compress;

			#### SET QUAL FILE NAME
			my $qualfile = $fastafile;
			$qualfile =~ s/fasta$/fasta.qual/;
			$qualfile =~ s/fasta\.gz$/fasta.qual.gz/;


			#### OPEN INPUT FILE AND SET RECORD SEPARATOR
			if( $fastafile =~ /\.gz$/ )
			{
				my $fasta_command = "zcat $fastafile |";
				my $qual_command = "zcat $qualfile |";
				open(FASTA, $fasta_command) or die "Can't open input file: $fastafile\n";
				open(QUAL, $qual_command) or die "Can't open input file: $qualfile\n";
			}
			else
			{
				open(FASTA, $fastafile) or die "Can't open input file: $fastafile\n";
				open(QUAL, $fastafile) or die "Can't open input file: $qualfile\n";
			}

			#### PRINT RECORDS IN THIS LIST FILE TO FASTA FILE:
			#### 1. OPEN LIST FILE AND READ BYTE POSITION ON EACH LINE
			#### 2. SEEK INSIDE INPUT FILE FOR BYTE POSITION
			#### 3. EXTRACT FASTA RECORD AND PRINT TO FASTA FILE
			open(LISTFILE, $listfile) or die "Can't open list file: $listfile\n";
			$/ = "\n";
			while ( <LISTFILE> )
			{
				next if $_ =~ /^\s*$/;

				my $fasta;
				seek(FASTA, $_, 0);
				read(FASTA, $fasta, $width + 3);

				my $qual;
				seek(QUAL, $_, 0);
				read(QUAL, $qual, $width + 3);

				$fasta =~ s/^>/@/;
				$qual =~ s/^>/+/;
				print FASTQ "$fasta$qual";				
			}
			close(LISTFILE);
			close(FASTA);
			close(QUAL);

		}	#### foreach my $infofile ( @$infofiles )

		close(FASTQ);
		print "Fastq file printed:\n\n$fastqfile\n";	

	}	#### FASTQ FILES

	return $fastqfiles;	
}		


sub printfastq
{
	my $args		=	shift;

	#my $fastafile	=	$args->{fastafile};
	#my $fastqfile	=	$args->{fastqfile};
	#my $qualfile	=	$args->{qualfile};
	#my $listfile	=	$args->{listfile};
	#my $width		=	$args->{width};
	#
	my $fastqfile =$args->{fastqfile};
	my $infofiles = $args->{infofiles};
	my $readlists = $args->{readlists};
	my $index = $args->{index};

	print "fastqfile: $fastqfile\n";
	print "infofiles: $infofiles\n";
	print "readlists: $readlists\n";
	print "index: $index\n";

	##### OPEN OUTPUT FASTA FILE
	#open(FASTQ, ">$fastqfile") or die "Can't open for writing FASTA file: $fastqfile\n";
	#	
	#for ( my $fileindex = 0; $fileindex < @$infofiles; $fileindex++ )
	#{
	#	my $infofile = $$infofiles[$fileindex];
	#	my $readlist = $$readlists[$fileindex];
	#	#print "infofile: \n";
	#	#print Dumper $infofile;
	#	#print "readlist:\n";
	#	#print Dumper $readlist;
	#	
	#	##### GET INFO FILE NAME AND THEN INPUT FILE NAME
	#	#my @keys = keys %$infofile;
	#	#my $infofilename = $keys[0];
	#	##print "infofilename: $infofilename\n";
	#	##print "infofile:\n";print Dumper $infofile;
	#	#my $width = $infofile->{$infofilename}->{totallength};
	#	##print "Width: $width\n";
	#	#
	#	##### GO THROUGH LIST FILE NUMBER index FOR EACH ORIGINAL FILE
	#	#my $listfiles = $readlist->{$infofilename};
	#	##print "listfiles: \n"; print Dumper $listfiles;
	#	#my $listfile = $$listfiles[$index]->{filename};
	#	#my $fastafile = inputfile($args, $infofilename);
	#	#my $qualfile = $fastafile;
	#	#$qualfile =~ s/fasta$/fasta.qual/;
	#	#$qualfile =~ s/fasta\.gz$/fasta.qual.gz/;
	#	#
	#	#print "listfile: $listfile\n";
	#	#print "fastafile: $fastafile\n";
	#	#print "qualfile: $qualfile\n";
	#	#
	#	###### OPEN INPUT FILE AND SET RECORD SEPARATOR
	#	##if( $fastafile =~ /\.gz$/ )
	#	##{
	#	##	my $fasta_command = "zcat $fastafile |";
	#	##	my $qual_command = "zcat $qualfile |";
	#	##	open(FASTA, $fasta_command) or die "Can't open input file: $fastafile\n";
	#	##	open(QUAL, $qual_command) or die "Can't open input file: $qualfile\n";
	#	##}
	#	##else
	#	##{
	#	##	open(FASTA, $fastafile) or die "Can't open input file: $fastafile\n";
	#	##	open(QUAL, $fastafile) or die "Can't open input file: $qualfile\n";
	#	##}
	#	##
	#	###### PRINT RECORDS IN THIS LIST FILE TO FASTA FILE:
	#	###### 1. OPEN LIST FILE AND READ BYTE POSITION ON EACH LINE
	#	###### 2. SEEK INSIDE INPUT FILE FOR BYTE POSITION
	#	###### 3. EXTRACT FASTA RECORD AND PRINT TO FASTA FILE
	#	##open(LISTFILE, $listfile) or die "Can't open list file: $listfile\n";
	#	##$/ = "\n";
	#	##while ( <LISTFILE> )
	#	##{
	#	##	#print "\$_: $_\n";
	#	##	next if $_ =~ /^\s*$/;
	#	##	
	#	##	my $fasta;
	#	##	seek(FASTA, $_, 0);
	#	##	read(FASTA, $fasta, $width + 3);
	#	##
	#	##	my $qual;
	#	##	seek(QUAL, $_, 0);
	#	##	read(QUAL, $qual, $width + 3);
	#	##
	#	##	$fasta =~ s/^>/@/;
	#	##	$qual =~ s/^>/+/;
	#	##	#print "$fasta$qual";
	#	##	print FASTQ "$fasta$qual";				
	#	##}
	#	##close(LISTFILE);
	#	##close(FASTA);
	#	##close(QUAL);
	#
	#}	#### foreach my $infofile ( @$infofiles )
	#
	#close(FASTQ);
	print "Fastq file printed:\n\n$fastqfile\n";	
	#last;

}




=head2

	SUBROUTINE		mergefiles

	PURPOSE

		MERGE A LIST OF FILES INTO A SINGLE FILE

=cut

sub mergefiles
{
	my $samples			=	shift;
	my $numberfiles		=	shift;
	my $type			=	shift;

	#### MERGE FILES
	for ( my $index = 0; $index < $numberfiles; $index++ )
	{
		my $singlenumber = $index + 1;
		my $suffix = "fasta";
		$suffix = "fasta.qual" if $type =~ /^qualfiles$/;
		$suffix = "fastq" if $type =~ /^fastqfiles$/;

		my $singlefile = "$outputdir/reads.$singlenumber-$numberfiles.$suffix";
		print "singlefile: $singlefile\n";

		my $smallfiles = [];
		foreach my $directory ( keys %$samples )
		{
			print "directory: $directory\n";
			my $files = $samples->{$directory}->{$type};
			my $file = $$files[$index];
			print "file: $file\n";
			push @$smallfiles, $file;
		}

		next if scalar(@$smallfiles) == 1;	

		print "Doing merge...\n";
		my $command = "cat @$smallfiles > $singlefile";
		print "$command\n";
		system($command);
		print "done\n";
	}
}

=head2

	SUBROUTINE		listfiles

	PURPOSE

		RETURN INPUT FILE GIVEN INFO FILE AND COMPRESS OPTION

=cut

sub listfiles
{

	my $inputfile		=	shift;
	my $numberfiles		=	shift;

	my $listfiles = [];
	for ( my $i = 1; $i < $numberfiles + 1; $i++ )
	{
		my $filename = "$inputfile.$i.list";
		my $listfile;
		$listfile->{filename} = $filename;
		push @$listfiles, $listfile;		
	}

	return $listfiles;
}


=head2

	SUBROUTINE		open_listfiles

	PURPOSE

		OPEN A FILE HANDLE FOR EACH LIST FILE AND ATTACH

		IT TO THE LIST FILE ENTRY

=cut

sub open_listfiles
{
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

sub close_listfiles
{
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

sub inputfile
{

	my $args			=	shift;
	my $infofile		=	shift;

	my $compress		=	$args->{compress};

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

sub infofiles
{

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
		print "infofile: $infofile\n";
		die "Can't find info file: $infofile\n" if not -f $infofile;

		open(FILE, $infofile) or die "Can't open info file: $infofile\n";
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

	SUBROUTINE		total

	PURPOSE

		GET THE TOTAL NUMBER OF READS IN ALL FILES IN THE INPUT DIRECTORY

=cut


sub total
{

	my $args			=	shift;
	my $infofiles		=	shift;

	my $inputdir		=	$args->{inputdir};
	my $compress		=	$args->{compress};
	my $dot				=	$args->{dot};

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

	print "total: $total\n";
	return $total;
}



=head2

	SUBROUTINE		files

	PURPOSE

		GET REQUIRED FILES IN THE INPUT DIRECTORY

=cut

sub files
{

	my $args			=	shift;

	my $inputdir		=	$args->{inputdir};
	my $compress		=	$args->{compress};

	#### GET FILES
	my $files = Util::files($inputdir);
	die "No files in input directory: $inputdir\n" if not defined $files or scalar(@$files) == 0;

	if ( defined $compress and $compress =~ /^gzip$/ )
	{
		$files = Util::by_suffix($files, "\.fasta\.gz");
	}
	else
	{
		$files = Util::by_suffix($files, "\.fasta");
	}
	die "No files in input directory: $inputdir\n" if not @$files or scalar(@$files) == 0;

	return $files;
}



sub symbolic2numeric
{
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

sub ascii
{
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

=cut

sub data
{
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






sub usage
{
	print `perldoc $0`;

	exit;
}


