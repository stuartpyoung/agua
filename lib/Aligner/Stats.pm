package Align::Stats;


=head2

	PACKAGE		Align::Stats

	PURPOSE

		1. GENERATE Align::Stats REPORTS IN PARALLEL ON A CLUSTER

			- GET TOTAL INPUT READ COUNT

			- GET TOTAL UNIQUE READ HIT COUNTS PER CHROMOSOME

			- GET TOTAL UNIQUE READ HIT COUNT FOR WHOLE GENOME

	NOTES

		1. COUNT INPUT READS PER SPLIT FILE TO GET TOTAL INPUT READS

		2. GET THE TOTAL NUMBER OF UNIQUE READ HITS FOR WHOLE GENOME

			i.e., THE VENN UNION OF UNIQUE READ HITS FOR ALL CHROMOSOMES

			STRATEGY 1: USE MYSQL TO 

				1. CONVERT SAM FILES TO *.tsv AND *-index.tsv FILES

				2. LOAD ALL chr*/accepted_hits.tsv FILES INTO sam TABLE,
					ONE SQLITE DATABASE PER REFERENCE CHROMOSOME

				3. LOAD ALL chr*/accepted_hits-index.tsv FILES INTO samindex TABLE
					IN A SINGLE SQLITE DATABASE FOR ALL REFERENCE CHROMOSOMES

				4. QUERY COUNT ALL DISTINCT READ IDS

				5. QUERY IDENTITIES AND HIT COUNTS OF READS WITH MOST HITS

				6. EXTRACT READ SEQUENCES OF 100 'MOST HITS' READS

			STRATEGY 2: DO PYRAMID UNIX MERGE SORT OF SAM FILES

				1. SORT ALL CHROMOSOME SAM FILES 

					sort +1 -2 inputfile1 

					OR:

					sort inputfile1 

				2. DO PYRAMID MERGE SORT (RETAIN DUPLICATES)

					sort -m inputfile1 inputfile2  > outfile 

				3. TRY WITH:

					- READ ID, CHROMOSOME AND POSITION FIELDS

					- ONLY READ ID


					UNIX SORT OPTIONS:

					-c     Check whether the given files are  already  sorted:
							if  they are not all sorted, print an error message
							and exit with a status of 1.

					-m     Merge the given files by sorting them as  a  group.
							Each  input  file  should  already  be individually sorted.
							sort normally merges input as it sorts. This option informs
							sort that the input is already sorted, thus sort runs much
							faster.

					-u      Suppress duplicates


				4. COUNT DUPLICATES WITH uniq

					uniq -c out.sam > out.hitcounts

					sort +1 -2 -n out.hitcounts

=cut

use strict;
use warnings;
use Carp;

#### INTERNAL MODULES
use Sampler;

#### HAS A
use Monitor::PBS;
use Monitor::LSF;
use Cluster;
#use LSF::Job;

#### EXTERNAL MODULES
use FindBin qw($Bin);
use Data::Dumper;
use File::Path;
use File::Copy;
use JSON;

require Exporter;
our @ISA = qw(Exporter Cluster);
#our @EXPORT_OK = qw();
our $AUTOLOAD;


#### SET SLOTS
our @DATA = qw(

INPUTDIR
OUTPUTDIR
FILENAME
FILETYPE
LABEL
KEY
REFSEQFILE
TRACKFILE
JBROWSE

CLUSTER
QUEUE
MAXJOBS
CPUS
SLEEP
CLEANUP
VERBOSE
TEMPDIR
DOT
VERBOSE

COMMAND
);
our $DATAHASH;
foreach my $key ( @DATA )	{	$DATAHASH->{lc($key)} = 1;	}





=head2

	SUBROUTINE		readHits

	PURPOSE

		MERGE THE 'READ' COLUMN OF ALL *-index.sam FILES 

	NOTES

		1. EXTRACT READS IDS ONLY FROM ALIGNMENT FILE
			AND PRINT TO .hitid FILE

		2. MERGE ALL .hitid FILES FOR ALL SPLIT INPUT FILES

			FOR EACH REFERENCE INTO SINGLE .hitid FILE

		3. SORT .hitid FILE

		4. COUNT DUPLICATE LINES WITH uniq -c AND OUTPUT

			INTO .hitct FILE

		5. GENERATE DISTRIBUTION OF HITS PER READ

		6. REPEAT 2 - 5 GOING FROM chromosome.hitid TO

			genome.hitct FILE

=cut

sub readHits
{
	my $self			=	shift;
	my $outputdir		=	shift;
	my $references		=	shift;
	my $splitfiles		=	shift;
	my $label			=	shift;



	my $hitfile_name = "out.hitid";
	my $countfile_name = "out.hitct";

	####1. EXTRACT HIT READS IDS ONLY FROM ALIGNMENT FILE
	####	AND PRINT TO .hitid FILE
	####
	$self->hitIds($outputdir, $references, $hitfile_name);

	####2. MERGE ALL .hitid FILES FOR ALL SPLIT INPUT FILES
	####	FOR EACH REFERENCE INTO SINGLE .hitid FILE
	####
	print "Align::Stats::readHits    Doing _cumulativeConcatMerge of $hitfile_name files\n";
	$self->_cumulativeConcatMerge($outputdir, $references, $splitfiles, $hitfile_name, $hitfile_name, "cumulativeMergeHits");
	print "Align::Stats::readHits    COMPLETED _cumulativeConcatMerge of $hitfile_name files\n";

exit;


	####3. SORT .hitid FILE
	####
	foreach my $reference ( $references )
	{
		my $hitfile = "$outputdir/$reference/$hitfile_name";
		my $command = "sort $hitfile";
		print "Align::Stats::readHits    sort command: $command\n";
		print `$command`;
	}

	####4. COUNT DUPLICATE LINES WITH uniq -c AND OUTPUT
	####	INTO .hitct FILE
	####	
	foreach my $reference ( $references )
	{
		my $hitfile = "$outputdir/$reference/$hitfile_name";
		my $countfile = "$outputdir/$reference/$countfile_name";
		my $command = "uniq -c $hitfile > $countfile";
		print "Align::Stats::readHits    count command: $command\n";
		print `$command`;
	}


	####5. GENERATE DISTRIBUTION OF HITS PER READ
	####

	####6. REPEAT 2 - 5 GOING FROM chromosome.hitid TO
	####	genome.hitct FILE
	####


	my $duplicatefile = "$outputdir/duplicate.txt";
	my $uniquefile = "$outputdir/unique.txt";
	File::Remove::rm($duplicatefile);
	File::Remove::rm($uniquefile);

	my $tsvindexfiles = [];
	foreach my $tsvfile ( @$tsvfiles )
	{
		my $tsvindexfile = $tsvfile;
		$tsvindexfile =~ s/\.tsv$/-index.tsv/;
		push @$tsvindexfiles, $tsvindexfile;
	}

	#### MERGE DUPLICATE READS INTO ONE FILE FOR ALL GENOME HITS
	my $commands = [];
	my $counter = 0;
	my $total = scalar(@$tsvfiles);
	foreach my $tsvfile ( @$tsvfiles )
	{
		$counter++;
		push @$commands, "echo 'Loading $counter of $total into duplicatefile: $duplicatefile'";
		push @$commands, "cut -f 1 $tsvfile >> $duplicatefile";
	}


	##### GET UNIQUE READS FOR ALL GENOME
	#my $uniq_command = qq{sort $duplicatefile | uniq -u > $uniquefile};
	#push @$commands, "echo '$uniq_command'";
	#push @$commands, $uniq_command;


	my $job = $self->setJob( $commands, "SAMINDEX", $outputdir );

	#### RUN COMMANDS IN SERIES ON ONE CLUSTER NODE
	$self->runJobs( [ $job ], $label );	




}




=head2

	SUBROUTINE		gatherFiles

	PURPOSE

		RETURN AN ARRAY OF sam FILES, CORRESPONDING tsv FILES AND A HASH

		ARRAY OF reference => sam FILES FOR A GIVEN LIST OF REFERENCE NAMES

		AND LIST OF SPLIT FILES 

=cut

sub gatherFiles
{
	my $self				=	shift;
	my $references			=	shift;
	my $splitfiles			=	shift;


	#### SEARCH THROUGH ARRAY	
	my $counter = 0;
	while ( not defined $$splitfiles[$counter] )
	{
		$counter++;
	}
	my ($basedir) = ${$$splitfiles[$counter]}[0] =~ /^(.+?)\/\d+\/[^\/]+$/;

	my $samfiles = [];
	my $tsvfiles = [];
	my $reference_samfiles = {};

	for my $reference ( @$references )
	{
		my $outputdir = "$basedir/$reference";

		my $files = $self->subfiles($basedir, $reference, $splitfiles, "accepted_hits.sam");
		(@$samfiles)  = (@$samfiles, @$files);

		my $tsvs = $self->samToTsvFilename($samfiles);
		(@$tsvfiles) = (@$tsvfiles, @$tsvs);

		#### ADD TO REFERENCE-HASHFILES HASH FOR LATER SQLITE TABLE LOADING
		$reference_samfiles->{$reference} = $files;
	}	

	return ($samfiles, $tsvfiles, $reference_samfiles);
}


=head2

	SUBROUTINE		splitfilesList

	PURPOSE

		PARSE SPLIT FILE NAMES FROM LIST FILE

=cut

sub splitfilesList
{
	my $self		=	shift;    
    my $filename    =   shift;


    my $splitfiles;
	open(FILE, $filename) or die "Can't open split filename file: $filename\n";
	my $counter = 0;
	while ( <FILE> )
	{
		my ($filenumber, $index, $outputfile) = $_ =~ /^(\S+)\s+(\S+)\s+(\S+)$/;
		$splitfiles->[$counter][$index] = $outputfile;
		$counter++;
	}

    return $splitfiles;
}



=head2

	SUBROUTINE		loadSamIndexTable

	PURPOSE

		RETURN AN ARRAY OF COMMANDS FOR LOADING ONE OR MORE

		TSV FILES INTO AN SQLITE TABLE

=cut

sub loadSamIndexTable
{
	my $self				=	shift;
	my $outputdir			=	shift;
	my $tsvfiles			=	shift;
	my $sqlite				=	shift;
	my $cluster				=	shift;
	my $label				=	shift;


	#### SET TABLE
	my $table 		=	"samindex";

	my $tsvindexfiles = [];
	foreach my $tsvfile ( @$tsvfiles )
	{
		my $tsvindexfile = $tsvfile;
		$tsvindexfile =~ s/\.tsv$/-index.tsv/;
		push @$tsvindexfiles, $tsvindexfile;
	}
	my $commands 	=	$self->tsvToSqliteCommands($tsvindexfiles, $table);
	print "Align::Stats::loadSamIndexTable    No. commands: ", scalar(@$commands), "\n";

	#### SET DB FILE FOR WHOLE GENOME
	my $dbfile		=	"$outputdir/stats-all.dbl";

	#### PRINT COMMAND FILE
	my $commandfile = $dbfile;
	$commandfile =~ s/\.dbl$/.cmd/;
	$self->printCommandFile($commandfile, $commands);

	#### CREATE SQLITE DATABASE
	$self->createSamIndexTable($dbfile);

	if ( not defined $cluster )
	{

	}
	else
	{

		#### SET LOAD COMMAND
		my $load_command = qq{$sqlite $dbfile < $commandfile};
		my $job = $self->setJob( [ $load_command ], $label, $outputdir );

		#### RUN JOB ON CLUSTER
		$self->runJobs( [ $job ], "SAMINDEX-$label" );	
	}
}


=head2

	SUBROUTINE		loadSamTables

	PURPOSE

		RETURN AN ARRAY OF COMMANDS FOR LOADING ONE OR MORE

		TSV FILES INTO AN SQLITE TABLE

=cut

sub loadSamTables
{
	my $self				=	shift;
	my $outputdir			=	shift;
	my $references			=	shift;
	my $reference_samfiles	=	shift;
	my $sqlite				=	shift;
	my $cluster				=	shift;
	my $label				=	shift;

	#### SET TABLE
	my $table 		=	"sam";

	my $load_commands = [];
	foreach my $reference ( @$references )
	{
		my $samfiles 	= 	$reference_samfiles->{$reference};
		my $tsvfiles	=	$self->samToTsvFilename($samfiles);
		my $commands 	=	$self->tsvToSqliteCommands($tsvfiles, $table);

		#### SET DB FILE
		my $dbfile		=	"$outputdir/$reference/stats.dbl";

		#### PRINT COMMAND FILE
		my $commandfile = $dbfile;
		$commandfile =~ s/\.dbl$/.cmd/;
		$self->printCommandFile($commandfile, $commands);

		#### SET LOAD COMMAND
		my $load_command = qq{$sqlite $dbfile < $commandfile};
		push @$load_commands, $load_command;

		#### CREATE SQLITE DATABASE
		$self->createSamTable($dbfile);
	}

	#### LOAD TSV FILES INTO DATABASE LOCALLY IN SERIES IF cluster NOT DEFINED
	if ( not defined $cluster )
	{
		for ( my $i = 0; $i < @$load_commands; $i++ )
		{
			my $load_command = $$load_commands[$i];
			print "$load_command\n";
			print `$load_command`;
		}
	}
	#### OTHERWISE, RUN IN PARALLEL ON CLUSTER
	else
	{
		my $jobs = [];
		for ( my $i = 0; $i < @$load_commands; $i++ )
		{
			my $load_command = $$load_commands[$i];
			my $job = $self->setJob( [ $load_command ], "$label-$i", $outputdir );
			push @$jobs, $job;
		}

		#### RUN JOB ON CLUSTER


#### DEBUG
		#$self->runJobs( $jobs, "SAM-$label" );	



	}
}



=head2

	SUBROUTINE		tsvToSqliteCommands

	PURPOSE

		RETURN AN ARRAY OF COMMANDS FOR LOADING ONE OR MORE

		TSV FILES INTO AN SQLITE TABLE

=cut

sub tsvToSqliteCommands
{
	my $self		=	shift;
	my $tsvfiles	=	shift;
	my $table		=	shift;



	my $commands = [];
	unshift @$commands, qq{.separator "\\t"\n};

	foreach my $tsvfile ( @$tsvfiles )
	{
		push @$commands, ".import $tsvfile $table\n";
	}

	return $commands;
}


=head2

	SUBROUTINE		printCommandFile

	PURPOSE

		PRINT ARRAY OF COMMANDS TO FILE
=cut

sub printCommandFile
{
	my $self		=	shift;
	my $commandfile	=	shift;
	my $commands	=	shift;

	#### CHECK INPUTS
	print "Align::Stats::printCommandFile    commandfile not defined\n" and exit if not defined $commandfile;
	print "Align::Stats::printCommandFile    commands not defined\n" and exit if not defined $commands;

	#### PRINT LOAD COMMANDS TO COMMAND FILE

	open(CMD, ">$commandfile") or die "Can't open command file: $commandfile\n";
	foreach my $command ( @$commands )
	{
		print CMD $command;
	}
}



=head2

	SUBROUTINE		createSamTable

	PURPOSE

		1. CREATE sam TABLE IN SQLITE DATABASE

		2. RETURN DATABASE HANDLE

	NOTES

		SAM Format Specification 0.1.2-draft (20090820)
		http://samtools.sourceforge.net/SAM1.pdf

		TAB-DELIMITED FORMAT:

		<QNAME> <FLAG> <RNAME> <POS> <MAPQ> <CIGAR> <MRNM> <MPOS> <ISIZE> <SEQ> <QUAL> \
		[<TAG>:<VTYPE>:<VALUE> [...]]

		EXAMPLE

		HWI-EAS185:1:3:1636:994#0       0       chr1    3035487 255     52M     *       0       0       CTTTGGGTGCTCTCACACTCCCGAGACTAGAGTGGGGGTCCCCGAAAGGGGG      aaa_a`a\_aaaaaaaaaZaaa_X`Xa`a``VJY\`\ZP`aBBBBBBBBBBB      NM:i:1
		HWI-EAS185:1:4:972:1117#0       0       chr1    3044403 1       52M     *       0       0       TTTACAAGGCCTAATGGTGACTCCTACAGTGGTTAACACCGACTACCATGGG      bbaaba[a[bbab^``^U[_a_ab``aaVPU`^_]_a__aY_a__aa^YZZS      NM:i:1
		HWI-EAS185:1:2:1737:216#0       16      chr1    3045067 0       50M     *       0       0       CACAAAGAGATAATGCGTTAGAATTTGAAGATTGTTGGATGTGCTTTTCT        \X]]Y_a^`^_]Y_`_`_`_a]^[^_`]_a`___U_aaaab`b[aa^\]]        NM:i:2

		FIELDS

			QNAME [^ \t\n\r]+ Query pair NAME if paired; or Query NAME if unpaired 2
			FLAG [0-9]+ [0,216-1] bitwise FLAG (Section 2.2.2)
			RNAME [^ \t\n\r@=]+ Reference sequence NAME 3
			POS [0-9]+ [0,229-1] 1-based leftmost POSition/coordinate of the clipped sequence
			MAPQ [0-9]+ [0,28-1] MAPping Quality (phred-scaled posterior probability that the mapping position of this read is incorrect) 4
			CIGAR ([0-9]+[MIDNSHP])+|\* extended CIGAR string
			MRNM [^ \t\n\r@]+ Mate Reference sequence NaMe; “=” if the same as <RNAME> 3
			MPOS [0-9]+ [0,229-1] 1-based leftmost Mate POSition of the clipped sequence
			ISIZE -?[0-9]+ [-229,229] inferred Insert SIZE 5
			SEQ [acgtnACGTN.=]+|\* query SEQuence; “=” for a match to the reference; n/N/. for ambiguity;
			cases are not maintained 6,7
			QUAL [!-~]+|\* [0,93] query QUALity; ASCII-33 gives the Phred base quality 6,7
			TAG [A-Z][A-Z0-9] TAG
			VTYPE [AifZH] Value TYPE
			VALUE [^\t\n\r]+ match <VTYPE> (space allowed)

=cut

sub createSamTable
{
	my $self		=	shift;
	my $dbfile		=	shift;



	#### REMOVE DB FILE IF EXISTS
	my $recursive = 0;
	File::Remove::remove(\$recursive, $dbfile) if -f $dbfile;

	#### CREATE DATABASE
	my $dbobject = DBase::SQLite->new(
		{
			dbfile	=>	$dbfile
		}
	);

	#### CREATE TABLE sam
    my $query = qq{
CREATE TABLE IF NOT EXISTS sam (
    qname   VARCHAR(100) NOT NULL,
    flag    VARCHAR(8) NOT NULL,
    rname   VARCHAR(100) NOT NULL,
    pos     INT(20) NOT NULL,
    mapq    INT(20) NOT NULL,
    cigar   VARCHAR(255) NOT NULL,
    mrnm    VARCHAR(100) NOT NULL,
    mpos    INT(20) NOT NULL,
    isize   INT(20) NOT NULL,
    seq     VARCHAR(255),
    qual    VARCHAR(255),
    tag     VARCHAR(255),
    PRIMARY KEY (qname, rname, pos, cigar, mrnm, mpos)
)
};    
    return $dbobject->do($query);    
}



=head2

	SUBROUTINE		createSamIndexTable

	PURPOSE

		1. CREATE samindex TABLE IN SQLITE DATABASE

		2. RETURN DATABASE HANDLE

	NOTES

		(SEE createSamTable)

=cut

sub createSamIndexTable
{
	my $self		=	shift;
	my $dbfile		=	shift;



	#### REMOVE DB FILE IF EXISTS
	my $recursive = 0;
	File::Remove::remove(\$recursive, $dbfile) if -f $dbfile;

	#### CREATE DATABASE
	my $dbobject = DBase::SQLite->new(
		{
			dbfile	=>	$dbfile
		}
	);

####	#### CREATE TABLE samindex
####    my $query = qq{
####CREATE TABLE IF NOT EXISTS samindex (
####    qname   VARCHAR(100) NOT NULL,
####    rname   VARCHAR(100) NOT NULL,
####    pos     INT(20) NOT NULL,
####    mrnm    VARCHAR(100) NOT NULL,
####    mpos    INT(20) NOT NULL,
####    PRIMARY KEY (qname, rname, pos, mrnm, mpos)
####)
####};    

    my $query = qq{
CREATE TABLE IF NOT EXISTS samindex (
    qname   VARCHAR(100) NOT NULL,
    rname   VARCHAR(100) NOT NULL,
    pos     INT(20) NOT NULL,
    PRIMARY KEY (qname, rname, pos)
)
};    

    return $dbobject->do($query);
}


=head2

	SUBROUTINE		samToTsvfiles

	PURPOSE

		1. RUN ALL JOBS CONCURRENTLY IN PARALLEL

		2. FOR EACH SAM FILE, CONCATENATE 'TAG' FIELDS AT END OF LINES,

		REMOVE DUPLICATES AND PRINT TO TWO DIFFERENT TSV FILES:

			FULL - CONTAINS ALL SAM FIELDS

			INDEX (SHORT) - ONLY READ (AND MATE) MAPPING INFO

	NOTES

		SINGLE 'TAG' ENTRY:

HWI-EAS185:1:2:1491:1093#0      0       chr1    4816687 0       52M     *       0       0       ACAAANCCACCTGCCTCTGCCTCCCGAGTGCTGGGATTAAAGGCGTGCGCCA    `a\`ZDZ``aa`^aa\a\Z`aZa`aZ\ZRZ`ZYR\W\RXJ\ZN`ZT\^S`^^\    NM:i:2

		MULTIPLE 'TAG' ENTRIES:

HWI-EAS185:1:2:1076:333#0       0       chr1    4822432 255     31M4619N19M     *       0       0     GAATCTGGAATTAAACAGGCAGCAGAAACCGTAAAAGCCTTGATAGATCA      a[aabaaa]_aaaa`baa_baaaY]VZaaa_]\_XUYaa^``^]_`S^a`      NM:i:0  XS:A:+  NS:i:0

		AFTER CONVERSION

HWI-EAS185:1:2:1076:333#0       0       chr1    4822432 255     31M4619N19M     *       0       0     GAATCTGGAATTAAACAGGCAGCAGAAACCGTAAAAGCCTTGATAGATCA      a[aabaaa]_aaaa`baa_baaaY]VZaaa_]\_XUYaa^``^]_`S^a`      NM:i:0,XS:A:+,NS:i:0

=cut

sub samToTsvfiles
{
	my $self		=	shift;
	my $outputdir	=	shift;
	my $samfiles	=	shift;
	my $tsvfiles	=	shift;
	my $cluster		=	shift;
	my $label		=	shift;


	print "Align::Stats::samToTsvfiles    Align::Stats::samToTsvfiles(outputdir, samfiles, label)\n";
	print "Align::Stats::samToTsvfiles    outputdir: $outputdir\n";
	print "Align::Stats::samToTsvfiles    samfiles: ";
	print join "\n", @$samfiles;
	print "\n";
	print "Align::Stats::samToTsvfiles    no. samfiles: ", scalar(@$samfiles), "\n";
	print "Align::Stats::samToTsvfiles    label: $label\n";

	#### IF CLUSTER NOT DEFINED, RUN IT LOCALLY
	if ( not defined $cluster )
	{

#### DEBUG		
#exit;

		print "Creating TSV files...\n";
		for ( my $i = 0; $i < @$samfiles; $i++ )
		{
			my $samfile = $$samfiles[$i];
			my $tsvfile = $$tsvfiles[$i];

			print ".";

			#### GENERATE TSV FILES TO LOAD INTO SQLITE
			$self->samToTsv($samfile, $tsvfile);
		}
		print "\n";

	}

	#### RUN JOBS IN PARALLEL ON CLUSTER
	else
	{
		#### GENERATE TSV FILES TO LOAD INTO SQLITE
		print "Creating TSV files in parallel...\n";

		#### SET EXECUTABLE SCRIPT TO CONVERT SAM TO TSV
		my $executable = "$Bin/samToTsv.pl";
		print "Align::Stats::samToTsvfiles    executable: $executable\n";

		my $jobs = [];
		for ( my $i = 0; $i < @$samfiles; $i++ )
		{
			my $samfile = $$samfiles[$i];
			my $tsvfile = $$tsvfiles[$i];
			my $command = "/usr/bin/perl $executable --samfile $samfile --tsvfile $tsvfile";

			my $job = $self->setJob( [$command], "$label-$i", $outputdir );
			push @$jobs, $job;
		}

		#### RUN COMMANDS IN PARALLEL
		$self->runJobs( $jobs, $label );	
	}
}

=head2

	SUBROUTINE		samToTsvFilenames

	PURPOSE

		PARSE SAM FILE NAMES TO TSV FILE NAMES

=cut


sub samToTsvFilename
{
	my $self		=	shift;
	my $samfiles	=	shift;

	my $tsvfiles = [];
	foreach my $samfile ( @$samfiles )
	{
		my $tsvfile = $samfile;
		$tsvfile =~ s/\.sam$/.tsv/;
		push @$tsvfiles, $tsvfile;
	}

	return $tsvfiles;
}




=head2

	SUBROUTINE		samToTsv

	PURPOSE

		GENERATE TSV FILES TO LOAD INTO SQLITE:

		   1. CONCATENATE 'TAG' FIELDS AT END OF SAM LINES

		   2. REMOVE DUPLICATES

		   3. SPLIT INTO TWO FILES:

				   FULL - CONTAINS ALL SAM FIELDS

				   INDEX (SHORT) - ONLY READ (AND MATE) MAPPING INFO

		   4. RETURN THE LIST OF FULL TSV FILES

	NOTES

		SINGLE 'TAG' ENTRY:

HWI-EAS185:1:2:1491:1093#0      0       chr1    4816687 0       52M     *       0       0       ACAAANCCACCTGCCTCTGCCTCCCGAGTGCTGGGATTAAAGGCGTGCGCCA    `a\`ZDZ``aa`^aa\a\Z`aZa`aZ\ZRZ`ZYR\W\RXJ\ZN`ZT\^S`^^\    NM:i:2

		MULTIPLE 'TAG' ENTRIES:

HWI-EAS185:1:2:1076:333#0       0       chr1    4822432 255     31M4619N19M     *       0       0     GAATCTGGAATTAAACAGGCAGCAGAAACCGTAAAAGCCTTGATAGATCA      a[aabaaa]_aaaa`baa_baaaY]VZaaa_]\_XUYaa^``^]_`S^a`      NM:i:0  XS:A:+  NS:i:0

		AFTER CONVERSION

HWI-EAS185:1:2:1076:333#0       0       chr1    4822432 255     31M4619N19M     *       0       0     GAATCTGGAATTAAACAGGCAGCAGAAACCGTAAAAGCCTTGATAGATCA      a[aabaaa]_aaaa`baa_baaaY]VZaaa_]\_XUYaa^``^]_`S^a`      NM:i:0,XS:A:+,NS:i:0

=cut

sub samToTsv
{
	my $self		=	shift;
	my $samfile		=	shift;
	my $tsvfile		=	shift;
	my $clean		=	shift;


	if ( not defined $tsvfile )
	{
		$tsvfile = $samfile;
		$tsvfile =~ s/\.sam$/.tsv/;
	}

	if ( -f $tsvfile and not -z $tsvfile and not $clean )
	{
		return;
	}

	#### SET TEMP FILE (MAY CONTAIN DUPLICATES)
	my $tempfile = "$tsvfile.tmp";

	#### SET INDEX FILE AND INDEX TEMP FILE (MAY CONTAIN DUPLICATES)
	my $indexfile	= $tsvfile;
	$indexfile =~ s/\.tsv$/-index.tsv/;
	my $tempindexfile = "$indexfile.tmp";

	#### OPEN SAMS
	open(SAM, $samfile) or die "Can't open samfile: $samfile\n";
	open(INDEX, ">$tempindexfile") or die "Can't open tempindexfile: $tempindexfile\n";
	open(TSV, ">$tempfile") or die "Can't open tempfile: $tempfile\n";
	while(<SAM>)
	{	
		my @elements = split "\t", $_;

		my $tag_index = 11;
		my $tag = '';
		while ( $tag_index < $#elements + 1 )
		{
			$tag .= $elements[$tag_index] . ",";
			$tag_index++;
		}
		$tag =~ s/,$//;
		$#elements = 10;
		push @elements, $tag;
		print TSV join "\t", @elements;

		print INDEX "$elements[0]\t$elements[2]\t$elements[3]\n";
	}
	close(SAM) or die "Can't close sam file: $samfile\n";
	close(TSV) or die "Can't close tsv file: $tempfile\n";
	close(INDEX) or die "Can't close index file: $tempindexfile\n";

	#### REMOVE DUPLICATES
	my $command = "uniq -u $tempfile > $tsvfile";
	`$command`;

	#### REMOVE DUPLICATES
	$command = "uniq -u $tempindexfile > $indexfile";
	`$command`;

	print "Tophat::samToTsv    tsvfile printed:\n$tsvfile\n";
	print "Tophat::samToTsv    indexfile printed:\n$indexfile\n";

	return $tsvfile;
}







=head2

	SUBROUTINE		splitReads

	PURPOSE

		1. IF SPLIT READS FILE IS NOT PRESENT OR EMPTY:

			-	COUNT TOTAL NUMBER OF READS PER INPUT SPLIT FILE

		2. OTHERWISE, GATHER THE ABOVE READ COUNTS FROM THE

			SPLIT READS FILE

		3. 	RETURN HASH OF FILENAMES VS. READ COUNTS

=cut

sub splitReads
{	
	my $self			=	shift;
	my $splitreadsfile	=	shift;
	my $splitfiles		=	shift;
	my $clean			=	shift;



	my $reads = {};
	my $total_reads = 0;
	if ( defined $clean or not -f $splitreadsfile )
	{
		open(OUT, ">$splitreadsfile") or die "Can't open splitreadsfile: $splitreadsfile\n";
		foreach my $filepair ( @$splitfiles )
		{
			foreach my $file ( @$filepair )
			{
				my $lines = $self->countLines($file);
				my $read_count = $lines / 4;
				$total_reads += $read_count;
				$reads->{splitfiles}->{$file} = $read_count;
				print OUT "$file\t$read_count\n";
			}
		}
		close(OUT) or die "Can't close splitreadsfile: $splitreadsfile\n";
		print "Align::Stats::splitReads    splitreadsfile printed:\n\n$splitreadsfile\n\n";
	}
	else
	{
		open(FILE, $splitreadsfile) or die "Can't open splitreadsfile: $splitreadsfile\n";
		while ( <FILE> )
		{
			next if $_ =~ /^\s*$/;
			my ($file, $read_count) = /^(\S+)\s+(\d+)/;

			$reads->{splitfiles}->{$file} = $read_count;
			$total_reads += $read_count;
		}
		close(FILE) or die "Can't close splitreadsfile: $splitreadsfile\n";
	}
	$reads->{total} = $total_reads;


	return $reads;
}




=head2

	SUBROUTINE		splitHits

	PURPOSE

		1. IF SPLIT HITS FILE IS NOT PRESENT OR EMPTY:

			-	COUNT THE HITS PER CHROMOSOME, BREAKDOWN BY

				INPUT SPLIT FILE

		2. OTHERWISE, GATHER THE ABOVE HIT COUNTS FROM THE

				SPLIT HITS FILE

		3. 	RETURN HASH OF FILENAMES VS. HIT COUNTS

=cut

sub splitHits
{	
	my $self			=	shift;
	my $splithitsfile	=	shift;
	my $splitfiles		=	shift;
	my $references		=	shift;
	my $clean			=	shift;



	my $hits = {};
	my $total_hits = 0;
	if ( defined $clean or not -f $splithitsfile )
	{

		open(OUT, ">$splithitsfile") or die "Can't open splithitsfile: $splithitsfile\n";

		my $samfiles;
		for my $reference ( @$references )
		{	
			#### GATHER STATS FROM ALL SPLITFILES
			my $total_hits = 0;
			foreach my $splitfile ( @$splitfiles )
			{
				#### SET *.sam FILE
				my ($basedir, $index) = $$splitfile[0] =~ /^(.+?)\/(\d+)\/([^\/]+)$/;
				my $outputdir = "$basedir/$reference";
				my $samfile = "$outputdir/$index/accepted_hits.sam";
				push @$samfiles, $samfile;

				my $hit_count = $self->countLines($samfile);
				$hits->{$reference}->{splitfiles}->{$$splitfile[0]} = $hit_count;
				print OUT "$reference\t$$splitfile[0]\t$hit_count\n";

				$total_hits += $hit_count;
			}

			#### SAVE TOTAL HITS COUNT
			$hits->{$reference}->{totalhits} = $total_hits;
		}	
		close(OUT) or die "Can't close splithitsfile: $splithitsfile\n";
		print "Align::Stats::readHits    splithitsfile printed:\n\n$splithitsfile\n\n";
	}
	else
	{
		open(FILE, $splithitsfile) or die "Can't open splithitsfile: $splithitsfile\n";
		while ( <FILE> )
		{
			next if $_ =~ /^\s*$/;
			my ($file, $hit_count) = /^\S+\s+(\S+)\s+(\d+)/;
			$hits->{splitfiles}->{$file} = $hit_count;
			$total_hits += $hit_count;
		}
		close(FILE) or die "Can't close splithitsfile: $splithitsfile\n";
	}


	return $hits;
}




=head2

	SUBROUTINE		countLines

	PURPOSE

		QUICKLY COUNT AND RETURN THE NUMBER OF LINES IN A FILE

=cut

sub countLines
{
	my $self		=	shift;
	my $file		=	shift;

	open(FILE, $file) or die "Can't open file: $file\n";
	my $counter = 0;
	while(<FILE>)
	{
		$counter++;
	}
	close(FILE);

	return $counter;
}



################################################################################
##################			HOUSEKEEPING SUBROUTINES			################
################################################################################


=head2

	SUBROUTINE		new

	PURPOSE

		CREATE THE NEW self OBJECT AND INITIALISE IT, FIRST WITH DEFAULT 

		ARGUMENTS, THEN WITH PROVIDED ARGUMENTS

=cut

sub new
{
 	my $class 		=	shift;
	my $arguments 	=	shift;

	my $self = {};
    bless $self, $class;

	#### SET DEFAULT TEMP DIR
	$self->{_tempdir} = "/tmp";

	#### INITIALISE THE OBJECT'S ELEMENTS
	$self->initialise($arguments);

	#### SET DEBUG IF verbose

    return $self;
}

=head2

	SUBROUTINE		initialise

	PURPOSE

		INITIALISE THIS OBJECT:

			1. LOAD THE ARGUMENTS

			2. CALL PARENT initialise TO LOAD ADDITIONAL INFO

=cut
sub initialise
{
    my $self		=	shift;
	my $arguments	=	shift;


    #### VALIDATE USER-PROVIDED ARGUMENTS
	($arguments) = $self->validate_arguments($arguments, $DATAHASH);	

	#### LOAD THE USER-PROVIDED ARGUMENTS
	foreach my $key ( keys %$arguments )
	{		
		#### LOAD THE KEY-VALUE PAIR
		$self->value($key, $arguments->{$key});
	}

	#### REQUIRED: CALL THE PARENT FOR ADDITIONAL INITIALISATION
	$self->SUPER::initialise();	
}



1;

__END__

sub mergeUniqueReads
{
	my $self				=	shift;
	my $outputdir			=	shift;
	my $tsvfiles			=	shift;
	my $sqlite				=	shift;
	my $cluster				=	shift;
	my $label				=	shift;



	my $duplicatefile = "$outputdir/duplicate.txt";
	my $uniquefile = "$outputdir/unique.txt";
	File::Remove::rm($duplicatefile);
	File::Remove::rm($uniquefile);

	my $tsvindexfiles = [];
	foreach my $tsvfile ( @$tsvfiles )
	{
		my $tsvindexfile = $tsvfile;
		$tsvindexfile =~ s/\.tsv$/-index.tsv/;
		push @$tsvindexfiles, $tsvindexfile;
	}

	#### MERGE DUPLICATE READS INTO ONE FILE FOR ALL GENOME HITS
	my $commands = [];
	my $counter = 0;
	my $total = scalar(@$tsvfiles);
	foreach my $tsvfile ( @$tsvfiles )
	{
		$counter++;
		push @$commands, "echo 'Loading $counter of $total into duplicatefile: $duplicatefile'";
		push @$commands, "cut -f 1 $tsvfile >> $duplicatefile";
	}


	##### GET UNIQUE READS FOR ALL GENOME
	#my $uniq_command = qq{sort $duplicatefile | uniq -u > $uniquefile};
	#push @$commands, "echo '$uniq_command'";
	#push @$commands, $uniq_command;


	my $job = $self->setJob( $commands, "SAMINDEX", $outputdir );

	#### RUN COMMANDS IN SERIES ON ONE CLUSTER NODE
	$self->runJobs( [ $job ], $label );	




}



=head2

	SUBROUTINE		readHits

	PURPOSE

		1. GENERATE Align::Stats REPORTS IN PARALLEL ON A CLUSTER

			- GET TOTAL INPUT READ COUNT

			- GET TOTAL UNIQUE READ HIT COUNTS PER CHROMOSOME

			- GET TOTAL UNIQUE READ HIT COUNT FOR WHOLE GENOME

=cut

sub readHits
{
	my $self			=	shift;
	my $outputdir		=	shift;
	my $references		=	shift;
	my $splitfile		=	shift;

	my $clean			=	$self->get_clean();
	my $sqlite			=	$self->get_sqlite();
	my $cluster			=	$self->get_cluster();


	print "Align::Stats::readHits    Align::Stats::readHits(outputdir, splitfiles, references, clean, sqlite, cluster)\n";
	print "Align::Stats::readHits    outputdir: $outputdir\n";
	print "Align::Stats::readHits    references: @$references\n";
	print "Align::Stats::readHits    clean: $clean\n" if defined $clean;
	print "Align::Stats::readHits    cluster: $cluster\n" if defined $cluster;

	#### SET SQLITE DEFAULT
	$sqlite = "sqlite3" if not defined $sqlite or not $sqlite;

	#### GET splitfile CONTAINING LIST OF SPLIT INPUT FILES
	$splitfile = "$outputdir/splitfiles" if not defined $splitfile or not $splitfile;
	print "Align::Stats::readHits    Can't find splitfile: $splitfile\n" and exit if not -f $splitfile;
	print "Align::Stats::readHits    splitfile is empty: $splitfile\n" and exit if -z $splitfile;

	#### GET LIST OF SPLIT INPUT FILES
	print "Align::Stats::readHits    Getting split files\n";
	my $splitfiles = $self->splitfilesList($splitfile);
	print "Align::Stats::readHits    No splitfiles produced. Quitting\n" and exit if not defined $splitfiles or scalar(@$splitfiles) == 0;

	#### HASH OF SPLIT FILE NAMES VS. COUNTS
	my $stats = {};

	#### SET SPLIT FILE reads FILE
	print "Align::Stats::readHits    Getting stats->{reads}\n";
	my $splitreadsfile = "$splitfile.reads";
	$stats->{reads} = $self->splitReads($splitreadsfile, $splitfiles, $clean);

	#### SET SPLIT FILE hits FILE IF NOT DEFINED
	#### GET THE HITS PER CHROMOSOME, BREAKDOWN BY INPUT SPLIT FILE
	print "Align::Stats::readHits    Getting stats->{hits}\n";
	my $splithitsfile = "$splitfile.hits";
	$stats->{hits} = $self->splitHits($splithitsfile, $splitfiles, $references, $clean);

	print "Align::Stats::readHits    Gathering samfiles, tsvfiles and reference_samfiles hash\n";
	my ($samfiles, $tsvfiles, $reference_samfiles) = $self->gatherFiles($references, $splitfiles);
	#exit;

	####	1. CONVERT SAM FILES TO *.tsv AND *-index.tsv FILES
	my ($label) = $outputdir =~ /\/([^\/]+)$/;
	#$self->samToTsvfiles($outputdir, $samfiles, $tsvfiles, $cluster, $label);	

	####	2. LOAD ALL chr*/accepted_hits.tsv FILES INTO sam TABLE,
	####		ONE SQLITE DATABASE PER REFERENCE CHROMOSOME
#	$self->loadSamTables($outputdir, $references, $reference_samfiles, $sqlite, $cluster, $label);


	####	3. LOAD ALL chr*/accepted_hits-index.tsv FILES INTO samindex TABLE
	####		IN A SINGLE SQLITE DATABASE FOR ALL REFERENCE CHROMOSOMES
	#$self->loadSamIndexTable($outputdir, $tsvfiles, $sqlite, $cluster, $label);


	#### 	4. QUERY COUNT ALL DISTINCT READ IDS
	print "Align::Stats::readHits    Doing readHits()\n";
	$self->readHits($outputdir, $tsvfiles, $sqlite, $cluster, $label);
	print "Align::Stats::readHits    Completed readHits()\n";


return;


	#### 	5. QUERY IDENTITIES AND HIT COUNTS OF READS WITH MOST HITS


	#### 	6. EXTRACT READ SEQUENCES OF 100 'MOST HITS' READS


	#### PRINT STATS
	my $outfile = "$outputdir/TOPHAT.stats";
	open(OUTFILE, ">$outfile") or die "Can't open output file: $outfile\n";

	close(OUTFILE) or die "Can't close output file: $outfile\n";
}
