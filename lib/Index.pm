package IndexRead;

#use Moose;
#use MooseX::FollowPBP;


=head2

	PACKAGE		IndexRead

	PURPOSE

		FAST RETRIEVAL OF READ HIT POSITIONS GIVEN THE READ NAME

	    APPLICATION     indexReads

    VERSION         0.01

    PURPOSE

        FAST RETRIEVAL OF READ HIT POSITIONS GIVEN THE READ NAME

    INPUT

        1. *.bam BAM-FORMAT FILE

    OUTPUT

        1. *.bam.coverage COVERAGE FILE OF WINDOW POSITIONS AND READ COVERAGES

		0	12
		100	23
		200	45
		300	32
		...

    NOTES

The tie() function attaches a hash name to an on-disk database file.

Perl supports several simple database formats through this mechanism. The one I chose is the "Berkeley DB" package, which supports B-Tree files on disk.
Because the C library implementing this file format is freely available, most Perl interpreters support it.
Plus, Berkeley DB supports any byte ordering, so you can move database files across platforms with impunity.
Berkeley DB supports C, C++, Java, and Perl APIs, for both UNIX and Windows 95/NT.


Building the Index Database

Perl's string-handling abilities make the basic indexing chore pretty easy. Listing One is indexer_basic.pl, a complete program to build the index.

The IndexWords subroutine first uses Perl's split function to break the $words string into an array of separate words. It then uses Perl's grep command to discard certain nonwords and remove duplicates from the list. Finally, it simply iterates over the words, appending the file key to the database entry for each word.


About B-Trees
Recall that a B-Tree file consists of pages stored on disk. Each page is either an internal page, which has keys and pointers to other pages, or a leaf page, which stores keys and their associated data.

To find a particular entry, you start with the top-level page, usually an internal page. You compare the desired key to each key on the page. Your key will fall between two of the keys on that page, which determines the next page number to check. You may walk through several internal pages before you arrive at a leaf page that contains the data you want. Figure 1 illustrates this process.

The depth of this tree depends on the size of the pages. If you increase the page size, you can reduce the depth of the tree, thereby reducing the number of disk accesses needed to retrieve an item, and speeding up your program. In my case, I found that a 32-KB page size was sufficient to produce a two-level tree; any item in my database can be found with only two disk accesses. Better yet, since the top page remains in the database cache, only the first database request requires more than a single disk access. This juggling of page size was especially important because I wanted to put the final database on a CD-ROM. CD-ROMs are exceptionally slow at seeking, which makes it especially important to minimize the number of pages that must be read, since each page read requires a seek.

Because a B-Tree stores items in sorted order, you can speed up multiple searches tremendously by doing your searches in sorted order. If your next search happens to land on the same page, you won't have to go to disk at all, since the required pages will already be in the cache.

Performance
To try out my indexing tools, I used as a testbed 10-years worth of DDJ articles in HTML format. Even with this large collection, search time was simply never an issue.
Building the index, however, presented some challenges.

As Table 1 shows, the indexer in Listing One has a serious speed problem. It required over three hours to index my sample corpus. However, it used only 14 minutes of CPU time, which suggests that the program is spending most of its time waiting on disk I/O.

Every time you add a word to the database, you have to find that word in the database first, then write it back with an additional key. Since the replacement is longer, the update is likely to require moving data to accommodate it. In addition, the words in a single file probably span much of the database, so if you're indexing a large number of files, you quickly get to the point where you're skimming through several megabytes of the database for every file you index.

To address this, the current version of my indexer defines %wordCache to be an ordinary in-memory hash, and adds new word entries to that. After indexing 500 files, the indexer calls FlushWordCache (Listing Two) to update the disk database from the in-memory hash. With this change, the total running time dropped from over three hours to a reasonable 20 minutes. The CPU time was practically unchanged.


    NB: IF THE APPLICATION CRASHES WHEN THE B-Tree HAS MANY PAGES (E.G., 1000+ PAGES), INCREASING THE PAGE SIZE MAY HELP.


=cut

use strict;
use warnings;
use Carp;

#### INTERNAL MODULES
use Sampler;

#### HAS A
use Cluster;

#### EXTERNAL MODULES
use Data::Dumper;
use File::Path;
use DB_File;
$DB_File::DB_BTREE->{cachesize} = 10_000_000; # 10meg cache
$DB_File::DB_BTREE->{psize} = 32*1024; # 32k pages
use Fcntl;      				#### Needed for DB_File
use Fcntl qw(O_RDWR O_CREAT);	#### import the O_CREAT and O_RDWR constants


require Exporter;
our @ISA = qw(Exporter Cluster);
#our @EXPORT_OK = qw();
our $AUTOLOAD;


#### SET SLOTS
our @DATA = qw(

INPUTFILES
DBFILE
OUTPUTFILE
READNAME
SAMTOOLS

MAXJOBS
CPUS
CLUSTER
QUEUE
WALLTIME
QSTAT
QSUB
SLEEP
CLEANUP
VERBOSE
TEMPDIR
DOT

BATCHSTATS
COMMAND
);
our $DATAHASH;
foreach my $key ( @DATA )	{	$DATAHASH->{lc($key)} = 1;	}


=head2

	SUBROUTINE		getValue

	PURPOSE

		RETURN AN ARRAY OF DUPLICATE VALUES FOR THE GIVEN KEY

	INPUTS

		1. DBTREE DB FILE HANDLE

		2. KEY

=cut

sub getValue
{
	my $self	=	shift;
	my $db		=	shift;
	my $key		=	shift;

	#$self->printOut($db);

	my $value;	
	#my $got = $db->get_dup($key);
	my $got = $db->get($key, $value);

	#my $status = $db->get_dup($key);

	return $value;
}

### PRINT KEY/VALUE PAIRS BY ITERATING THROUGH BTREE USING seq
sub printOut
{
	my $self	=	shift;
	my $db		=	shift;
	my $limit	=	shift;

	print "IndexRead::printOut    IndexRead::printOut(db, limit)\n";
	print "IndexRead::printOut    db: $db\n";
	print "IndexRead::printOut    limit: $limit\n" if defined $limit;

	my $counter = 0;
	my $status;
	my $key = 0;
	my $value = 0;
	for ($status = $db->seq($key, $value, R_FIRST) ;
		 $status == 0 ;
		 $status = $db->seq($key, $value, R_NEXT) )
	{
		$counter++;
		print "$key\t$value\n";
		last if defined $limit and $counter >= $limit;
	}

}

=head2

	SUBROUTINE		compare

	PURPOSE

		COMPARE ONLY THE KEYS OF TWO INPUT DB FILES

	INPUTS

		1. A DB FILE CONTAINING KEYS FOUND ONLY IN THE QUERY DB

		2. A DB FILE CONTAINING KEYS FOUND IN BOTH QUERY AND TARGET DB

	OUTPUTS

		1. LOCATIONS OF QUERY AND TARGET DFILES

		2. LOCATIONS OF OUTPUT DB FILES (QUERYONLY, QUERYBOTH)

=cut

sub indexVenn
{
	my $self		=	shift;	
	my $queryfile	=	shift;
	my $targetfile	=	shift;
	my $queryonlyfile	=	shift;
	my $querybothfile	=	shift;

	print "IndexRead::compare    queryfile: $queryfile\n";
	print "IndexRead::compare    targetfile: $targetfile\n";

	my $queryonly_counter = 0;
	my $queryboth_counter = 0;

	#### CHECK INPUTS
	print "IndexRead::compare   queryfile not defined.\n" and exit if not defined $queryfile;
	print "IndexRead::compare   targetfile not defined.\n" and exit if not defined $targetfile;
	print "IndexRead::compare   queryonlyfile not defined.\n" and exit if not defined $queryonlyfile;
	print "IndexRead::compare   querybothfile not defined.\n" and exit if not defined $querybothfile;

	# ENABLE DUPLICATE RECORDS
	$DB_BTREE->{'flags'} = R_DUP;

	sub Compare
    {
		my ($key1, $key2) = @_ ;
		my ($key1_number1, $key1_number2, $key1_number3) = $key1 =~ /^\D+(\d+)\.(\d+).+\/(\d)$/;
		my ($key2_number1, $key2_number2, $key2_number3) = $key2 =~ /^\D+(\d+)\.(\d+).+\/(\d)/;
		return -1 if $key1_number1 < $key2_number1;
		return 1 if $key1_number1 > $key2_number1;
		return -1 if $key1_number2 < $key2_number2;
		return 1 if $key1_number2 > $key2_number2;

		if ( defined $key1_number3 and defined $key2_number3 )
		{
			return -1 if $key1_number3 < $key2_number3;
			return 1 if $key1_number3 > $key2_number3;
		}

        return 0;
    }

    #### COMPARISON SUBROUTINE
    $DB_BTREE->{'compare'} = \&Compare;


	#### TIE TARGET HASH TO DB FILE	
	my $targethash;
	my $targetdb = tie %$targethash, "DB_File", $targetfile, O_RDWR|O_CREAT, 0666, $DB_BTREE or die "Cannot open targetfile: $targetfile: $!\n";

	#### TIE QUERY HASH TO INPUT DB FILE	
	my $queryhash;
	my $querydb = tie %$queryhash, "DB_File", $queryfile, O_RDWR|O_CREAT, 0666, $DB_BTREE or die "Cannot open queryfile: $queryfile: $!\n";

	#$self->printOut($querydb);

	#$self->printOut($targetdb);


	### CHECK FOR MATCHING KEYS IN TARGET DATABASE
	my $total = 0;
	my $matched = 0;

	my ($key, $value, $status);
#	$key = $value = 0 ;
#	# valid flag values are:
#	# R_CURSOR, R_FIRST, R_LAST, R_NEXT and R_PREV.
    for ($status = $querydb->seq($key, $value, R_FIRST) ;
         $status == 0 ;
         $status = $querydb->seq($key, $value, R_NEXT) )
	{
		$total++;

#	foreach my $key (keys %$queryhash)
#    {
		#my $value = $queryhash->{$key};

		#### CHECK FOR UNIQUE KEYS ONLY
		#my $found = $targetdb->get_dup($key);

		#### CHECK FOR UNIQUE KEY-VALUE PAIRS
		my $found = $targetdb->find_dup($key, $value);
		print "key: $key found: $found\n";

		if ( defined $found )	
		{
			$matched++;
		}
		else
		{
		}
	}
	print "IndexRead::compare    matched: $matched\n" ;

	undef $querydb;
	untie $queryhash;

	undef $targetdb;
	untie $targethash;
}





=head2

	SUBROUTINE		getReads

	PURPOSE

		RETURN THE TOTAL NUMBER OF UNIQUE READS IN THE DB FILE

=cut

sub getReads
{
	my $self		=	shift;	
	my $start		=	shift;
	my $stop		=	shift;

	my $dbfile		=	$self->get_dbfile();

	$start = 0 if not defined $start;


	#### CHECK INPUTS
	print "IndexRead::getReads   dbfile not defined. Exiting\n"
		and exit if not defined $dbfile;

#    #### MEMORY USAGE BEFORE    
#	warn "memory usage BEFORE hash:\n";
#    warn qx{ ps -o rss,vsz $$ }, "\n";

    #### SET TIMER
	use Time::HiRes qw[gettimeofday tv_interval];
	my $time1=[gettimeofday()];	
	my $milliseconds  = 0;

	#### TIE HASH TO EXISTING DB FILE
	my $hash;
	my $db = tie %$hash, "DB_File", $dbfile, O_RDWR|O_CREAT, 0666, $DB_BTREE 
		or die "Cannot open $dbfile: $!\n";


	### PRINT KEY/VALUE PAIRS BY ITERATING THROUGH BTREE USING seq
	my $total = 0;
	my ($key, $value, $status);
	$key = $value = 0 ;
    for ($status = $db->seq($key, $value, R_FIRST) ;
         $status == 0 ;
         $status = $db->seq($key, $value, R_NEXT) )
	{
		$total++;
		next if $total < $start;
		print "$key\t$value\n";
	}


	#### DELETE DB OBJECT REFERENCE
	undef $db;

	#### UNTIE hash FROM DB FILE
	untie(%$hash);

	return $total;
}


=head2

	SUBROUTINE		mergeIndexes

	PURPOSE

		RETURN THE TOTAL NUMBER OF UNIQUE READS IN THE DB FILE

	INPUTS

		1. LIST OF INPUT DB FILE LOCATIONS

		2. LOCATION OF OUTPUT DB FILE

	OUTPUTS

		1. LOCATION OF OUTPUT DB FILE

=cut

sub mergeIndexes
{
	my $self		=	shift;	

	my $inputfiles	=	$self->get_inputfiles();
	my $outputfile	=	$self->get_outputfile();


	#### CHECK INPUTS
	print "IndexRead::mergeIndexes   outputfile not defined. Exiting\n"
		and exit if not defined $outputfile;

	# ENABLE DUPLICATE RECORDS
	$DB_BTREE->{'flags'} = R_DUP;

	#### TIE HASH TO EMPTY DB FILE	
	my $hash;
	my $db = tie %$hash, "DB_File", $outputfile, O_RDWR|O_CREAT, 0666, $DB_BTREE 
	or die "Cannot open $outputfile: $!\n";

	##### DEBUG
	#my $max = 100;
	#my $counter = 0;
	#my ($key, $value, $status);
	#$key = $value = 0 ;
	#for ($status = $db->seq($key, $value, R_FIRST) ;
	#	 $status == 0 ;
	#	 $status = $db->seq($key, $value, R_NEXT) )
	#{
	#	$counter++;		
	#	last if $counter >= $max;
	#	my @entries = $db->get_dup($key);
	#	foreach my $entry ( @entries )
	#	{
	#	}
	#}


	##### SET TIMER
	#use Time::HiRes qw[gettimeofday tv_interval];
	#my $time1=[gettimeofday()];	

	foreach my $inputfile ( @$inputfiles )
	{
		print "\n\nIndexRead::mergeIndexes    Loading inputfile: $inputfile\n";

		##### MEMORY USAGE BEFORE    
		#warn "IndexRead::mergeIndexes    memory usage BEFORE counter:\n";
		#warn qx{ ps -o rss,vsz $$ }, "\n";

		#### TIE HASH TO INPUT FILE	
		print "IndexRead::mergeIndexes    tying hash\n";
		my $filehash;
		my $subdb = tie %$filehash, "DB_File", $inputfile, O_RDWR|O_CREAT, 0666, $DB_BTREE 
		or die "Cannot open inputfile: $inputfile: $!\n";
		print "IndexRead::mergeIndexes    hash tying completed\n";



		#### MERGE INPUTFILE WITH HASH
		### PRINT KEY/VALUE PAIRS BY ITERATING THROUGH BTREE USING seq
		print "IndexRead::mergeIndexes    Merging database entries\n";
		my $counter = 0;
		my $dot = 10000;
		my $entries = 0;
		my ($key, $value, $status);
		$key = $value = 0 ;
		for ($status = $subdb->seq($key, $value, R_FIRST) ;
			 $status == 0 ;
			 $status = $subdb->seq($key, $value, R_NEXT) )
		{
			$entries++;
			next if exists $hash->{$key};
			$hash->{$key} = $value;
			$counter++;
			print "$counter\n" if $counter % $dot == 0;
		}
		print "IndexRead::mergeIndexes    original entries: $entries\n";
		print "IndexRead::mergeIndexes    added counter: $counter\n";


		##### MERGE INPUTFILE WITH HASH
		####
		#### NB:
		#### INSERT NEW READ ENTRY IF THE HIT LOCATION DIFFERS FROM
		#### ALL EXISTING LOCATIONS
		#### 
		#### PRINT KEY/VALUE PAIRS BY ITERATING THROUGH BTREE USING seq
		#my $counter = 0;
		#my $dot = 100000;
		#my $entries = 0;
		#my ($key, $value, $status);
		#$key = $value = 0 ;
		#for ($status = $subdb->seq($key, $value, R_FIRST) ;
		#	 $status == 0 ;
		#	 $status = $subdb->seq($key, $value, R_NEXT) )
		#{
		#	##next if exists $existing->{$key};
		#	#$existing->{$key} = 1;
		#	$counter++;
		#	
		#	my $entries;
		#	my $existing;
		#	@$entries = $subdb->get_dup($key);
		#	%$existing = $db->get_dup($key, 1);
		#	
		#	
		#	foreach my $entry ( @$entries )
		#	{
		#		#print "    skipped $key\t$entry\n" if exists $existing->{$entry};
		#		next if exists $existing->{$entry};
		#		#print "loaded $key\t$entry\n";
		#		$hash->{$key} = $entry;
		#		$entries++;
		#	}
		#}
		#



		##### MEMORY USAGE AFTER
		#warn "IndexRead::mergeIndexes    memory usage AFTER counter:\n";
		#warn qx{ ps -o rss,vsz $$ }, "\n";

		##### TIME AFTER
		#my $milliseconds = tv_interval($time1)*1000;

		##### RESET TIMER
		#$time1=[gettimeofday()];	

		##### DEBUG
		#$max = 100;
		#$counter = 0;
		#$key = $value = 0 ;
		#for ($status = $db->seq($key, $value, R_FIRST) ;
		#	 $status == 0 ;
		#	 $status = $db->seq($key, $value, R_NEXT) )
		#{
		#	$counter++;		
		#	last if $counter >= $max;
		#	my @entries = $db->get_dup($key);
		#	foreach my $entry ( @entries )
		#	{
		#	}
		#}


	}

	##### DEBUG
	#$max = 100;
	#$counter = 0;
	#$key = $value = 0 ;
	#for ($status = $db->seq($key, $value, R_FIRST) ;
	#	 $status == 0 ;
	#	 $status = $db->seq($key, $value, R_NEXT) )
	#{
	#	$counter++;		
	#	last if $counter >= $max;
	#	my @entries = $db->get_dup($key);
	#		foreach my $entry ( @entries )
	#		{
	#		}
	#}

	#### DELETE DB OBJECT REFERENCE
	undef $db;

	#### UNTIE hash FROM DB FILE
	untie(%$hash);
}





=head2

	SUBROUTINE		numberReads

	PURPOSE

		RETURN THE TOTAL NUMBER OF UNIQUE READS IN THE DB FILE

=cut

sub numberReads
{
	my $self		=	shift;	

	my $dbfile		=	$self->get_dbfile();

	#### CHECK INPUTS
	print "IndexRead::numberReads   dbfile not defined. Exiting\n"
		and exit if not defined $dbfile;

    #### MEMORY USAGE BEFORE    
	warn "memory usage BEFORE tie:\n";
    warn qx{ ps -o rss,vsz $$ }, "\n";

    #### SET TIMER
	use Time::HiRes qw[gettimeofday tv_interval];
	my $time1=[gettimeofday()];	
	my $milliseconds  = 0;

	# ENABLE DUPLICATE RECORDS
	$DB_BTREE->{'flags'} = R_DUP;

	no warnings;
	sub Compare
    {
		my ($key1, $key2) = @_ ;
		my ($key1_number1, $key1_number2, $key1_number3) = $key1 =~ /^\D+(\d+)\.(\d+).+\/(\d)$/;
		my ($key2_number1, $key2_number2, $key2_number3) = $key2 =~ /^\D+(\d+)\.(\d+).+\/(\d)/;
		return -1 if $key1_number1 < $key2_number1;
		return 1 if $key1_number1 > $key2_number1;
		return -1 if $key1_number2 < $key2_number2;
		return 1 if $key1_number2 > $key2_number2;

		if ( defined $key1_number3 and defined $key2_number3 )
		{
			return -1 if $key1_number3 < $key2_number3;
			return 1 if $key1_number3 > $key2_number3;
		}

        return 0;
    }
	use warnings;

    #### COMPARISON SUBROUTINE
    $DB_BTREE->{'compare'} = \&Compare;

	#### TIE HASH TO EXISTING DB FILE
	my $hash;
	my $db = tie %$hash, "DB_File", $dbfile, O_RDWR|O_CREAT, 0666, $DB_BTREE 
		or die "Cannot open $dbfile: $!\n";

    #### GET TIME
	$milliseconds = tv_interval($time1)*1000;
	$time1=[gettimeofday()];
	printf "Elapsed time (milliseconds): $milliseconds\n";

    #### MEMORY USAGE BEFORE    
	warn "memory usage AFTER tie:\n";
    warn qx{ ps -o rss,vsz $$ }, "\n";

	#$self->printOut($db);

	### PRINT KEY/VALUE PAIRS BY ITERATING THROUGH BTREE USING seq
	my $existing;
	my $total = 0;
	my ($key, $value, $status);
	$key = $value = 0 ;
    for ($status = $db->seq($key, $value, R_FIRST) ;
         $status == 0 ;
         $status = $db->seq($key, $value, R_NEXT) )
	{
		$total++;
	}

    #### GET TIME
	$milliseconds = tv_interval($time1)*1000;
	$time1=[gettimeofday()];
	printf "Elapsed time (milliseconds): $milliseconds\n";

    #### MEMORY USAGE BEFORE    
	warn "memory usage AFTER count:\n";
    warn qx{ ps -o rss,vsz $$ }, "\n";

	#### DELETE DB OBJECT REFERENCE
	undef $db;

	#### UNTIE hash FROM DB FILE
	untie(%$hash);

	return $total;
}




=head2

	SUBROUTINE		numberHits

	PURPOSE

		RETURN THE TOTAL NUMBER OF UNIQUE READS IN THE DB FILE

=cut

sub numberHits
{
	my $self		=	shift;	

	my $dbfile		=	$self->get_dbfile();


	#### CHECK INPUTS
	print "IndexRead::numberHits   dbfile not defined. Exiting\n"
		and exit if not defined $dbfile;

    #### MEMORY USAGE BEFORE    
	warn "memory usage BEFORE tie:\n";
    warn qx{ ps -o rss,vsz $$ }, "\n";

    #### SET TIMER
	use Time::HiRes qw[gettimeofday tv_interval];
	my $time1=[gettimeofday()];	
	my $milliseconds  = 0;

	## ENABLE DUPLICATE RECORDS
	#$DB_BTREE->{'flags'} = R_DUP;

#	no warnings;
#	sub Compare
#    {
#		my ($key1, $key2) = @_ ;
#		my ($key1_number1, $key1_number2, $key1_number3) = $key1 =~ /^\D+(\d+)\.(\d+).+\/(\d)$/;
#		my ($key2_number1, $key2_number2, $key2_number3) = $key2 =~ /^\D+(\d+)\.(\d+).+\/(\d)/;
#		return -1 if $key1_number1 < $key2_number1;
#		return 1 if $key1_number1 > $key2_number1;
#		return -1 if $key1_number2 < $key2_number2;
#		return 1 if $key1_number2 > $key2_number2;
#		
#		if ( defined $key1_number3 and defined $key2_number3 )
#		{
#			#print "Compare mate pair: $key1_number3 vs $key2_number3\n";
#			return -1 if $key1_number3 < $key2_number3;
#			return 1 if $key1_number3 > $key2_number3;
#		}
#		
#        return 0;
#    }
#	use warnings;
#
#    #### COMPARISON SUBROUTINE
#    $DB_BTREE->{'compare'} = \&Compare;

	#### TIE HASH TO EXISTING DB FILE
	my $hash;
	my $db = tie %$hash, "DB_File", $dbfile, O_RDWR|O_CREAT, 0666, $DB_BTREE 
		or die "Cannot open $dbfile: $!\n";

    #### GET TIME
	$milliseconds = tv_interval($time1)*1000;
	$time1=[gettimeofday()];
	printf "Elapsed time (milliseconds): $milliseconds\n";

    #### MEMORY USAGE BEFORE    
	warn "memory usage AFTER tie:\n";
    warn qx{ ps -o rss,vsz $$ }, "\n";

	$self->printOut($db);

	my $total = 0;
	my %seen;
	foreach my $key ( keys %$hash )
	{
		next if $seen{$key};
		$seen{$key} = 1;
		$total++;
	}
	print "total: $total\n";

    #### GET TIME
	$milliseconds = tv_interval($time1)*1000;
	$time1=[gettimeofday()];
	printf "Elapsed time (milliseconds): $milliseconds\n";

    #### MEMORY USAGE BEFORE    
	warn "memory usage AFTER count:\n";
    warn qx{ ps -o rss,vsz $$ }, "\n";

	#### DELETE DB OBJECT REFERENCE
	undef $db;

	#### UNTIE hash FROM DB FILE
	untie(%$hash);

	return $total;
}



=head2

	SUBROUTINE		lookupRead

	PURPOSE

		RETURN ALL HIT LOCATIONS FOR THE INPUT READ NAME

=cut

sub lookupRead
{
	my $self		=	shift;	
	my $readname	=	shift;
	my $mode		=	shift;

	#### SET DEFAULT MODE TO LIST
	$mode = "array" if not defined $mode;

	my $dbfile		=	$self->get_dbfile();


	#### CHECK INPUTS
	print "IndexRead::lookupRead   readname not defined. Exiting\n"
		and exit if not defined $readname;
	print "IndexRead::lookupRead   dbfile not defined. Exiting\n"
		and exit if not defined $dbfile;

	# ENABLE DUPLICATE RECORDS
	$DB_BTREE->{'flags'} = R_DUP;

	my ($key, $value, $db, $status);

	my $hash;
	$db = tie %$hash, "DB_File", $dbfile, O_RDWR|O_CREAT, 0666, $DB_BTREE 
	or die "Cannot open $dbfile: $!\n";

#	### PRINT KEY/VALUE PAIRS BY ITERATING THROUGH BTREE USING seq
#	$key = $value = 0 ;
#    for ($status = $db->seq($key, $value, R_FIRST) ;
#         $status == 0 ;
#         $status = $db->seq($key, $value, R_NEXT) )
#	{
#	}

	#### get_dup()
	#### 
	#### 	NB: CAREFUL WITH THE SYNTAX
	####
	####	$count = $db->get_dup($key) ;		<== OK
	####	@list  = $db->get_dup($key) ;		<== OK
	####	%list  = $db->get_dup($key, 1) ;  	<== OK
	####	%list  = $db->get_dup($key) ;  		<== BAD - RETURNS A TRUNCATED ARRAY!
	####	@list  = $db->get_dup($key, 1) ;	<== BAD - RETURNS HASH AS ARRAY!

	my $output;	
	if ( $mode eq "count" )
	{
		$output = $db->get_dup($readname, 1);
	}
	elsif ( $mode eq "array" )
	{
		@$output = $db->get_dup($readname);
		print "\n";print Dumper ;
	}
	else
	{
		%$output = $db->get_dup($readname, 1);
	}

	#### DELETE DB OBJECT REFERENCE
	undef $db;

	#### UNTIE hash FROM DB FILE
	untie(%$hash);

	return $output;
}




sub compare
{
	print "compare: \@_:\n";
	print join "\n", @_;
	print "\n";
	exit;

	my ($key1, $key2) = @_ ;
	my ($key1_number1, $key1_number2, $key1_number3) = $key1 =~ /^\D+(\d+)\.(\d+).+\/(\d)$/;
	my ($key2_number1, $key2_number2, $key2_number3) = $key2 =~ /^\D+(\d+)\.(\d+).+\/(\d)/;
	return -1 if $key1_number1 < $key2_number1;
	return 1 if $key1_number1 > $key2_number1;
	return -1 if $key1_number2 < $key2_number2;
	return 1 if $key1_number2 > $key2_number2;

	if ( defined $key1_number3 and defined $key2_number3 )
	{
		return -1 if $key1_number3 < $key2_number3;
		return 1 if $key1_number3 > $key2_number3;
	}

	return 0;
}



=head2

	SUBROUTINE		buildIndex

	PURPOSE

		INDEX MULTIPLE SAM-FORMAT FILES IN THE SAME DATABASE FILE

=cut


sub buildIndex
{
	my $self	=	shift;

	my $inputfiles	=	$self->get_inputfiles();
	my $dbfile		=	$self->get_dbfile();

#
#my $pair = $self->pairOrder(99);
#exit;

#### LATER: USE THIS AS ARGUMENT
### $DB_HASH->{'cachesize'} = 10000 ;


	#### GET MAX CACHE AND dot OR SET TO DEFAULT
	my $dot			=	$self->get_dot() || 100000;
	my $max_cache	=	$self->get_max_cache() || 500;


	#### CHECK INPUTS
	print "IndexRead::buildIndex   inputfiles not defined. Exiting\n" and exit if not defined $inputfiles;
	print "IndexRead::buildIndex   dbfile not defined. Exiting\n" and exit if not defined $dbfile;

	#### HOLD RECENT INDEX DATA IN-MEMORY
	$self->{_wordCache} = {};

	#### ENABLE DUPLICATE RECORDS
    $DB_BTREE->{'flags'} = R_DUP ;

	no warnings;
	sub Compare
    {
		my ($key1, $key2) = @_ ;
		my ($key1_number1, $key1_number2, $key1_number3) = $key1 =~ /^\D+(\d+)\.(\d+).+\/(\d)$/;
		my ($key2_number1, $key2_number2, $key2_number3) = $key2 =~ /^\D+(\d+)\.(\d+).+\/(\d)/;
		return -1 if $key1_number1 < $key2_number1;
		return 1 if $key1_number1 > $key2_number1;
		return -1 if $key1_number2 < $key2_number2;
		return 1 if $key1_number2 > $key2_number2;

		if ( defined $key1_number3 and defined $key2_number3 )
		{
			return -1 if $key1_number3 < $key2_number3;
			return 1 if $key1_number3 > $key2_number3;
		}

        return 0;
    }
	use warnings;

    #### COMPARISON SUBROUTINE
    $DB_BTREE->{'compare'} = \&compare;


	#### BACK UP PREVIOUS DB FILE
	rename $dbfile, "$dbfile.old";

	#### TIE THE HASH TO A DB FILE
	my $hash;
	my $db = tie %$hash, "DB_File", $dbfile, O_RDWR|O_CREAT, 0666, $DB_BTREE 
		or die "Cannot open $dbfile: $!\n";
	$self->{_hash} = $hash;
	$self->{_db} = $db;

	#### SET COUNTER
	$self->{_counter} = 0;


	#### LOAD INPUT FILES INTO DB FILE USING CACHE
	foreach my $inputfile ( @$inputfiles )
	{
		my $counter  = 0;
		my $total_reads  = 0;
		if ( $inputfile =~ /\.sam\.gz$/ )
		{
			open(FILE, "zcat $inputfile|" ) or die "Can't zcat inputfile: $inputfile\n";
		}
		elsif ( $inputfile =~ /\.bam$/ )
		{
			my $samtools = $self->get_samtools();
			print "IndexRead::buildIndex    samtools not defined. Exiting\n" and exit if not defined $samtools;
			open(FILE, "$samtools/samtools view $inputfile |") or die "Can't samtools view inputfile: $inputfile\n";
		}
		else
		{
			open(FILE, $inputfile) or die "Can't open inputfile: $inputfile\n";
			$/ = "\n";
		}

		while ( <FILE> )
		{
			next if $_ =~ /^@/;
			$counter++;
			$total_reads++;
			print "$total_reads\n" if $total_reads % $dot == 0;

			my @elements = split "\t", $_;
			my $id 			=	$elements[0];
			my $flag		=	$elements[1];
			my $chromosome	=	$elements[2];
			my $position	=	$elements[3];

			my $pair_order = $self->pairOrder($flag);

			$self->{_wordCache}{"$id/$pair_order||$chromosome,$position,$flag"} = 1;

			#### If we've processed 500 entries into %self->{_wordCache}, sync to disk
			if ( $counter >= 500)
			{	
				$self->FlushWordCache();
				$counter = 0;
			}
		}
		close(FILE) or die "Can't close inputfile: $inputfile\n";

    }


	#### FLUSH REMAINING IN-MEMORY INDEX TO DISK
	$self->FlushWordCache();


#	### PRINT KEY/VALUE PAIRS BY ITERATING THROUGH BTREE USING seq
#	my $existing;
#	my $total = 0;
#	my ($key, $value, $status);
#	$key = $value = 0 ;
#    for ($status = $db->seq($key, $value, R_FIRST) ;
#         $status == 0 ;
#         $status = $db->seq($key, $value, R_NEXT) )
#	{
#		my @entries = $db->get_dup($key);
#		next if exists $existing->{$key};
#		$existing->{$key} = 1;
#		$total++;
#		
#		#last if $total >= 10;
#		#print "IndexRead::buildIndex   $key -> $value\n"
#	}


	#### SYNC DB TO DISC -- NOT NECESSARY?
	#$db->sync();

	#### DELETE REFERENCE TO DATABASE OBJECT TO AVOID THIS ERROR:
	#### 'untie attempted while 1 inner references still exist'
	undef $db;

	#my @keys = keys %$hash;


	#### UNTIE HASH FROM DATABASE
	#### no warnings TO AVOID THIS ERROR:
	#### 'untie attempted while 1 inner references still exist'
	no warnings;
	untie(%$hash);
	use warnings;


	return $self->{_counter};
}


=head2

	SUBROUTINE		FlushWordCache

	PURPOSE

		FLUSH TEMPORARY IN-MEMORY INDEX CACHE TO HASH 

=cut

sub FlushWordCache
{
	my $self	=	shift;

    my( $id, $entry, $count, $wordcount);	


	sub sort_ids  {
		my ($aa) = $a =~ /^(.+?)\|\|/;
		my ($bb) = $b =~ /^(.+?)\|\|/;
		$a cmp $b;
	}

	#### LOAD IN SORTED ORDER FOR FASTER LOADING TO DB FILE
	my $sorted_ids = {};
	foreach my $key ( sort sort_ids keys %{$self->{_wordCache}} )
	{
		my ($id, $entry) = split "\\|\\|", $key;
		$self->{_hash}->{$id} = $entry;
		$self->{_counter}++;

    }

#	#### iterate through the btree using seq
#    #### and print each key/value pair.
#    my $key = 0;
#	my $value = 0 ;
#	my $status;
#	my $db = $self->{_db};
#    for ($status = $db->seq($key, $value, R_FIRST) ;
#         $status == 0 ;
#         $status = $db->seq($key, $value, R_NEXT) )
#    {
#	}


#	foreach ( keys %{$self->{_hash}} )
#    {
#	}


	#### EMPTY CACHE
    $self->{_wordCache} = {}; 
}

=head2

	SUBROUTINE		pairOrder

	PURPOSE

		DEDUCE PAIR ORDER NUMBER FROM HEXADECIMAL FLAG

	INPUT

		HEXADECIMAL FLAG (2nd COLUMN) OF SAM FILE ENTRY

	OUTPUT

		1 OR 2 READ PAIR ORDER NUMBER IF SPECIFIED IN FLAG

		0 IF NOT ENOUGH INFORMATION

	NOTES

		0x0001 the read is paired in sequencing, no matter whether it is mapped in a pair
		0x0002 the read is mapped in a proper pair (depends on the protocol, normally inferred during alignment) 1
		0x0004 the query sequence itself is unmapped
		0x0008 the mate is unmapped 1
		0x0010 strand of the query (0 for forward; 1 for reverse strand)
		0x0020 strand of the mate 1
		0x0040 the read is the first read in a pair 1,2
		0x0080 the read is the second read in a pair 1,2
		0x0100 the alignment is not primary (a read having split hits may have multiple primary alignment records)
		0x0200 the read fails platform/vendor quality checks
		0x0400 the read is either a PCR duplicate or an optical duplicate

		Notes:
		1. Flag 0x02, 0x08, 0x20, 0x40 and 0x80 are only meaningful when flag 0x01 is present.
		2. If in a read pair the information on which read is the first in the pair is lost in the upstream analysis, flag 0x01 should
		be present and 0x40 and 0x80 are both zero.

=cut

sub pairOrder
{
	my $self	=	shift;
	my $flag	=	shift;

	#### NB: FLAG IS PARSED AUTOMATICALLY FROM STRING INTO DECIMAL.
	#### SET PAIR ORDER TO 1 OR 2, RESPECTIVELY, IF 0x40 OR 0x80 ARE SET
	my $pair_order = 0;
	if ( ($flag % 100) > 40 )
	{
		$pair_order = 1;

		if ( ($flag % 100) > 80 )
		{
			$pair_order = 2;
		}
	}

	return $pair_order;
}

=head2

	SUBROUTINE		compactFile

	PURPOSE

		COMPACT INDEX BY COPYING INTO A NEW DB FILE

=cut

sub compactFile
{
	my $self	=	shift;
	my $dbfile	=	shift;
	my $newfile	=	shift;

	##### Compact index by copying into a fresh file #####

	#rename "index.db","indexold.db"; # Move index.db to indexold.db
	#
	#tie(%self->{_db}->old,'DB_File',"indexold.db", O_RDWR, 0644, $DB_File::DB_BTREE)
	#	|| die "Can't open/create DB file: $!";
	#
	#tie(%self->{_db}->,'DB_File',"index.db",
	#	O_RDWR | O_CREAT, 0644, $DB_File::DB_BTREE)
	#	|| die "Can't open/create DB file: $!";
	#
	#my($topcount,$count,$key,$val);
	#while(($key,$val) = each %self->{_db}->old) {
	#	$topcount++;
	#	if(++$count >= 1000) {
	#	&Progress("Copying index.db: $topcount");
	#	$count = 0;
	#	}
	#	$self->{_db}->{$key} = $val;
	#}
	#untie(%self->{_db}->);
	#untie(%self->{_db}->old);	
}





=head2

	SUBROUTINE		tieFile

	PURPOSE

		SET THE HANDLE FOR THE TIED DATABASE FILE AND HASH

=cut

sub tieFile
{
	my $self	=	shift;
	my $hash	=	shift;
	my $dbfile	=	shift;

	# Attach %$hash to an on-disk DB file
	rename $dbfile, "$dbfile.bkp" if -f $dbfile;
	unlink($dbfile) if -f $dbfile;
	my ($db) = tie(%$hash, 'DB_File', $dbfile,
		O_RDWR | O_CREAT, 0644, $DB_BTREE)
		or die "Can't tie dbfile: $dbfile: $!";

	$self->set_hash($hash);

	return $db;
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

	#### INITIALISE THE OBJECT'S ELEMENTS
	$self->initialise($arguments);

	#### SET REFERENCE
	if ( defined $self->{_referencefile} )
	{
		my ($reference) = $self->{_referencefile} =~ /([^\/]+)$/;
		$reference =~ s/\.[^\.]+$//;

		$self->{_reference} = $reference;
	}

	#### SET DEBUG IF verbose

    return $self;
}






=head2

	SUBROUTINE		initialise

	PURPOSE

		INITIALISE THE self OBJECT:

			1. LOAD THE DATABASE, USER AND PASSWORD FROM THE ARGUMENTS

			2. FILL OUT %VARIABLES% IN XML AND LOAD XML

			3. LOAD THE ARGUMENTS

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

	#### SET TIME
	$self->{_starttime} = time();

	#### REQUIRED: CALL THE PARENT FOR ADDITIONAL INITIALISATION
	$self->SUPER::initialise();		
}



=head2

	SUBROUTINE		AUTOLOAD

	PURPOSE

		AUTOMATICALLY DO 'set_' OR 'get_' FUNCTIONS IF THE

		SUBROUTINES ARE NOT DEFINED.

=cut

sub AUTOLOAD {
    my ($self, $newvalue) = @_;

    my ($operation, $attribute) = ($AUTOLOAD =~ /(get|set)(_\w+)$/);

    # Is this a legal method name?
    unless($operation && $attribute) {
        croak "Method name $AUTOLOAD is not in the recognized form (get|set)_attribute\n";
    }
    unless( exists $self->{$attribute} )
	{
        #croak "No such attribute '$attribute' exists in the class ", ref($self);
		#return undef;
    }

    # Turn off strict references to enable "magic" AUTOLOAD speedup
    no strict 'refs';

    # AUTOLOAD accessors
    if($operation eq 'get') {

		return $self->{$attribute};

        # define subroutine
        *{$AUTOLOAD} = sub { shift->{$attribute} };

    # AUTOLOAD mutators
    }elsif($operation eq 'set') {


        ## define subroutine
        #*{$AUTOLOAD} = sub { shift->{$attribute} = shift; };

        # set the new attribute value
        $self->{$attribute} = $newvalue;
    }

    # Turn strict references back on
    use strict 'refs';

    # return the attribute value
    return $self->{$attribute};
}


# When an object is no longer being used, this will be automatically called
# and will adjust the count of existing objects
sub DESTROY {
    my($self) = @_;

	#if ( defined $self->{_databasehandle} )
	#{
	#	my $dbh =  $self->{_databasehandle};
	#	$dbh->disconnect();
	#}

#    my($self) = @_;
}





1;

