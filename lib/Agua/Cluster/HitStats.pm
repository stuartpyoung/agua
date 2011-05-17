package Agua::Cluster::HitStats;
use Moose::Role;
use Moose::Util::TypeConstraints;

use Data::Dumper;

#has 'submit'	=> ( isa => 'Int', is => 'rw', default => 0 );
#has 'cluster'	=> ( isa => 'Str', is => 'rw', default => '' );
#has 'queue'		=> ( isa => 'Str', is => 'rw', default => '' );
#has 'walltime'	=> ( isa => 'Str', is => 'rw', default => '' );
#has 'cpus'		=> ( isa => 'Str', is => 'rw', default => '' );
#has 'qstat'		=> ( isa => 'Str', is => 'rw', default => '' );
#has 'qsub'		=> ( isa => 'Str', is => 'rw', default => '' );
#has 'maxjobs'	=> ( isa => 'Str', is => 'rw', default => '' );
#has 'sleep'		=> ( isa => 'Str', is => 'rw', default => '' );
#has 'cleanup'	=> ( isa => 'Str', is => 'rw', default => '' );
#has 'dot'		=> ( isa => 'Str', is => 'rw', default => '' );
#has 'verbose'	=> ( isa => 'Str', is => 'rw', default => '' );
#


################################################################################
##########################       HIT STATS METHODS        ######################
################################################################################
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

sub readHits {
	my $self			=	shift;
	my $outputdir		=	shift;
	my $references		=	shift;
	my $splitfiles		=	shift;
	my $label			=	shift;
	my $samfile			=	shift;

	print "Agua::Cluster::HitStats::readHits    samfile not defined\n" and exit if not defined $samfile;

	my $hitfile = "out.hitid";
	my $countfile = "out.hitct";
	my $count_filepath = "$outputdir/$countfile";
	my $sorted_filepath = "$outputdir/$countfile.sorted";
	my $bin_filepath = "$outputdir/$countfile.bins";

	#### 1. EXTRACT HIT READS IDS FROM CHROMOSOME SAM FILE
	####	AND PRINT TO CHROMOSOME .hitid FILE
	####
	$self->hitIdFiles($outputdir, $references, $samfile, $hitfile);
	print "Agua::Cluster::HitStats::readHits    Completed hitIdFiles\n";

	#### 2. SORT .hitid FILES
	####
	print "Agua::Cluster::HitStats::readHits    Doing sortHitFiles()\n";
	$self->sortHitFiles($outputdir, $references, $hitfile, $hitfile);
	print "Agua::Cluster::HitStats::readHits    Completed sortHitFiles()\n";

	######## 3. MERGE .hitid FILES FOR ALL CHROMOSOMES INTO
	########		SINGLE .hitid FILE FOR WHOLE GENOME
	########
	#####print "Agua::Cluster::HitStats::readHits    Doing mergeHits()\n";
	#####$self->mergeHits($outputdir, $references, $hitfile, $hitfile, $label);
	#####print "Agua::Cluster::HitStats::readHits    Completed mergeHits()\n";

	### 3. MERGE .hitid FILES FOR ALL CHROMOSOMES INTO
	###		SINGLE .hitid FILE FOR WHOLE GENOME
	###
	print "Agua::Cluster::HitStats::readHits    Doing pyramidMergeHits()\n";
	$self->pyramidMergeHits($outputdir, $references, $hitfile, "$hitfile", $label);
	print "Agua::Cluster::HitStats::readHits    Completed pyramidMergeHits()\n";


	######### 3. SORT .hitid FILE
	#########
	#####my $sort_command = "sort $outputdir/$hitfile -o $outputdir/$hitfile";
	#####print "Agua::Cluster::HitStats::readHits    sort command: $sort_command\n";
	#####my $sort_job = $self->setJob( [$sort_command], "sortHitId", $outputdir );
	#####$self->runJobs( [ $sort_job ], "sortHitId" );	
	#####print "Agua::Cluster::HitStats::readHits    Completed sortHitId\n";


	#### 4. COUNT DUPLICATE LINES WITH uniq -c AND OUTPUT
	####	INTO .hitct FILE
	####	
	my $uniq_command = "uniq -c $outputdir/$hitfile > $outputdir/$countfile";
	print "Agua::Cluster::HitStats::readHits    uniq command: $uniq_command\n";
	print "Agua::Cluster::HitStats::readHits    uniq command: $uniq_command\n";
	my $uniq_job = $self->setJob( [ $uniq_command ], "uniqHitId", $outputdir );
	$self->runJobs( [ $uniq_job ], "uniqHitId" );	
	print "Agua::Cluster::HitStats::readHits    Completed uniqHitId\n";



	##### 5. SORT HIT COUNTS FILE BY HIT COUNTS
	#####
	my $count_command = qq{sort -t "\t" -b -k 1,1n $count_filepath -o $sorted_filepath};
	print "Agua::Cluster::HitStats::readHits    count command: $count_command\n";
	my $count_job = $self->setJob( [ $count_command ], "countHitId", $outputdir );
	$self->runJobs( [ $count_job ], "countHitId" );	
	print "Agua::Cluster::HitStats::readHits    Completed countHitId\n";


	#### 6. BIN HIT COUNTS
	####
	my $bin_command = qq{sed -e 's/[ ]*\\([0-9]*\\) [A-Z0-9\\.\\:\\#\\/]*/\\1/' < $sorted_filepath | sed -e 's/[ ]*//g' | uniq -c > $bin_filepath};
	print "Agua::Cluster::HitStats::readHits    bin command: $bin_command\n";

	my $bin_job = $self->setJob( [ $bin_command ], "binHitCounts", $outputdir );
	$self->runJobs( [ $bin_job ], "binHitCounts" );	
	print "Agua::Cluster::HitStats::readHits    Completed binHitCounts\n";


}



=head2

	SUBROUTINE		sortHitFiles

	PURPOSE

		SORT out.hitid FILE FOR EACH REFERENCE CHROMOSOME

=cut

sub sortHitFiles {
	my $self			=	shift;
	my $outputdir		=	shift;
	my $references		=	shift;
	my $infile			=	shift;
	my $outfile			=	shift;

	print "Agua::Cluster::HitStats::sortHitFiles    Agua::Cluster::sortHitFiles(outputdir, references, splitfiles, infile, outfile, toplines)\n";
	print "Agua::Cluster::HitStats::sortHitFiles    outputdir: $outputdir\n";
	print "Agua::Cluster::HitStats::sortHitFiles    No. references: ", scalar(@$references), "\n";
	print "Agua::Cluster::HitStats::sortHitFiles    references: ";
	print join "\n", @$references;
	print "\n";
	print "Agua::Cluster::HitStats::sortHitFiles    infile: $infile\n";
	print "Agua::Cluster::HitStats::sortHitFiles    outfile: $outfile\n";

	#### SORT SINGLE SAM FILE
	my $jobs = [];
	foreach my $reference ( @$references )
	{
		my $inputfile = "$outputdir/$reference/$infile";
		my $outputfile = "$outputdir/$reference/$outfile";

		#### SET EXECUTABLE SCRIPT TO CONVERT SAM TO TSV
		my $command = qq{sort $inputfile -o $outputfile; };

		print "Agua::Cluster::HitStats::sortHitFiles    command: $command\n";
		my $job = $self->setJob( [$command], "sortHitFiles-$reference", "$outputdir/$reference" );
		push @$jobs, $job;
	}
	print "Doing runJobs, no. jobs: " , scalar(@$jobs), "\n";

	#### RUN COMMANDS IN PARALLEL
	$self->runJobs( $jobs, "sortHitFiles" );	
	print "Agua::Cluster::HitStats::sortHitFiles    Completed.\n";
}

=head2

	SUBROUTINE		hitIdFiles

	PURPOSE

		1. EXTRACT READS IDS ONLY FROM ALIGNMENT FILE
			AND PRINT TO .hitid FILE

=cut

sub hitIdFiles {
	my $self			=	shift;
	my $outputdir		=	shift;
	my $references		=	shift;
	my $samfile			=	shift;
	my $hitfile_name	=	shift;
	my $samfile_name	=	shift;



	$samfile_name = "out.sam" if not defined $samfile_name;

	my $jobs = [];
	foreach my $reference ( @$references )
	{
		print "Creating hit files in parallel for reference: $reference\n";

		#### SET EXECUTABLE SCRIPT TO CONVERT SAM TO TSV
		my $executable = "$Bin/samHits.pl";
		my $inputfile = "$outputdir/$reference/$samfile";
		my $outputfile = "$outputdir/$reference/$hitfile_name";		
		my $command = "/usr/bin/perl $executable --inputfile $inputfile | cut -f 1 > $outputfile";
		my $job = $self->setJob( [$command], "hitIdFiles-$reference", "$outputdir/$reference" );
		push @$jobs, $job;
	}

	#### RUN COMMANDS IN PARALLEL
	$self->runJobs( $jobs, "hitIdFiles" );	
}

=head2

	SUBROUTINE		mergeHits

=cut

sub mergeHits {
	my $self			=	shift;
	my $outputdir		=	shift;
	my $references		=	shift;
	my $hitfile_name	=	shift;
	my $outfile_name	=	shift;
	my $label			=	shift;



	#### DO MERGE OF SORTED FILES
	my $outputfile = "$outputdir/$outfile_name";
	my $commands = [];
	push @$commands, "rm -fr $outputfile";

	my $command = "sort -m ";
	foreach my $reference ( @$references )
	{
		$command .= " $outputdir/$reference/$hitfile_name ";
	}
	$command .= " -o $outputfile";
	push @$commands, $command;

	my $job = $self->setJob($commands, $label, $outputdir);

	#### RUN MERGE JOBS
	$self->runJobs([ $job ], "mergeHits-$label" );
}

=head2

	SUBROUTINE		pyramidMergeHits

=cut

sub pyramidMergeHits {
	my $self			=	shift;
	my $outputdir		=	shift;
	my $references		=	shift;
	my $hitfile_name	=	shift;
	my $outfile_name	=	shift;
	my $label			=	shift;


	#### DO MERGE OF SORTED FILES
	my $outputfile = "$outputdir/$outfile_name";

	my $mergesub = sub {
		my $firstfile	=	shift;
		my $secondfile	=	shift;

		print "Agua::Cluster::HitStats::pyramidMergeHits::mergesub    firstfile: $firstfile\n";
		print "Agua::Cluster::HitStats::pyramidMergeHits::mergesub    secondfile: $secondfile\n";

		#### KEEP THIS PART OF THE LOGIC HERE JUST IN CASE THE
		#### MERGE FUNCTION EXPECTS A PARTICULAR FILE SUFFIX
		my $outfile = $firstfile;
		if ( $outfile =~ /\.merge\.(\d+)$/ )
		{
			my $number = $1;
			$number++;
			$outfile =~ s/\d+$/$number/;	 
		}
		else
		{
			$outfile .= ".merge.1"
		}

		#### MERGE SORTED FILES
		my $command = "sort -m $firstfile $secondfile -o $outfile\n";

		return ($command, $outfile);
	};

	#### REMOVE OUTPUT FILE IF EXISTS
	print "Doing rm -fr $outputfile\n" if -f $outputfile;
	print `rm -fr $outputfile` if -f $outputfile;
	print "Agua::Cluster::HitStats::pyramidMergeHits    Can't delete outputfile: $outputfile\n" if -f $outputfile;

	#### SET CURRENT TIME (START TIMER)
	my $current_time =  time();

	#### GENERATE INPUT FILE LIST
	my $inputfiles = [];
	foreach my $reference ( @$references )
	{
		push @$inputfiles, "$outputdir/$reference/$hitfile_name";
	}

	my $alljobs = [];
	my $running = 1;
	my $counter = 0;
	while ( $running )
	{
		print "Agua::Cluster::HitStats::pyramidMergeHits    running\n";
		$counter++;

		my $outputfile = "$outputdir/$outfile_name";
		my $label = "_pyramidMerge-$counter";
		print "Agua::Cluster::HitStats::pyramidMergeHits    label: $label\n";

		my $jobs = [];
		($jobs, $inputfiles) = $self->mergeJobs(
			{
				inputfiles 	=> 	$inputfiles,
				outputfile 	=> 	$outputfile,
				mergesub	=>	$mergesub,
				label		=>	$label
			}
		);
		print "Agua::Cluster::HitStats::pyramidMergeHits    No. jobs: ", scalar(@$jobs), "\n";
		#exit;

		#### ADD TO ALL JOBS
		(@$alljobs) = (@$alljobs, @$jobs) if defined $jobs and @$jobs;

		if ( not @$alljobs )
		{
			$running = 0;
			last;
		}

		#### RUN JOBS
		$self->runJobs($alljobs, $label);

		#### EMPTY ALLJOBS
		$alljobs = [];
	}

	#### SET EMPTY USAGE STATISTICS
	my $usage = [];

	#### GET DURATION (STOP TIMER)
	my $duration = Timer::runtime( $current_time, time() );

	#### PUSH ONTO GLOBAL USAGE STATS
	$self->addUsageStatistic($label, $duration, $usage);

	#### PRINT DURATION
	print "Agua::Cluster::HitStats::_pyramidMerge    Completed $label at ", Timer::current_datetime(), ", duration: $duration\n";
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

sub splitReads {	
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
		print "Agua::Cluster::HitStats::splitReads    splitreadsfile printed:\n\n$splitreadsfile\n\n";
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

sub splitHits {	
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
		print "Agua::Cluster::HitStats::readHits    splithitsfile printed:\n\n$splithitsfile\n\n";
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

sub countLines {
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




1;
