package Agua::Cluster::Merge;
use Moose::Role;
use Moose::Util::TypeConstraints;

use Data::Dumper;


##########################       MERGE METHODS        ##########################

#### METHODS sortSubdirSam AND subdirSamHIts REQUIRE Jobs::getIndex


=head2

	SUBROUTINE		pyramidMergeBam

	PURPOSE

		1. MERGE ALL FILES PER REFERENCE IN A PYRAMID OF DECREASING

			(BY HALF) BATCHES OF MERGES UNTIL THE FINAL OUTPUT FILE

			IS PRODUCED.

		2. EACH STAGE IS COMPLETED WHEN ALL OF THE MERGES FOR ALL

			REFERENCES ARE COMPLETE.

		3. EACH PAIRWISE MERGE OPERATION IS CARRIED OUT ON A SEPARATE

			EXECUTION HOST

		4. THIS METHOD ASSUMES:

			1. THERE IS NO ORDER TO THE FILES

			2. ALL FILES MUST BE MERGED INTO A SINGLE FILE

			3. THE PROVIDED SUBROUTINE MERGES THE FILES

		INPUTS

=cut

sub pyramidMergeBam {
	my $self			=	shift;
	my $outputdir		=	shift;
	my $references		=	shift;
	my $splitfiles		=	shift;
	my $infile			=	shift;
	my $outfile			=	shift;

print "Agua::Cluster::Merge::pyramidMergeBam    infile: $infile\n";print Dumper ;


	#### SET DEFAULT infile AND outfile
	$infile = "out.bam" if not defined $infile;
	$outfile = "out.bam" if not defined $outfile;

	#### GET SAMTOOLS
	my $samtools = $self->samtools();

	#### LOAD UP WITH INITIAL BAM FILES
	my $reference_bamfiles;
	foreach my $reference ( @$references )
	{
		#### GET SAM FILES
		my $bamfiles = $self->subfiles($outputdir, $reference, $splitfiles, $infile);
		$reference_bamfiles->{$reference} = $bamfiles;
	}

	my $mergesub = sub {
		my $firstfile	=	shift;
		my $secondfile	=	shift;

		print "Agua::Cluster::Merge::pyramidMergeBam::mergesub    outfile: $outfile\n";
		print "Agua::Cluster::Merge::pyramidMergeBam::mergesub    firstfile: $firstfile\n";
		print "Agua::Cluster::Merge::pyramidMergeBam::mergesub    secondfile: $secondfile\n";

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

		#### MERGE BAM FILES
		my $command = "$samtools/samtools merge $outfile $firstfile $secondfile;\n";

		return ($command, $outfile);
	};


	#### RUN PYRAMID MERGE ON ALL REFERENCES IN PARALLEL
	my $label = "pyramidMergeBam";
	$self->_pyramidMerge($outputdir, $references, $splitfiles, $infile, $outfile, $reference_bamfiles, $label, $mergesub);

	print "Agua::Cluster::Merge::pyramidMergeBam    Completed.\n";
}


=head2

	SUBROUTINE		_pyramidMerge

	PURPOSE

		MERGE BY REFERENCE FILE IN PARALLEL USING THE SUPPLIED

		MERGE FUNCTION

		1. FOR EACH ROUND OF pyramidMerge:

				- COLLECT JOBS FOR EACH REFERENCE FROM mergeJobs

				- mergeJobs CALLS pairwiseMergeJobs

				- RUN ALL JOBS AND RETURN REMAINING NUMBER OF

					FILES FOR EACH REFERENCE IN HASH

		2. IF ONLY ONE FILE REMAINING IN INPUT FILES, mergeJobs

			WILL RETURN NULL

		3. IF THE reference-jobs HASH VALUES FOR ALL REFERENCES IS NULL,

			STOP RUNNING pyramidMerge AND RETURN

		4. THIS METHOD ASSUMES:

			1. THERE IS NO ORDER TO THE FILES

			2. ALL FILES MUST BE MERGED INTO A SINGLE FILE

			3. THE PROVIDED SUBROUTINE MERGES THE FILES

=cut

sub _pyramidMerge {
	my $self			=	shift;
	my $outputdir		=	shift;
	my $references		=	shift;
	my $splitfiles		=	shift;
	my $infile			=	shift;
	my $outfile			=	shift;
	my $reference_infiles	=	shift;
	my $label			=	shift;
	my $mergesub		=	shift;

	print "Agua::Cluster::Merge::_pyramidMerge    Agua::Cluster::_pyramidMerge(outputdir, references, splitfiles, infile, outfile, reference_infiles, label, mergesub)\n";
	print "Agua::Cluster::Merge::_pyramidMerge    outputdir: $outputdir\n";
	print "Agua::Cluster::Merge::_pyramidMerge    No. references: ", scalar(@$references), "\n";
	print "Agua::Cluster::Merge::_pyramidMerge    references: ";
	print join "\n", @$references;
	print "\n";

	#### SET CURRENT TIME (START TIMER)
	my $current_time =  time();

	my $alljobs = [];
	my $running = 1;
	my $counter = 0;
	while ( $running )
	{
		$counter++;

		for my $reference ( @$references )
		{
			print "reference_infiles->{$reference} not defined\n" and next if not defined $reference_infiles->{$reference};
			next if not @{$reference_infiles->{$reference}};

			my $outputfile = "$outputdir/$reference/$outfile";
			my $label = "_pyramidMerge-$reference-$counter";
			print "Agua::Cluster::Merge::_pyramidMerge    label: $label\n";

			my $jobs = [];
			my $inputfiles = $reference_infiles->{$reference};
			($jobs, $reference_infiles->{$reference}) = $self->mergeJobs(
				{
					inputfiles 	=> 	$inputfiles,
					outputfile 	=> 	$outputfile,
					mergesub	=>	$mergesub,
					label		=>	$label
				}
			);


		print "Agua::Cluster::Merge::_pyramidMerge    No. jobs: ", scalar(@$jobs), "\n";
		#exit;


			#### ADD TO ALL JOBS
			(@$alljobs) = (@$alljobs, @$jobs) if defined $jobs and @$jobs;
		}

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
	print "Agua::Cluster::Merge::_pyramidMerge    Completed $label at ", Timer::current_datetime(), ", duration: $duration\n";

	print "Agua::Cluster::Merge::_pyramidMerge completed\n";	
}




=head2

	SUBROUTINE		mergeJobs

	PURPOSE

		RETURN JOBS TO MERGE FILES ON A CLUSTER.

		THIS METHOD ASSUMES:

			1. THERE IS NO ORDER TO THE FILES

			2. ALL FILES MUST BE MERGED INTO A SINGLE FILE

			3. THE PROVIDED SUBROUTINE MERGES THE FILES

=cut

sub mergeJobs {
	my $self			=	shift;
	my $args			=	shift;

	my $inputfiles		=	$args->{inputfiles};
	my $outputfile		=	$args->{outputfile};
	my $mergesub		=	$args->{mergesub};
	my $label			=	$args->{label};

	print "Agua::Cluster::Merge::mergeJobs    Agua::Cluster::mergeJobs(inputfiles, outputfile, mergesub)\n";
	print "Agua::Cluster::Merge::mergeJobs    No. inputfiles: ", scalar(@$inputfiles), "\n";
	print "Agua::Cluster::Merge::mergeJobs    outputfile: $outputfile\n";


	#### IF ONLY ONE FILE LEFT, MOVE THE LAST MERGED INPUT FILE TO THE OUTPUT FILE
	if ( scalar(@$inputfiles) == 1 )
	{
		print "Agua::Cluster::Merge::mergeJobs    Completed merging. Doing mv to outputfile\n";
		my $move = "mv -f $$inputfiles[0] $outputfile";
		print "Agua::Cluster::Merge::mergeJobs    move: $move\n";
		`$move`;

		#### RETURN NO JOBS AND NO INPUT FILES
		return ([], []);
	}

	#### MOVE THE LAST MERGED INPUT FILE TO THE OUTPUT FILE
	else
	{
		my ($outdir) = $$inputfiles[0] =~ /^(.+)\/[^\/]+$/; 
		print "Agua::Cluster::Merge::mergeJobs    label: $label\n";
		print "Agua::Cluster::Merge::mergeJobs    outdir: $outdir\n";

		#### SET JOBS
		my $jobs;
		($jobs, $inputfiles) = $self->pairwiseMergeJobs($inputfiles, $mergesub, $label, $outdir);
		print "Agua::Cluster::Merge::mergeJobs    No. REMAINING inputfiles: ", scalar(@$inputfiles), "\n";

		#exit;

		#### RETURN JOBS AND INPUT FILES
		return ($jobs, $inputfiles);
	}

}



=head2

	SUBROUTINE		pairwiseMergeJobs

	PURPOSE

		RETURN AN ARRAY OF MERGE COMMANDS FOR PAIRWISE MERGE 

		OF ALL FILES IN INPUT ARRAY OF SUBFILES

=cut

sub pairwiseMergeJobs {
	my $self		=	shift;
	my $inputfiles	=	shift;
	my $mergesub	=	shift;
	my $label		=	shift;
	my $outdir		=	shift;


	#### DO split MERGE OF INPUT FILES
	my $jobs = [];
	my $outputfiles = [];
	my $totalfiles = scalar(@$inputfiles);
	for ( my $i = 0; $i < @$inputfiles; $i+=2 )
	{
		#### SKIP IF ONLY ONE FILE LEFT
		push @$outputfiles, $$inputfiles[$i] and last if $i == scalar(@$inputfiles) - 1;

		#### MERGE NOTE
		my $subfile = $$inputfiles[$i];
		my $next_index = $i + 1;
		my $note = "echo 'Merging subfiles $i and $next_index of $totalfiles'";

		my $commands = [];
		push @$commands, $note;

		#### MERGE COMMAND
		my ($merge_command, $mergedfile) = &$mergesub($$inputfiles[$i], $$inputfiles[$i + 1]);
		push @$commands, $merge_command;

		#### SET LABEL
		my $this_label = $label . "-$i";

		#### SET JOB
		my $job = $self->setJob($commands, $this_label, $outdir);
		push @$jobs, $job;

		#### ADD MERGED FILE TO FILE TO BE MERGED IN NEXT ROUND
		push @$outputfiles, $mergedfile;
	}

	return ($jobs, $outputfiles);
}	





=head2

	SUBROUTINE		cumulativeMergeSam

	PURPOSE

		1. SUCCESSIVELY MERGE ALL OF THE SAM FILES FOR

			EACH SPLIT INPUT FILE FOR ALL PROVIDED

			REFERENCES

		2. CUMULATIVELY ADD TO A TEMP FILE (USED TO SHOW ONGOING)

		3. ONCE COMPLETE, RENAME TEMP FILE TO OUTPUT FILE

=cut

sub cumulativeMergeSam {
	my $self			=	shift;
	my $outputdir		=	shift;
	my $references		=	shift;
	my $splitfiles		=	shift;
	my $infile			=	shift;
	my $outfile			=	shift;

	#### CHECK INPUTS
	print "Agua::Cluster::Merge::cumulativeMergeSam    outputdir not defined. Exiting.\n" and exit if not defined $outputdir;
	print "Agua::Cluster::Merge::cumulativeMergeSam    references not defined. Exiting.\n" and exit if not defined $references;
	print "Agua::Cluster::Merge::cumulativeMergeSam    splitfiles not defined. Exiting.\n" and exit if not defined $splitfiles;
	print "Agua::Cluster::Merge::cumulativeMergeSam    infile not defined. Exiting.\n" and exit if not defined $infile;
	print "Agua::Cluster::Merge::cumulativeMergeSam    outfile not defined. Exiting.\n" and exit if not defined $outfile;

	print "Agua::Cluster::Merge::cumulativeMergeSam    Agua::Cluster::cumulativeMergeSam(outputdir, references, splitfiles, infile, outfile, toplines)\n";


	print "Agua::Cluster::Merge::cumulativeMergeSam    outputdir: $outputdir\n";
	print "Agua::Cluster::Merge::cumulativeMergeSam    No. references: ", scalar(@$references), "\n";

	print "Agua::Cluster::Merge::cumulativeMergeSam    infile: $infile\n";
	print "Agua::Cluster::Merge::cumulativeMergeSam    outfile: $outfile\n";

	#####################################################################################
	######## STRATEGY 1. MERGE SUBDIR out.sam FILES FIRST, THEN SORT SINGLE SAM FILE
	#####################################################################################

	#### MERGE ALL SPLIT *.sam FILES INTO A SINGLE CHROMOSOME *.sam FILE
	$self->mergeSam($outputdir, $references, $splitfiles, $infile, $outfile);

	#### SORT SINGLE CHROMOSOME SAM FILE
	$self->sortSam($outputdir, $references, $outfile, $outfile);


	#####################################################################################
	##### STRATEGY 2. SORT SPLIT out.sam FILES FIRST, THEN MERGE AND SORT SINGLE SAM FILE
	#####################################################################################

	##### SORT ALL SPLIT out.sam FILES
	#$self->sortSubdirSam($outputdir, $references, $splitfiles, $infile, $outfile);

	##### MERGE ALL SPLIT out.sam FILES
	#$self->mergeSam($outputdir, $references, $splitfiles, $outfile, $outfile);

	##### SORT SINGLE SAM FILE
	#$self->sortSam($outputdir, $references, $infile, $outfile);
}

=head2

	SUBROUTINE		mergeSam

	PURPOSE

		1. SUCCESSIVELY MERGE ALL OF THE SAM FILES FOR EACH SPLIT

			INPUT FILE FOR ALL PROVIDED REFERENCES

=cut

sub mergeSam {
	my $self			=	shift;
	my $outputdir		=	shift;
	my $references		=	shift;
	my $splitfiles		=	shift;
	my $infile			=	shift;
	my $outfile			=	shift;

	print "Agua::Cluster::Merge::mergeSam    Agua::Cluster::mergeSam(outputdir, references, splitfiles, infile, outfile, toplines)\n";
	print "Agua::Cluster::Merge::mergeSam    outputdir: $outputdir\n";
	print "Agua::Cluster::Merge::mergeSam    No. references: ", scalar(@$references), "\n";
	print "Agua::Cluster::Merge::mergeSam    infile: $infile\n";
	print "Agua::Cluster::Merge::mergeSam    outfile: $outfile\n";


	#### STRATEGY 1. MERGE SPLIT out.sam FILES FIRST, THEN SORT SINGLE SAM FILE
	my $jobs = [];
	foreach my $reference ( @$references )
	{
		#### GET SAM FILES
		my $inputfiles = $self->subfiles($outputdir, $reference, $splitfiles, $infile);
		#my $no_header = $self->noSamHeader($$inputfiles[0]);

		my $outputfile = "$outputdir/$reference/$outfile";
		print "Agua::Cluster::Merge::mergeSam    outputfile is a directory: $outputfile\n" and exit if -d $outputfile;
		print `rm -fr $outputfile`;
		print "Agua::Cluster::Merge::mergeSam    Can't remove outputfile: $outputfile\n" and exit if -f $outputfile;

		my $commands = [];
		foreach my $inputfile ( @$inputfiles )
		{
			my $command = "cat $inputfile >> $outputfile";
#exit;
			push @$commands, $command;
			push @$commands, "echo 'Merging file: ' $inputfile ";
		}
		my $job = $self->setJob( $commands, "mergeSam-$reference", "$outputdir/$reference" );
		#### OVERRIDE checkfile
		$job->{checkfile} = $outputfile;

		push @$jobs, $job;
	}

	#### RUN COMMANDS IN PARALLEL
	$self->runJobs( $jobs, "mergeSam" );	
	print "Agua::Cluster::Merge::mergeSam    Completed.\n";
}


=head2

	SUBROUTINE		sortSubdirSam

	PURPOSE

		SORT SAM FILES IN EACH SPLIT SUBDIR

=cut

sub sortSubdirSam {
	my $self			=	shift;
	my $outputdir		=	shift;
	my $references		=	shift;
	my $splitfiles		=	shift;
	my $infile			=	shift;
	my $outfile			=	shift;

	print "Agua::Cluster::Merge::sortSubdirSam    Agua::Cluster::sortSubdirSam(outputdir, references, splitfiles, infile, outfile, toplines)\n";
	print "Agua::Cluster::Merge::sortSubdirSam    outputdir: $outputdir\n";
	print "Agua::Cluster::Merge::sortSubdirSam    No. references: ", scalar(@$references), "\n";
	print "Agua::Cluster::Merge::sortSubdirSam    references: ";
	print join "\n", @$references;
	print "\n";
	print "Agua::Cluster::Merge::sortSubdirSam    infile: $infile\n";
	print "Agua::Cluster::Merge::sortSubdirSam    outfile: $outfile\n";

	#### GET CLUSTER
	my $cluster = $self->cluster();

	#### SET INDEX PATTERN FOR BATCH JOB
	my $index = $self->getIndex();

	#### DO BATCH JOBS
	my $jobs = [];
	foreach my $reference ( @$references )
	{
		my $inputfile = "$outputdir/$reference/$index/$infile";
		my $outputfile = "$outputdir/$reference/$index/$outfile";
		my $outdir = "$outputdir/$reference";
		print "Agua::Cluster::Merge::sortSubdirSam    inputfile: $inputfile\n";
		print "Agua::Cluster::Merge::sortSubdirSam    outputfile: $outputfile\n";
		print "Agua::Cluster::Merge::sortSubdirSam    outdir: $outdir\n";

		#### REMOVE OUTPUT FILE IF EXISTS
		print "Agua::Cluster::Merge::sortSubdirSam    outputfile is a directory: $outputfile\n" and exit if -d $outputfile;
		print `rm -fr $outputfile`;
		print "Agua::Cluster::Merge::sortSubdirSam    Can't remove outputfile: $outputfile\n" and exit if -f $outputfile;

		#### SET EXECUTABLE SCRIPT TO CONVERT SAM TO TSV
		my $command = qq{sort -t "\t" -k 3,3 -k 4,4n $inputfile -o $outputfile};
		my $job = $self->setBatchJob( [$command], "sortSubdirSam-$reference", $outputdir );

		push @$jobs, $job;
	}

	#### RUN COMMANDS IN PARALLEL
	$self->runJobs( $jobs, "sortSubdirSam" );	
	print "Agua::Cluster::Merge::sortSubdirSam    Completed.\n";
}

=head2

	SUBROUTINE		samHits

	PURPOSE

		1. FILTER HITS AND MISSES INTO SEPARATE SAM FILES

		2. DO ALL SAM FILES AT THE CHROMOSOME/REFERENCE LEVEL

=cut


sub samHits {
	my $self			=	shift;
	my $outputdir		=	shift;
	my $references		=	shift;
	my $splitfiles		=	shift;
	my $infile			=	shift;
	my $hitfile			=	shift;
	my $missfile		=	shift;



	#### CHECK INPUTS
	print "Agua::Cluster::Merge::samHits    infile not defined\n" and exit if not defined $infile;
	print "Agua::Cluster::Merge::samHits    hitfile not defined\n" and exit if not defined $hitfile;


	#### GET TASKS ( = NO. SPLITFILES)
	my $tasks = scalar(@$splitfiles);

	#### GET CLUSTER
	my $cluster = $self->cluster();

	#### DO BATCH JOBS
	my $jobs = [];
	foreach my $reference ( @$references )
	{
		my $inputfile = "$outputdir/$reference/$infile";
		my $outputfile = "$outputdir/$reference//$hitfile" if defined $hitfile;
		my $missed = "$outputdir/$reference/$missfile" if defined $missfile;
		my $outdir = "$outputdir/$reference";

		#### REMOVE OUTPUT FILE IF EXISTS
		print "Agua::Cluster::Merge::samHits    outputfile is a directory: $outputfile\n" and exit if -d $outputfile;
		print `rm -fr $outputfile`;
		print "Agua::Cluster::Merge::samHits    Can't remove outputfile: $outputfile\n" and exit if -f $outputfile;

		#### SET EXECUTABLE SCRIPT TO CONVERT SAM TO TSV
		use FindBin qw($Bin);
		my $executable = "$Bin/samHits.pl";
		print "Agua::Cluster::Merge::samHits    executable: $executable\n";
		my $command = qq{/usr/bin/perl $executable --inputfile $inputfile};
		$command .= qq{ --outputfile $outputfile} if defined $hitfile;
		$command .= qq{ --missfile $missed} if defined $missfile;

		my $job = $self->setJob([$command], "samHits-$reference", $outputdir, $tasks);
		push @$jobs, $job;
	}

	#### RUN COMMANDS IN PARALLEL
	$self->runJobs( $jobs, "samHits" );	
	print "Agua::Cluster::Merge::samHits    Completed.\n";
}



=head2

	SUBROUTINE		subdirSamHits

	PURPOSE

		1. FILTER HITS AND MISSES INTO SEPARATE SAM FILES

		2. DO ALL SAM OUTPUT FILES AT THE SUBDIR LEVEL

=cut


sub subdirSamHits {
	my $self			=	shift;
	my $outputdir		=	shift;
	my $references		=	shift;
	my $splitfiles		=	shift;
	my $infile			=	shift;
	my $hitfile			=	shift;
	my $missfile		=	shift;



	#### CHECK INPUTS - infile AND hitfile REQUIRED
	print "Agua::Cluster::Merge::subdirSamHits    infile not defined\n" and exit if not defined $infile;
	print "Agua::Cluster::Merge::subdirSamHits    hitfile not defined\n" and exit if not defined $hitfile;


	#### GET TASKS ( = NO. SPLITFILES)
	my $tasks = scalar(@$splitfiles);

	#### GET CLUSTER
	my $cluster = $self->cluster();

	#### SET INDEX PATTERN FOR BATCH JOB
	my $index = $self->getIndex();

	#### DO BATCH JOBS
	my $jobs = [];
	foreach my $reference ( @$references )
	{
		my $inputfile = "$outputdir/$reference/$index/$infile";
		my $outputfile = "$outputdir/$reference/$index/$hitfile";
		my $missed = "$outputdir/$reference/$index/$missfile" if defined $missfile;
		my $outdir = "$outputdir/$reference";

		#### REMOVE OUTPUT FILE IF EXISTS
		print "Agua::Cluster::Merge::subdirSamHits    outputfile is a directory: $outputfile\n" and exit if -d $outputfile;
		print `rm -fr $outputfile`;
		print "Agua::Cluster::Merge::subdirSamHits    Can't remove outputfile: $outputfile\n" and exit if -f $outputfile;

		#### SET EXECUTABLE SCRIPT TO CONVERT SAM TO TSV
		#my $executable = "/nethome/bioinfo/apps/agua/0.5/bin/apps/samHits.pl";	
		use FindBin qw($Bin);
		my $executable = "$Bin/samHits.pl";
		print "Agua::Cluster::Merge::subdirSamHits    executable: $executable\n";
		my $command = qq{/usr/bin/perl $executable --inputfile $inputfile};
		$command .= qq{ --outputfile $outputfile};
		$command .= qq{ --missfile $missed} if defined $missfile;


		my $job = $self->setBatchJob([$command], "subdirSamHits-$reference", "$outputdir/$reference", $tasks);

		#### OVERRIDE checkfile
		$job->{checkfile} = $outputfile;
#exit;
		push @$jobs, $job;
	}


#exit;

	#### RUN COMMANDS IN PARALLEL
	$self->runJobs( $jobs, "subdirSamHits" );	
	print "Agua::Cluster::Merge::subdirSamHits    Completed.\n";
}


=head2

	SUBROUTINE		noSamHeader

	PURPOSE

		1. RETURN 1 IF FILE CONFORMS TO THE SAM FORMAT, 0 OTHERWISE

		2. SAM REQUIREMENTS:

			- HAS AT LEAST 11 FIELDS

			- 4TH AND 5TH FIELDS STRICTLY NUMERIC

			- 10TH FIELD IS SEQUENCE (I.E., NON-NUMERIC)

	NOTES

		Sequence Alignment/Map (SAM) Format
		Version 0.1.2-draft (20090820)
		http://samtools.sourceforge.net/SAM1.pdf

		The alignment section consists of multiple TAB-delimited lines with each line describing an alignment. Each line is:
		<QNAME> <FLAG> <RNAME> <POS> <MAPQ> <CIGAR> <MRNM> <MPOS> <ISIZE> <SEQ> <QUAL> \
		[<TAG>:<VTYPE>:<VALUE> [...]]

		EXAMPLE

			HWI-EAS185:1:1:24:584#0/1
			0
			chrY
			2536991
			255
			52M
			*
			0
			0
			AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
			\XXXX\^\\RNNNNMX\XYMN\X`\^JYX]W\\\[XJW[VW\X]BBBBBBBB
			XA:i:0
			MD:Z:52
			NM:i:0

=cut

sub noSamHeader {
	my $self		=	shift;
	my $file		=	shift;


	#### NB: DON'T SKIP HEADER
	open(FILE, $file) or die "Agua::Cluster::Merge::noSamHeader Can't open file: $file\n";
	my $line = <FILE>;
	my @elements = split " ", $line;
	return 0 if $#elements < 10;
	return 0 if $elements[3] !~ /^\d+$/;
	return 0 if $elements[4] !~ /^\d+$/;
	return 0 if $elements[9] !~ /^[A-Z]+$/i;

	return 1;	
}


=head2

	SUBROUTINE		removeTop

	PURPOSE

		REMOVE AN ARBITRARY NUMBER OF THE TOP ROWS OF A FILE

=cut

sub removeTop {
	my $self		=	shift;
	my $inputfile	=	shift;
	my $number_lines=	shift;
	my $outputfile	=	shift;
	my $topfile		=	shift;



	#### SET CURRENT TIME (START TIMER)
	my $current_time =  time();

	my $tempfile = $outputfile;
	$tempfile = $inputfile . "-removeTop" if not defined $outputfile;
	$topfile = $inputfile . ".top" if not defined $topfile;

	#### SAVE RECORD SEPARATOR
	my $oldsep = $/;
	$/ = "\n";

	open(FILE, $inputfile) or die "Can't open inputfile: $inputfile\n";
	open(TOPFILE, ">$topfile") or die "Can't open topfile: $topfile\n";
	for ( my $i = 0; $i < $number_lines; $i++ )
	{
		my $line = <FILE>;
		print TOPFILE $line;
	}
	close(TOPFILE);

	open(TEMPFILE, ">$tempfile") or die "Can't open tempfile: $tempfile\n";
	while ( <FILE> )
	{
		print TEMPFILE $_;
	}
	close(FILE) or die "Can't close inputfile: $inputfile\n";
	close(TEMPFILE) or die "Can't close tempfile: $tempfile\n";

	print `mv $tempfile $inputfile` if not defined $outputfile;

	#### RESTORE RECORD SEPARATOR
	$/ = $oldsep;

	my $label = "removeTop";

	#### SET EMPTY USAGE STATISTICS
	my $usage = [];

	#### GET DURATION (STOP TIMER)
	my $duration = Timer::runtime( $current_time, time() );

	#### PUSH ONTO GLOBAL USAGE STATS
	$self->addUsageStatistic($label, $duration, $usage);
}

=head2

	SUBROUTINE		isSamfile

	PURPOSE

		1. RETURN 1 IF FILE CONFORMS TO THE SAM FORMAT, 0 OTHERWISE

		2. SAM REQUIREMENTS:

			- HAS AT LEAST 11 FIELDS

			- 4TH AND 5TH FIELDS STRICTLY NUMERIC

			- 10TH FIELD IS SEQUENCE (I.E., NON-NUMERIC)

	NOTES

		Sequence Alignment/Map (SAM) Format
		Version 0.1.2-draft (20090820)
		http://samtools.sourceforge.net/SAM1.pdf

		The alignment section consists of multiple TAB-delimited lines with each line describing an alignment. Each line is:
		<QNAME> <FLAG> <RNAME> <POS> <MAPQ> <CIGAR> <MRNM> <MPOS> <ISIZE> <SEQ> <QUAL> \
		[<TAG>:<VTYPE>:<VALUE> [...]]

		EXAMPLE

			HWI-EAS185:1:1:24:584#0/1
			0
			chrY
			2536991
			255
			52M
			*
			0
			0
			AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
			\XXXX\^\\RNNNNMX\XYMN\X`\^JYX]W\\\[XJW[VW\X]BBBBBBBB
			XA:i:0
			MD:Z:52
			NM:i:0

=cut

sub isSamfile {
	my $self		=	shift;
	my $file		=	shift;
	my $skip_header	=	shift;


	open(FILE, $file) or die "Agua::Cluster::Merge::isSamfile Can't open file: $file\n";
	if ( defined $skip_header )
	{
		for ( 0..2 )
		{
			print "Agua::Cluster::Merge::isSamfile    Skipping line\n";

			<FILE>;
		}
	}

	my $line = <FILE>;
	my @elements = split " ", $line;
	print "Agua::Cluster::Merge::isSamfile    No. elements: ", $#elements + 1, "\n";
	return 0 if $#elements < 10;
	return 0 if $elements[3] !~ /^\d+$/;
	return 0 if $elements[4] !~ /^\d+$/;
	return 0 if $elements[9] !~ /^[A-Z]+$/i;

	return 1;	
}


=head2

	SUBROUTINE		_cumulativeConcatMerge

	PURPOSE

		1. SUCCESSIVELY MERGE ALL OF THE SAM FILES FOR

			EACH SPLIT INPUT FILE FOR ALL PROVIDED

			REFERENCES

		2. CUMULATIVELY ADD TO A TEMP FILE (USED TO SHOW ONGOING)

		3. ONCE COMPLETE, RENAME TEMP FILE TO OUTPUT FILE

=cut

sub _cumulativeConcatMerge {
	my $self			=	shift;
	my $outputdir		=	shift;
	my $references		=	shift;
	my $splitfiles		=	shift;
	my $infile			=	shift;
	my $outfile			=	shift;
	my $label			=	shift;

	print "Agua::Cluster::Merge::_cumulativeConcatMerge    Agua::Cluster::_cumulativeConcatMerge(outputdir, references, splitfiles, infile, outfile, label)\n";
	print "Agua::Cluster::Merge::_cumulativeConcatMerge    outputdir: $outputdir\n";
	print "Agua::Cluster::Merge::_cumulativeConcatMerge    No. references: ", scalar(@$references), "\n";
	print "Agua::Cluster::Merge::_cumulativeConcatMerge    references: ";
	print join "\n", @$references;
	print "\n";

	my $jobs = [];
	foreach my $reference ( @$references )
	{
		#### GET SAM FILES
		my $inputfiles = $self->subfiles($outputdir, $reference, $splitfiles, $infile);

		#### SET BAM OUTPUT FILE 
		my $refererence_outfile = "$outputdir/$reference/$outfile";

		#### SET OUTPUT DIR
		my ($outdir) = "$outputdir/$reference";

		my $label = "_cumulativeConcatMerge-$reference";

		#### GET CUMULATIVE MERGE COMMANDS
		my $commands = [];
		push @$commands, "rm -fr $refererence_outfile";
		foreach my $inputfile ( @$inputfiles )
		{
			push @$commands, "cat $inputfile >> $refererence_outfile";
		}

		my $job = $self->setJob($commands, $label, $outdir);
		push @$jobs, $job;
	}

	#### RUN MERGE JOBS
	$self->runJobs($jobs, $label );
}




=head2

	SUBROUTINE		cumulativeMergeBam

	PURPOSE

		1. SUCCESSIVELY MERGE ALL OF THE BAM FILES FOR

			EACH SPLIT INPUT FILE FOR ALL PROVIDED

			REFERENCES

		2. CUMULATIVELY ADD TO A TEMP FILE (USED TO SHOW ONGOING)

		3. ONCE COMPLETE, RENAME TEMP FILE TO OUTPUT FILE

=cut

sub cumulativeMergeBam {
	my $self			=	shift;
	my $outputdir		=	shift;
	my $references	=	shift;
	my $splitfiles		=	shift;
	my $infile			=	shift;
	my $outfile			=	shift;

	#### GET SAMTOOLS
	my $samtools = $self->samtools();

	my $jobs = [];
	foreach my $reference ( @$references )
	{
		#### GET SAM FILES
		my $bamfiles = $self->subfiles($outputdir, $reference, $splitfiles, $infile);

		#### SET BAM FILE 
		my $bamfile = "$outputdir/$reference/$outfile";

		#### NOTE THAT MERGE REQUIRES BOTH INPUT FILES TO MERGE TO A DISTINCT
		#### OUTPUT FILE, I.E., NOT LIKE A STRAIGHT CONCAT: infile >> outfile
		my $mergesub = sub {
			my $outfile		=	shift;
			my $tempfile	=	shift;
			my $subfile		=	shift;

			#### MERGE BAM FILES
			return "$samtools/samtools merge $outfile $tempfile $subfile";
		};

		#### GET CUMULATIVE MERGE COMMANDS
		my $commands = $self->cumulativeMergeCommands($bamfiles, $bamfile, $mergesub);
		print "Tophat::run    merge commands:\n";
		print join "\n", @$commands;
		print "\n\n";

		my $label = "cumulativeMergeBam-$reference";
		my ($outdir) = $bamfile =~ /^(.+)\/[^\/]+$/; 
		print "Agua::Cluster::Merge::run    label: $label\n";
		print "Agua::Cluster::Merge::run    outdir: $outdir\n";

		my $job = $self->setJob($commands, $label, $outdir);
		push @$jobs, $job;
	}

	#### RUN MERGE JOBS
	$self->runJobs($jobs, 'cumulativeMergeBam');
}



=head2

	SUBROUTINE		cumulativeMergeCommands

	PURPOSE

		RETURN AN ARRAY OF MERGE COMMANDS TO MERGE THE PROVIDED

		ARRAY OF SUBFILES INTO A SINGLE OUTPUT FILE


qq{time $maq/maq filemerge $outfile $tempfile $subfile}; 	

=cut

sub cumulativeMergeCommands {
	my $self		=	shift;
	my $subfiles	=	shift;
	my $outfile		=	shift;
	my $mergesub	=	shift;



	#### SET temp.file FILE FOR INCREMENTAL MERGE
	my $tempfile = "$outfile.temp";

	#### COPY FIRST SUBFILE TO TEMP FILE
	my $first_subfile = splice( @$subfiles, 0, 1);

	my $commands = [];
	#push @$commands, "echo 'Copying first subfile ($first_subfile) to temp file: $tempfile'";
	push @$commands, "time cp $first_subfile $tempfile\n";

	#### DO CUMULATIVE MERGE OF REMAINING SUBFILE FILES INTO TEMP FILE
	my $total_subfiles = scalar(@$subfiles);
	for ( my $i = 0; $i < @$subfiles; $i++ )
	{
		my $subfile = $$subfiles[$i];
		#push @$commands, "echo `ls -al $tempfile`";
		#push @$commands, "echo 'Doing subfile $i of $total_subfiles'";

		#### MERGE COMMAND
		my $merge_command = &$mergesub($outfile, $tempfile, $subfile);

		push @$commands, $merge_command;
		#push @$commands, "echo `ls -al $outfile`";	
		push @$commands, "mv $outfile $tempfile";
	}
	push @$commands, "time mv $tempfile $outfile";

	return $commands;
}	






=head2

	SUBROUTINE		merge

	PURPOSE

		1. MERGE ALL FILES PER REFERENCE IN A SERIES OF PAIRWISE

			MERGE BATCHES

		2. RUN ALL PAIRWISE MERGES ON DIFFERENT EXECUTION HOSTS

			IN PARALLEL

		THIS METHOD ASSUMES:

			1. THERE IS NO ORDER TO THE FILES

			2. ALL FILES MUST BE MERGED INTO A SINGLE FILE

			3. THE PROVIDED SUBROUTINE MERGES THE FILES

=cut

sub merge {
	my $self			=	shift;
	my $args			=	shift;

	my $inputfiles		=	$args->{inputfiles};
	my $outputfile		=	$args->{outputfile};
	my $mergesub		=	$args->{mergesub};

	print "Agua::Cluster::Merge::merge    Agua::Cluster::merge(inputfiles, outputfile, mergesub)\n";
	print "Agua::Cluster::Merge::merge    No. inputfiles: ", scalar(@$inputfiles), "\n";
	print "Agua::Cluster::Merge::merge    outputfile: $outputfile\n";

	my $counter = 0;
	while ( scalar(@$inputfiles) > 1 )
	{
		my $jobs; 

		my $label = "merge-$counter";
		$counter++;
		my ($outdir) = $$inputfiles[0] =~ /^(.+)\/[^\/]+$/; 
		print "Agua::Cluster::Merge::merge    label: $label\n";
		print "Agua::Cluster::Merge::merge    outdir: $outdir\n";

		($jobs, $inputfiles) = $self->pairwiseMergeJobs($inputfiles, $mergesub, $label, $outdir);
		print "Agua::Cluster::Merge::merge    No. inputfiles: ", scalar(@$inputfiles), "\n";
		print "Agua::Cluster::Merge::merge    Last inputfile: ", $$inputfiles[scalar(@$inputfiles) - 1], "\n";
#exit;

		## RUN JOBS
		$self->runJobs($jobs, "merge-$counter");
	}

	#### MOVE THE LAST MERGED INPUT FILE TO THE OUTPUT FILE
	my $move = "mv -f $$inputfiles[0] $outputfile";
	print "Agua::Cluster::Merge::merge    move: $move\n";
	print `$move`;
}


1;
