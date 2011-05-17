package MAQ;


=head2

	PACKAGE		MAQ

	PURPOSE

		WRAPPER SCRIPT FOR RUNNING MAQ ASSEMBLY AND SNP PREDICTION

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
#use MD5;


require Exporter;
our @ISA = qw(Exporter Cluster);
#our @EXPORT_OK = qw();
our $AUTOLOAD;


#### SET SLOTS
our @DATA = qw(

INPUTFILES
MATEFILES
REFERENCEDIR
OUTPUTDIR
TEMPDIR
MAXLINES
CLEAN
SPLITFILE
SOLEXA

PARAMS
LABEL
MAQ
SAMTOOLS
SAMTOOLSINDEX
SPECIES
VERBOSE
CONVERT

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

	SUBROUTINE		run

	PURPOSE

		1. CHECK INPUT MATE PAIR FILES

		2. DO MAQ ALIGNMENT AND SNP CALLS

			- LOCALLY (IN SERIES)

 			- ON CLUSTER (IN PARALLEL)

=cut

sub run {
	my $self			=	shift;


	#### FILES AND DIRS	
	my $referencedir 	=	$self->get_referencedir();
	my $outputdir 		=	$self->get_outputdir();
	my $inputfiles 		=	$self->get_inputfiles();
	my $matefiles 		=	$self->get_matefiles();

	#### CHECK INPUT FILES EXIST AND NOT EMPTY
	my @infiles = split ",", $inputfiles;
	my @mates = split ",", $matefiles if defined $matefiles;
	checkFiles(\@infiles);
	checkFiles(\@mates);

	#### GET LABEL, SPLIT FILES, CONVERT
	my $label 			=	$self->get_label();
	my $splitfile 		=	$self->get_splitfile();
	my $convert 		=	$self->get_convert();

	#### SET DEFAULT SPLITFILE IF NOT DEFINED
	$splitfile = "$outputdir/splitfile.txt" if not defined $splitfile;

#### INPUTS 

	###############################################
	###########    GET REFERENCEFILES    ##########
	###############################################
	my $referencefiles = $self->listReferenceFiles($referencedir, "\*.bfa");

	###############################################
	###########   SET REFERENCE NAMES    ##########
	###############################################
	my $references = $self->getFilenames($referencefiles);
	foreach my $reference ( @$references )	{	$reference =~ s/\..{1,5}$//	}
	print "MAQ::run    references: @$references\n";

	##############################################
	############	SPLIT INPUT FILES   ###########
	##############################################
	print "MAQ::run    Doing doSplitfiles()  ", Timer::current_datetime(), "\n";
	my $splitfiles = $self->doSplitfiles($splitfile, $label);
	print "MAQ::run    After doSplitfiles()  ", Timer::current_datetime(), "\n";

	##############################################
	####     CONVERT SOLEXA TO SANGER FASTQ   #####
	###############################################
	print "MAQ::run    Doing solexaToSanger  ", Timer::current_datetime(), "\n" if defined $convert;
	$self->solexaToSanger($splitfiles) if defined $convert;
	print "MAQ::run    After solexaToSanger  ", Timer::current_datetime(), "\n" if defined $convert;

	###############################################
	############# CONVERT FASTQ TO BFQ  ###########
	###############################################
	print "MAQ::run    Doing fastqToBfq      ", Timer::current_datetime(), "\n";
	$self->fastqToBfq($splitfiles, $outputdir, $label);
	print "MAQ::run    After fastqToBfq      ", Timer::current_datetime(), "\n";



#### ALIGNMENT


	###############################################
	#############     MAQ ALIGNMENT     ###########
	###############################################
	print "MAQ::run    Doing doBatchAlignment()   ", Timer::current_datetime(), "\n";
	$self->doBatchAlignment($splitfiles, $referencedir, $outputdir, $label);
	print "MAQ::run    After doBatchAlignment()   ", Timer::current_datetime(), "\n";




#### CONVERT TO CHROMOSOME SAM 


####################################################################
#### STRATEGY 1
#### CONVERT SUBDIR-LEVEL out.map TO out.sam ---> mergeSam
####################################################################

	################################################
	######     CONVERT out.map TO out.sam     ######
	################################################
	print "MAQ::run    Doing mapToSam        ", Timer::current_datetime(), "\n";
	$self->subdirMapToSam($outputdir, $references, $splitfiles, "out.map", "hit.sam");
	print "MAQ::run    After mapToSam        ", Timer::current_datetime(), "\n";

	################################################
	#######        FILTER SAM HITS          ######
	################################################
	#$self->subdirSamHits($outputdir, $references, $splitfiles, "out.sam", "hit.sam", "miss.sam");


	################################################
	######        CUMULATIVE MERGE SAM        ######
	################################################
	print "MAQ::run    Doing cumulativeMergeSam        ", Timer::current_datetime(), "\n";
	$self->cumulativeMergeSam($outputdir, $references, $splitfiles, "hit.sam", "hit.sam");
	print "MAQ::run    After cumulativeMergeSam        ", Timer::current_datetime(), "\n";

	#################################################
	########         PYRAMID MERGE SAM         ######
	#################################################
	#$self->pyramidMergeSam($outputdir, $references, $splitfiles, "hit.sam", "hit.sam");



####################################################################
#### STRATEGY 2
#### 1. MERGE SUBDIR out.map TO CHROMOSOME out.map
#### 2. CONVERT CHROMOSOME out.map TO out.sam
#####################################################################
#
#	####################################################
#	##############    CUMULATIVE MERGE MAP      ########
#	####################################################
#	#####print "MAQ::run    Doing maqSnps     ", Timer::current_datetime(), "\n";
#	#####$self->cumulativeMergeMap($splitfiles, $references, $outputdir);
#	#####print "MAQ::run    After maqSnps     ", Timer::current_datetime(), "\n";
#	####
#
#	###################################################
#	#############      PYRAMID MERGE MAP       ########
#	###############################################
#	$self->pyramidMergeMap($splitfiles, $references, $outputdir);
#
#	################################################
#	######     CONVERT out.map TO out.sam     ######
#	################################################
#	$self->mapToSam($splitfiles, $referencedir, $outputdir);
#
#	###############################################
#	######        FILTER SAM HITS          ######
#	###############################################
#	$self->samHits($outputdir, $references, $splitfiles, "out.sam", "hit.sam", "miss.sam");
#



#### READ HIT STATS


	################################################
	############        READ HITS         ##########
	################################################
	#$self->readHits($outputdir, $references, $splitfiles, "eland-readHits", "hit.sam");



#### SNP CALLING


	##############################################
	####     PREDICT SNPs WITH SAMTOOLS     ######
	##############################################
	print "MAQ::run    Doing samtoolSnps      ", Timer::current_datetime(), "\n";
	$self->samtoolSnps($splitfiles, $references, $outputdir, "hit.sam");
	print "MAQ::run    After samtoolSnps      ", Timer::current_datetime(), "\n";


	####################################################
	##############    PREDICT SNPS WITH MAQ     ########
	####################################################
	#####print "MAQ::run    Doing maqSnps     ", Timer::current_datetime(), "\n";
	#####$self->maqSnps($splitfiles, $referencedir, $outputdir);
	#####print "MAQ::run    After maqSnps     ", Timer::current_datetime(), "\n";


	###############################################
	#####          PRINT USAGE STATS         ######
	###############################################
	print "MAQ::run    Doing printUsage       ", Timer::current_datetime(), "\n";
	my $usagefile = "$outputdir/MAQ-USAGE.txt";
	$self->printUsage($usagefile);
	print "MAQ::run    After printUsage       ", Timer::current_datetime(), "\n";


	###############################################
	#####          MOVE TO FINALDIR          ######
	###############################################
	my $finaldir 		=	$self->get_finaldir();
	if ( defined $finaldir )
	{
		File::Path::mkpath($finaldir) if not -d $finaldir;
		print "Skipping move because can't create finaldir: $finaldir\n" if not -d $finaldir;
		$self->moveToDir($outputdir, $finaldir) if -d $finaldir;
	}

}	# run

=head2

	SUBROUTINE		doBatchAlignment

	PURPOSE

		RUNNING 

=cut

sub doBatchAlignment {
	my $self			=	shift;
	my $splitfiles		=	shift;
	my $referencedir	=	shift;
	my $outputdir		=	shift;
	my $label			=	shift;


	#### GET REFERENCE FILES
	my $referencefiles = $self->listReferenceFiles($referencedir, "\*.bfa");

	#### CHANGE SPLITFILE SUFFIXES TO .bfq
	foreach my $splitfile ( @$splitfiles )
	{
		$$splitfile[0] =~ s/\.fastq$/.bfq/;
		$$splitfile[1] =~ s/\.fastq$/.bfq/ if defined $$splitfile[1];
	}

#### DEBUG
@$referencefiles = reverse @$referencefiles;

#
##print "splitfiles:\n";
##print Dumper $splitfiles;
##exit;
#
#	#### SET NUMBER OF SPLITFILES FOR GENERATING BATCH JOB LATER
#	my $number_splitfiles = scalar(@$splitfiles);
#
#
##exit;
#
#	#### DO ALIGNMENTS AGAINST ALL REFERENCE FILES
#	my $jobs = [];
#	for my $referencefile ( @$referencefiles )
#	{
#		next if $referencefile =~ /^[\.]+$/;
#	
#		$self->set_referencefile($referencefile);
#		
#		#### SET REFERENCE BINARY FILE
#		my $referencebinary = $referencefile;
#		$referencebinary =~ s/\.[^\.]+?$/.bfa/;
#	
#		#### SET REFERENCE
#		my ($reference) = $referencefile =~ /([^\/]+)\.bfa$/i;
#
#		#### SET OUTPUT DIR FOR THIS REFERENCE
#		my $outdir = "$outputdir/$reference";
#
#		#### DO ALIGNMENTS USING BATCH COMMAND
#		my $batch_command = $self->batchCommand($outputdir, $referencebinary, $splitfiles);
#		
#		#### SET LABEL
#		my $label = "maq-$reference";
#		
#		#### SET JOB
#		my $job = $self->setBatchJob( [$batch_command], $label, $outdir, $number_splitfiles);
#		push @$jobs, $job;
#
#
###### DEBUG
##last;
#
#
#	}	

	#### COLLECT ALL JOBS FOR EACH INPUT FILE AGAINST ALL REFERENCE FILES
	my $jobs = $self->generateBatchJobs($outputdir, $referencefiles, $splitfiles, $label);


	print "MAQ::doBatchAlignment    No. jobs: ", scalar(@$jobs), "\n";

	#### RUN ALIGNMENT JOBS
	$self->runJobs($jobs, "doBatchAlignment");
}


=head2

	SUBROUTINE		batchCommand

	PURPOSE

		CREATE BATCH COMMAND WITH PLACEHOLDERS FOR TASK NUMBER

	INPUTS

		1. BASE OUTPUT DIRECTORY CONTAINING SUBDIR ALIGNMENTS

		2. PATH TO BINARY REFERENCE .bfa FILE

		3. ARRAY OF INPUT FILES, WHERE EACH ELEMENT IS A SHORT

			ARRAY OF ONE OR TWO FILES: READ AND ITS MATE IF AVAILABLE

=cut

sub batchCommand {
	my $self			=	shift;
	my $outputdir		=	shift;
	my $referencebinary	=	shift;
	my $splitfiles		=	shift;


	#### GET LABEL FOR LATER USE TO GENERATE INPUTFILE NAMES
	my $label 			=	$self->get_label();	
	my $maq 			=	$self->get_maq();	
	my $params 			=	$self->get_params();
	$params = "" if not defined $params;

	##### CLUSTER
	my $cluster 			=	$self->get_cluster();

	#### GET THE BASE DIRECTORY OF THE SPLIT FILES - ALWAYS
	#### TWO DIRECTORIES DOWN FROM THE SPLIT FILE
	my $splitfile = $$splitfiles[0][0];
	print "splitfile: $splitfile\n";
	my ($basedir) = $splitfile =~ /^(.+?)\/\d+\/([^\/]+)$/;
	print "basedir: $basedir\n";

	#### SET BINARY FASTQ .bfq INPUT FILE SUFFIX
	my $suffix = ".bfq";

	#### SET SCHEDULAR JOB INDEX IN INPUT BINARY FILES
	my $index;
	$index = "\$LSB_JOBINDEX" if $cluster eq "LSF";
	$index = "\$PBS_TASKNUM" if $cluster eq "PBS";	
	my $inputbinaries;
	my $firstmate = $label . "_1";
	my $secondmate = $label . "_2";
	push @$inputbinaries, "$basedir/$index/$firstmate.$index$suffix";
	push @$inputbinaries, "$basedir/$index/$secondmate.$index$suffix" if defined $$splitfiles[0][1];	

	#### GET REFERENCE NAME
	my ($reference) = $referencebinary =~ /^.+?\/([^\/]+)\.bfa$/;

	#### SET OUTPUTDIR FILE
	my $outdir = "$outputdir/$reference/$index";

	#### SET UNIQUE OUTERR FILE FOR EACH REFERENCE ALIGNMENT
	my $outerrfile = "$outdir/maq-$reference-outerr.txt";

	#### SET UNIQUE OUTERR FILE FOR EACH REFERENCE ALIGNMENT
	my $mapfile = "$outdir/out.map";

	#### OUTPUT TO outputdir ACROSS NFS
    my $command = qq{
mkdir -p $outdir;
cd $outdir;
time $maq/maq match $params $mapfile $referencebinary @$inputbinaries  &> $outerrfile
};

	#### CHANGE COMMAND IF TEMP DIR DEFINED
	my $tempdir = $self->get_tempdir();
	if ( defined $tempdir and $tempdir )
	{

		#### SET TEMP OUTPUT DIR
		my $temp_outputdir = $tempdir . "/" . $outdir;

		#### SET UNIQUE OUTERR FILE FOR EACH REFERENCE ALIGNMENT
		my $outerrfile = "$temp_outputdir/maqBatch-$reference-outerr.txt";

		#### OUTPUT TO temp_outputdir ON LOCAL EXECUTION HOST
		my $mapfile = "$temp_outputdir/out.map";

		#### OUTPUT TO /tmp ON EXECUTION HOST AND MOVE AFTER COMPLETED
		$command = qq{
mkdir -p $temp_outputdir
time $maq/maq match $params $mapfile $referencebinary @$inputbinaries  &> $outerrfile
mv $temp_outputdir/* $outdir
rm -fr $temp_outputdir
};
	}

	return $command;

}	#	batchCommand







=head2

	SUBROUTINE		faToBfa

	PURPOSE

		CONVERT .fa FASTA REFERENCE FILES TO .bfa BINARY FASTA FORMAT

=cut

sub faToBfa {
	my $self			=	shift;
	my $inputdir		=	shift;
	my $outputdir		=	shift;



	#### SANITY CHECK
	print "MAQ::faToBfa    inputdir not defined: $inputdir\n" and exit if not defined $inputdir;
	print "MAQ::faToBfa    outputdir not defined: $outputdir\n" and exit if not defined $outputdir;
	print "MAQ::faToBfa    inputdir is a file: $inputdir\n" and exit if -f $inputdir;
	print "MAQ::faToBfa    outputdir is a file: $outputdir\n" and exit if -f $outputdir;
	print "MAQ::faToBfa    Can't find inputdir: $inputdir\n" and exit if not -d $inputdir;
	print "MAQ::faToBfa    Can't find outputdir: $outputdir\n" and exit if not -d $outputdir;

	#### GET REFERENCE FILES
	my $referencefiles = $self->listReferenceFiles($inputdir, "\*.fa");
	print "referencefiles: @$referencefiles\n";

	###############################################
	############# MAKE BINARY REF FILES ###########
	###############################################

	my $current_time = time();	
	for my $referencefile ( @$referencefiles )
	{
		next if $referencefile =~ /^[\.]+$/;

		my ($reference) = $referencefile =~ /^.+?\/([^\/]+)$/;

		#### SET REFERENCE BINARY FILE
		my $referencebinary = "$outputdir/$reference";
		$referencebinary =~ s/\.fa$/.bfa/;

		#### CONVERT REFERENCE FILE INTO BINARY FILE
		my $convert_command = $self->convertReference($referencefile, $referencebinary);
		print `$convert_command` if defined $convert_command;
	}

	#### GET DURATION OF LOCAL JOB AND ADD AS USAGE STATISTIC
	my $duration = Timer::runtime( $current_time, time() );
	$self->addUsageStatistic("faToBfa", $duration);
}

=head2

	SUBROUTINE		convertReference

	PURPOSE

		CONVERT REFERENCE FASTA INTO MAQ BINARY FASTA FORMAT

=cut

sub convertReference {
	my $self			=	shift;
	my $referencefile	=	shift;
	my $referencebinary	=	shift;

	my $maq		=	$self->get_maq();
	print "MAQ::convertReference    Maq not defined\n" and exit if not defined $maq or not $maq;

	#### SKIP IF ALREADY EXISTS
	print "MAQ::convertReference    Skipping conversion because binary reference file already exists: $referencebinary\n" and next if -f $referencebinary and not -z $referencebinary;

	#### OTHERWISE, DO CONVERSION
	my $command = "time $maq/maq fasta2bfa $referencefile $referencebinary" if not -f $referencebinary or -z $referencebinary;

	return $command;
}




=head2

	SUBROUTINE		fastqToBfq

	PURPOSE

		CONVERT A LIST OF FASTQ FILES (AND MATES) TO BFQ FORMAT

=cut

sub fastqToBfq {
	my $self		=	shift;
	my $splitfiles	=	shift;
	my $outputdir	=	shift;
	my $label		=	shift;



	#### MAQ, CONVERT, AND CLEAN
	my $maq				=	$self->get_maq();
	my $convert			=	$self->get_convert();
	my $clean			=	$self->get_clean();

	#### INITIALISE JOBS
	my $jobs = [];

	my $counter = 0;
	foreach my $splitfile ( @$splitfiles )
	{
		$counter++;

		#my $commands = $self->fastqToBfqCommands($splitfile);

		my $commands = [];
		foreach my $inputfile ( @$splitfile )
		{
			my $bfqfile = $inputfile;
			$bfqfile =~ s/\.[^\.]+?$/.bfq/;



			push @$commands, "#### Converting .fastq file to .bfq file...";
			push @$commands, "time $maq/maq fastq2bfq $inputfile $bfqfile";		
		}


		next if not defined $commands or scalar(@$commands) == 0;

		#### SET LABEL
		my $this_label = "fstqToBfq-$counter";

		#### SET OUTPUT DIRECTORY
		my ($outdir) = $$splitfile[0] =~ /^(.+?)\/[^\/]+$/;

		#### SET JOB
		my $job = $self->setJob($commands, $this_label, $outdir);

		push @$jobs, $job;
	}		

	print "MAQ::fastqToBfq    Skipping fastqToBfq file conversion because no jobs to run\n" and return if scalar(@$jobs) == 0;

	#### RUN ALIGNMENT JOBS
	print "MAQ::    RUNNING ", scalar(@$jobs), " jobs...\n";
	$self->runJobs($jobs, 'fastqToBfq');
}


=head2

	SUBROUTINE		solexaToSanger

	PURPOSE

		CONVERT SOLEXA FASTQ FILES TO SANGER FASTQ FILES

=cut

sub solexaToSanger {
	my $self		=	shift;
	my $splitfiles	=	shift;


	my $jobs = [];
	my $counter = 0;
	foreach my $splitfile ( @$splitfiles )
	{
		$counter++;


		my $infile = $$splitfile[0];

		my $commandsHash = $self->solToSangerCommands($splitfile);
		my $commands = $commandsHash->{commands};
		my $files = $commandsHash->{files};
		next if not defined $commands or scalar(@$commands) == 0;

		#### SET LABEL
		my $label = "solexaToSanger-$counter";

		#### SET OUTPUT DIRECTORY
		my ($outputdir) = $infile =~ /^(.+?)\/[^\/]+$/;

		#### SET JOB
		my $job = $self->setJob($commands, $label, $outputdir);
		push @$jobs, $job;
	}


	#### RUN ALIGNMENT JOBS
	$self->runJobs($jobs, 'SolexaToSanger');
}

=head2

	SUBROUTINE		solToSangerCommands

	PURPOSE

		RETURN COMMANDS TO CONVERT FROM SOLEXA- TO SANGER-FORMAT 

		FASTQ FILES AND THE LIST OF .fastq FILES TO BE CREATED

	INPUTS

		1. ARRAY OF INPUT .txt OR .fastq FILES

	OUTPUTS

		1. ARRAY OF COMMANDS

		2. ARRAY OF .bfq FILE NAMES

=cut

sub solToSangerCommands {
	my $self			=	shift;
	my $inputfiles		=	shift;


	#### MAQ, CONVERT, AND CLEAN
	my $maq				=	$self->get_maq();
	my $convert			=	$self->get_convert();
	my $clean			=	$self->get_clean();


	#### GET TEMP DIR
	my $tempdir = $self->get_tempdir();

	my $commands = [];
	my $fastqfiles = [];
	if ( defined $convert and $$inputfiles[0] !~ /\.bfq$/ )
	{
		foreach my $inputfile ( @$inputfiles )
		{
            #### OUTPUT TO OUTPUTDIR ACROSS NFS
			#### SET CONVERSION OUTPUT FILE SUFFIX TO .fastq			
			my $fastqfile = $inputfile;
			if ( $inputfile =~ /fastq$/ )
			{
				$fastqfile =~ s/fastq$/sanger.fastq/;
			}
			elsif ( $inputfile =~ /txt$/ )
			{
				$fastqfile =~ s/\.txt$/.fastq/;
			}
			else
			{
				print "Quitting - Input file has incorrect suffix:	$inputfile\n";
				exit;
			}

			push @$commands, "echo 'Converting solexa sequence file to Sanger fastq file'";

			if ( $convert =~ /^post-1.3$/ )
			{
				push @$commands, "time $maq/maq ill2sanger $inputfile $fastqfile";
			}
			elsif ( $convert =~ /^pre-1.3$/ )
			{
				push @$commands, "time $maq/maq sol2sanger $inputfile $fastqfile";
			}
			else
			{
				print "MAQ::alignmentCommands    Conversion type not supported: $convert (must be 'post-1.3' or 'pre1.3')\n" and exit;
			}

			#### SAVE .fastq FILE NAME
			push @$fastqfiles, $fastqfile;
		}
	}

	return { commands => $commands, files => $fastqfiles };	
}


=head2

	SUBROUTINE		getReferenceFiles

	PURPOSE

		RETURN A LIST OF REFERENCE FILES (FULL PATHS TO FILES)

=cut

sub getReferenceFiles {
	my $self		=	shift;
	my $referencedir=	shift;


	return $self->listReferenceFiles($referencedir, "\*\.bfa");
}



=head2

	SUBROUTINE		getReferences

	PURPOSE

		RETURN A LIST OF REFERENCE NAMES (N.B.: NOT THE PATHS TO THE FILES)

=cut

sub getReferences {
	my $self		=	shift;
	my $referencedir=	shift;


	my $referencefiles = $self->getReferenceFiles($referencedir);
	my $references = $self->getFilenames($referencefiles);
	foreach my $reference ( @$references )	{	$reference =~ s/\.bfa$//;	}

	return $references;
}




=head2

	SUBROUTINE		subdirMapToSam

	PURPOSE

		CONVERT ALL .map FILES IN chr*/<NUMBER> SUBDIRECTORIES

		INTO .sam FILES

			maq2sam-long out.map > out.sam

		INPUTS

			1. OUTPUT DIRECTORY USED TO PRINT CHROMOSOME-SPECIFIC

				SAM FILES TO chr* SUB-DIRECTORIES

			2. LIST OF REFERENCE FILE NAMES

			3. SPLIT INPUT FILES LIST

			4. NAME OF MAP INPUTFILE (E.G., "out.map")

			4. NAME OF SAM OUTPUTFILE (E.G., "accepted_hits.sam")

=cut


sub subdirMapToSam {
	my $self			=	shift;
	my $outputdir		=	shift;
	my $references 		=	shift;
	my $splitfiles		=	shift;
	my $infile			=	shift;
	my $outfile			=	shift;

	print "Cluster::subdirMapToSam    Cluster::subdirMapToSam(outputdir, references, splitfiles, inputfile)\n";
	print "Cluster::subdirMapToSam    outputdir: $outputdir\n";
	print "Cluster::subdirMapToSam    infile: $infile\n";
	print "Cluster::subdirMapToSam    outfile: $outfile\n";
	print "Cluster::subdirMapToSam    references: @$references\n";

	#### GET REQUIRED VARIABLES
	my $samtools = $self->get_samtools();
	my $maq = $self->get_maq();
	print "MAQ::subdirMapToSam    samtools not defined\n" and return if not defined $samtools;

	#### NB: maq2sam COMMAND ASSUMES out.map FILES GENERATED BY maq-0.7.x.
	#### maq2sam-short IS FOR maq-0.6.x AND EARLIER
	#### PARSE THIS NUMBER OUT FROM MAQ LOCATION WITH FORMAT: ../maq/1.2.3
	my $maq2sam = "maq2sam-long";
	my ($maq_version) = $maq =~ /(\d)\.\d$/;
	$maq2sam = "maq2sam-short" if $maq_version < 7;

	#### SET INDEX
	my $cluster = $self->get_cluster();
	my $index;
	$index = "\$LSB_JOBINDEX" if $cluster eq "LSF";
	$index = "\$PBS_TASKNUM" if $cluster eq "PBS";

	#### SET NUMBER OF SPLITFILES FOR GENERATING BATCH JOB LATER
	my $number_splitfiles = scalar(@$splitfiles);

	my $jobs = [];
	foreach my $reference ( @$references )
	{
		#### SET JOB TO CONVERT ALL MAP FILES FOR THIS REFERENCE
		#### INTO SAM FILES 
		my $outdir = "$outputdir/$reference";
		my $mapfile = "$outdir/$index/$infile";
		my $samfile = "$outdir/$index/$outfile";
		my $command = "$samtools/$maq2sam $mapfile > $samfile";
		my $label = "subdirMapToSam-$reference";
		my $job = $self->setBatchJob( [$command], $label, $outdir, $number_splitfiles);

		push @$jobs, $job;
	}

	#### RUN CONVERSION JOBS
	print "Cluster::subdirMapToSam    DOING runJobs for " , scalar(@$jobs), " jobs\n";
	$self->runJobs($jobs, 'subdirMapToSam');
	print "Cluster::subdirMapToSam    Completed subdirMapToSam\n";
}



=head2

	SUBROUTINE		mapToSam

	PURPOSE

		RUNNING 

=cut

sub mapToSam {	
	my $self			=	shift;
	my $splitfiles		=	shift;
	my $referencedir	=	shift;
	my $outputdir		=	shift;
	my $mapfile			=	shift;
	my $samfile			=	shift;


	#### GET REFERENCE FILES
	my $referencefiles = $self->listReferenceFiles($referencedir, "\*.bfa");

	#### CONVERT ALL .map FILES
	my $jobs = [];
	for my $referencefile ( @$referencefiles )
	{
		next if $referencefile =~ /^[\.]+$/;

		$referencefile = "$referencedir/$referencefile";

		#### SET REFERENCE BINARY FILE
		my $referencebinary = $referencefile;
		$referencebinary =~ s/\.[^\.]+?$/.bfa/;

		#### SET REFERENCE
		my ($reference) = $referencefile =~ /^.+?\/([^\/]+)\.bfa$/i;

		#### DO ALIGNMENTS FOR ALL SPLITFILES 
		my $counter = 0;
		foreach my $splitfile ( @$splitfiles )
		{
			$counter++;

			#### SET OUTPUT DIR TO SUBDIRECTORY CONTAINING SPLITFILE
			my ($basedir, $index) = $$splitfile[0] =~ /^(.+?)\/(\d+)\/[^\/]+$/;
			my $outdir = "$basedir/$reference/$index";

			#### SET *.map FILE
			my $mapfile = "$outdir/$mapfile";

			#### SET .sam FILE	
			my $samfile = "$outdir/$samfile";

			#### CONVERT .map FILE TO .sam FILE
			my $commands;
			my $sam_command = $self->mapToSamCommand($mapfile, $samfile);
			push @$commands, $sam_command;

			#### SET LABEL
			my $label = "$reference-$counter-mapToSam";

			#### SET JOB
			my $job = $self->setJob($commands, $label, $outdir);
			push @$jobs, $job;
		}
	}	

	print "MAQ::mapToSam    No. mapToSam jobs: ", scalar(@$jobs), "\n";

	#### RUN ALIGNMENT JOBS
	$self->runJobs($jobs, "mapToSam");
}

=head2

	SUBROUTINE		mapToSamCommand

	PURPOSE

		RETURN THE COMMAND TO CONVERT A *.map FILE TO A .sam FILE

=cut

sub mapToSamCommand {
	my $self		=	shift;
	my $maqfile		=	shift;
	my $samfile		=	shift;


	#### SANITY CHECK
	print "MAQ::mapToSamCommand    maqfile not defined: $maqfile\n" and exit if not defined $maqfile;
	print "MAQ::mapToSamCommand    samfile not defined: $samfile\n" and exit if not defined $samfile;

	#### GET REQUIRED VARIABLES
	my $samtools = $self->get_samtools();
	my $maq = $self->get_maq();
	print "MAQ::mapToSamCommand    samtools not defined\n" and return if not defined $samtools;

	#### NB: maq2sam COMMAND ASSUMES out.map FILES GENERATED BY maq-0.7.x.
	#### maq2sam-short IS FOR maq-0.6.x AND EARLIER
	my $maq2sam = "maq2sam-long";
	my ($maq_version) = $maq =~ /(\d)\.\d$/;
	$maq2sam = "maq2sam-short" if $maq_version < 7;

    my $command = "$samtools/$maq2sam $maqfile > $samfile";

	return $command;	
}



=head2

	SUBROUTINE		maqSnps

	PURPOSE

		RUNNING 

=cut

sub maqSnps {
	my $self			=	shift;
	my $splitfiles		=	shift;
	my $referencedir	=	shift;
	my $outputdir		=	shift;


	#### SET TIMER
	my $current_time = time();

	#### GET REFERENCE FILES
	my $referencefiles = $self->listReferenceFiles($referencedir, "\*.bfa");

	#### GET COMMANDS TO MERGE MAPS AND PREDICT SNPS AND INDELS
	my $snp_jobs = [];
	for my $referencefile ( @$referencefiles )
	{
		next if $referencefile =~ /^[\.]+$/;

		#### SET REFERENCE BINARY FILE
		my $referencebinary = $referencefile;
		$referencebinary =~ s/\.[^\.]+?$/.bfa/;

		#### SET REFERENCE
		my ($reference) = $referencefile =~ /^.+\/([^\/]+)\.[^\.]{2,6}$/;

		#### SET OUTPUT DIRECTORY
		my $outdir = "$outputdir/$reference";

		#### GET SNP PREDICTION COMMANDS
		my $commands = $self->snpCommands($splitfiles, $referencefile, $referencebinary, $outdir);

		#### SET LABEL
		my $label = "maqSnps-" . $reference;

		#### GET LABEL, COMMANDS, ETC.
		my $job = $self->setJob($commands, $label, $outdir);	
		push @$snp_jobs, $job;
	}
	print "MAQ::    length(jobs): ", scalar(@$snp_jobs), "\n";

	#### RUN out.map MERGE AND SNP/INDEL CALLING JOBS
	$self->runJobs($snp_jobs, "maqSnps");
}


=head2

	SUBROUTINE		snpCommands

	PURPOSE

		RETURN AN ARRAY OF COMMANDS TO ANALYSE SNPS AND INDELS

=cut

sub snpCommands {
	my $self			=	shift;
	my $splitfiles		=	shift;
	my $referencefile 	=	shift;
	my $referencebinary	=	shift;
	my $outputdir		=	shift;


	#### EXECUTABLES
	my $maq 			=	$self->get_maq();

	#### STORE COMMANDS
	my $commands = [];
	push @$commands, "echo 'Changing to output directory: $outputdir'";
	push @$commands, "cd $outputdir";

	#######################
	##### DO SNPS
	#######################

	# 3. Build the mapping assembly
	push @$commands, "time $maq/maq assemble consensus.cns $referencebinary out.map 2> assemble.log";

	# 4. Extract consensus sequences and qualities
	push @$commands, "time $maq/maq cns2fq consensus.cns > cns.fq";

	# 5. Extract list of SNPs 
	push @$commands, "time $maq/maq cns2snp consensus.cns > cns.snp";


	#######################
	##### DO INDELS
	#######################

	#2. rmdup       remove pairs with identical outer coordinates (PE)
	push @$commands, "time $maq/maq rmdup out.rmdup out.map";

	#3. indelpe     indel calling (PAIRED READS ONLY)
	push @$commands, "time $maq/maq indelpe $referencebinary out.rmdup > out.indelpe";

	#4. indelsoa    state-of-art homozygous indel detectionll
	push @$commands, "time $maq/maq indelsoa $referencebinary out.map > out.indelsoa";

	#5. filter indels
	push @$commands, "awk '\$5+\$6-\$4 >= 3 \&\& \$4 <= 1' out.indelsoa > out.indelsoa.filter";

	#6. SNPfilter    filter SNP predictions
	push @$commands, "time $maq/scripts/maq.pl SNPfilter -d 1 -s out.indelsoa -F out.indelpe cns.snp &> out.SNPfilter";

#exit;

	return $commands;	
} 




=head2

	SUBROUTINE		pyramidMergeMap

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

sub pyramidMergeMap {
	my $self			=	shift;
	my $outputdir		=	shift;
	my $references		=	shift;
	my $splitfiles		=	shift;
	my $infile			=	shift;
	my $outfile			=	shift;



	#### SET DEFAULT infile AND outfile
	$infile = "out.map" if not defined $infile;
	$outfile = "out.map" if not defined $outfile;

	#### GET SAMTOOLS
	my $maq = $self->get_maq();

	#### LOAD UP WITH INITIAL BAM FILES
	my $reference_mapfiles;
	foreach my $reference ( @$references )
	{
		#### GET SAM FILES
		my $mapfiles = $self->subfiles($outputdir, $reference, $splitfiles, $infile);
		$reference_mapfiles->{$reference} = $mapfiles;
	}

	my $mergesub = sub {
		my $firstfile	=	shift;
		my $secondfile	=	shift;


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


		#### MERGE MAP FILES
		my $command .= "time $maq/maq mapmerge $outfile $firstfile $secondfile;\n";

		return ($command, $outfile);
	};


	#### RUN PYRAMID MERGE ON ALL REFERENCES IN PARALLEL
	my $label = "pyramidMergeMap";
	$self->_pyramidMerge($outputdir, $references, $splitfiles, $infile, $outfile, $reference_mapfiles, $label, $mergesub);

	print "MAQ::pyramidMergeMap    Completed.\n";
}


=head2

	SUBROUTINE		cumulativeMergeMap

	PURPOSE

		CUMULATIVELY MERGE ALL INPUT SPLIT FILE-BASED out.map FILES INTO A SINGLE

		out.map FILE FOR EACH REFERENCE

=cut

sub cumulativeMergeMap {
	my $self			=	shift;
	my $splitfiles		=	shift;
	my $references		=	shift;
	my $outputdir		=	shift;


	#### GET MAQ
	my $maq = $self->get_maq();

	#### SET TIMER
	my $current_time = time();

	#### GET COMMANDS TO MERGE MAPS AND PREDICT SNPS AND INDELS
	my $jobs = [];
	for my $reference ( @$references )
	{		
		#### GET SAM FILES
		my $infile = "out.map";
		my $mapfiles = $self->subfiles($outputdir, $reference, $splitfiles, $infile);

		#### CHECK FOR MISSED MAP FILES (FAILED ALIGNMENTS)
		$mapfiles = $self->missedMaps($mapfiles);

		#### SET OUTPUT DIRECTORY
		my $outdir = "$outputdir/$reference";

		#### SET OUTPUT FILE
		my $outfile = "$outdir/out.map";

		#### SET MERGE SUBROUTINE
		my $mergesub = sub {
			my $outfile		=	shift;
			my $tempfile	=	shift;
			my $subfile		=	shift;

			return "time $maq/maq mapmerge $outfile $tempfile $subfile";
		};

		#### GET SNP PREDICTION COMMANDS
		my $commands = $self->cumulativeMergeCommands($mapfiles, $outfile, $mergesub);

		#### SET LABEL
		my $label = "cumulativeMergeMap-" . $reference;

		#### GET LABEL, COMMANDS, ETC.
		my $job = $self->setJob($commands, $label, $outdir);	
		push @$jobs, $job;
	}
	print "MAQ::cumulativeMergeMap    length(jobs): ", scalar(@$jobs), "\n";

	#### RUN out.map MERGE AND SNP/INDEL CALLING JOBS
	$self->runJobs($jobs, "cumulativeMergeMap");
} 



=head2

	SUBROUTINE		missedMaps

	PURPOSE

		CHECK TO SEE IF ANY out.map FILES ARE MISSING

=cut

sub missedMaps {
	my $self		=	shift;
	my $mapfiles	=	shift;

	#### CHECK THAT ALL FILES EXIST
	my $missedmaps = [];
	for ( my $i = 0; $i < @$mapfiles; $i++ )
	{
		my $mapfile = $$mapfiles[$i];
		if ( not -f $mapfile )
		{
			push @$missedmaps, splice (@$mapfiles, $i, 1);
			$i--;
		}
	}

	#### PRINT ANY MISSING out.map FILES
	if ( scalar(@$missedmaps) > 0 )
	{
		print "MAQ::missedMaps    Missed maps:\n";
		print join "\n", @$missedmaps;
		print "\n\n";
		print "MAQ::missedMaps    Exiting\n";
		exit;
	}

	return $mapfiles;	
}


=head2

	SUBROUTINE		checkFiles

	PURPOSE

		CHECK INPUT FILES ARE ACCESSIBLE AND NON-EMPTY

=cut

sub checkFiles {
	my $self		=	shift;
	my $files		=	shift;

	#### SANITY CHECK
	foreach my $file ( @$files )
	{
		if ( $file !~ /\.(txt|fastq|bfq)$/ )
		{
			print "MAQ::::alignmentCommands     file must end in .txt, .fastq or .bfq\n";
			print "MAQ::::alignmentCommands    file: $file\n";
			exit;
		}
		if ( -z $file )
		{
			print "MAQ::::alignmentCommands     file is empty. Quitting.\n";
			exit;
		}
	}
}



=head2

	SUBROUTINE		fastqToBfqCommands

	PURPOSE

		RETURN COMMANDS TO CONVERT .fastq FILES TO .bfq FORMAT AND

		THE LIST OF .bfq FILES TO BE CREATED

	INPUTS

		1. ARRAY OF INPUT .fastq FILES

	OUTPUTS

		1. ARRAY OF COMMANDS

		2. ARRAY OF .bfq FILE NAMES

	NOTES

		IF THE convert FLAG IS SET, .txt OR .fastq INPUT FILES WILL 

		BE CONVERTED FROM SOLEXA TO SANGER FASTQ SEQUENCES. IN THE

		CASE OF .bfq INPUT FILES, CONVERSION WILL BE SKIPPED.

=cut

sub fastqToBfqCommands {
	my $self			=	shift;
	my $inputfiles		=	shift;



	#### MAQ, CONVERT, AND CLEAN
	my $maq				=	$self->get_maq();
	my $convert			=	$self->get_convert();
	my $clean			=	$self->get_clean();

	my $commands = [];
	foreach my $inputfile ( @$inputfiles )
	{
		my $bfqfile = $inputfile;
		$bfqfile =~ s/\.[^\.]+?$/.bfq/;

		push @$commands, "#### Converting .fastq file to .bfq file...";
		push @$commands, "time $maq/maq fastq2bfq $inputfile $bfqfile";		
	}

	return $commands;
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


1;


__END__

#=head2
#
#	SUBROUTINE		mapcheck
#	
#	PURPOSE
#			
#		GET STATISTICS ABOUT MAQ ALIGNMENTS FOR EACH INPUT FILE CHUNK
#
#		 NB: WATCH THIS: MOVED TO align SUBROUTINE EARLIER BECAUSE GOT 
#		 SEGMENTATION FAULT WITH MAP FILE OF 331 MB
#
#=cut
#
#sub mapcheck
#{
#	my $self			=	shift;
#	
#
#	my $splitfiles		=	shift;
#	my $referencedir	=	shift;
#	my $outputdir		=	shift;
#
#
#	#### SET TIMER
#	my $current_time = time();
#
#	#### GET REFERENCE FILES
#	my $referencefiles = $self->listReferenceFiles($referencedir, "\*.bfa");
#	
#	#### GET COMMANDS TO MERGE MAPS AND PREDICT SNPS AND INDELS
#	my $snp_jobs = [];
#	for my $referencefile ( @$referencefiles )
#	{
#		my $command =  "time $maq/maq mapcheck $referencefile out.map > mapcheck.txt";
#		my $job = $self->setJob( [$command], )
#	}
#}

####=head2
####
####	SUBROUTINE		matchCommands
####	
####	PURPOSE
####	
####		RETURN COMMANDS TO USE MAQ match TO ALIGN SEQUENCES AGAINST REFERENCE
####		
####=cut
####
####sub matchCommands
####{
####	#print "\@_:\n";
####	#print Dumper @_;
####
####	my $self			=	shift;
####	my $inputfiles		=	shift;
####	my $outputdir		=	shift;
####	my $referencefile	=	shift;
####	my $referencebinary	=	shift;
####	my $mapfile			=	shift;
####
####
####	#### SET REFERENCE
####	my ($reference)		= $referencefile =~ /([^\/]+?)\.[^\.]{2,6}$/;
####	#print "MAQ::matchCommands    reference: $reference\n";
####	#exit;
####
####	#### GET MAQ BINARY
####	my $maq				=	$self->get_maq();
####
####	#### INPUT FILES MUST BE .bfq FORMAT, CHECK THEY HAVE THE RIGHT SUFFIX	
####	foreach my $inputfile ( @$inputfiles )
####	{
####		if ( $inputfile !~ /\.bfq$/ )
####		{
####			print "MAQ::matchCommands    input file must end in .bfq\n";
####			print "MAQ::matchCommands    inputfile: $inputfile\n";
####			exit;
####		}
####	}
####	
####	#### CHECK FILES AND DIRS
####	print "MAQ::matchCommands    Output dir is a file: $outputdir\n" and exit if -f $outputdir;
####	File::Path::mkpath($outputdir) if not -d $outputdir;
####	die "MAQ::matchCommands    Couldn't create output dir: $outputdir\n" if not -d $outputdir;
####
####	#### ALIGNMENT COMMANDS
####	my $commands;
####
####	# 3. SET INPUT BINARIES
####	my $inputbinaries = join " ", @$inputfiles;
####	
####	#### SET UNIQUE OUTERR FILE FOR EACH REFERENCE ALIGNMENT
####	my $outerrfile = "$outputdir/$reference-outerr.txt";
####
####	# 4. COMMAND TO ALIGN THE READS AGAINST THE REFERENCE
####	push @$commands, "#### Doing maq 'match' alignment...";
####	
####	#### GET TEMP DIR
####	my $tempdir = $self->get_tempdir();
####	#print "MAQ::matchCommands    self->{_tempdir}: $self->{_tempdir}\n";
####	#print "MAQ::matchCommands    tempdir: $tempdir\n";
####
####	#### OUTPUT TO outputdir ACROSS NFS
####    my $command = "time $maq/maq match $mapfile $referencebinary $inputbinaries  &> $outerrfile";
####	push @$commands, $command;
####
####	#### OUTPUT TO /tmp ON EXECUTION HOST AND MOVE AFTER COMPLETED
#####    my $command = "time $maq/maq match $tempdir/$reference-out.map $referencebinary $inputbinaries  &> $outerrfile";
#####	push @$commands, $command;
#####	print "command: $command\n";
#####    $command = "time mv $tempdir/$reference-out.map $outputdir/$reference-out.map";
####	#push @$commands, $command;
####
####	return $commands;
####}	
####	
####
####
####
####=head2
####
####	SUBROUTINE		doAlignments
####	
####	PURPOSE
####	
####		ALIGN ALL SEQUENCE READ SPLITFILES USING maq AGAINST ALL REFERENCE SEQUENCES
####
####=cut
####
####sub doAlignments
####{
####	my $self		=	shift;
####	my $splitfiles	=	shift;
####	my $referencefile	=	shift;
####
####
####	
####	#### INPUTS
####	my $dot = $self->get_dot();
####	my $outputdir = $self->get_outputdir();
####
####	#### SET REFERENCE
####	my ($reference) = $referencefile =~ /([^\/]+)$/;
####	$reference =~ s/(\.fa|\.fasta)$//i;
####
####	#### MONITOR
####	my $monitor = $self->get_monitor();
####	
####	#### CLUSTER
####	my $cluster = $self->get_cluster();
####	my $jobs 	= $self->get_jobs();
####	my $qsub 	= $self->get_qsub();
####	my $sleep	= $self->get_sleep();
####	my $qstat 	= $self->get_qstat();
####	my $queue 	= $self->get_queue();	
####
####
####	#### SET REFERENCE BINARY FOR LATER
####	my $referencebinary = $referencefile;
####	$referencebinary =~ s/\.[^\.]+?$/.bfa/;
####	
####	#### DO THE MAQ match ALIGNMENT TO PRODUCE THE out.map FILES
####	my $align_commands = [];
####	for ( my $index = 0; $index < @$splitfiles; $index++ )
####	{
####		#### SET ALIGNMENT COMMANDS
####		my ($outputdir) = $splitfiles->[$index][0] =~ /^(.+?)\/[^\/]+$/;
####		my $subcommands = $self->alignmentCommands($$splitfiles[$index], $outputdir, $referencefile, $referencebinary);
####		die "No commands returned for splitfiles: $splitfiles->[$index][0], $splitfiles->[$index][1]\n" if not defined $subcommands;
####		push @$align_commands, $subcommands;
####	}
####	
####
####	#### RUN maq match ON EVERY FILE
####	my $counter = 0;
####	my $jobids;
####	my $scriptfiles;
####	foreach my $commands ( @$align_commands )
####	{	
####		$counter++;
####		#print "MAQ::doAlignments    DOING alignment $counter\n" if $counter % $dot == 0;
####		
####		#### CREATE DIRECTORIES
####		if ( not -d "$outputdir/$counter" )
####		{
####			File::Path::mkpath("$outputdir/$counter");
####		}
####	
####		#### MOVE TO OUTPUT SUBDIR
####		chdir("$outputdir/$counter") or die "Can't move to output subdir: $outputdir/$counter\n";
####	
####		#### SET LABEL
####		my $label = "$reference-$counter";
####		#print "MAQ::doAlignments    label: $label\n";	
####		
####		#### PRINT SHELL SCRIPT	
####		my $scriptfile = "$outputdir/$counter/$label.sh";
####		my $usagefile = "$outputdir/$counter/$label-usage.txt";
####		my $stdoutfile = "$outputdir/$counter/$label-stdout.txt";
####		my $stderrfile = "$outputdir/$counter/$label-stderr.txt";
####		$self->printScriptfile($scriptfile, $commands, $label, $usagefile, $stdoutfile, $stderrfile);
####	
####		#### USE THE RETURNED CLUSTER JOB ID AS THE PID
####		#### (submitJob INHERITED FROM Monitor::PBS)
####		my $jobid = $monitor->submitJob(
####			{
####				scriptfile  => $scriptfile,
####				queue       => $queue,
####				qsub       => $qsub,
####				stdoutfile  => $stdoutfile,
####				stderrfile  => $stderrfile
####			}
####		);
####
####
####
####		#### SAVE PID FOR CHECKING 
####		push @$jobids, $jobid;
####		
####		#### CHECK TO MAKE SURE WE HAVEN'T REACHED THE LIMIT OF MAX CONCURRENT JOBS
####		while ( scalar(@$jobids) >= $jobs )
####		{
####			sleep($sleep);
####			$jobids = $monitor->remainingJobs($jobids, $qstat);   
####		}
####	
####		#### CLEAN UP
####		#`rm -fr $scriptfile`;
####		push @$scriptfiles, $scriptfile;
####	}
####	
####	#### WAIT TIL ALL JOBS ARE FINISHED
####	while ( defined $jobids and scalar(@$jobids) > 0 )
####	{
####		sleep($sleep);
####		$jobids = $self->get_monitor()->remainingJobs($jobids);   
####	}
####	
####	##### CLEAN UP CLUSTER STDOUT AND STDERR FILES
####	#print "MAQ::doAlignments    Cleaning up scriptfiles...\n";
####	#foreach my $scriptfile ( @$scriptfiles )
####	#{
####	#	print "MAQ::doAlignments    Removing scriptfile: $scriptfile\n";
####	#	`rm -fr $scriptfile*`;
####	#}
####
####	print "MAQ::doAlignments    END OF MAQ::doAlignments\n";
####}
####
####
####
####


