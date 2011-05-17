package Aligner::NOVOALIGN;


=head2

		PACKAGE		Aligner::NOVOALIGN

		VERSION		0.01

		PURPOSE

	        WRAPPER SCRIPT FOR RUNNING NOVOALIGN ALIGNMENT

		HISTORY
					0.01 BASIC VERSION
=cut

use strict;
use warnings;
use Carp;

require Exporter;
our @ISA = qw(Exporter Cluster);
#our @EXPORT_OK = qw();
our $AUTOLOAD;

#### INTERNAL MODULES
use Cluster;
use Util;
use Sampler;

#### EXTERNAL MODULES
use Data::Dumper;
use File::Path;
use File::Remove;
use FindBin qw($Bin);
#use MD5;

#### SET SLOTS
our @DATA = qw(

THRESHOLD
MIN
MAX
PAIRED
REPLICATES
INPUTFILES
MATEFILES
DISTANCE
DEVIATION
PARAMS
OUTPUTDIR
SPLITFILE
MAXLINES
CHUNKS
CLEAN
KEEP

REFERENCEDIR
SPECIES
LABEL

NOVOALIGN
SAMTOOLS
SAMTOOLSINDEX

CLUSTER
QUEUE
WALLTIME
QSTAT
MAXJOBS
CPUS
SLEEP
CLEANUP
VERBOSE
TEMPDIR
DOT
QSUB

BATCHSTATS
COMMAND
);
our $DATAHASH;
foreach my $key ( @DATA )	{	$DATAHASH->{lc($key)} = 1;	}


=head2

	SUBROUTINE		run

	PURPOSE

		CREATE .sh FILE

=cut

sub run {
	my $self		=	shift;	



	#### GET CLUSTER
	my $cluster = $self->get_cluster();

	#### FILES AND DIRS	
	my $referencedir 	= 	$self->get_referencedir();
	my $outputdir 		=	$self->get_outputdir();
	my $inputfiles 		=	$self->get_inputfiles();
	my $matefiles 		=	$self->get_matefiles();

	#### GET LABEL, SPLITFILE, CHUNKS
	my $label 			= $self->get_label();
	my $splitfile 		= $self->get_splitfile();
	my $chunks 			= $self->get_chunks();

	#### SET DEFAULT SPLITFILE IF NOT DEFINED
	$splitfile = "$outputdir/splitfile.txt" if not defined $splitfile;


#### INPUTS

	###############################################
	###########    GET REFERENCEFILES    ##########
	###############################################
	print "Aligner::NOVOALIGN::run    Doing listReferenceFiles()  ", Timer::current_datetime(), "\n";
	my $referencefiles = $self->listReferenceFiles($referencedir);
	print "Aligner::NOVOALIGN::run    After listReferenceFiles()  ", Timer::current_datetime(), "\n";
	print "Aligner::NOVOALIGN::::run    No. referencefiles: ", scalar(@$referencefiles), "\n";

	###############################################
	###########   SET REFERENCE NAMES    ##########
	###############################################
	my $references = [];
	foreach my $referencefile ( @$referencefiles )
	{
		my ($reference) = $referencefile =~ /^.+?\/([^\/]+)$/;
		$reference =~ s/\.idx$//;
		push @$references, $reference if defined $reference;
		print "Reference not defined for referencefile: $referencefile\n" if not defined $reference;
	}

	##############################################
	############   SPLIT INPUT FILES   ###########
	##############################################
	my $splitfiles = $self->doSplitfiles($splitfile, $label);

	##############################################
	############	   SET CHUNKS      ###########
	##############################################
	$splitfiles = $self->splitfileChunks($splitfiles, $chunks) if defined $chunks;

	###############################################
	#############  CONVERT REFERENCES   ###########
	###############################################
	#$self->convertReferences($referencedir, $referencedir);


#### ALIGNMENT

	##############################################
	##########       RUN Aligner::NOVOALIGN         ##########
	##############################################
	print "Aligner::NOVOALIGN::run    Doing doAlignment()   ", Timer::current_datetime(), "\n";
	$self->doBatchAlignment($outputdir, $referencefiles, $splitfiles, $label);
	print "Aligner::NOVOALIGN::run    After doAlignment()   ", Timer::current_datetime(), "\n";


#### CONVERT TO CHROMOSOME SAM 

	###############################################
	##########    CHROMOSOME SAM FILES    #########
	###############################################
	#### OBJECTIVE: CREATE ONE .sam FILE PER CHROMOSOME
	#### FROM THE SMALL .sam FILES CREATED FOR EACH
	#### INPUT SPLIT FILE
	####
	#### STRATEGY 1: USE BAM FILE INTERMEDIARY
	####	
	#### 	1. CONVERT chr*/<NUMBER> SUBDIR SAM FILES INTO BAM FILES
	#### 	2. MERGE BAM FILES INTO A SINGLE BAM FILE AND SORT IT
	#### 	3. CONVERT SINGLE BAM FILE INTO A SAM FILE
	####
	#### STRATEGY 2: SORT MERGED SAM FILE
	####	
	#### 	1. MERGE SAM FILES INTO SINGLE SAM FILE
	#### 	2. SORT SAM FILES WITH sort
	#### 	

	################################################
	#####  ******      STRATEGY 1       ******  ####
	#####         BAM FILE INTERMEDIARY         ####
	################################################	
	#
	################################################
	###########    CONVERT SAM TO BAM     ##########
	################################################
	#$self->subdirSamToBam($outputdir, $references, $splitfiles, "accepted_hits.sam");
	#
	###############################################
	##########   CUMULATIVE MERGE BAM    ##########
	###############################################
	#$self->cumulativeMergeBam($outputdir, $references, $splitfiles, "accepted_hits.bam", "accepted_hits.bam");

	###############################################
	##########     PYRAMID MERGE BAM     ##########
	###############################################
	#$self->pyramidMergeBam($outputdir, $references, $splitfiles, "accepted_hits.bam", "out.bam");
	#
	################################################
	###########         SORT BAM         ##########
	################################################
	#$self->sortBam($outputdir, $references, "out.bam", "out.sam");
	#
	#
	################################################
	###########         BAM TO SAM        ##########
	################################################
	#$self->bamToSam($outputdir, $references);
	#


		###############################################
		####  ******      STRATEGY 2       ******  ####
		####         SORT MERGED SAM FILE          ####
		###############################################	
		##
		##################################################
		#############         MERGE SAM         ##########
		##################################################
		###print "Aligner::NOVOALIGN::run    Doing mergeSam     ", Timer::current_datetime(), "\n";
		###$self->pyramidMergeSam($outputdir, $references, $splitfiles, "accepted_hits.sam", "out.sam");
		###print "Aligner::NOVOALIGN::run    After mergeSam     ", Timer::current_datetime(), "\n";

		################################################
		###########         MERGE SAM         ##########
		################################################
		#$self->cumulativeMergeSam($outputdir, $references, $splitfiles, "accepted_hits.sam", "out.sam");

		################################################
		###########         SORT SAM         ###########
		################################################
		#$self->sortSam($outputdir, $references, "out.sam");

		###############################################
		##########       SAM TO BAM         ###########
		###############################################
		#$self->samToBam($outputdir, $references, "out.sam", "out.bam");
		#


	###############################################
	######        FILTER SAM HITS          ######
	###############################################
	print "Aligner::NOVOALIGN::runSubdirSamHits    Doing subdirSamHits        ", Timer::current_datetime(), "\n";
	$self->subdirSamHits($outputdir, $references, $splitfiles, "out.sam", "hit.sam", "miss.sam");
	print "Aligner::NOVOALIGN::runSubdirSamHits    After subdirSamHits        ", Timer::current_datetime(), "\n";


	################################################
	######        CUMULATIVE MERGE SAM        ######
	################################################
	print "Aligner::NOVOALIGN::run    Doing cumulativeMergeSam        ", Timer::current_datetime(), "\n";
	$self->cumulativeMergeSam($outputdir, $references, $splitfiles, "hit.sam", "hit.sam");
	print "Aligner::NOVOALIGN::run    After cumulativeMergeSam        ", Timer::current_datetime(), "\n";

	##############################################
	#########         SORT SAM         ###########
	##############################################
	print "Aligner::NOVOALIGN::run    Doing samToBam     ", Timer::current_datetime(), "\n";
	$self->samToBam($outputdir, $references, "hit.sam", "hit.bam");
	print "Aligner::NOVOALIGN::run    After samToBam     ", Timer::current_datetime(), "\n";



#### READ HIT STATS


	################################################
	############        READ HITS         ##########
	################################################
	#$self->readHits($outputdir, $references, $splitfiles, "eland-readHits", "hit.sam");


#### SNP CALLING


	###############################################
	#####     PREDICT SNPs WITH SAMTOOLS     ######
	###############################################
	#my $samtools_index = $self->get_samtoolsindex();
	#$self->samtoolSnps($splitfiles, $references, $outputdir, "hit.sam");


	################################################
	##########          USAGE STATS      ###########
	################################################
	#$self->printUsage("$outputdir/novoalign-USAGE.txt");

}	#	run




=head2

	SUBROUTINE		convertReferences

	PURPOSE

		CONVERT .fa REFERENCE FASTA FILES INTO .idx NOVOALIGN INDEX-FORMAT FILES

=cut

sub convertReferences {
	my $self			=	shift;
	my $inputdir		=	shift;
	my $outputdir		=	shift;


	$outputdir = $inputdir if not defined $outputdir;

	#### SANITY CHECK
	print "Aligner::NOVOALIGN::convertReferences    inputdir not defined: $inputdir\n" and exit if not defined $inputdir;
	print "Aligner::NOVOALIGN::convertReferences    outputdir not defined: $outputdir\n" and exit if not defined $outputdir;
	print "Aligner::NOVOALIGN::convertReferences    inputdir is a file: $inputdir\n" and exit if -f $inputdir;
	print "Aligner::NOVOALIGN::convertReferences    outputdir is a file: $outputdir\n" and exit if -f $outputdir;
	print "Aligner::NOVOALIGN::convertReferences    Can't find inputdir: $inputdir\n" and exit if not -d $inputdir;
	print "Aligner::NOVOALIGN::convertReferences    Can't find outputdir: $outputdir\n" and exit if not -d $outputdir;

	#### GET REFERENCE FILES
	my $fastafiles = $self->listFiles($inputdir, "\*.fa");

	#### GET NOVOALIGN
	my $novoalign		=	$self->get_novoalign();
	print "MAQ::convertReference    novoalign not defined\n" and exit if not defined $novoalign or not $novoalign;


	###############################################
	############# MAKE BINARY REF FILES ###########
	###############################################
	my $current_time = time();	
	my $jobs = [];
	for my $fastafile ( @$fastafiles )
	{
		next if $fastafile =~ /^[\.]+$/;

		my ($reference) = $fastafile =~ /^.+?\/([^\/]+)$/;

		#### SET REFERENCE BINARY FILE
		my $referencebinary = "$outputdir/$reference";
		$referencebinary =~ s/\.fa$/.idx/;

		#### SKIP IF ALREADY EXISTS
		print "Aligner::NOVOALIGN::convertReference    Skipping conversion because binary reference file already exists: $referencebinary\n" and next if -f $referencebinary and not -z $referencebinary;


		#### CONVERT REFERENCE FILE INTO BINARY FILE
		my $command;
		$command = "time $novoalign/novoindex -k 14 -s 2 $referencebinary $fastafile" if not -f $referencebinary or -z $referencebinary;
		print "MAQ::convertReference    command: $command\n" if defined $command;
		next if not defined $command;

		my $job = $self->setJob( [$command], "convertReference", $outputdir );
		push @$jobs, $job;
	}

	#### SKIP IF NO JOBS
	return if scalar(@$jobs) == 0;

	$self->runJobs( $jobs, "convertReference" );	

	#### GET DURATION OF LOCAL JOB AND ADD AS USAGE STATISTIC
	my $duration = Timer::runtime( $current_time, time() );
	$self->addUsageStatistic("convertReference", $duration);
}










=head2

	SUBROUTINE		doBatchAlignment

	PURPOSE

		1. RUN NOVOALIGN AGAINST ALL REFERENCE FILES

		2. DO ALL REFERENCES IN PARALLEL AS A BATCH JOBS

=cut

sub doBatchAlignment {
	my $self			=	shift;
	my $outputdir		=	shift;
	my $referencefiles 	=	shift;
	my $splitfiles		=	shift;
	my $label			=	shift;	


	#### COLLECT ALL JOBS FOR EACH INPUT FILE AGAINST ALL REFERENCE FILES
	my $jobs = $self->generateBatchJobs($outputdir, $referencefiles, $splitfiles, $label);

	##### RUN ALIGNMENT JOBS
	$self->runJobs($jobs, "Aligner::NOVOALIGN");
}

=head2

	SUBROUTINE		batchCommand

	PURPOSE

		CREATE .sh FILE

=cut

sub batchCommand {
	my $self			=	shift;
	my $outputdir		=	shift;
	my $referencepath	=	shift;
	my $splitfiles		=	shift;



	#### USER INPUTS
	my $threshold		= 	$self->get_threshold();
	my $paired 			= 	$self->get_paired();
	my $distance 		=	$self->get_distance();
	my $deviation 		=	$self->get_deviation();
	my $params 			=	$self->get_params();  #### OVERRIDE PARAMS IN DEFAULT COMMAND
	my $label 			=	$self->get_label();	#### USED TO GENERATE INPUTFILE NAMES	
	my $keep 			=	$self->get_keep();

	#### EXECUTABLES
	my $novoalign 		=	$self->get_novoalign();

	##### CLUSTER AND CPUS
	my $cluster 		=	$self->get_cluster();
	my $cpus 			=	$self->get_cpus();

	print "Aligner::NOVOALIGN::batchCommand    distance not defined. Exiting\n" if not defined $distance and $paired;
	print "Aligner::NOVOALIGN::batchCommand    deviation not defined. Exiting\n" if not defined $deviation and $paired;
	print "Aligner::NOVOALIGN::batchCommand    label not defined. Exiting\n" if not defined $label;


	#### SET PAIRED USING MATEFILES, I.E., DURING NORMAL RUN, NOT CALLED BY check
	$paired = 1 if defined $self->get_matefiles();

	#### GET THE BASE DIRECTORY OF THE SPLIT FILES - ALWAYS TWO DIRECTORIES DOWN
	my $splitfile = $$splitfiles[0][0];
	my ($basedir) = $splitfile =~ /^(.+?)\/\d+\/([^\/]+)$/;

	#### GET SUFFIX OF SPLIT FILE IF EXISTS
	my ($suffix) = $self->fileSuffix($splitfile);
	$suffix = '' if not defined $suffix;

	#### SET INDEX PATTERN FOR BATCH JOB
	my $index;
	$index = "\$LSB_JOBINDEX" if $cluster eq "LSF";
	$index = "\$PBS_TASKNUM" if $cluster eq "PBS";

	#### SET OUTPUT FILES
	my ($reference) = $referencepath =~ /([^\/]+)$/;
	$reference =~ s/\.idx//;

	my $outdir = "$outputdir/$index";
	my $outputfile = "$outdir/out.sam";
	my $alignedfile = "$outdir/aligned.txt";
	my $unalignedfile = "$outdir/unaligned.txt";

	#### CREATE OUTPUT DIR IF NOT EXISTS
	File::Path::mkpath($outdir) if not -d $outdir;
	print "Can't create outdir: $outdir\n" if not -d $outdir;

	#### SET INPUT AND MATE FILES, E.G.:
	#### /scratch/syoung/base/pipeline/bixby/run1/ln/$LSB_JOBINDEX/ln_1.$LSB_JOBINDEX.txt
	#### /scratch/syoung/base/pipeline/bixby/run1/ln/$LSB_JOBINDEX/ln_2.$LSB_JOBINDEX.txt	
	my $firstmate = $label . "_1";
	my $inputfile = "$basedir/$index/$firstmate.$index$suffix";

	#### DO SECOND MATE IF matefiles DEFINED
	my $secondmate = $label . "_2" if $paired;
	my $matefile = "$basedir/$index/$secondmate.$index$suffix" if $paired;

	my $command; 

	#### CHECK TEMPDIR EXISTS AND IS NOT A FILE
	my $tempdir = $self->get_tempdir();
	if ( defined $tempdir and $tempdir and not -d $tempdir )
	{
		print "Aligner::NOVOALIGN::batchCommand    tempdir directory not found: $tempdir\n" if not -d $tempdir;
		print "Aligner::NOVOALIGN::batchCommand    tempdir is a file: $tempdir\n" if -f $tempdir;
	}

	#### SET TEMP-RELATED DIRS
	my $old_outdir	=	$outdir;
	my $temp_outdir 	=	"$tempdir/$outdir" if defined $tempdir;

	#### USE TEMPDIR IF DEFINED
	if ( defined $tempdir and $tempdir and -d $tempdir )
	{

		#### SET TEMP OUTPUT DIR
		$outdir = $temp_outdir;
	}

	#### SET NOVOALIGN COMMAND
	my $novoalign_command = qq{time $novoalign/novoalign \\
-o SAM \\\n};

	#### ADD PARAMS
	$novoalign_command .= qq{$params \\\n} if defined $params;

	#### ADD DISTANCE & DEVIATION IF PAIRED
	$novoalign_command .= qq{-i $distance $deviation \\\n} if $paired;

	#### ADD REFERENCE FILE
	$novoalign_command .= qq{-d $referencepath \\\n};

	#### ADD INPUT FILES
	$novoalign_command .= qq{-f $inputfile $threshold \\\n} if not $paired;
	$novoalign_command .= qq{-f $inputfile $matefile $threshold\\\n} if $paired;

	#### OUTPUT FILE
	$novoalign_command =~ s/\\+\n$//g;
	$novoalign_command .= qq{ > $outputfile\n};

print "Aligner::NOVOALIGN    $novoalign_command\n";
exit;


	#### SET SCRIPT
	$command = qq{

export PATH=$novoalign:\$PATH
export PATH=$outdir:\$PATH
mkdir -p $outdir
cd $outdir

$novoalign_command

};

	if ( defined $tempdir )
	{
		$command .= qq{
mv $outdir/* $old_outdir
rm -fr $outdir

};
	}

	return $command;

}	#	batchCommand







=head2

	SUBROUTINE		listReferenceFiles

	PURPOSE

		RETRIEVE A LIST OF REFERENCE FILES FROM THE REFERENCE DIRECTORY

=cut

sub listReferenceFiles {
	my $self		=	shift;
	my $reference	=	shift;


	my $referencefiles = $self->listFiles($reference, "*\.idx");
	print "Novoalign::listReferenceFiles    No reference files in directory: $reference\n" and exit if not defined $referencefiles or scalar(@$referencefiles) == 0;

	#### SORT BY NUMBER
	@$referencefiles = Util::sort_naturally(\@$referencefiles);

	#### DEBUG
	@$referencefiles = reverse @$referencefiles;

	return $referencefiles;
}


=head2

	SUBROUTINE		getReferenceFiles

	PURPOSE

		RETURN A LIST OF REFERENCE FILES (FULL PATHS TO FILES)

=cut

sub getReferenceFiles {
	my $self		=	shift;
	my $referencedir=	shift;

	return $self->listReferenceFiles($referencedir, "\*\.idx");
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
	foreach my $reference ( @$references )	{	$reference =~ s/\.idx$//;	}

	return $references;
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

	#$self->SUPER::new();	

	#### INITIALISE THE OBJECT'S ELEMENTS
	$self->initialise($arguments);

	#### REQUIRED: CALL THE PARENT FOR ADDITIONAL INITIALISATION
	$self->SUPER::initialise();	

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

