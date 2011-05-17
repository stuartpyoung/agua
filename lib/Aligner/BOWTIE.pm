package BOWTIE;


=head2

		PACKAGE		BOWTIE

		VERSION		0.02

		PURPOSE

	        WRAPPER SCRIPT FOR RUNNING BOWTIE SNP PREDICTION

		HISTORY
					0.02 ADDED CHUNK-BY-CHROMOSOME IF referencedir SPECIFIED
						AND RUN CUFFLINKS COMMAND
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
use Sampler;
use Util;

#### EXTERNAL MODULES
use Data::Dumper;
use File::Path;
use File::Remove;
use FindBin qw($Bin);
#use MD5;

#### SET SLOTS
our @DATA = qw(

MIN
MAX
PAIRED
REPLICATES
INPUTFILES
MATEFILES
SEQUENCEDIR
DISTANCE
OUTPUTDIR
REFERENCEDIR
SPLITFILE
MAXLINES
CHUNKS
CLEAN

PARAMS
LABEL
GTF
KEEP

SPECIES
BOWTIE
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


	print "BOWTIE::run    BOWTIE::run()\n";

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
	my $chunks = $self->get_chunks();

	#### SET DEFAULT SPLITFILE IF NOT DEFINED
	$splitfile = "$outputdir/splitfile.txt" if not defined $splitfile;

	###############################################
	###########    GET REFERENCEFILES    ##########
	###############################################
	print "BOWTIE::run    Doing listReferenceFiles()  ", Timer::current_datetime(), "\n";
	my $referencefiles = $self->listReferenceFiles($referencedir);
	print "BOWTIE::run    After listReferenceFiles()  ", Timer::current_datetime(), "\n";
	print "BOWTIE::run    No. referencefiles: ", scalar(@$referencefiles), "\n";
	print "BOWTIE::run    referencefiles: \n";
	print join "\n", @$referencefiles;
	print "\n";
#exit;

	###############################################
	###########   SET REFERENCE NAMES    ##########
	###############################################
	my $references = [];
	foreach my $referencefile ( @$referencefiles )
	{
		my ($reference) = $referencefile =~ /^.+?\/([^\/]+)$/;
		push @$references, $reference if defined $reference;
		print "Reference not defined for referencefile: $referencefile\n" if not defined $reference;
	}
	print "BOWTIE::::run    references: @$references\n";

	##############################################
	############	SPLIT INPUT FILES   ###########
	##############################################
	my $splitfiles = $self->doSplitfiles($splitfile, $label);

	##############################################
	############	   SET CHUNKS      ###########
	##############################################
	$splitfiles = $self->splitfileChunks($splitfiles, $chunks) if defined $chunks;


#### DEBUG
#### DEBUG
#### DEBUG
#### DEBUG

	##############################################
	##########       RUN BOWTIE         ##########
	##############################################
	print "BOWTIE::run    Doing doAlignment()   ", Timer::current_datetime(), "\n";
	$self->doBatchAlignment($outputdir, $referencefiles, $splitfiles, $label);
	print "BOWTIE::run    After doAlignment()   ", Timer::current_datetime(), "\n";

#### DEBUG
#### DEBUG
#### DEBUG
#### DEBUG




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

	###############################################
	####  ******      STRATEGY 1       ******  ####
	####         BAM FILE INTERMEDIARY         ####
	###############################################	

	################################################
	###########    CONVERT SAM TO BAM     ##########
	################################################
	#$self->subdirSamToBam($outputdir, $references, $splitfiles, "accepted_hits.sam");

	################################################
	###########   CUMULATIVE MERGE BAM    ##########
	################################################
	##print "BOWTIE::run    Doing cumulativeMergeBam     ", Timer::current_datetime(), "\n";
	##$self->cumulativeMergeBam($outputdir, $references, $splitfiles, "accepted_hits.bam", "accepted_hits.bam");
	##print "BOWTIE::run    After cumulativeMergeBam     ", Timer::current_datetime(), "\n";
	#
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

		################################################
		###########         MERGE SAM         ##########
		################################################
		#$self->pyramidMergeSam($outputdir, $references, $splitfiles, "accepted_hits.sam", "out.sam");

		################################################
		###########         SORT SAM         ###########
		################################################
		#$self->sortSam($outputdir, $references, "hit.sam");

		###############################################
		##########         SORT SAM         ###########
		###############################################
		#$self->samToBam($outputdir, $references, "out.sam", "out.bam");
		#
		#


### DEBUG
	###############################################
	######        FILTER SAM HITS          ######
	###############################################
	print "BOWTIE::runSubdirSamHits    Doing subdirSamHits        ", Timer::current_datetime(), "\n";
	$self->subdirSamHits($outputdir, $references, $splitfiles, "out.sam", "hit.sam", "miss.sam");
	print "BOWTIE::runSubdirSamHits    After subdirSamHits        ", Timer::current_datetime(), "\n";

	################################################
	######        CUMULATIVE MERGE SAM        ######
	################################################
	print "BOWTIE::run    Doing cumulativeMergeSam        ", Timer::current_datetime(), "\n";
	$self->cumulativeMergeSam($outputdir, $references, $splitfiles, "hit.sam", "hit.sam");
	print "BOWTIE::run    After cumulativeMergeSam        ", Timer::current_datetime(), "\n";
### DEBUG


	##############################################
	#########         SORT SAM         ###########
	##############################################
	print "BOWTIE::run    Doing samToBam     ", Timer::current_datetime(), "\n";
	$self->samToBam($outputdir, $references, "hit.sam", "hit.bam");
	print "BOWTIE::run    After samToBam     ", Timer::current_datetime(), "\n";




#### READ HIT STATS


	################################################
	############        READ HITS         ##########
	################################################
	#$self->readHits($outputdir, $references, $splitfiles, "BOWTIE-readHits", "hit.sam");


#### SNP CALLING

	###############################################
	#####     PREDICT SNPs WITH SAMTOOLS     ######
	###############################################
	#my $samtools_index = $self->get_samtoolsindex();
	#$self->samtoolSnps($splitfiles, $references, $outputdir, "hit.sam");


	################################################
	##########          USAGE STATS      ###########
	################################################
	#$self->printUsage("$outputdir/bowtie-USAGE.txt");
}	#	run


=head2

	SUBROUTINE		getReferenceFiles

	PURPOSE

		RETURN A LIST OF REFERENCE FILES (FULL PATHS TO FILES)

=cut

sub getReferenceFiles {
	my $self		=	shift;
	my $referencedir=	shift;

	return $self->listReferenceFiles($referencedir, "\*\.rev\.1\.ebwt");
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
	foreach my $reference ( @$references )	{	$reference =~ s/\.rev\.1\.ebwt$//;	}

	return $references;
}




=head2

	SUBROUTINE		generateBatchJobs

	PURPOSE

		GENERATE LIST OF BATCH JOBS TO RUN ELAND AGAINST ALL REFERENCES

		USING SUBFILES OF INPUT FILES

=cut

sub generateBatchJobs {
	my $self			=	shift;
	my $outputdir		=	shift;
	my $referencefiles 	=	shift;
	my $splitfiles		=	shift;	


	#### COLLECT ALL JOBS FOR EACH INPUT FILE AGAINST ALL REFERENCE FILES
	my $jobs = [];
	my $number_splitfiles = scalar(@$splitfiles);

	my $reference_counter = 0;
	foreach my $referencefile ( @$referencefiles )
	{
		$reference_counter++;

		#### CREATE A BATCH JOB
		my ($reference) = $referencefile =~ /^.+?\/([^\/]+)$/;

		#### SET OUTPUT DIR
		my $outdir = "$outputdir/$reference";
		File::Path::mkpath($outdir) if not -d $outdir;

		my $command = $self->batchCommand($outdir, $referencefile, $splitfiles);

		#### SET LABEL
		my $label = "bowtie-$reference";

		#### SET *** BATCH *** JOB 
		my $job = $self->setBatchJob([$command], $label, $outdir, $number_splitfiles);

		push @$jobs, $job;
	}

	return $jobs;
}






=head2

	SUBROUTINE		doAlignment

	PURPOSE

		RUN BOWTIE AGAINST ALL REFERENCE FILES

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
	$self->runJobs($jobs, "BOWTIE");
}



=head2

	SUBROUTINE		listReferenceFiles

	PURPOSE

		RETRIEVE A LIST OF REFERENCE FILES FROM THE REFERENCE DIRECTORY

=cut

sub listReferenceFiles {
	my $self		=	shift;
	my $reference	=	shift;


	my $referencefiles = $self->listFiles($reference, "\*rev.1.ebwt");
	print "Bowtie::listReferenceFiles    No reference files in directory: $reference\n" and exit if not defined $referencefiles or scalar(@$referencefiles) == 0;

	#### TRUNCATE REFERENCE FILES TO CREATE CORRECT STUB IDENTIFIER
	foreach my $referencefile ( @$referencefiles )	{ $referencefile =~ s/\.rev\.1\.ebwt$//; }

	#### SORT BY NUMBER
	@$referencefiles = Util::sort_naturally(\@$referencefiles);

	#### DEBUG
	@$referencefiles = reverse @$referencefiles;

	return $referencefiles;
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
	my $paired 			= 	$self->get_paired();
	my $distance 		=	$self->get_distance();
	my $params 			=	$self->get_params();  #### OVERRIDE PARAMS IN DEFAULT COMMAND
	my $label 			=	$self->get_label();	#### USED TO GENERATE INPUTFILE NAMES	
	my $keep 			=	$self->get_keep();

	$paired = 1 if defined $self->get_matefiles();

	print "BOWTIE::batchCommand    distance not defined. Exiting\n" if not defined $distance;
	print "BOWTIE::batchCommand    label not defined. Exiting\n" if not defined $label;

	#### EXECUTABLES
	my $bowtie 			=	$self->get_bowtie();

	##### CLUSTER AND CPUS
	my $cluster 		=	$self->get_cluster();
	my $cpus 			=	$self->get_cpus();

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
	$index = "\$TASKNUM" if $cluster eq "SGE";

	#### SET OUTPUT FILES
	my ($reference) = $referencepath =~ /([^\/]+)$/;
	$reference =~ s/\.bfa//;
	my $outdir = "$outputdir/$index";

	#### SET TEMP OUTPUT DIR
	my $tempdir = $self->get_tempdir();
	if ( defined $tempdir and $tempdir )
	{
		print "BOWTIE::batchCommand    Cant' find tempdir: $tempdir\n"
			and exit if not -d $tempdir;
		print "BOWTIE::batchCommand    tempdir is a file: $tempdir\n"
			and exit if -f $tempdir;
	}
	$outdir = "$tempdir/$outdir" if defined $tempdir and $tempdir;

	my $outputfile = "$outdir/out.sam";
	my $alignedfile = "$outdir/aligned.txt";
	my $unalignedfile = "$outdir/unaligned.txt";

	#### SET INPUT AND MATE FILES, E.G.:
	#### /scratch/syoung/base/pipeline/bixby/run1/ln/$LSB_JOBINDEX/ln_1.$LSB_JOBINDEX.txt
	#### /scratch/syoung/base/pipeline/bixby/run1/ln/$LSB_JOBINDEX/ln_2.$LSB_JOBINDEX.txt	
	my $firstmate = $label . "_1";
	my $inputfile = "$basedir/$index/$firstmate.$index$suffix";

	#### DO SECOND MATE IF matefiles DEFINED
	my $secondmate = $label . "_2" if $paired;
	my $matefile = "$basedir/$index/$secondmate.$index$suffix" if $paired;

	my $command; 

	##### ADD OUTPUT DIR TO PATH 
	#$command = qq{\nexport PATH=$bowtie:$referencepath/\$LSB_JOBINDEX:\$PATH\n} if $cluster eq "LSF";
	#$command = qq{\nexport PATH=$bowtie:$referencepath/\$PBS_TASKNUM:\$PATH\n} if $cluster eq "PBS";
	#$command = qq{\nexport PATH=$bowtie:$referencepath/\$TASK_ID:\$PATH\n} if $cluster eq "SGE";

	$command .= qq{mkdir -p $outdir
cd $outdir\n};

	#### SET BOWTIE COMMAND
	my $bowtie_command = qq{time $bowtie/bowtie \\
--sam \\
--rf \\
--threads $cpus \\\n};
####--verbose \\

	#### DISTANCE
	$bowtie_command .= qq{-X $distance \\\n} if defined $distance;

	#### PARAMS
	$bowtie_command .= qq{$params \\\n} if defined $params;

	#### OUTPUT FILES
	$bowtie_command .= qq{--al $alignedfile \\
--un $unalignedfile \\\n};

	#### REFERENCE FILE
	$bowtie_command .= qq{$referencepath \\\n};

	#### INPUT FILES
	$bowtie_command .= qq{-1 $inputfile \\
-2 $matefile \\\n} if $paired;
	$bowtie_command .= qq{$inputfile \\\n} if not $paired;

	#### OUTPUT FILE
	$bowtie_command .= qq{$outputfile};


	#### SET BOWTIE COMMAND
	$command .= $bowtie_command;

	#### DO MOVE IF tempdir IS DEFINED
	$command .= qq{
mkdir -p $outputdir/$index
mv $outdir/* $outputdir/$index
rm -fr $outdir\n} if defined $tempdir and $tempdir;

	print "BOWTIE::batchCommand    command: $command\n";

	return $command;

}	#	batchCommand


=head2

	SUBROUTINE		indexReferences

	PURPOSE

		CONVERT .fa REFERENCE FILES INTO BOWTIE *ebwt BINARY REFERENCE FILES

=cut

sub indexReferences {
	my $self		=	shift;
	my $inputdir	=	shift;

	#### CHECK INPUTS
die "BOWTIE::indexReferences    inputdir not defined (Use --help for usage)\n" if not defined $inputdir;

	my $bowtie 	=	$self->get_bowtie();
	print "BOWTIE::indexReferences    bowtie is not defined. Exiting\n" if not defined $bowtie;

	chdir($inputdir) or die "BOWTIE::indexReferences    Can't change to inputdir directory: $inputdir\n";

	my @files = <*\.fa>;
	foreach my $file ( @files )
	{
		#### CREATE STUB IDENTIFIER
		my ($stub) = $file =~ /^(.+)\.fa$/;

		#### IGNORE IF INDEX FILES ALREADY EXISTS 
		my $indexfile = $stub . ".rev.1.ebwt";
		next if -f $indexfile and not -z $indexfile;

		my $command = "time $bowtie/bowtie-build  $file $stub";
		print "command: $command\n";
		`$command`;
	}
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

	print "BOWTIE::new    BOWTIE::new(arguments)\n";

	my $self = {};
    bless $self, $class;

	#$self->SUPER::new();	

	#### INITIALISE THE OBJECT'S ELEMENTS
	$self->initialise($arguments);

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

