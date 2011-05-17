package Expression::TOPHAT;


=head2

		PACKAGE		Expression::TOPHAT

		VERSION		0.02

		PURPOSE

	        WRAPPER SCRIPT FOR RUNNING Expression::TOPHAT SNP PREDICTION

		HISTORY
					0.02 ADDED CHUNK-BY-CHROMOSOME IF referencedir SPECIFIED
						AND RUN CUFFLINKS COMMAND
					0.01 BASIC VERSION
=cut

use strict;
use warnings;
use Carp;

require Exporter;
our @ISA = qw(Exporter Agua::Cluster);
#our @EXPORT_OK = qw();
our $AUTOLOAD;

#### USE LIBRARY
use FindBin qw($Bin);
use lib "$Bin/../../lib";	
use lib "$Bin/../../lib/external";	

#### INTERNAL MODULES
use Cluster;
#use Util;
#use Sampler;

#### EXTERNAL MODULES
use Data::Dumper;
use File::Path;
use File::Remove;
use FindBin qw($Bin);
#use MD5;

#### SET SLOTS
our @DATA = qw(

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
BATCHSTATS

SPECIES
PARAMS
LABEL
GTF
KEEP

TOPHAT
BOWTIE
CUFFLINKS
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

COMMAND
);
our $DATAHASH;
foreach my $key ( @DATA )	{	$DATAHASH->{lc($key)} = 1;	}


=head2

	SUBROUTINE		run

	PURPOSE

		CREATE .sh FILE

=cut

sub run
{
	my $self		=	shift;	


	print "Expression::TOPHAT::run    Expression::TOPHAT::run()\n";

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
	print "Expression::TOPHAT::run    Doing getReferencefiles()  ", Timer::current_datetime(), "\n";
	my $referencefiles = $self->getReferencefiles($referencedir);
	print "Expression::TOPHAT::run    After getReferencefiles()  ", Timer::current_datetime(), "\n";
	print "Expression::TOPHAT::::run    No. referencefiles: ", scalar(@$referencefiles), "\n";

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
	print "Expression::TOPHAT::::run    references: @$references\n";

	##############################################
	############	SPLIT INPUT FILES   ###########
	##############################################
	my $splitfiles = $self->doSplitfiles($splitfile, $label);

	##############################################
	############	   SET CHUNKS      ###########
	##############################################
	$splitfiles = $self->splitfileChunks($splitfiles, $chunks) if defined $chunks;

	###############################################
	###########       RUN Expression::TOPHAT         ##########
	###############################################
	print "Expression::TOPHAT::run    Doing runTophat()   ", Timer::current_datetime(), "\n";
	$self->runTophat($outputdir, $referencefiles, $splitfiles);
	print "Expression::TOPHAT::run    After runTophat()   ", Timer::current_datetime(), "\n";

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
	print "Expression::TOPHAT::run    Doing STRATEGY 1: USE BAM FILE INTERMEDIARY\n";

		###############################################
		##########    CONVERT SAM TO BAM     ##########
		###############################################
		print "Expression::TOPHAT::run    Doing subdirSamToBam     ", Timer::current_datetime(), "\n";
		$self->subdirSamToBam($outputdir, $references, $splitfiles, "accepted_hits.sam");
		print "Expression::TOPHAT::run    After subdirSamToBam     ", Timer::current_datetime(), "\n";

		####################################################
		###############   CUMULATIVE MERGE BAM    ##########
		####################################################
		######print "Expression::TOPHAT::run    Doing cumulativeMergeBam     ", Timer::current_datetime(), "\n";
		######$self->cumulativeMergeBam($outputdir, $references, $splitfiles, "accepted_hits.bam", "accepted_hits.bam");
		######print "Expression::TOPHAT::run    After cumulativeMergeBam     ", Timer::current_datetime(), "\n";

		##############################################
		#########     PYRAMID MERGE BAM     ##########
		##############################################
		print "Expression::TOPHAT::run    Doing cumulativeMergeBam     ", Timer::current_datetime(), "\n";
		$self->pyramidMergeBam($outputdir, $references, $splitfiles, "accepted_hits.bam", "out.bam");
		print "Expression::TOPHAT::run    After cumulativeMergeBam     ", Timer::current_datetime(), "\n";


		###############################################
		##########         SORT BAM         ##########
		###############################################
		print "Expression::TOPHAT::run    Doing sortBam     ", Timer::current_datetime(), "\n";
		$self->sortBam($outputdir, $references, "out.bam", "out.bam");
		print "Expression::TOPHAT::run    After sortBam     ", Timer::current_datetime(), "\n";


		##############################################
		#########         BAM TO SAM        ##########
		##############################################
		print "Expression::TOPHAT::run    Doing bamToSam     ", Timer::current_datetime(), "\n";
		$self->bamToSam($outputdir, $references, "out.bam", "out.sam");
		print "Expression::TOPHAT::run    After bamToSam     ", Timer::current_datetime(), "\n";

		print "Expression::TOPHAT::run    Completed STRATEGY 1: USE BAM FILE INTERMEDIARY\n";


	###################################################
	########  ******      STRATEGY 2       ******  ####
	########         SORT MERGED SAM FILE          ####
	###################################################	
	#####print "Expression::TOPHAT::run    Doing STRATEGY 2: SORT MERGED SAM FILE\n";
	####	
	####	################################################
	####	###########   CUMULATIVE MERGE SAM    ##########
	####	################################################
	####	#print "Expression::TOPHAT::run    Doing mergeSam     ", Timer::current_datetime(), "\n";
	####	#$self->cumulativeMergeSam($outputdir, $references, $splitfiles, "accepted_hits.sam", "out.sam");
	####	#print "Expression::TOPHAT::run    After mergeSam     ", Timer::current_datetime(), "\n";
	####
	####	################################################
	####	###########         SORT SAM         ###########
	####	################################################
	####	#print "Expression::TOPHAT::run    Doing sortSam     ", Timer::current_datetime(), "\n";
	####	#$self->sortSam($outputdir, $references, "out.sam");
	####	#print "Expression::TOPHAT::run    After sortSam     ", Timer::current_datetime(), "\n";
	####
	####	###############################################
	####	##########         SORT SAM         ###########
	####	###############################################
	####	#print "Expression::TOPHAT::run    Doing samToBam     ", Timer::current_datetime(), "\n";
	####	#$self->samToBam($outputdir, $references, "out.sam", "out.bam");
	####	#print "Expression::TOPHAT::run    After samToBam     ", Timer::current_datetime(), "\n";
	####	#
	####	#print "Expression::TOPHAT::run    Completed STRATEGY 2: SORT MERGED SAM FILE\n";
	####	#

	################################################
	###########         CUFFLINKS         ##########
	################################################
	print "Expression::TOPHAT::run    Doing runCufflinks     ", Timer::current_datetime(), "\n";
	$self->runCufflinks($outputdir, $references, "out.sam");
	print "Expression::TOPHAT::run    After runCufflinks     ", Timer::current_datetime(), "\n";

	###############################################
	#########          USAGE STATS      ###########
	###############################################
	print "Expression::TOPHAT::run    Doing printUsage       ", Timer::current_datetime(), "\n";
	$self->printUsage("$outputdir/tophat-USAGE.txt");
	print "Expression::TOPHAT::run    After printUsage       ", Timer::current_datetime(), "\n";

}	#	run




=head2

	SUBROUTINE		runTophat

	PURPOSE

		RUN Expression::TOPHAT AGAINST ALL REFERENCE FILES

=cut

sub runTophat
{
	my $self			=	shift;
	my $outputdir		=	shift;
	my $referencefiles 	=	shift;
	my $splitfiles		=	shift;	


	print "Expression::TOPHAT::runTophat    Expression::TOPHAT::runTophat(outputfile, referencefile, splitfiles)\n";


	#### COLLECT ALL JOBS FOR EACH INPUT FILE AGAINST ALL REFERENCE FILES
	my $jobs = [];
	my $number_splitfiles = scalar(@$splitfiles);

	my $reference_counter = 0;
	foreach my $referencefile ( @$referencefiles )
	{
		$reference_counter++;

		#### CREATE A BATCH JOB
		my ($reference) = $referencefile =~ /^.+?\/([^\/]+)$/;

		my $batch_command = $self->tophatBatchCommand("$outputdir/$reference", $referencefile, $splitfiles);

		#### SET LABEL AND OUTPUT DIRECTORY
		my $label = "tophatBatchAlignment-$reference";
		my $outdir = "$outputdir/$reference";

		#### SET *** BATCH *** JOB 
		my $job = $self->setBatchJob([$batch_command], $label, $outdir, $number_splitfiles);
		push @$jobs, $job;
	}
	print "Tophat::runTophat    length(jobs): ", scalar(@$jobs), "\n";


	##### RUN ALIGNMENT JOBS
	$self->runJobs($jobs, "Expression::TOPHAT");
}
=head2

	SUBROUTINE		getReferencefiles

	PURPOSE

		RETRIEVE A LIST OF REFERENCE FILES FROM THE REFERENCE DIRECTORY

=cut

sub getReferencefiles
{
	my $self		=	shift;
	my $reference	=	shift;


	my $referencefiles = $self->listFiles($reference, "\*rev.1.ebwt");
	print "Tophat::getReferencefiles    No reference files in directory: $reference\n" and exit if not defined $referencefiles or scalar(@$referencefiles) == 0;

	#### TRUNCATE REFERENCE FILES TO CREATE CORRECT STUB IDENTIFIER
	foreach my $referencefile ( @$referencefiles )	{ $referencefile =~ s/\.rev\.1\.ebwt$//; }

	#### SORT BY NUMBER
	@$referencefiles = Util::sort_naturally(\@$referencefiles);

	#### DEBUG
	@$referencefiles = reverse @$referencefiles;

	return $referencefiles;
}


=head2

	SUBROUTINE		tophatBatchCommand

	PURPOSE

		CREATE .sh FILE

=cut

sub tophatBatchCommand
{
	my $self			=	shift;
	my $outputdir		=	shift;
	my $referencepath	=	shift;
	my $splitfiles		=	shift;


	#### USER INPUTS
	my $distance 		=	$self->get_distance();
	my $params 			=	$self->get_params();
	my $label 			=	$self->get_label();	#### USED TO GENERATE INPUTFILE NAMES

	my $keep 			=	$self->get_keep();

	#### EXECUTABLES
	my $tophat 			=	$self->get_tophat();
	my $bowtie 			=	$self->get_bowtie();

	##### CLUSTER
	my $cluster 		=	$self->get_cluster();
	my $cpus 			=	$self->get_cpus();

	#### GET THE BASE DIRECTORY OF THE SPLIT FILES - ALWAYS ONE DIRECTORY DOWN
	#### basedir/1/sequence.1.fastq
	my $splitfile = $$splitfiles[0][0];
	my ($basedir) = $splitfile =~ /^(.+?)\/\d+\/([^\/]+)$/;

	#### GET SUFFIX OF SPLIT FILE IF EXISTS
	my ($suffix) = $self->fileSuffix($splitfile);
	$suffix = '' if not defined $suffix;

	my $matefiles = $self->get_matefiles();

	#### SET INPUT AND MATE FILES, E.G.:
	#### /scratch/syoung/base/pipeline/bixby/run1/ln/$LSB_JOBINDEX/ln_1.$LSB_JOBINDEX.txt
	#### /scratch/syoung/base/pipeline/bixby/run1/ln/$LSB_JOBINDEX/ln_2.$LSB_JOBINDEX.txt	
	my $index;
	$index = "\$LSB_JOBINDEX" if $cluster eq "LSF";
	$index = "\$PBS_TASKNUM" if $cluster eq "PBS";
	$index = "\$SGE_TASK_ID" if $cluster eq "SGE";
	my $firstmate = $label . "_1";
	my $secondmate = $label . "_2";
	my $inputfile = "$basedir/$index/$firstmate.$index$suffix";
	my $matefile = "$basedir/$index/$secondmate.$index$suffix";

	#### SET OUTPUT DIR AS mainfolder/chromosome/\$LSB_INDEX
	my $outdir = "$outputdir/$index";

	#### HOLDER FOR COMPLETE COMMAND:
	#### Expression::TOPHAT COMMAND + ENVIRONMENT VARIABLES, ETC.
	my $command; 

	my $tempdir = $self->get_tempdir();
	if ( defined $tempdir and $tempdir and not -d $tempdir )
	{
		print "Expression::TOPHAT::tophatBatchCommand    tempdir directory not found: $tempdir\n" if not -d $tempdir;
		print "Expression::TOPHAT::tophatBatchCommand    tempdir is a file: $tempdir\n" if -f $tempdir;
	}

	if ( defined $tempdir and $tempdir and -d $tempdir )
	{

		#### SET TEMP OUTPUT DIR
		my $temp_outdir = $tempdir . $outdir;

		#### SET Expression::TOPHAT COMMAND
		my $tophat_command = qq{time $tophat/tophat \\
--num-threads $cpus \\\n};

		#### ADD PARAMS
		$tophat_command .= qq{$params \\\n} if defined $params;

		#### KEEP INTERMEDIATE FILES
		$tophat_command .= qq{--keep-tmp \\\n} if defined $keep;

		#### SPECIFY OUTPUT DIR, REFERENCE AND INPUT FILES
		$tophat_command .= qq{--output-dir $temp_outdir/\$LSB_JOBINDEX \\
--mate-inner-dist $distance \\
$referencepath \\
$inputfile };

		$tophat_command .= qq{\\\n$matefile } if defined $matefiles and $matefiles;


		#### SET SCRIPT
		$command = qq{
export PATH=$bowtie:\$PATH
export PATH=$tophat:\$PATH
export PATH=$referencepath/\\$index:\$PATH


cd $outdir

echo
echo "outdir: " $outdir
echo
echo "inputfile: " $inputfile

echo "mkdir -p $temp_outdir"
mkdir -p $temp_outdir

echo "$tophat_command"
$tophat_command

echo "mv $temp_outdir/* $outdir"
mv $temp_outdir/* $outdir

echo "rm -fr $temp_outdir"
rm -fr $temp_outdir

exit;
};

	}
	else
	{
		my $tophat_command = qq{time $tophat/tophat \\
--num-threads $cpus \\\n};

		#### ADD PARAMS
		$tophat_command .= qq{$params \\\n} if defined $params;

		#### KEEP INTERMEDIATE FILES
		$tophat_command .= qq{--keep-tmp \\\n} if defined $keep;

		#### OUTPUT DIR AND MATE PAIR DISTANCE
		$tophat_command .= qq{--output-dir $outdir \\\n};
		$tophat_command .= qq{--mate-inner-dist $distance \\\n} if defined $distance;

		#### REFERENCE AND INPUT FILES
$tophat_command .= qq{$referencepath \\
$inputfile };

		$tophat_command .= qq{\\\n$matefile } if defined $matefiles and $matefiles;

		$command = qq{
export PATH=$bowtie:\$PATH
export PATH=$tophat:\$PATH
export PATH=$referencepath/$index:\$PATH

cd $outdir

echo
echo "outdir: " $outdir
echo
echo "inputfile: " $inputfile
echo
echo "$tophat_command"
echo
$tophat_command

exit;
};

	}

	return $command;

}	#	tophatBatchCommand




=head2

	SUBROUTINE		tophatCommand

	PURPOSE

		CREATE .sh FILE

=cut


sub tophatCommand
{
	my $self			=	shift;
	my $inputfiles		=	shift;
	my $matefiles		=	shift;
	my $outputdir		=	shift;
	my $reference		=	shift;



	#### USER INPUTS
	my $distance 		=	$self->get_distance();
	my $params 			=	$self->get_params();
	my $label 			=	$self->get_label();
	my $keep 			=	$self->get_keep();

	#### EXECUTABLES
	my $tophat 			=	$self->get_tophat();
	my $bowtie 			=	$self->get_bowtie();

	##### CLUSTER
	#my $qstat 			=	$self->get_qstat();
	#my $jobs 			=	$self->get_jobs();
	my $cpus 			=	$self->get_cpus();
	#my $sleep 			=	$self->get_sleep();
	#my $qsub 			=	$self->get_qsub();
	#my $queue 			=	$self->get_queue();

	my $command = qq{
export PATH=$bowtie:\$PATH
export PATH=$tophat:\$PATH
export PATH=$outputdir:\$PATH

cd $outputdir

time $tophat/tophat \\
--num-threads $cpus \\\n};

	#### ADD PARAMS
	$command .= qq{$params \\\n} if defined $params;

	#### KEEP INTERMEDIATE FILES
	$command .= qq{--keep-tmp \\\n} if defined $keep;

	#### SPECIFY OUTPUT DIR, REFERENCE AND INPUT FILES
	$command .= qq{--output-dir $outputdir \\
--mate-inner-dist $distance \\
$reference \\
$inputfiles };

	$command .= qq{\\\n$matefiles } if defined $matefiles and $matefiles;

	return $command;

}	#	tophatCommand






=head2

	SUBROUTINE		runCufflinks

	PURPOSE

		RUN CUFFLINKS ON THE MERGED SAM FILE FOR EACH REFERENCE

		TO GET EXPRESSION CALCULATIONS.

=cut

sub runCufflinks
{
	my $self			=	shift;
	my $outputdir		=	shift;
	my $references	=	shift;
	my $infile			=	shift;

	print "Expression::TOPHAT::runCufflinks    Expression::TOPHAT::runCufflinks    (outputdir, references, infile)\n";


	my $jobs = [];
	foreach my $reference ( @$references )
	{
		#### SET INPUT FILE
		my $samfile = "$outputdir/$reference/$infile";

		#### SET COMMAND
		my $commands = [];
		my $command = $self->cufflinksCommand($samfile);
		push @$commands, $command;

		#### SET JOB		
		my $label = "cufflinks-$reference";
		my $outdir = "$outputdir/$reference";
		my $job = $self->setJob([$command], $label, $outdir);
		push @$jobs, $job;
	}

	#### RUN CUFFLINKS
	$self->runJobs($jobs, 'CUFFLINKS');
	print "FINISHED RUNNING CUFFLINKS\n";
}



=head2

	SUBROUTINE		cufflinksCommand

	PURPOSE

		CREATE .sh FILE

=cut


sub cufflinksCommand
{
	my $self			=	shift;
	my $samfile			=	shift;



	my ($outdir) = $samfile =~ /^(.+?)\/[^\/]+$/;

	#### USER INPUTS
	my $distance 		=	$self->get_distance();
	my $label 			=	$self->get_label();
	my $gtf 			=	$self->get_gtf();

	#### EXECUTABLES
	my $cufflinks 		=	$self->get_cufflinks();

	##### CLUSTER
	my $cpus 			=	$self->get_cpus();

	my $command = qq{
cd $outdir

time $cufflinks/cufflinks \\
--num-threads $cpus \\};

	$command .= qq{--inner-dist-mean $distance \\} if defined $distance;

	#### ADD GTF IF DEFINED
	$command .= qq{--GTF $gtf \\\n} if defined $gtf;

	##### ADD ADDITIONAL PARAMS
	#$command .= qq{$params \\\n} if defined $params;

	#### END COMMAND
	$command .= qq{$samfile\n};


	return $command;

}	#	cufflinksCommand



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

	print "Expression::TOPHAT::new    Expression::TOPHAT::new(arguments)\n";

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

		INITIALISE THE self OBJECT:

			1. LOAD THE DATABASE, USER AND PASSWORD FROM THE ARGUMENTS

			2. FILL OUT %VARIABLES% IN XML AND LOAD XML

			3. LOAD THE ARGUMENTS

=cut


sub initialise
{
    my $self		=	shift;
	my $arguments	=	shift;

print "Expression::TOPHAT::initialise    Expression::TOPHAT::initialise(arguments)\n";

	#### SET TIME
	$self->{_starttime} = time();

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

