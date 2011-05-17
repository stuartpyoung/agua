#use Getopt::Simple;
use MooseX::Declare;
#use Moose::Util::TypeConstraints;

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

#### EXTERNAL MODULES
use FindBin qw($Bin);

#### USE LIB
use lib "$Bin/..";

#### USES ROLES
use Agua::Cluster::Checker;
use Agua::Cluster::Jobs;

use strict;
use warnings;
use Carp;

our $VERSION = 0.01;

class Expression::TOPHAT with (Agua::Cluster::Jobs, Agua::Cluster::Checker, Agua::Cluster::Util, Agua::Cluster::Convert, Agua::Cluster::Merge, Agua::Cluster::Sort, Agua::Cluster::Usage) {

# BOOLEAN
has 'clean'		=> ( isa => 'Bool|Undef', is => 'rw', default => '' );

# INTS
has 'walltime'	=> ( isa => 'Int|Undef', is => 'rw', default => 24 );

# STRINGS
has 'referencedir'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'outputdir'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'replicates'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'label'			=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'splitfile'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'cluster'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );

has 'tophat'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'bowtie'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'cufflinks'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'gtf'			=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'keep'			=> ( isa => 'Str|Undef', is => 'rw', default => '' );

has 'check'			=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'inputfiles'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'matefiles'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'distance'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'params'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'sequencedir'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'rundir'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'tempdir'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'maxlines'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'chunks'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'keep'			=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'readhits'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'referencedir'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'species'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'subdirs'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'inputtype'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'pairparams'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'quality'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'seedlength'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'samtools'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'samtoolsindex'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'lane'			=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'firstlane'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'convert'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );

# OBJECTS
has 'splitfiles'=> ( isa => 'ArrayRef|Undef', is => 'rw' );
has 'conf' 	=> (
	is =>	'rw',
	'isa' => 'Conf::Agua',
	default	=>	sub { Conf::Agua->new(	backup	=>	1, separator => "\t"	);	}
);

#### INTERNAL MODULES
use Util;

#### EXTERNAL MODULES
use Data::Dumper;
use File::Path;
use File::Remove;

#### DEBUG

####/////}

sub BUILD {
	my $self	=	shift;
}
=head2

	SUBROUTINE		run

	PURPOSE

		CREATE .sh FILE

=cut

sub run {
	my $self		=	shift;	

	print "Expression::TOPHAT::run    Expression::TOPHAT::run()\n";
#exit;

	#### GET CLUSTER
	my $cluster = $self->cluster();

	#### FILES AND DIRS	
	my $referencedir 	= 	$self->referencedir();
	my $outputdir 		=	$self->outputdir();
	my $inputfiles 		=	$self->inputfiles();
	my $matefiles 		=	$self->matefiles();

	#### GET LABEL, SPLITFILE, CHUNKS
	my $label 			= $self->label();
	my $splitfile 		= $self->splitfile();
	my $chunks = $self->chunks();

	#### SET DEFAULT SPLITFILE IF NOT DEFINED
	$splitfile = "$outputdir/splitfile.txt" if not defined $splitfile;

	###############################################
	###########    GET REFERENCEFILES    ##########
	###############################################
	print "Expression::TOPHAT::run    Doing getReferencefiles()  ", Timer::current_datetime(), "\n";
	my $referencefiles = $self->getReferencefiles($referencedir);
	print "Expression::TOPHAT::run    After getReferencefiles()  ", Timer::current_datetime(), "\n";
	print "Expression::TOPHAT::run    No. referencefiles: ", scalar(@$referencefiles), "\n";

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
	print "Expression::TOPHAT::run    references: @$references\n";

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

}	#	run




=head2

	SUBROUTINE		runTophat

	PURPOSE

		RUN Expression::TOPHAT AGAINST ALL REFERENCE FILES

=cut

sub runTophat {
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

		my $batch_command = $self->batchCommand("$outputdir/$reference", $referencefile, $splitfiles);

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

sub getReferencefiles {
	my $self		=	shift;
	my $reference	=	shift;


	my $referencefiles = $self->listReferenceFiles($reference, "\*rev.1.ebwt");
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
	my $distance 		=	$self->distance();
	my $params 			=	$self->params();
	my $label 			=	$self->label();	#### USED TO GENERATE INPUTFILE NAMES
	my $keep 			=	$self->keep();

	#### EXECUTABLES
	my $tophat 			=	$self->tophat();
	my $bowtie 			=	$self->bowtie();

	##### CLUSTER
	my $cluster 		=	$self->cluster();
	my $cpus 			=	$self->cpus();

	#### GET THE BASE DIRECTORY OF THE SPLIT FILES - ALWAYS ONE DIRECTORY DOWN
	#### basedir/1/sequence.1.fastq
	my $splitfile = $$splitfiles[0][0];
	my ($basedir) = $splitfile =~ /^(.+?)\/\d+\/([^\/]+)$/;

	#### GET SUFFIX OF SPLIT FILE IF EXISTS
	my ($suffix) = $self->fileSuffix($splitfile);
	$suffix = '' if not defined $suffix;

	my $matefiles = $self->matefiles();

	#### SET INPUT AND MATE FILES, E.G.:
	#### /scratch/syoung/base/pipeline/bixby/run1/ln/$LSB_JOBINDEX/ln_1.$LSB_JOBINDEX.txt
	#### /scratch/syoung/base/pipeline/bixby/run1/ln/$LSB_JOBINDEX/ln_2.$LSB_JOBINDEX.txt	
	my $index = $self->getIndex();
	$index =~ s/^\\//;

	#$index = "\$LSB_JOBINDEX" if $cluster eq "LSF";
	#$index = "\$PBS_TASKNUM" if $cluster eq "PBS";
	#$index = "\$SGE_TASK_ID" if $cluster eq "SGE";
	my $firstmate = $label . "_1";
	my $secondmate = $label . "_2";
	my $inputfile = "$basedir/$index/$firstmate.$index$suffix";
	my $matefile = "$basedir/$index/$secondmate.$index$suffix";

	#### SET OUTPUT DIR AS mainfolder/chromosome/\$LSB_INDEX
	my $outdir = "$outputdir/$index";

	#### HOLDER FOR COMPLETE COMMAND:
	#### TOPHAT COMMAND + ENVIRONMENT VARIABLES, ETC.
	my $command; 

	my $tempdir = $self->tempdir();
	if ( defined $tempdir and $tempdir and not -d $tempdir )
	{
		print "Expression::TOPHAT::batchCommand    tempdir directory not found: $tempdir\n" if not -d $tempdir;
		print "Expression::TOPHAT::batchCommand    tempdir is a file: $tempdir\n" if -f $tempdir;
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

mkdir -p $outdir
cd $outdir

mkdir -p $temp_outdir

$tophat_command

mv $temp_outdir/* $outdir

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

mkdir -p $outdir
cd $outdir

$tophat_command};

	}

	return $command;

}	#	batchCommand




=head2

	SUBROUTINE		tophatCommand

	PURPOSE

		CREATE .sh FILE

=cut


sub tophatCommand {
	my $self			=	shift;
	my $inputfiles		=	shift;
	my $matefiles		=	shift;
	my $outputdir		=	shift;
	my $reference		=	shift;



	#### USER INPUTS
	my $distance 		=	$self->distance();
	my $params 			=	$self->params();
	my $label 			=	$self->label();
	my $keep 			=	$self->keep();

	#### EXECUTABLES
	my $tophat 			=	$self->tophat();
	my $bowtie 			=	$self->bowtie();

	##### CLUSTER
	my $cpus 			=	$self->cpus();

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

sub runCufflinks {
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


sub cufflinksCommand {
	my $self			=	shift;
	my $samfile			=	shift;



	my ($outdir) = $samfile =~ /^(.+?)\/[^\/]+$/;

	#### USER INPUTS
	my $distance 		=	$self->distance();
	my $label 			=	$self->label();
	my $gtf 			=	$self->gtf();

	#### EXECUTABLES
	my $cufflinks 		=	$self->cufflinks();

	##### CLUSTER
	my $cpus 			=	$self->cpus();

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




}	#### Expression::TOPHAT

