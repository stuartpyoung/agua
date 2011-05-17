use Getopt::Simple;
use MooseX::Declare;
#use Moose::Util::TypeConstraints;

=head2

		PACKAGE		ELAND

		PURPOSE

	        WRAPPER SCRIPT FOR RUNNING ELAND SEQUENCE ALIGNMENT

=cut#### EXTERNAL MODULES
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

#class Aligner::ELAND {
class Aligner::ELAND with (Agua::Cluster::Checker, Agua::Cluster::Jobs) {
	has 'referencedir'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
	has 'outputdir'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
	has 'replicates'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
	has 'label'			=> ( isa => 'Str|Undef', is => 'rw', default => '' );
	has 'splitfiles'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
	has 'cluster'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );

	has 'conf'		=> ( isa => 'HashRef', is => 'rw', default => '' );
	has 'check'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
	has 'inputfiles'=> ( isa => 'Str|Undef', is => 'rw', default => '' );
	has 'matefiles'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
	has 'distance'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
	has 'params'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
	has 'sequencedir'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
	has 'rundir'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
	has 'maxlines'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
	has 'chunks'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
	has 'clean'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
	has 'keep'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
	has 'readhits'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
	has 'referencedir'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
	has 'species'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
	has 'subdirs'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
	has 'casava'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
	has 'inputtype'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
	has 'pairparams'=> ( isa => 'Str|Undef', is => 'rw', default => '' );
	has 'quality'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
	has 'seedlength'=> ( isa => 'Str|Undef', is => 'rw', default => '' );
	has 'samtools'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
	has 'samtoolsindex'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
	has 'lane'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
	has 'firstlane'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
	has 'convert'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );

	##### INTERNAL MODULES
	#use Agua::Cluster::Util;

	#### EXTERNAL MODULES
	use Data::Dumper;
	use File::Path;
	use File::Remove;

	#### DEBUG


	#/////}


=head

	SUBROUTINE		run

	PURPOSE

		DO BATCH ALIGNMENT FOLLOWED BY SNP CALLS

=cut

sub run {
	my $self		=	shift;


	#### INPUTS
	my $referencedir 	=	$self->get_referencedir();
	my $outputdir 		=	$self->get_outputdir();
	my $inputfiles 		=	$self->get_inputfiles();
	my $matefiles 		=	$self->get_matefiles();
	my $readhits 		=	$self->get_readhits();

	#### CHECK INPUT FILES EXIST AND NOT EMPTY
	my @infiles = split ",", $inputfiles;
	my @mates = split ",", $matefiles if defined $matefiles;
	$self->checkFiles(\@infiles);
	$self->checkFiles(\@mates);

	#### SPLIT FILES
	my $maxlines 		=	$self->get_maxlines();
	my $clean 			=	$self->get_clean();
	my $label 			=	$self->get_label();
	my $splitfile 		=	$self->get_splitfile();

	#### SET DEFAULT SPLITFILE IF NOT DEFINED
	$splitfile = "$outputdir/splitfile.txt" if not defined $splitfile;

#### INPUTS

	#### GET REFERENCES -- LATER CHANGE THIS: USES chr*.vld FILES TO GET NAMES OF chr* DIRS
	my $referencefiles = $self->listReferenceFiles($referencedir, "\*\.vld");
	my $references = $self->getReferences($referencedir);

	##############################################
	############	SPLIT INPUT FILES   ###########
	##############################################
	print "ELAND::run    Doing doSplitfiles()   ", Timer::current_datetime(), "\n";
	my $splitfiles = $self->doSplitfiles($splitfile, $label);
	print "ELAND::run    After doSplitfiles()   ", Timer::current_datetime(), "\n";


#### ALIGNMENT

	#############################################
	###########	ELAND ALIGNMENT    ###########
	#############################################
	print "ELAND::run    Doing doAlignment()    ", Timer::current_datetime(), "\n";
	$self->doBatchAlignment($outputdir, $referencefiles, $splitfiles, $label);
	print "ELAND::run    After doAlignment()    ", Timer::current_datetime(), "\n";


#### CONVERT TO CHROMOSOME SAM 

	################################################
	######   CONVERT *_export.txt TO out.sam    ######
	###################################ll#############
	print "ELAND::run    Doing exportToSam      ", Timer::current_datetime(), "\n";
	$self->exportToSam($splitfiles, $referencedir, $outputdir, $references, $matefiles);
	print "ELAND::run    After exportToSam      ", Timer::current_datetime(), "\n";

	################################################
	#######        FILTER SAM HITS          ######print Dumper ;
	################################################
	#$self->chromosomeHits($outputdir, $references, $splitfiles, "out.sam", "hit.sam", "miss.sam");

	################################################
	#######        FILTER SAM HITS          ######
	################################################
	print "ELAND::run    Doing subdirSamHits        ", Timer::current_datetime(), "\n";
	$self->subdirSamHits($outputdir, $references, $splitfiles, "out.sam", "hit.sam", "miss.sam");
	print "ELAND::run    After subdirSamHits        ", Timer::current_datetime(), "\n";

	################################################
	######        CUMULATIVE MERGE SAM        ######
	################################################
	print "ELAND::run    Doing cumulativeMergeSam        ", Timer::current_datetime(), "\n";
	$self->cumulativeMergeSam($outputdir, $references, $splitfiles, "hit.sam", "hit.sam");
	print "ELAND::run    After cumulativeMergeSam        ", Timer::current_datetime(), "\n";

	##############################################
	#########         SORT SAM         ###########
	##############################################
	print "ELAND::run    Doing samToBam     ", Timer::current_datetime(), "\n";
	$self->samToBam($outputdir, $references, "hit.sam", "hit.bam");
	print "ELAND::run    After samToBam     ", Timer::current_datetime(), "\n";

#### READ HIT STATS

	###############################################
	###########        READ HITS         ##########
	###############################################
	#if ( defined $readhits )
	#{
	#	$self->readHits($outputdir, $references, $splitfiles, "eland-readHits", "hit.sam");
	#}

#### SNP CALLING

	###############################################
	#####     PREDICT SNPs WITH SAMTOOLS     ######
	###############################################
	#$self->samtoolSnps($splitfiles, $references, $outputdir, "hit.sam");


#### USAGE STATS

	###############################################
	#####          PRINT USAGE STATS         ######
	###############################################
	#my $usagefile = "$outputdir/ELAND-USAGE.txt";
	#$self->printUsage($usagefile);
}
=head2

	SUBROUTINE		doBatchAlignment

	PURPOSE

		RUN ELAND AGAINST ALL REFERENCE FILES USING CLUSTER BATCH JOB

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
	$self->runJobs($jobs, "eland");
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
	my $matefiles 		= 	$self->get_matefiles();
	my $params 			=	$self->get_params();  #### OVERRIDE PARAMS IN DEFAULT COMMAND
	my $label 			=	$self->get_label();	#### USED TO GENERATE INPUTFILE NAMES	
	my $inputtype 		=	$self->get_inputtype();

	#### CHECK INPUTS
	print "ELAND::batchCommand    inputtype not defined. Exiting\n" and exit if not defined $inputtype;
	print "ELAND::batchCommand    label not defined. Exiting\n" if not defined $label;

	#### EXECUTABLES
	my $casava 			= 	$self->get_casava();
	print "ELAND::batchCommand    casava not defined. Exiting\n" and exit if not defined $casava;

	##### CLUSTER
	my $cluster 		=	$self->get_cluster();

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
	my ($referencedir, $reference) = $referencepath =~ /^(.+?)\/([^\/]+)$/;
	$reference =~ s/\.vld$//;
	$reference =~ s/\.fa$//;

	#### SET INPUT AND MATE FILES, E.G.:
	#### /scratch/syoung/base/pipeline/bixby/run1/ln/$LSB_JOBINDEX/ln_1.$LSB_JOBINDEX.txt
	#### /scratch/syoung/base/pipeline/bixby/run1/ln/$LSB_JOBINDEX/ln_2.$LSB_JOBINDEX.txt	
	my $firstmate = $label . "_1";
	my $inputfile = "$basedir/$index/$firstmate.$index$suffix";

	#### DO SECOND MATE IF matefiles DEFINED
	my $secondmate = $label . "_2" if defined $matefiles;
	my $matefile = "$basedir/$index/$secondmate.$index$suffix" if defined $matefiles;

	#### GET OPTIONAL VARIABLES
	my $seedlength = $self->get_seedlength();
	my $quality = $self->get_quality();
	my $pairparams = $self->get_pairparams();

	#### ELAND COMMAND
	my $eland_command = qq{$casava/ELAND_standalone.pl \\\n};
	$eland_command .= qq{--input-type $inputtype \\\n};
	$eland_command .= qq{--eland-genome $referencedir/$reference \\\n};
	$eland_command .= qq{--input-file $inputfile \\\n};
	$eland_command .= qq{--input-file $matefile \\\n} if defined $matefile;	
	$eland_command .= qq{--seedlength $seedlength \\\n} if defined $seedlength;
	$eland_command .= qq{--base-quality $quality \\\n} if defined $quality;
	$eland_command .= qq{--pair-params $pairparams \\\n} if defined $pairparams;
	$eland_command .= "\n";

	#### CHECK TEMPDIR EXISTS AND IS NOT A FILE
	my $tempdir = $self->get_tempdir();
	if ( defined $tempdir and $tempdir and not -d $tempdir )
	{
		print "ELAND::batchCommand    tempdir directory not found: $tempdir\n" if not -d $tempdir;
		print "ELAND::batchCommand    tempdir is a file: $tempdir\n" if -f $tempdir;
	}

	#### SET OUTPUT DIRECTORY
	my $outdir = "$basedir/$reference/$index";

	#### SET TEMP-RELATED DIRS
	my $old_outdir	=	$outdir;
	my $temp_outdir =	"$tempdir/$outdir" if defined $tempdir;

	#### USE TEMPDIR IF DEFINED
	if ( defined $tempdir and $tempdir and -d $tempdir )
	{

		#### SET TEMP OUTPUT DIR
		$outdir = $temp_outdir;
	}

	#### SET SCRIPT
	my $command = qq{
export PATH=$casava:\$PATH
export PATH=$referencedir/\$LSB_JOBINDEX:\$PATH
mkdir -p $outdir
cd $outdir
$eland_command
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

	SUBROUTINE		convertReferences

	PURPOSE

		CONVERT .fa REFERENCE FASTA FILES INTO .idx ELAND INDEX FORMAT FILES

=cut

sub convertReferences {
	my $self			=	shift;
	my $inputdir		=	shift;
	my $outputdir		=	shift;


	$outputdir = $inputdir if not defined $outputdir;

	#### SANITY CHECK
	print "ELAND::convertReferences    inputdir not defined: $inputdir\n" and exit if not defined $inputdir;
	print "ELAND::convertReferences    outputdir not defined: $outputdir\n" and exit if not defined $outputdir;
	print "ELAND::convertReferences    inputdir is a file: $inputdir\n" and exit if -f $inputdir;
	print "ELAND::convertReferences    outputdir is a file: $outputdir\n" and exit if -f $outputdir;
	print "ELAND::convertReferences    Can't find inputdir: $inputdir\n" and exit if not -d $inputdir;
	print "ELAND::convertReferences    Can't find outputdir: $outputdir\n" and exit if not -d $outputdir;

	#### GET REFERENCE FILES
	my $fastafiles = $self->listFiles($inputdir, "\*.fa");

	#### GET CASAVA
	my $casava		=	$self->get_casava();
	print "ELAND::convertReferences    casava not defined\n" and exit if not defined $casava or not $casava;

	###############################################
	############# MAKE BINARY REF FILES ###########
	###############################################
	my $current_time = time();	
	for my $fastafile ( @$fastafiles )
	{
		next if $fastafile =~ /^[\.]+$/;

		my ($reference) = $fastafile =~ /^.+?\/([^\/]+)$/;

		#### SET REFERENCE BINARY FILE
		my $referencebinary = "$outputdir/$reference";
		$referencebinary .= ".vld";

		#### CONVERT REFERENCE FILE INTO BINARY FILE
		my $command;
		$command = "time $casava/squashGenome $outputdir $fastafile" if not -f $referencebinary or -z $referencebinary;
		print "ELAND::convertReferences    command: $command\n" if defined $command;
		print `$command` if defined $command;
		print "ELAND::convertReferences    Binary reference files exist. Skipping squash\n" if not defined $command;
	}

	#### GET DURATION OF LOCAL JOB AND ADD AS USAGE STATISTIC
	my $duration = Timer::runtime( $current_time, time() );
	$self->addUsageStatistic("convertReferences", $duration);
}





=head2

	SUBROUTINE		exportToSam

	PURPOSE

		1. CONVERT *_export.txt ELAND OUTPUT FILE TO out.sam SAMTOOLS FILE

		2. DO FOR *_export.txt FILES PRODUCED BY ALIGNMENTS OF ALL SPLITFILES

			AGAINST ALL REFERENCE SUBDIRECTORIES

	NOTES

		NB: THIS USES --qlogodds OPTION FOR export2sam.p WHICH ASSUMES THAT

			THE INPUT FILES ARE IN solexa FORMAT

=cut

sub exportToSam {
	my $self			=	shift;
	my $splitfiles 		=	shift;
	my $referencedir 	=	shift;
	my $outputdir 		=	shift;
	my $references 		=	shift;
	my $matefiles		=	shift;


	#### GET REQUIRED VARIABLES
	my $samtools = $self->get_samtools();
	print "ELAND::exportToSam    samtools not defined. Exiting\n" and exit if not defined $samtools;

	#### COLLECT ALL JOBS FOR EACH INPUT FILE AGAINST ALL REFERENCE FILES
	my $jobs = [];
	foreach my $reference ( @$references )
	{
		my $refdir = "$referencedir/$reference";
		next if $refdir =~ /^[\.]+$/;

		my ($basedir, $task) = $$splitfiles[0][0] =~ /^(.+?)\/(\d+)\/[^\/]+$/;

		##### CLUSTER
		my $cluster 		=	$self->get_cluster();

		#### SET INDEX PATTERN FOR BATCH JOB
		my $index;
		$index = "\$LSB_JOBINDEX" if $cluster eq "LSF";
		$index = "\$PBS_TASKNUM" if $cluster eq "PBS";

		my $outdir = "$basedir/$reference";
		File::Path::mkpath($outdir) if not -d $outdir;
		print "Could not create output dir: $outdir\n" if not -d $outdir;

		#### OUTPUT SAM FILE
		my $samfile = "$outdir/$index/out.sam";

		#### INPUT EXPORT FILES
		my $exportfile_1 = "$outdir/$index/reanalysis_export.txt";
		my $exportfile_2;
		if ( defined $matefiles )
		{
			$exportfile_1 = "$outdir/$index/reanalysis_1_export.txt";
			$exportfile_2 = "$outdir/$index/reanalysis_2_export.txt";
		}

		#### SET LOCATION OF export2sam
		my $export2sam = "$samtools/export2sam.pl";
		$samtools =~ /(\d+)\.(\d+)\.(\d+)(\/)?$/;
		my $hundreds= 	$1 || 0;
		my $tens	=	$2 || 0;
		my $ones	=	$3 || 0;
		my ($samtools_version) = $hundreds * 100 + $tens * 10 + $ones;

		#### SET COMMAND
		my $command = qq{$samtools/export2sam.pl $exportfile_1 $exportfile_2 > $samfile\n};

		if ( $samtools_version >= 18 )
		{
			####	export2sam.pl converts GERALD export files to SAM format.
			####	
			####	Usage: export2sam.pl --read1=FILENAME [ options ] | --version | --help
			####	
			####	  --read1=FILENAME  read1 export file (mandatory)
			####	  --read2=FILENAME  read2 export file
			####	  --nofilter        include reads that failed the pipeline/RTA purity filter
			####	  --qlogodds        assume export file(s) use logodds quality values as reported
			####						  by pipeline prior to v1.3 (default: phred quality values)

			$command = "$samtools/misc/export2sam.pl --read1=$exportfile_1 ";
			$command .= " --read2=$exportfile_2" if defined $exportfile_2;
			$command .= " --qlogodds ";
			$command .= " > $samfile\n";
		}	

		#### SET LABEL AND TASKS
		my $label = "exportToSam-$reference";
		my $tasks = scalar(@$splitfiles);

		#### SET JOB
		my $job = $self->setBatchJob( [ $command ], $label, $outdir, $tasks);
#exit;

		push @$jobs, $job;	
	}

	#### RUN ALIGNMENT JOBS
	$self->runJobs($jobs, "exportToSam");	
}

=head2

	SUBROUTINE		getReferences

	PURPOSE

		RETURN A LIST OF REFERENCE NAMES (N.B.: NOT THE PATHS TO THE FILES)

=cut

sub getReferences {
	my $self		=	shift;
	my $referencedir=	shift;


	my $referencefiles = $self->listReferenceFiles($referencedir, "\*\.vld");
	my $references = $self->getFilenames($referencefiles);
	foreach my $reference ( @$references )
	{
		$reference =~ s/\.vld$//;
		$reference =~ s/\.fa$//;
	}

	return $references;
}






}

