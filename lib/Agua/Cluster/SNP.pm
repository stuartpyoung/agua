package Agua::Cluster::SNP;
use Moose::Role;
use Moose::Util::TypeConstraints;

use Data::Dumper;


################################################################################
##########################         SNP METHODS        ##########################
################################################################################
=head2

	SUBROUTINE		samtoolSnps

	PURPOSE

		USE SAMTOOLS TO PREDICT SNPs

	NOTES

		THE (REQUIRED) maxdepth ARGUMENT IS USED AS THE VALUE FOR THE -D OPTION OF varFilter

			http://samtools.sourceforge.net/samtools.shtml
			The -D option of varFilter controls the maximum read depth, which should be adjusted to about twice the average read depth. One may consider to add -C50 to mpileup if mapping quality is overestimated for reads containing excessive mismatches. Applying this option usually helps BWA-short but may not other mappers.

		THE SNP CALLING PART OF THIS SUBROUTINE CONSISTS ESSENTIALLY OF TWO COMMANDS:

			1. GET VARIANTS

				samtools pileup -vcf ref.fa aln.bam | tee raw.txt | samtools.pl varFilter -D100 > flt.txt

				1. Get the raw variant calls on the X chromosome. Assume X chromosome is named X in both hsRef.fa and aln.bam.

				   samtools view -u aln.bam X | samtools pileup -vcf hsRef.fa - > var-X.raw.txt

				   A region is specified with chr, chr:begin or chr:begin-end.

				2. Filter the raw variant calls.

				   samtools.pl varFilter -D100 var-X.raw.txt > var-X.flt.txt

				   Please remember to set the -D to set a maximum read depth according to the average read depth. In whole genome shotgun (WGS) resequencing, SNPs with excessively high read depth are usually caused by structural variations or alignment artifacts and should not be trusted.

			2. QUALITY FILTER VARIANTS	

				Acquire final variant calls by setting a quality threshold.

			   awk '($3=="*"&&$6>=50)||($3!="*"&&$6>=20)' var-X.flt.txt > var-X.final.txt

			   Here 50 is the quality threshold for indels and 20 for substitutions.


	#### 1. MERGE SPLITFILE SAM FILES INTO ONE SAM FILE

	#### 2. CONVERT SAM TO BAM FILE

	#### 3. CONVERT BAM TO SNP FILE

=cut

sub samtoolSnps {
	my $self			=	shift;
	my $splitfiles		=	shift;
	my $references		=	shift;
	my $outputdir		=	shift;
	my $samfile			=	shift;
	my $maxdepth		=	shift;

	my $bamfile = $samfile;
	$bamfile =~ s/\.sam$/.bam/;


	#########################################################################################
	############ STRATEGY 1. MERGE SUBDIR out.sam FILES FIRST, THEN SORT SINGLE SAM FILE
	#########################################################################################
	#### THIS METHOD IS FASTER THAN STRATEGY 2

		### 1. CONVERT SAM TO BAM FILE
		print "Agua::Cluster::SNP::samtoolSnps    Doing samToBam      ", Timer::current_datetime(), "\n";
		$self->samToBam($outputdir, $references, $samfile, $bamfile);
		print "Agua::Cluster::SNP::samtoolSnps    After samToBam      ", Timer::current_datetime(), "\n";

		### 2. SORT BAM FILE
		### Sort alignments by leftmost coordinates. File <out.prefix>.bam will be created.
		### This command may also create temporary files <out.prefix>.%d.bam when the
		### whole alignment cannot be fitted into memory (controlled by option -m).
		###	DO THIS TO AVOID THIS ERROR:
		###	[bam_sort_load] fail to load BAM sort.
		###	[main_samview] random alignment retrieval only works for sorted BAM files.
		print "Agua::Cluster::SNP::samtoolSnps    Doing sortBam      ", Timer::current_datetime(), "\n";
		$self->sortBam($outputdir, $references, $bamfile, $bamfile);
		print "Agua::Cluster::SNP::samtoolSnps    After sortBam      ", Timer::current_datetime(), "\n";

		### 3. INDEX BAM FILE
		### Index sorted alignment for fast random access.
		### Index file <aln.bam>.bai will be created.
		###	DO THIS TO AVOID THIS ERROR:
		###	[bam_index_load] fail to load BAM index.
		###	[main_samview] random alignment retrieval only works for indexed BAM files.
		print "Agua::Cluster::SNP::samtoolSnps    Doing indexBam      ", Timer::current_datetime(), "\n";
		$self->indexBam($outputdir, $references);
		print "Agua::Cluster::SNP::samtoolSnps    After indexBam      ", Timer::current_datetime(), "\n";

		### 4. CONVERT BAM TO SNP FILE
		print "Agua::Cluster::SNP::samtoolSnps    Doing bamToSnp      ", Timer::current_datetime(), "\n";
		$self->bamToSnp($outputdir, $references, $bamfile, $maxdepth);
		print "Agua::Cluster::SNP::samtoolSnps    After bamToSnp      ", Timer::current_datetime(), "\n";


	#####################################################################################
	##### STRATEGY 2. GENERATE SNPS FROM SAM FILE
	#####################################################################################
	#### THIS METHOD IS SLOWER THAN STRATEGY 1 BUT DOES HAVE THE ADVANTAGE OF EASILY 
	#### VIEWED INTERMEDIATE FILES, I.E., *.sam AND *.sam.gz FILES

	######## 1. PRINT head AND tail OF SAM FILE TO *.sam.head AND *.sam.tail FILES
	####print "Agua::Cluster::SNP::samtoolSnps    Doing headTail      ", Timer::current_datetime(), "\n";
	####$self->headTail($outputdir, $references, $samfile);
	####print "Agua::Cluster::SNP::samtoolSnps    After headTail      ", Timer::current_datetime(), "\n";
	####
	######## 2. ZIP UP SAM FILE (IN ORDER TO USE SAMTOOLS TO CALL SNPS)
	####print "Agua::Cluster::SNP::samtoolSnps    Doing gzipFile      ", Timer::current_datetime(), "\n";
	####$self->gzipFile($outputdir, $references, $samfile);
	####print "Agua::Cluster::SNP::samtoolSnps    After gzipFile      ", Timer::current_datetime(), "\n";
	####
	######## 3. CONVERT SAM TO SNP FILE
	####print "Agua::Cluster::SNP::samtoolSnps    Doing samToSnp      ", Timer::current_datetime(), "\n";
	####$self->samToSnp($outputdir, $references, "$samfile.gz");
	####print "Agua::Cluster::SNP::samtoolSnps    After samToSnp      ", Timer::current_datetime(), "\n";
}


=head2

	SUBROUTINE		headTail

	PURPOSE

		1. PRINT HEAD AND TAIL OF A SAM FILE TO A *.sam.headtail FILE

		2. DO SO FOR ALL CHROMOSOME SAM FILES

=cut

sub headTail {
	my $self			=	shift;
	my $outputdir		=	shift;
	my $references 		=	shift;
	my $inputfile		=	shift;

	print "Agua::Cluster::SNP::headTail    Agua::Cluster::headTail(outputdir, references, inputfile)\n";

	my $commands = [];
	foreach my $reference ( @$references )
	{
		push @$commands, "head $outputdir/$reference/$inputfile > $outputdir/$reference/$inputfile.head";
		push @$commands, "tail $outputdir/$reference/$inputfile > $outputdir/$reference/$inputfile.tail";
	}

	my $label = "headTail";		
	my $job = $self->setJob( $commands, $label, $outputdir);
	my $jobs = [];
	push @$jobs, $job;

	#### RUN JOBS
	$self->runJobs($jobs, "headTail");
}



=head2

	SUBROUTINE		gzipFile

	PURPOSE

		1. ZIP UP A *.sam FILE

		2. DO SO FOR ALL CHROMOSOME SAM FILES

=cut

sub gzipFile {
	my $self			=	shift;
	my $outputdir		=	shift;
	my $references 		=	shift;
	my $inputfile		=	shift;

	print "Agua::Cluster::SNP::gzipFile    Agua::Cluster::gzipFile(outputdir, references, inputfile)\n";

	my $jobs = [];
	foreach my $reference ( @$references )
	{
		my $command = "gzip -f $outputdir/$reference/$inputfile; sleep 5; echo 'gzipFile completed'";
		my $label = "gzipFile";		
		my $job = $self->setJob( [ $command ], $label, "$outputdir/$reference");
		push @$jobs, $job;
	}

	#### RUN JOBS
	$self->runJobs($jobs, "gzipFile");
}


=head2

	SUBROUTINE		samToSnp

	PURPOSE

		PREDICT SNPs FROM SAM FILE ALIGNMENTS AGAINST REFERENCES

	NOTES

		BASICALLY RUN TWO COMMANDS:

			1. GET VARIANTS

				samtools pileup -vcf ref.fa aln.bam | tee raw.txt | samtools.pl varFilter -D100 > flt.txt

				1. Get the raw variant calls on the X chromosome. Assume X chromosome is named X in both hsRef.fa and aln.bam.

				   samtools view -u aln.bam X | samtools pileup -vcf hsRef.fa - > var-X.raw.txt

				   A region is specified with chr, chr:begin or chr:begin-end.

				2. Filter the raw variant calls.

				   samtools.pl varFilter -D100 var-X.raw.txt > var-X.flt.txt

				   Please remember to set the -D to set a maximum read depth according to the average read depth. In whole genome shotgun (WGS) resequencing, SNPs with excessively high read depth are usually caused by structural variations or alignment artifacts and should not be trusted.

			2. QUALITY FILTER VARIANTS	

				Acquire final variant calls by setting a quality threshold.

			   awk '($3=="*"&&$6>=50)||($3!="*"&&$6>=20)' var-X.flt.txt > var-X.final.txt

			   Here 50 is the quality threshold for indels and 20 for substitutions.

		INPUTS

			1. OUTPUT DIRECTORY

			2. LIST OF FASTA REFERENCE FILES

			3. ***SORTED*** SAM FILE IN EACH OUTPUTDIR/<REFERENCE> SUB-DIRECTORY

		OUTPUTS

			1. PILEUP-FORMAT SAMTOOLS SNP PREDICTIONS FILE

=cut

sub samToSnp {
	my $self			=	shift;
	my $outputdir		=	shift;
	my $references 		=	shift;
	my $inputfile		=	shift;
	my $maxdepth		=	shift;

	print "Agua::Cluster::SNP::samToSnp    Agua::Cluster::samToSnp(outputdir, references, inputfile)\n";
	print "Agua::Cluster::SNP::samToSnp    outputdir: $outputdir\n";
	print "Agua::Cluster::SNP::samToSnp    No. references: ", scalar(@$references), "\n";
	print "Agua::Cluster::SNP::samToSnp    inputfile: $inputfile\n";

	#### CHECK INPUTS
	print "Agua::Cluster::SNP::samToSnp    outputdir not defined\n" and exit if not defined $outputdir;
	print "Agua::Cluster::SNP::samToSnp    references not defined\n" and exit if not defined $references;
	print "Agua::Cluster::SNP::samToSnp    inputfile not defined\n" and exit if not defined $inputfile;
	print "Agua::Cluster::SNP::samToSnp    maxdepth not defined\n" and exit if not defined $maxdepth;

	#### GET REQUIRED VARIABLES
	my $samtools = $self->samtools();
	my $samtools_index = $self->samtoolsindex();

	my $jobs = [];
	foreach my $reference ( @$references )
	{
		#### GET SAM FILE
		my $samfile = "$outputdir/$reference/$inputfile";

		#### STORE COMMANDS
		my $commands = [];

		#### SET SNP FILE
		my ($filestub) = $samfile =~ /^(.+)\.sam(\.gz)?$/;
		my $rawfile = "$filestub.raw";
		my $filterfile = "$filestub.filter";
		my $snpfile = "$filestub.snp";

		#### GET RAW VARIANT CALLS
		#### AND FILTER VARIANT CALLS
		my $referencefile = $reference . ".fa";
		my $indexfile = $reference . ".fai";
		my $command = "$samtools/samtools pileup -vcf $samtools_index/$referencefile -t $samtools_index/$indexfile $samfile | tee $rawfile | $samtools/misc/samtools.pl varFilter -D$maxdepth> $filterfile";
		push @$commands, $command;
		print "Agua::Cluster::SNP::samToSnp    $command\n";

		#### QUALITY FILTER VARIANTS
		$command = qq{awk '(\$3=="*"&&\$6>=50)||(\$3!="*"&&\$6>=20)' $filterfile > $snpfile};
		push @$commands, $command;

		my $label = "samToSnp-$reference";
		my $outdir = "$outputdir/$reference";

		my $job = $self->setJob( $commands, $label, $outdir);
		push @$jobs, $job;
	}

	#### RUN JOBS
	$self->runJobs($jobs, "samToSnp");
}




=head2

	SUBROUTINE		bamToSnp

	PURPOSE

		PREDICT SNPs FROM BAM FILE ALIGNMENTS AGAINST REFERENCES

	NOTES

		BASICALLY RUN TWO COMMANDS:

			1. GET VARIANTS

				samtools pileup -vcf ref.fa aln.bam | tee raw.txt | samtools.pl varFilter -D100 > flt.txt

				1. Get the raw variant calls on the X chromosome. Assume X chromosome is named X in both hsRef.fa and aln.bam.

				   samtools view -u aln.bam X | samtools pileup -vcf hsRef.fa - > var-X.raw.txt

				   A region is specified with chr, chr:begin or chr:begin-end.

				2. Filter the raw variant calls.

				   samtools.pl varFilter -D100 var-X.raw.txt > var-X.flt.txt

				   Please remember to set the -D to set a maximum read depth according to the average read depth. In whole genome shotgun (WGS) resequencing, SNPs with excessively high read depth are usually caused by structural variations or alignment artifacts and should not be trusted.

			2. QUALITY FILTER VARIANTS	

				Acquire final variant calls by setting a quality threshold.

			   awk '($3=="*"&&$6>=50)||($3!="*"&&$6>=20)' var-X.flt.txt > var-X.final.txt

			   Here 50 is the quality threshold for indels and 20 for substitutions.

		INPUTS

			1. OUTPUT DIRECTORY

			2. LIST OF FASTA REFERENCE FILES

=cut

sub bamToSnp {
	my $self			=	shift;
	my $outputdir		=	shift;
	my $references 		=	shift;
	my $bamfile			=	shift;
	my $maxdepth		= 	shift;

	#### CHECK INPUTS
	print "Agua::Cluster::SNP::samToSnp    outputdir not defined\n" and exit if not defined $outputdir;
	print "Agua::Cluster::SNP::samToSnp    references not defined\n" and exit if not defined $references;
	print "Agua::Cluster::SNP::samToSnp    maxdepth not defined\n" and exit if not defined $maxdepth;
	print "Agua::Cluster::SNP::bamToSnp    bamfile not defined. Exiting\n" and exit if not defined $bamfile;

	#### GET REQUIRED VARIABLES
	my $samtools = $self->samtools();
	my $samtools_index = $self->samtoolsindex();

	my $jobs = [];
	foreach my $reference ( @$references )
	{
		#### GET BAM FILE
		my $bam = "$outputdir/$reference/$bamfile";

		#### STORE COMMANDS
		my $commands = [];

		#### SET SNP FILE
		my ($filestub) = $bam =~ /^(.+)\.bam$/;
		my $rawfile = "$filestub.raw";
		my $filterfile = "$filestub.filter";
		my $snpfile = "$filestub.snp";

		my $referencefile = $reference . ".fa";

		#### GET RAW VARIANT CALLS
		#### AND FILTER VARIANT CALLS
		# USE pileup TO GET RAW SNP PREDICTIONS
		# THEN FILTER RAW VARIANTS WITH varFilter
		#
		# http://sourceforge.net/apps/mediawiki/samtools/index.php?title=SAM_protocol#Support_Protocol_2:_Local_Realignment
		# varFilter - Filter the raw variant calls
		# Please remember to set the -D to set a maximum read depth according to the average read depth
		my $command = "$samtools/samtools pileup -vcf $samtools_index/$referencefile $bam | tee $rawfile | $samtools/samtools.pl varFilter -D$maxdepth > $filterfile";
		push @$commands, $command;
print "Agua::Cluster::SNP::bamToSnp    $command\n";

		#### QUALITY FILTER VARIANTS
		# 
		# Acquire final variant calls by setting a quality threshold
		# Here 50 is the quality threshold for indels and 20 for substitutions.
		$command = qq{awk '(\$3=="*"&&\$6>=50)||(\$3!="*"&&\$6>=20)' $filterfile > $snpfile};
		push @$commands, $command;

		my $label = "bamToSnp-$reference";
		my $outdir = "$outputdir/$reference";

		my $job = $self->setJob( $commands, $label, $outdir);
		push @$jobs, $job;
	}

	#### RUN JOBS
	$self->runJobs($jobs, "bamToSnp");
}




1;
