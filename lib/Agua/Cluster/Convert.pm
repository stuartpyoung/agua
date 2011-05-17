package Agua::Cluster::Convert;
use Moose::Role;
use Moose::Util::TypeConstraints;

use Data::Dumper;


################################################################################
##########################      CONVERT METHODS        #########################
################################################################################
=head2

	SUBROUTINE		simpleFastqHeader

	PURPOSE

		CONVERT FASTQ HEADER: REMOVE MACHINE ID AND ADD READ '#/NUMBER' 

=cut

sub simpleFastqHeader {
	my $self		=	shift;
	my $inputfiles	=	shift;
	my $matefiles	=	shift;



	my $infiles;
	push @$infiles, split ",", $inputfiles;

	my $outfiles;
	my $matenumber = 1;
	#my $matenumber;
	#$matenumber = 1 if defined $matefiles;
	foreach my $infile ( @$infiles )
	{
		my $outfile	=	$infile;
		$outfile =~ s/(\.[^\.]{1,5})$/.simpleHeader$1/;
		push @$outfiles, $outfile;
		$self->simpleHeader($infile, $outfile, $matenumber);
	}
	$self->set_inputfiles(join ",", @$outfiles);
	print "Agua::Cluster::Convert::simpleFastqHeader    NEW inputfiles: ", $self->inputfiles(), "\n";

	return if not defined $matefiles;

	#### DO MATEFILES
	$matenumber = 2;
	$infiles = undef;	
	push @$infiles, split ",", $matefiles;

	$outfiles = undef;
	foreach my $infile ( @$infiles )
	{
		my $outfile	=	$infile;
		$outfile =~ s/(\.[^\.]{1,5})$/.simpleHeader$1/;
		push @$outfiles, $outfile;
		$self->simpleHeader($infile, $outfile, $matenumber);
	}
	$self->set_matefiles(join ",", @$outfiles);
	print "Agua::Cluster::Convert::simpleFastqHeader    NEW matefiles: ", $self->matefiles(), "\n";	
}



=head2

	SUBROUTINE		simpleHeader

	PURPOSE

		CONVERT FASTQ HEADER: REMOVE MACHINE ID AND ADD READ '#/NUMBER' 

=cut


sub simpleHeader {
	my $self		=	shift;
	my $infile		=	shift;
	my $outfile		=	shift;
	my $matenumber	=	shift;


	if ( not -f $outfile or -z $outfile )
	{
		my $command = "sed -e 's/ [A-Za-z0-9\\.\\_\\-]*:\\([0-9\\:]*[0-9\\:]*[0-9\\:]*[0-9\\:]*\\)/:\\1#0\\/$matenumber/' < $infile  | sed -e 's/[ ]*length=[0-9]*[ ]*//' > $outfile";
		print "Agua::Cluster::Convert::simpleHeader    command: $command\n";
		print `$command`;
	}
	else
	{
		print "Agua::Cluster::Convert::simpleHeader    Skipping conversion. outfile already exists and is non-empty: $outfile\n";
	}
}

=head2

	SUBROUTINE		subdirSamToBam

	PURPOSE

		CONVERT ALL SAM FILES IN chr*/<NUMBER> SUBDIRECTORIES

		INTO BAM FILES

			samtools view -bt ref_list.txt -o aln.bam aln.sam.gz

		INPUTS

			1. OUTPUT DIRECTORY USED TO PRINT CHROMOSOME-SPECIFIC

				SAM FILES TO chr* SUB-DIRECTORIES

			2. LIST OF REFERENCE FILE NAMES

			3. SPLIT INPUT FILES LIST

			4. NAME OF SAM INPUTFILE (E.G., "accepted_hits.sam")

=cut

sub subdirSamToBam {
	my $self			=	shift;
	my $outputdir		=	shift;
	my $references 		=	shift;
	my $splitfiles		=	shift;
	my $inputfile		=	shift;

	#### SET SAMTOOLS INDEX
	my $samtools = $self->samtools();
	my $samtools_index = $self->samtoolsindex();

	my $jobs = [];
	foreach my $reference ( @$references )
	{
		#### GET SAM FILES
		my $samfiles = $self->subfiles($outputdir, $reference, $splitfiles, $inputfile);

		#### 	1. CONVERT SAM FILES INTO BAM FILES
		# samtools view -bt ref_list.txt -o aln.bam aln.sam.gz
		my $bamfiles = [];
		my $counter = 0;
		my $commands = [];
		foreach my $samfile ( @$samfiles )
		{
			$counter++;

			my $bamfile = $samfile;
			$bamfile =~ s/sam$/bam/;
			my $command = "$samtools/samtools view -bt $samtools_index/$reference.fai -o $bamfile $samfile";
			push @$commands, $command;
		}

		my $label = "samToBam-$reference";
		my $outdir = "$outputdir/$reference";
		my $job = $self->setJob($commands, $label, $outdir);
		push @$jobs, $job;
	}

	#### RUN CONVERSION JOBS
	print "Agua::Cluster::Convert::subdirSamToBam    DOING runJobs for " , scalar(@$jobs), " jobs\n";
	$self->runJobs($jobs, 'subdirSamToBam');
	print "Agua::Cluster::Convert::subdirSamToBam    Completed samToBam\n";
}



=head2

	SUBROUTINE		bamToSam

	PURPOSE

		CONVERT EACH CHROMOSOME BAM FILE INTO A SAM FILE

=cut

sub bamToSam {
	my $self			=	shift;
	my $outputdir		=	shift;
	my $references		=	shift;
	my $infile			=	shift;
	my $outfile			=	shift;


	#### SET DEFAULT FILE NAMES IF NOT DEFINED
	$infile = "accepted_hits.bam" if not defined $infile;
	$outfile = "accepted_hits.sam" if not defined $outfile;

	#### GET SAMTOOLS
	my $samtools = $self->samtools();

	my $jobs = [];
	# CONVERT BAM TO SAM: samtools view -o out.sam in.bam 
	foreach my $reference ( @$references )
	{
		my $bamfile = "$outputdir/$reference/$infile";
		my $samfile = "$outputdir/$reference/$outfile";
		my $command = "$samtools/samtools view -o $samfile $bamfile";

		my $label = "bamToSam-$reference";
		my ($outputdir) = $bamfile =~ /^(.+)\/[^\/]+$/; 
		print "Agua::Cluster::Convert::bamToSam    label: $label\n";
		print "Agua::Cluster::Convert::bamToSam    outputdir: $outputdir\n";

		my $job = $self->setJob([$command], $label, $outputdir);
		push @$jobs, $job;
	}

	#### RUN CONVERSION JOBS
	print "Agua::Cluster::Convert::bamToSam    DOING bamToSam conversion...\n";
	$self->runJobs($jobs, "bamToSam");
	print "Agua::Cluster::Convert::bamToSam    FINISHED bamToSam conversion.\n";
}


=head2

	SUBROUTINE		samToBam

	PURPOSE

		CONVERT SAM FILES INTO (UNSORTED) BAM FILES

			samtools view -bt ref_list.txt -o aln.bam aln.sam.gz

=cut

sub samToBam {
	my $self			=	shift;
	my $outputdir		=	shift;
	my $references 		=	shift;
	my $infile			=	shift;
	my $outfile			=	shift;

	#### SET DEFAULT FILE NAMES IF NOT DEFINED
	$infile = "out.sam" if not defined $infile;
	$outfile = "out.bam" if not defined $outfile;

	#### SET SAMTOOLS INDEX
	my $samtools = $self->samtools();
	my $samtools_index = $self->samtoolsindex();

	my $jobs = [];
	foreach my $reference ( @$references )
	{
		#### GET SAM FILES
		my $samfile = "$outputdir/$reference/$infile";
		my $bamfile = "$outputdir/$reference/$outfile";

		my $command = "$samtools/samtools view -bt $samtools_index/$reference.fai -o $bamfile $samfile";
		my $label = "samToBam-$reference";
		my $outdir = "$outputdir/$reference";

		my $job = $self->setJob( [ $command ], $label, $outdir);
		push @$jobs, $job;
	}

	#### RUN JOBS
	$self->runJobs($jobs, "samToBam");
}

1;
