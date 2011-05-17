package SNP;


=head2

	PACKAGE		SNP

	PURPOSE

		SNP CALLING USING MULTIPLE *.sam HIT FILES FOR EACH REFERENCE CHROMOSOME

	VERSION		0.02

	HISTORY

		0.02 ADDED realignBam TO REALIGN READS AROUND INDELS AND MISMATCHES
		0.01 BASIC *bam FILE MERGING, SNP CALLING AND BINNING TO SPREAD WORKLOAD

=cut

use strict;
use warnings;
use Carp;

#### INTERNAL MODULES
use Sampler;

#### HAS A
use Monitor::PBS;
use Monitor::LSF;
use Cluster;
use UCSCBin;
#use LSF::Job;

#### EXTERNAL MODULES
use FindBin qw($Bin);
use Data::Dumper;
use File::Path;
#use MD5;

require Exporter;
our @ISA = qw(Exporter Cluster);
#our @EXPORT_OK = qw();
our $AUTOLOAD;


#### SET SLOTS
our @DATA = qw(

CHUNKSIZE
INPUTDIRS
OUTPUTDIR
BAMFILE
QUERYDIR
TARGETDIR
QUERYLABEL
TARGETLABEL
WINDOW
INPUTFILE
INPUTFILES
OUTPUTFILE
FILENAME

BINLEVEL
SUFFIX
QUERYINDEX
TARGETINDEX

BINDIR
MAXDEPTH	
ZIPPED
CLEAN
START
LABEL

PARAMS
REFERENCEDIR
SPECIES
GATK
JAVA
SAMTOOLS
SAMTOOLSINDEX

COMMAND
CLUSTER
QUEUE
WALLTIME
CPUS
QSTAT
QSUB
MAXJOBS
SLEEP
CLEANUP
VERBOSE
TEMPDIR
DOT
);
our $DATAHASH;
foreach my $key ( @DATA )	{	$DATAHASH->{lc($key)} = 1;	}

=head2

	SUBROUTINE		generateChunkInputfiles

	PURPOSE

		COLLATE ALL INPUT FILE GENERATION JOBS FOR ALL

		CHROMOSOMES AND INPUTDIRS THEN RUN TOGETHER

=cut
sub generateChunkInputfiles {
	my $self		=	shift;

	##### FILES AND DIRS	
	my $inputdirs 		=	$self->get_inputdirs();
	my $outputdir 		=	$self->get_outputdir();
	my $binlevel		=	$self->get_binlevel();
	my $bindir			=	$self->get_bindir();
	my $maxdepth 		=	$self->get_maxdepth();
	my $chunksize 		=	$self->get_chunksize();	
	my $filename 		=	$self->get_filename();
	my $samtools 		=	$self->get_samtools();
	my $samtools_index 	=	$self->get_samtoolsindex();
	my $start 			=	$self->get_start();


	#### GET CHROMOSOMES
	my $inputdir = $$inputdirs[0];
	my $chromosomes = $self->getChromosomes($inputdir);

	#### REPEAT FOR ALL *bam BIN FILES FOR EACH CHROMOSOME IN ALL INPUTDIRS
	my $outfile_chunks = {};
	my $jobs = [];
	my $counter = 0;
	$counter = $start - 1 if defined $start;
	foreach my $inputdir ( @$inputdirs )
	{
		$counter++;

		my ($dirname) = $inputdir =~ /([^\/]+)$/;

		#### DO ALL CHROMOSOMES
		foreach my $chromosome ( @$chromosomes )
		{
			#### INPUTFILE AND OUTDIR
			my $bamfile = "$inputdir/$chromosome/$filename";
			print "SNP::generateChunkInputfiles    Can't find bamfile: $bamfile\n" and exit if not -f $bamfile;

			#### CREATE BINNER AND GET CHROMOSOME SIZE
			my $binner = UCSCBin->new({	samtools => $samtools, outputdir => $outputdir	});
			my $chromosome_size = $binner->getChromosomeSize($bamfile);

			#### GET ARRAY OF BIN HASHES
			my $bins = $binner->getBins($binlevel, $chromosome_size);

			#### SET BINS BY HIT RANGE	
			my $rangefile = "$inputdir/$chromosome/$filename.range";
			$bins = $binner->setBinsByRange($bins, $rangefile);

			#### DO ALL BINS
			foreach my $bin ( @$bins )
			{
				#### SET OUTPUT DIR
				my $outdir = "$outputdir/$dirname/$chromosome";
				File::Path::mkpath($outdir) if not -d $outdir;
				print "SNP::generateChunkInputfiles     Can't create outdir: $outdir\n" if not -d $outdir;

				#### GET PREVIOUS AND INTERMEDIARY FILES
				my $binfile = $binner->setBinfile($filename, $outputdir, $binlevel, $bin);
				my ($filestub) = $binfile =~ /([^\/]+)$/;
				$filestub =~ s/\.bam$//;

				#### SET INPUT FILE			
				my $bamfile = "$outputdir/$chromosome/$filestub-$dirname.bam";
				print "SNP::generateChunkInputfiles    Skipping missing bamfile: $bamfile\n"
					and next if not -f $bamfile;

				#### SET SUBDIR FOR SAV AND STDOUT FILES AS BIN NUMBER
				#### I.E., POSITION ALONG THE CHROMOSOME
				my $binnumber = $bin->{number};
				my $binstart = $bin->{start};
				my $binstop = $bin->{stop};
				if ( not defined $binnumber )
				{
					exit;
				}
				my $chunkdir = "$outputdir/$chromosome/indir$dirname/bin$binnumber";

				#### PROCESS INPUT FILE INTO CHUNKS
				my ($chunkjobs, $chunkfiles) = $self->printChunkfileJobs(
					$chunksize,
					$bamfile,
					$chunkdir,
					{
						filestub	=>	$filestub,
						outputdir	=>	$outputdir,
						chromosome	=>	$chromosome,
						filename	=>	$filename,
						samtools	=>	$samtools,
						samtools_index	=>	$samtools_index,
						binnumber	=>	$binnumber,
						binstart	=>	$binstart,
						binstop		=>	$binstop,
						outdir		=>	$outdir,
						dirname		=>	$dirname
					}
				);
				@$jobs = (@$jobs, @$chunkjobs);

				#### SET OUTPUT FILE			
				my $snpfile = "$outputdir/$chromosome/$filestub.snp";
				$outfile_chunks->{$snpfile} = $chunkfiles;

				$counter++;
			}

		}	####	chromosomes

	}	####	inputdirs


	print "SNP::generateChunkInputfiles    No. jobs: ", scalar(@$jobs), "\n";
#exit;


	#### PRINT CHUNK FILE
	my $chunkfile = "$outputdir/binsnptosav.binlevel$binlevel.chunk$chunksize.json";
	use JSON;
	my $jsonparser = new JSON;
	my $json = $jsonparser->objToJson($outfile_chunks, {pretty => 1, indent => 2});
	print "SNP::binBamToSnp    json not defined!! Skipping print to chunkfile: $chunkfile\n";
	open(OUT, ">$chunkfile") or die "Can't open chunkfile: $chunkfile\n";
	print OUT $json;
	close(OUT) or die "Can't close chunkfile: $chunkfile\n";

	#### RUN JOBS
	print "\nSNP::generateChunkInputfiles    Running jobs\n";
	$self->runJobs( $jobs, "generateChunkInputfiles");
}

=head2

	SUBROUTINE		binBamToSnp

	PURPOSE

		CALL SNPS IN BINNED *bam FILE SUB-SECTIONS

=cut
sub binBamToSnp {
	my $self		=	shift;

	##### FILES AND DIRS	
	my $inputdirs 		=	$self->get_inputdirs();
	my $outputdir 		=	$self->get_outputdir();
	my $binlevel		=	$self->get_binlevel();
	my $bindir			=	$self->get_bindir();
	my $maxdepth 		=	$self->get_maxdepth();
	my $chunksize 		=	$self->get_chunksize();	
	my $filename 		=	$self->get_filename();
	my $samtools 		=	$self->get_samtools();
	my $samtools_index 	=	$self->get_samtoolsindex();
	my $start 			=	$self->get_start();


	#### GET CHROMOSOMES
	my $inputdir = $$inputdirs[0];
	my $chromosomes = $self->getChromosomes($inputdir);

	#### PRINT ALL RANGE FILES AT ONCE IF FIRST ONE IS MISSING	
	$self->printRangefiles();

	#### SORT ALL BAM FILES	
	$self->sortBam();

	#### INDEX ALL BAM FILES	
	$self->indexBam();

	#### PRINT CHUNK INPUT FILES FOR ALL CHROMOSOMES AND INPUTDIRS
	$self->generateChunkInputfiles();

	#### REPEAT FOR ALL *bam BIN FILES FOR EACH CHROMOSOME IN ALL INPUTDIRS
	my $outfile_chunks = {};
	my $jobs = [];
	my $counter = 0;
	$counter = $start - 1 if defined $start;
	foreach my $inputdir ( @$inputdirs )
	{
		$counter++;

		my ($dirname) = $inputdir =~ /([^\/]+)$/;

		#### DO ALL CHROMOSOMES
		foreach my $chromosome ( @$chromosomes )
		{
			#### INPUTFILE AND OUTDIR
			my $bamfile = "$inputdir/$chromosome/$filename";

			#### CHECK INPUT FILE		
			print "SNP::binBamToSnp    Can't find bamfile: $bamfile\n" and exit if not -f $bamfile;

			#### CREATE BINNER AND GET CHROMOSOME SIZE
			my $binner = UCSCBin->new({	samtools => $samtools, outputdir => $outputdir	});
			my $chromosome_size = $binner->getChromosomeSize($bamfile);

			#### GET ARRAY OF BIN HASHES
			my $bins = $binner->getBins($binlevel, $chromosome_size);

			#### SET BINS BY HIT RANGE	
			my $rangefile = "$inputdir/$chromosome/$filename.range";
			$bins = $binner->setBinsByRange($bins, $rangefile);

			#### DO ALL BINS
			foreach my $bin ( @$bins )
			{
				#### SET OUTPUT DIR
				my $outdir = "$outputdir/$dirname/$chromosome";
				File::Path::mkpath($outdir) if not -d $outdir;
				print "SNP::binBamToSnp     Can't create outdir: $outdir\n" if not -d $outdir;

				#### GET PREVIOUS AND INTERMEDIARY FILES
				my $binfile = $binner->setBinfile($filename, $outputdir, $binlevel, $bin);
				my ($filestub) = $binfile =~ /([^\/]+)$/;
				$filestub =~ s/\.bam$//;

				#### SET INPUT FILE			
				my $bamfile = "$outputdir/$chromosome/$filestub-$dirname.bam";
				print "SNP::binBamToSnp    Skipping missing bamfile: $bamfile\n"
					and next if not -f $bamfile;

				#### SET SUBDIR FOR SAV AND STDOUT FILES AS BIN NUMBER
				#### I.E., POSITION ALONG THE CHROMOSOME
				my $binnumber = $bin->{number};
				my $binstart = $bin->{start};
				my $binstop = $bin->{stop};
				if ( not defined $binnumber )
				{
					exit;
				}
				my $chunkdir = "$outputdir/$chromosome/indir$dirname/bin$binnumber";

				my $args = {
					filestub	=>	$filestub,
					outputdir	=>	$outputdir,
					chromosome	=>	$chromosome,
					maxdepth	=>	$maxdepth,
					filename	=>	$filename,
					samtools	=>	$samtools,
					samtools_index	=>	$samtools_index,
					binnumber	=>	$binnumber,
					binstart	=>	$binstart,
					binstop		=>	$binstop,
					outdir		=>	$outdir,
					dirname		=>	$dirname
				};

				#### PROCESS INPUT FILE INTO CHUNKS
				my $inputfiles = $self->getChunkfiles($chunksize, $bamfile, $chunkdir, $args);
				print "SNP::binBamToSnp    No. inputfiles: ", scalar(@$inputfiles), "\n";
				print "SNP::binBamToSnp    inputfiles:\n";
				print "SNP::binBamToSnp    inputfiles[0]: $$inputfiles[0]\n";
#exit;

				my $chunks = scalar(@$inputfiles);
				print "SNP::binBamToSnp    chunks: $chunks\n";

				#### PROCESS INPUT FILE INTO CHUNKS
				my ($chunkjobs, $chunkfiles) = $self->setChunkJobs(
					$chunksize,
					$bamfile,
					$chunkdir,
					$args
				);
				print "chunkjobs not defined. Skipping binnumber: $binnumber\n"
					and next if not defined $chunkjobs;
				print "chunkfiles not defined. Skipping binnumber: $binnumber\n"
					and next if not defined $chunkfiles;



				@$jobs = (@$jobs, @$chunkjobs);

				#### SET OUTPUT FILE			
				my $snpfile = "$outputdir/$chromosome/$filestub.snp";
				$outfile_chunks->{$snpfile} = $chunkfiles;


			}

		}	####	chromosomes

	}	####	inputdirs

	#### RUN JOBS
	print "\nSNP::binBamToSnp    Running jobs\n";
	$self->runJobs( $jobs, "binBamToSnp");

	#### PRINT CHUNK FILE
	my $chunkfile = "$outputdir/binsnptosav.binlevel$binlevel.chunk$chunksize.json";
	use JSON;
	my $jsonparser = new JSON;
	my $json = $jsonparser->objToJson($outfile_chunks, {pretty => 1, indent => 2});
	print "SNP::binBamToSnp    json not defined!! Skipping print to chunkfile: $chunkfile\n";
	open(OUT, ">$chunkfile") or die "Can't open chunkfile: $chunkfile\n";
	print OUT $json;
	close(OUT) or die "Can't close chunkfile: $chunkfile\n";

	#### MERGE ALL CHUNKS INTO SINGLE SNP FILE PER CHROMOSOME PER BIN
	$self->mergeChunks($outfile_chunks);

}

=head2

	SUBROUTINE		printChunkfileJobs

	PURPOSE

		SPLIT INPUT FILE INTO EQUAL CHUNKS

=cut
sub printChunkfileJobs {
	my $self		=	shift;
	my $chunksize	=	shift;
	my $inputfile	=	shift;
	my $outputdir	=	shift;
	my $args		=	shift;


	#### SET CHUNKFILE
	my ($filestub, $suffix) = $inputfile =~ /\/([^\/]+)\.([^\.]+)$/;

	#### GET SAMTOOLS	
	my $samtools = $self->get_samtools();
	print "UCSCBin::binCommand    samtools not defined\n" and exit if not defined $samtools;

	my $chromosome = $args->{chromosome};
	my $binstart = $args->{binstart};
	my $binstop = $args->{binstop};

	my $jobs = [];
	my $chunkfiles = [];
	my $chunknumber = 1;
	my $chunkfile = "$outputdir/$filestub.chunk$chunknumber.$suffix";
	my $start = 0;
	my $stop = 0;
	while ( $stop < $binstop )	
	{
		$start = $binstart + ($chunknumber * $chunksize);
		$stop = $binstart + (($chunknumber + 1) * $chunksize);
		$stop = $binstop if $stop > $binstop;
		my $region = "$chromosome:$start-$stop";
		my $commands = [];
		my $command = "$samtools/samtools view -h -b $inputfile $region > $chunkfile";
		push @$chunkfiles, $chunkfile;

		my $label = "$filestub.chunk$chunknumber";
		my $job = $self->setJob( [ $command ], $label, $outputdir);
		push @$jobs, $job;

		$chunknumber++;
		$chunkfile = "$outputdir/$filestub.chunk$chunknumber.$suffix";
	}

	##### RUN JOBS
	#$self->runJobs( $jobs, "printChunkfileJobs");

	return ($jobs, $chunkfiles);
}





=head2

	SUBROUTINE		binBamToSnp

	PURPOSE

		CALL SNPS IN BINNED *bam FILE SUB-SECTIONS


sub binBamToSnp {
	my $self		=	shift;
	my $inputdirs	=	shift;
	my $outputdir	=	shift;
	my $filename	=	shift;
	my $binlevel	=	shift;
	my $bindir		=	shift;
	my $maxdepth	=	shift;

	##### FILES AND DIRS	
	$inputdirs 		=	$self->get_inputdirs() if not defined $inputdirs;
	$outputdir 		=	$self->get_outputdir() if not defined $outputdir;
	$filename		=	$self->get_filename() if not defined $filename;
	$binlevel		=	$self->get_binlevel() if not defined $binlevel;
	$bindir			=	$self->get_bindir() if not defined $bindir;
	$maxdepth 		=	$self->get_maxdepth() if not defined $maxdepth;


	#### GET CHROMOSOMES
	my $inputdir = $$inputdirs[0];
	my $chromosomes = $self->getChromosomes($inputdir);

	######## PRINT ALL RANGE FILES AT ONCE IF FIRST ONE IS MISSING	
	####$self->printRangefiles();

	#### GET REQUIRED VARIABLES
	my $start 	=	$self->get_start();
	my $samtools = $self->get_samtools();
	my $samtools_index = $self->get_samtoolsindex();

	#### REPEAT FOR ALL *bam BIN FILES FOR EACH CHROMOSOME IN ALL INPUTDIRS
	my $jobs = [];
	my $counter = 0;
	$counter = $start - 1 if defined $start;
	foreach my $inputdir ( @$inputdirs )
	{
		$counter++;

		my ($dirname) = $inputdir =~ /([^\/]+)$/;

		#### DO ALL CHROMOSOMES
		foreach my $chromosome ( @$chromosomes )
		{
			#### INPUTFILE AND OUTDIR
			my $bamfile = "$inputdir/$chromosome/$filename";
			#my $outputdir = "$inputdir/$chromosome";

			#### CHECK INPUT FILE		
			print "SNP::binBamToSnp    Can't find bamfile: $bamfile\n" and exit if not -f $bamfile;

			##### CREATE CHROMOSOME OUTDIR
			#File::Path::mkpath($outputdir) if not -d $outputdir;

			#### CREATE BINNER AND GET CHROMOSOME SIZE
			my $binner = UCSCBin->new({	samtools => $samtools, outputdir => $outputdir	});
			my $chromosome_size = $binner->getChromosomeSize($bamfile);

			#### GET ARRAY OF BIN HASHES
			my $bins = $binner->getBins($binlevel, $chromosome_size);

			#### SET BINS BY HIT RANGE	
			my $rangefile = "$inputdir/$chromosome/$filename.range";
			$bins = $binner->setBinsByRange($bins, $rangefile);

			#### DO ALL BINS
			foreach my $bin ( @$bins )
			{
				my $commands = [];
				my $binfile = $binner->setBinfile($filename, $outputdir, $binlevel, $bin);

				#### GET PREVIOUS AND INTERMEDIARY FILES
				my ($filestub) = $binfile =~ /([^\/]+)$/;
				$filestub =~ s/\.bam$//;

				my $sortedfile_stub = "$outputdir/$filestub-$dirname";
				my $mergedfile = "$outputdir/$filestub-$dirname.bam";
				print "SNP::binBamToSnp    Skipping missing mergedfile: $mergedfile\n"
					and next if not -f $mergedfile;

				#### SET REFERENCE AND OUTPUT FILES
				my $rawfile = "$sortedfile_stub.snp.raw";
				my $filterfile = "$sortedfile_stub.snp.filter";
				my $snpfile = "$sortedfile_stub.snp";
				my $referencefile = $chromosome . ".fa";

				#### CHANGE TO OUTPUT DIR
				my $changedir = "cd $outputdir";
				push @$commands, $changedir;

				#### GET RAW VARIANT CALLS AND FILTER VARIANT CALLS
				my $predict = "$samtools/samtools pileup -vcf $samtools_index/$referencefile $mergedfile | tee $rawfile | $samtools/samtools.pl varFilter -D$maxdepth > $filterfile";
				push @$commands, $predict;

				#### QUALITY FILTER VARIANTS
				my $filter = qq{awk '(\$3=="*"&&\$6>=50)||(\$3!="*"&&\$6>=20)' $filterfile > $snpfile};
				push @$commands, $filter;
				my $main_command = join ";\n", @$commands;

				#### SET LABEL
				my $binnumber = $bin->{number};
				print "\nSNP::binBamToSnp    bin->{number} not defined\n" if not defined $bin->{number};
				print Dumper $bin and exit if not defined $bin->{number};
				my $label = "binBamToSnp-bin$binnumber-dir$dirname-$chromosome";

				##### RUN JOBS
				my $job = $self->setJob( [ $main_command ], $label, $outputdir);
				push @$jobs, $job;
			}

		}	####	chromosomes

	}	####	inputdirs

	#### RUN JOBS
	print "\nSNP::binBamToSnp    Running jobs\n";
	$self->runJobs( $jobs, "binBamToSnp");
}

=cut


=head2

	SUBROUTINE		getChunkfiles

	PURPOSE

		SPLIT INPUT FILE INTO EQUAL CHUNKS

=cut
sub getChunkfiles {
	my $self		=	shift;
	my $chunksize	=	shift;
	my $inputfile	=	shift;
	my $outputdir	=	shift;
	my $args		=	shift;


	#### SET CHUNKFILE
	my ($filestub, $suffix) = $inputfile =~ /\/([^\/]+)\.([^\.]+)$/;

	#### GET SAMTOOLS	
	my $samtools = $self->get_samtools();
	print "UCSCBin::binCommand    samtools not defined\n" and exit if not defined $samtools;

	my $chromosome = $args->{chromosome};
	my $binstart = $args->{binstart};
	my $binstop = $args->{binstop};

	my $jobs = [];
	my $chunkfiles = [];
	my $chunknumber = 1;
	my $chunkfile = "$outputdir/$filestub.chunk$chunknumber.$suffix";
	my $start = 0;
	my $stop = 0;
	while ( $stop < $binstop )	
	{
		$start = $binstart + ($chunknumber * $chunksize);
		$stop = $binstart + (($chunknumber + 1) * $chunksize);
		$stop = $binstop if $stop > $binstop;
		#my $region = "$chromosome:$start-$stop";
		#my $commands = [];
		#my $command = "$samtools/samtools view -h -b $inputfile $region > $chunkfile";
		push @$chunkfiles, $chunkfile;

		#my $label = "$filestub.chunk$chunknumber";
		#my $job = $self->setJob( [ $command ], $label, $outputdir);
		#push @$jobs, $job;

		$chunknumber++;
		$chunkfile = "$outputdir/$filestub.chunk$chunknumber.$suffix";
	}
#exit;


	##### RUN JOBS
	#$self->runJobs( $jobs, "getChunkfiles");

	return $chunkfiles;
}





=head2

	SUBROUTINE		getNumberChunks

	PURPOSE

		RETURN A LIST OF CHUNK FILES TO BE GENERATED TO SPLIT THE

		SNP FILE INTO chunksize-SIZED PORTIONS

sub getNumberChunks {
	my $self		=	shift;
	my $chunksize	=	shift;
	my $args		=	shift;

	my $bin		=	shift;


	my $region = $bin->{stop} - $bin->{start};	
	my $numberchunks = int($region/$chunksize);
	$numberchunks++ if $region % $chunksize != 0;

exit;

	return $numberchunks;
}
=cut

=head2

	SUBROUTINE		chunkCommand

	PURPOSE

		SET COMMAND FOR INDIVIDUAL CHUNK

=cut

sub chunkCommand {
	my $self		=	shift;
	my $chunknumber	=	shift;
	my $inputfile	=	shift;
	my $chunkdir	=	shift;
	my $args		=	shift;


	my $filename	=	$args->{filename};
	my $maxdepth	=	$args->{maxdepth};
	my $samtools	=	$args->{samtools};
	my $samtools_index	=	$args->{samtools_index};
	my $filestub	=	$args->{filestub};
	my $outputdir	=	$args->{outputdir};
	my $chromosome	=	$args->{chromosome};
	my $dirname 	=	$args->{dirname};
	my $binnumber	=	$args->{binnumber};
	my $outdir		=	$args->{outdir};


	#### SET OUTPUT FILE
	my $snpfile = "$chunkdir/$filestub-$dirname.chunk$chunknumber.snp";

	#### SET INTERMEDIARY FILES
	my $sortedfile_stub = "$chunkdir/$filestub-$dirname.chunk$chunknumber";
	my $bamfile = "$chunkdir/$filestub-$dirname.chunk$chunknumber.bam";
	print "SNP::chunkCommand    Skipping missing bamfile: $bamfile\n"
		and return if not -f $bamfile;

	#### SET REFERENCE AND OUTPUT FILES
	my $rawfile = "$sortedfile_stub.snp.raw";
	my $filterfile = "$sortedfile_stub.snp.filter";
	my $referencefile = $chromosome . ".fa";

	#### COLLECT COMMANDS
	my $commands = [];

	#### CHANGE TO OUTPUT DIR
	my $changedir = "cd $chunkdir";
	push @$commands, $changedir;

	my ($bamfile_stub) = $bamfile =~ /([^\/]+)$/;
	my $sort = "time $samtools/samtools sort $bamfile $bamfile_stub";
	push @$commands, $sort;

	my $index = "time $samtools/samtools index $bamfile";
	push @$commands, $index;

	#### GET RAW VARIANT CALLS AND FILTER VARIANT CALLS
	my $predict = "$samtools/samtools pileup -vcf $samtools_index/$referencefile $bamfile | tee $rawfile | $samtools/samtools.pl varFilter -D$maxdepth > $filterfile";


	push @$commands, $predict;

	#### QUALITY FILTER VARIANTS
	my $filter = qq{awk '(\$3=="*"&&\$6>=50)||(\$3!="*"&&\$6>=20)' $filterfile > $snpfile};
	push @$commands, $filter;

	##### JOIN COMMANDS INTO STRING
	#my $main_command = join ";\n", @$commands;

	#### SET LABEL
	my $label = "chunkCommand-$chromosome-indir$dirname-bin$binnumber-chunk$chunknumber";

	###### RUN JOBS
	#my $job = $self->setJob( [ $main_command ], $label, $outdir);

	##### RUN JOBS
	my $job = $self->setJob($commands, $label, $outdir);

	#### SET CHECKFILE
	$job->{checkfile} = $snpfile;




	return $job, $snpfile;
}



=head2

	SUBROUTINE		indexBam

	PURPOSE

		PRINT ALL RANGE FILES AT ONCE IF FIRST ONE IS MISSING

=cut

sub indexBam {
	my $self	=	shift;

	#### GET REQUIRED VARIABLES
	my $inputdirs	=	$self->get_inputdirs();
	my $outputdir	=	$self->get_outputdir();
	my $referencedir=	$self->get_referencedir();
	my $binlevel	=	$self->get_binlevel();
	my $filename	=	$self->get_filename();
	my $params		=	$self->get_params();
	my $samtools	=	$self->get_samtools();


	#### GET ALL REFERENCE CHROMOSOMES
	my $inputdir = $$inputdirs[0];
	print "SNP::indexBam    inputdir: $inputdir\n";
	my $counter = 0;
	my $chromosomes = $self->getChromosomes($inputdir);	
	foreach my $inputdir ( @$inputdirs )
	{
		$counter++;
		my ($dirname) = $inputdir =~ /([^\/]+)$/;

		#### COLLATE JOBS
		my $jobs = [];

		#### DO ALL CHROMOSOMES AT THE SAME START POINT
		foreach my $chromosome ( @$chromosomes )
		{
			#### INPUTFILE AND OUTDIR
			my $inputfile = "$inputdir/$chromosome/$filename";
			my $outdir = "$outputdir/$chromosome";

			#### CHECK INPUT FILE		
			print "SNP::indexBam    Can't find inputfile: $inputfile\n" and exit if not -f $inputfile;

			#### CREATE CHROMOSOME OUTDIR
			File::Path::mkpath($outdir) if not -d $outdir;
			print "Can't create outdir: $outdir\n" if not -d $outdir;

			#### SET BIN FILES DIRECTORY
			my $bindir = "$inputdir/$chromosome/bins";

			#### CREATE BINNER AND GET CHROMOSOME SIZE
			my $binner = UCSCBin->new({	samtools => $samtools, outputdir => $bindir	});
			my $chromosome_size = $binner->getChromosomeSize($inputfile);

			#### GET ARRAY OF BIN HASHES
			my $bins = $binner->getBins($binlevel, $chromosome_size);

			#### DO ALL BINS
			foreach my $bin ( @$bins )
			{
				#### GET PREVIOUS AND INTERMEDIARY FILES
				my $binfile = $binner->setBinfile($filename, $outputdir, $binlevel, $bin);
				my ($filestub) = $binfile =~ /([^\/]+)$/;
				$filestub =~ s/\.bam$//;

				#### SET LABEL
				my $binnumber = $bin->{number};
				print "\nSNP::indexBam    bin->{number} not defined\n" if not defined $bin->{number};
				print Dumper $bin and exit if not defined $bin->{number};
				my $label = "indexBam-$chromosome-num$binnumber-$counter";

				#### SET INPUT FILE			
				my $bamfile = "$outputdir/$chromosome/$filestub-$dirname.bam";
				print "SNP::binBamToSnp    Skipping missing bamfile: $bamfile\n"
					and next if not -f $bamfile;

				#### SORT COMMAND
				my $command = "time $samtools/samtools index $bamfile";
				my $job = $self->setJob( [$command], $label, $outdir);
				push @$jobs, $job;
			}
		}

		#### RUN JOBS
		print "\nSNP::indexBam    Running jobs\n";
		$self->runJobs( $jobs, "indexBam-$dirname");
	}

}



=head2

	SUBROUTINE		sortBam

	PURPOSE

		PRINT ALL RANGE FILES AT ONCE IF FIRST ONE IS MISSING

=cut

sub sortBam {
	my $self	=	shift;

	#### GET REQUIRED VARIABLES
	my $inputdirs	=	$self->get_inputdirs();
	my $outputdir	=	$self->get_outputdir();
	my $referencedir=	$self->get_referencedir();
	my $binlevel	=	$self->get_binlevel();
	my $filename	=	$self->get_filename();
	my $params		=	$self->get_params();
	my $samtools	=	$self->get_samtools();


	#### GET ALL REFERENCE CHROMOSOMES
	my $inputdir = $$inputdirs[0];
	print "SNP::sortBam    inputdir: $inputdir\n";
	my $counter = 0;
	my $chromosomes = $self->getChromosomes($inputdir);	
	foreach my $inputdir ( @$inputdirs )
	{
		$counter++;
		my ($dirname) = $inputdir =~ /([^\/]+)$/;

		#### COLLATE JOBS
		my $jobs = [];

		#### DO ALL CHROMOSOMES AT THE SAME START POINT
		foreach my $chromosome ( @$chromosomes )
		{
			#### INPUTFILE AND OUTDIR
			my $inputfile = "$inputdir/$chromosome/$filename";
			my $outdir = "$outputdir/$chromosome";

			#### CHECK INPUT FILE		
			print "SNP::sortBam    Can't find inputfile: $inputfile\n" and exit if not -f $inputfile;

			#### CREATE CHROMOSOME OUTDIR
			File::Path::mkpath($outdir) if not -d $outdir;
			print "Can't create outdir: $outdir\n" if not -d $outdir;

			#### SET BIN FILES DIRECTORY
			my $bindir = "$inputdir/$chromosome/bins";

			#### CREATE BINNER AND GET CHROMOSOME SIZE
			my $binner = UCSCBin->new({	samtools => $samtools, outputdir => $bindir	});
			my $chromosome_size = $binner->getChromosomeSize($inputfile);

			#### GET ARRAY OF BIN HASHES
			my $bins = $binner->getBins($binlevel, $chromosome_size);

			#### DO ALL BINS
			foreach my $bin ( @$bins )
			{
				#### GET PREVIOUS AND INTERMEDIARY FILES
				my $binfile = $binner->setBinfile($filename, $outputdir, $binlevel, $bin);
				my ($filestub) = $binfile =~ /([^\/]+)$/;
				$filestub =~ s/\.bam$//;

				#### SET LABEL
				my $binnumber = $bin->{number};
				print "\nSNP::sortBam    bin->{number} not defined\n" if not defined $bin->{number};
				print Dumper $bin and exit if not defined $bin->{number};
				my $label = "sortBam-$chromosome-num$binnumber-$counter";

				#### SET INPUT FILE			
				my $bamfile = "$outputdir/$chromosome/$filestub-$dirname.bam";
				print "SNP::binBamToSnp    Skipping missing bamfile: $bamfile\n"
					and next if not -f $bamfile;

				#### SORT COMMAND
				my ($bamfile_stub) = $bamfile =~ /([^\/]+)$/;
				my $command = "time $samtools/samtools sort $bamfile $bamfile_stub";
				my $job = $self->setJob( [$command], $label, $outdir);
				push @$jobs, $job;
			}
		}

		#### RUN JOBS
		print "\nSNP::sortBam    Running jobs\n";
		$self->runJobs( $jobs, "sortBam-$dirname");
	}

}


=head2

	SUBROUTINE		snpCounts

    PURPOSE

		FOR EVERY REFERENCE CHROMOSOME IN EACH INPUT DIRECTORY:

            1. COUNT NUMBER OF READ HITS IN *bam FILES USING samtools idxstats

            2. COUNT NUMBER OF RAW (*snp) SNPs

            3. COUNT NUMBER OF ANNOTATED (*sav) SNPs

            4. PRINT TO OUTPUT FILE

=cut



sub snpCounts {
	my $self	=	shift;

	my $inputdirs	=	$self->get_inputdirs();
	my $outputdir	=	$self->get_outputdir();
	my $outputfile	=	$self->get_outputfile();
	my $window		=	$self->get_window();
	my $referencedir=	$self->get_referencedir();
	my $binlevel	=	$self->get_binlevel();
	my $filename	=	$self->get_filename();
	my $samtools	=	$self->get_samtools();


	#### PRINT samtools idxstats TO FILE *bam.idxstats
	$self->printBamIdxstats();

	#### PRINT LINE COUNTS TO FILE
	$self->printLineCounts("snp");

	#### PRINT LINE COUNTS TO FILE
	$self->printLineCounts("sav");

	#### GET ALL REFERENCE CHROMOSOMES
	my $inputdir = $$inputdirs[0];
	my $chromosomes = $self->getChromosomes($inputdir);

	### GATHER COUNTS
	my $counts = [];
	foreach my $inputdir ( @$inputdirs )
	{
		foreach my $chromosome ( @$chromosomes )
		{	
			#### CREATE CHROMOSOME OUTDIR
			my $bamfile = "$inputdir/$chromosome/$filename";
			my ($dirname) = $inputdir =~ /([^\/]+)$/;
			my $outdir = "$outputdir/$chromosome";
			File::Path::mkpath($outdir) if not -d $outdir;
			print "Can't create outdir: $outdir\n" if not -d $outdir;

			#### SET BIN FILES DIRECTORY
			my $bindir = "$inputdir/$chromosome/bins";

			#### CREATE BINNER
			my $binner = UCSCBin->new({	samtools => $samtools, outputdir => $bindir	});

			#### GET CHROMOSOME SIZE
			my $chromosome_size = $binner->getChromosomeSize($bamfile);	

			#### GET ARRAY OF BINS
			my $bins = $binner->getBins($binlevel, $chromosome_size);

			#### GET RANGE FOR BINS
			my $rangefile = "$inputdir/$chromosome/hit.bam.range";

			#### SET BINS BY HIT RANGE	
			$bins = $binner->setBinsByRange($bins, $rangefile);

			#### ADD UP THE COUNTS FROM ALL OF THE BINS
			my $bincounts = [];
			foreach my $bin ( @$bins )
			{
				#### CHECK INPUT BAM FILE
				my $binfile = $binner->setBinfile($filename, $bindir, $binlevel, $bin);
				my $idxstatsfile = "$binfile.idxstats";
				my ($filestub) = $binfile =~ /([^\/]+)$/;
				$filestub =~ s/\.bam$//;
				my $snpfile = "$outdir/$filestub-$dirname.snp.lines";
				my $savfile = "$outdir/$filestub-$dirname.sav.lines";

				print "SNP::snpCounts    Can't find binfile: $binfile\n" and next if not -f $binfile;
				print "SNP::idxstatsCounts    Can't find idxstatsfile: $idxstatsfile\n" and next if not -f $idxstatsfile;
				print "SNP::snpCounts    Can't find snpfile: $snpfile\n" and next if not -f $snpfile;
				print "SNP::snpCounts    Can't find savfile: $savfile\n" and next if not -f $savfile;

				#### DEBUG

				my $hash;
				my $idstats = $self->fileContents($idxstatsfile);
				my ($hits, $misses) = $idstats =~ /^\S+\s+\S+\s+(\S+)\s+(\S+)/;
				$hash->{hits} = $hits;
				$hash->{misses} = $misses;
				$hash->{snps} = $self->fileContents($snpfile);
				$hash->{savs} = $self->fileContents($savfile);
				push @$bincounts, $hash;
			}

			#### COLLATE COUNTS FROM ALL BINS
			my ($hits, $misses, $snps, $savs);
			$hits = $misses = $snps = $savs = 0;
			foreach my $bincount ( @$bincounts )
			{
				$hits += $bincount->{hits};
				$misses += $bincount->{misses};
				$snps += $bincount->{snps};
				$savs += $bincount->{savs};
			}

			push @$counts, {
				inputdir 	=>	$inputdir,
				chromosome	=>	$chromosome,
				filename	=>	$filename,
				hits 		=>	$hits,
				misses 		=>	$misses,
				snps		=>	$snps,
				savs		=>	$savs
			}
		}
	}


	#### PRINT OUTPUT FILE
	my $report = '';
	my $total_hits = 0;
	foreach my $count ( @$counts )
	{
		$total_hits += $count->{hits};

		$report .= $count->{inputdir} 	. "\t";
		$report .= $count->{chromosome} . "\t";
		$report .= $count->{filename} 	. "\t";
		$report .= $count->{hits} 		. "\t";
		$report .= $total_hits	 		. "\t";
		$report .= $count->{misses} 	. "\t";
		$report .= $count->{snps} 		. "\t";
		$report .= $count->{savs} 		. "\n";
	}

	if ( defined $outputfile )
	{
		open(OUT, ">$outputfile") or die "Can't open outputfile: $outputfile\n";
		print OUT $report;
		close(OUT) or die "Can't close outputfile: $outputfile\n";
	}
	else
	{
		print "$report\n";
	}
}


=head2

	SUBROUTINE		bamToSnp

	PURPOSE

		CONVERT *bam FILE TO *snp FILE USING SAMTOOLS

=cut

sub bamToSnp {
	my $self			=	shift;


	#### FILES AND DIRS	
	my $inputdirs 		=	$self->get_inputdirs();
	my $outputdir 		=	$self->get_outputdir();
	my $maxdepth 		=	$self->get_maxdepth();
	my $binlevel		=	$self->get_binlevel();
	my $bindir			=	$self->get_bindir();
	my $filename		=	$self->get_filename();
	my $zipped			=	$self->get_zipped();
	my $stdout 			=	$self->get_stdout();

	print "SNP::bamToSnp    binlevel: $binlevel\n";

	#### PRINT TO STDOUT IF DEFINED stdout
	if ( defined $stdout )
	{
		print "bamToSnp.pl    Printing STDOUT to stdoutfile:\n\n$stdout\n\n";

		my ($stdout_path) = $stdout =~ /^(.+)\/[^\/]+$/;
		File::Path::mkpath($stdout_path) if not -d $stdout_path;
		open(STDOUT, ">$stdout") or die "Can't open STDOUT file: $stdout\n" if defined $stdout;
	}

	#### STARTING WITH THE FIRST DIRECTORY AND GOING IN ORDER
	if ( not defined $binlevel )
	{
		$self->bamSnps($inputdirs, $outputdir, $filename, $maxdepth, $zipped);
	}
	else
	{
		$self->binBamToSnp($inputdirs, $outputdir, $filename, $binlevel, $bindir, $maxdepth);
	}


}

=head2

	SUBROUTINE		bamSnps

	PURPOSE

		FOR A GIVEN CHROMOSOME chr0, FOR EACH inputdir/chr0 SUBDIRECTORY:

			1. CONVERT THE hit.sam FILE INTO A hit.bam FILE

			2. CONCAT THE FILE WITH THE EXISTING merged.bam FILE (IF ONE

				EXISTS) INTO A NEW merged.bam FILE

			3. COPY THE merged.bam FILE TO merged-N.bam

			4. CALL SNPS USING THE merged-N.bam FILE

			5. RETAIN THE merged.bam FILE FOR USE IN THE NEXT ITERATION

=cut

sub bamSnps {
	my $self		=	shift;

	my $inputdirs	=	shift;
	my $outputdir	=	shift;
	my $filename	=	shift;
	my $maxdepth	=	shift;
	my $zipped		=	shift;

	#### FILES AND DIRS	
	$inputdirs 		=	$self->get_inputdirs() if not defined $inputdirs;
	$outputdir 		=	$self->get_outputdir() if not defined $outputdir;
	$filename		=	$self->get_filename() if not defined $filename;
	$maxdepth		=	$self->get_maxdepth() if not defined $maxdepth;
	$zipped			=	$self->get_zipped() if not defined $zipped;

	#### GET CHROMOSOMES
	my $inputdir = $$inputdirs[0];
	my $chromosomes = $self->getChromosomes($inputdir);



	#### GET REQUIRED VARIABLES
	my $start 	=	$self->get_start();
	my $samtools = $self->get_samtools();
	my $samtools_index = $self->get_samtoolsindex();

	#### SET SLEEP IN BETWEEN JOBS
	my $sleep = $self->get_sleep();
	$sleep = 1500 if not defined $sleep;


	#### REPEAT FOR ALL INPUT DIRECTORIES	
	#### IF start IS DEFINED, RUN BEGINNING FROM merged.bam-<start> FILE
	my $jobs = [];
	my $counter = 0;
	$counter = $start - 1 if defined $start;
	foreach my $inputdir ( @$inputdirs )
	{
		$counter++;

		my ($dirname) = $inputdir =~ /([^\/]+)$/;

		foreach my $chromosome ( @$chromosomes )
		{
			#### INPUTFILE AND OUTDIR
			my $inputfile = "$inputdir/$chromosome/$filename";
			my $outdir = "$outputdir/$chromosome";

			#### CHECK INPUT FILE		
			print "SNP::bamSnps    Can't find inputfile: $inputfile\n" and exit if not -f $inputfile;

			#### CREATE CHROMOSOME OUTDIR
			File::Path::mkpath($outdir) if not -d $outdir;
			print "Can't create outdir: $outdir\n" if not -d $outdir;

			#### STORE COMMANDS HERE
			my $commands = [];

			#### UNZIP SAMFILE IF ZIPPED
			if ( $inputfile !~ /\.bam$/ )
			{
				if ( defined $zipped )
				{
					my $zippedfile = $inputfile;
					$inputfile =~ s/\.gz$// if $zipped eq "gzip";
					$inputfile =~ s/\.zip$// if $zipped eq "zip";
					$inputfile =~ s/\.bz2$// if $zipped eq "bz2";

					print "SNP::bamSnps    Can't find zippedfile: $zippedfile\n" and exit if not -f $zippedfile;

					my $unzip = "time zcat $zippedfile > $inputfile" if $zipped eq "gzip" or $zipped eq "zip";
					$unzip = "time bzip2 -dc $zippedfile > $inputfile" if $zipped eq "bz2";
					push @$commands, $unzip;
				}

				#### 1. CONVERT SAM TO BAM
				my $bamfile = $inputfile;
				$bamfile =~ s/\.sam$/.bam/;
				my $convert = "time $samtools/samtools view -bt $samtools_index/$chromosome.fai -o $bamfile $inputfile";
				push @$commands, $convert;
			}

			#### SET FILENAMES
			my $bamfile = $inputfile;
			$bamfile =~ s/\.sam$/.bam/;
			my $previous = $counter - 1;
			my $previousfile = "$outdir/merged-$previous.bam";
			my $tempfile = "$outdir/merged-$dirname.bam.temp";
			my $sortedfile_stub = "$outdir/merged-$dirname";
			my $mergedfile = "$outdir/merged-$dirname.bam";


			#### CALL SNPS USING THE mergedfile
			#### 
			my $rawfile = "$sortedfile_stub.snp.raw";
			my $filterfile = "$sortedfile_stub.snp.filter";
			my $snpfile = "$sortedfile_stub.snp";
			my $referencefile = $chromosome . ".fa";
			my $indexfile = $chromosome . ".fai";

			#### CHANGE TO OUTPUT DIR
			my $changedir = "cd $outdir";
			push @$commands, $changedir;

			#### GET RAW VARIANT CALLS AND FILTER VARIANT CALLS
			####
			####	RAW VARIANT CALLS USING samtools pileup
			####
			####			Usage:  samtools pileup [options] <in.bam>|<in.sam>
			####
			####Option: -s        simple (yet incomplete) pileup format
			####		-S        the input is in SAM
			####		-A        use the MAQ model for SNP calling
			####		-2        output the 2nd best call and quality
			####		-i        only show lines/consensus with indels
			####		-m INT    filtering reads with bits in INT [1796]
			####		-M INT    cap mapping quality at INT [60]
			####		-d INT    limit maximum depth for indels [unlimited]
			####		-t FILE   list of reference sequences (force -S)
			####		-l FILE   list of sites at which pileup is output
			####		-f FILE   reference sequence in the FASTA format
			####
			####		-c        output the SOAPsnp consensus sequence
			####		-v        print variants only (for -c)
			####		-g        output in the GLFv3 format (suppressing -c/-i/-s)
			####		-T FLOAT  theta in maq consensus calling model (for -c/-g) [0.850000]
			####		-N INT    number of haplotypes in the sample (for -c/-g) [2]
			####		-r FLOAT  prior of a difference between two haplotypes (for -c/-g) [0.001000]
			####		-G FLOAT  prior of an indel between two haplotypes (for -c/-g) [0.000150]
			####		-I INT    phred prob. of an indel in sequencing/prep. (for -c/-g) [40]

			####
			####	SNP Filtering
			####	Potential SNPs are filtered using the following SAMTOOLS default filtering criteria:
			####	
			####	Minimum mapping quality for SNPs [25]
			####	Minimum mapping quality for gaps [10]
			####	Minimum read depth [3]
			####	Maximum read depth [100]
			####	Min indel score for nearby SNP filtering [25]
			####	SNP within INT bp around a gap to be filtered [10]
			####	Window size for filtering dense SNPs [10]
			####	Max number of SNPs in a window [2]
			####	Window size for filtering adjacent gaps [30]


			my $predict = "$samtools/samtools pileup -vcf $samtools_index/$referencefile $mergedfile | tee $rawfile | $samtools/samtools.pl varFilter -D$maxdepth > $filterfile";
			push @$commands, $predict;

			#### QUALITY FILTER VARIANTS
			my $filter = qq{awk '(\$3=="*"&&\$6>=50)||(\$3!="*"&&\$6>=20)' $filterfile > $snpfile};
			push @$commands, $filter;

			#### SLEEP
			push @$commands, "sleep $sleep";	
			my $main_command = join ";\n", @$commands;

			#### SET LABEL
			my $label = "cumulativeSnp-$chromosome-$counter";

			##### RUN JOBS
			print "\nSNP::bamSnps    Running jobs:\n";
			my $job = $self->setJob( [ $main_command ], $label, $outdir);
			push @$jobs, $job;

		}	#### chromosomes

	} #### inputdirs

	#### RUN JOBS
	$self->runJobs( $jobs, "cumulativeSnps");	
}







=head2

	SUBROUTINE		fileContents

	PURPOSE

		1. COUNT THE NUMBER OF LINES IN EACH FILE AND PRINT TO FILE

		2. CREATE ONE FILE FOR EACH INPUT DIRECTORY AND CHROMOSOME

=cut

sub fileContents {
	my $self		=	shift;
	my $inputfile	=	shift;

	open(FILE, $inputfile) or die "SNP::fileContents    Can't open inputfile: $inputfile\n";
	my $contents = <FILE>;
	close(FILE) or die "SNP::fileContents    Can't close inputfile: $inputfile\n";
	$contents =~ s/\s+$//;

	return $contents;
}



=head2

	SUBROUTINE		printLineCounts

	PURPOSE

		1. COUNT THE NUMBER OF LINES IN EACH FILE AND PRINT TO FILE

		2. CREATE ONE FILE FOR EACH INPUT DIRECTORY AND CHROMOSOME

=cut

sub printLineCounts {
	my $self	=	shift;
	my $suffix	=	shift;

	my $inputdirs	=	$self->get_inputdirs();
	my $outputdir	=	$self->get_outputdir();
	my $window		=	$self->get_window();
	my $referencedir=	$self->get_referencedir();
	my $binlevel	=	$self->get_binlevel();
	my $filename	=	$self->get_filename();
	my $samtools	=	$self->get_samtools();


	#### GET ALL REFERENCE CHROMOSOMES
	my $inputdir = $$inputdirs[0];
	my $chromosomes = $self->getChromosomes($inputdir);

	#### SET EXECUTABLE
	my $executable = "$Bin/../utils/countLines.pl";

	#### PRINT ALL RANGE FILES AT ONCE IF FIRST ONE IS MISSING	
	$self->printRangefiles();

	#### REPEAT FOR ALL *bam BIN FILES FOR EACH INPUT DIR AND CHROMOSOME
	my $jobs = [];
	my $counter = 0;
	foreach my $inputdir ( @$inputdirs )
	{
		$counter++;

		my ($dirname) = $inputdir =~ /([^\/]+)$/;

		#### DO ALL CHROMOSOMES AT THE SAME START POINT
		foreach my $chromosome ( @$chromosomes )
		{
			#### CREATE CHROMOSOME OUTDIR
			my $outdir = "$outputdir/$chromosome";
			File::Path::mkpath($outdir) if not -d $outdir;
			print "Can't create outdir: $outdir\n" if not -d $outdir;

			#### CHECK INPUT BAM FILE
			my $bamfile = "$inputdir/$chromosome/$filename";
			print "SNP::printLineCounts    Can't find bamfile: $bamfile\n" and exit if not -f $bamfile;

			#### SET BIN FILES DIRECTORY
			my $bindir = "$inputdir/$chromosome/bins";

			#### CREATE BINNER AND GET CHROMOSOME SIZE
			my $binner = UCSCBin->new({	samtools => $samtools, outputdir => $bindir	});
			my $chromosome_size = $binner->getChromosomeSize($bamfile);

			#### GET ARRAY OF BIN HASHES
			my $bins = $binner->getBins($binlevel, $chromosome_size);

			#### SET BINS BY HIT RANGE	
			my $rangefile = "$inputdir/$chromosome/hit.bam.range";
			$bins = $binner->setBinsByRange($bins, $rangefile);

			#### DO ALL BINS
			foreach my $bin ( @$bins )
			{
				my $commands = [];

				#### CHECK INPUT BINNED BAM FILE
				my $binfile = $binner->setBinfile($filename, $bindir, $binlevel, $bin);
				print "SNP::snpCounts    Can't find binfile: $binfile\n" and exit if not -f $binfile;

				#### SET INPUT AND OUTPUT FILES
				my ($filestub) = $binfile =~ /([^\/]+)$/;
				$filestub =~ s/\.bam$//;
				my $infile = "$outdir/$filestub-$dirname.$suffix";
				my $outfile = "$infile.lines";

				#### SET COMMAND
				my $command = qq{$executable \\\n};
				$command .= qq{ --inputfile $infile \\\n};
				$command .= qq{ --outputfile $outfile\n};

				#### SET LABEL
				my $label = "printLineCounts-$dirname-$chromosome-$bin->{number}";

				##### SET OUTDIR
				#my $outdir = "$inputdir/$chromosome/idxstats";

				##### COLLECT JOB
				my $job = $self->setJob( [ $command ], $label, $outputdir);
				push @$jobs, $job;
			}

		}	####	chromosomes

	}	####	inputdirs

	#### RUN JOBS
	print "\nSNP::printLineCounts    Running " , scalar(@$jobs), " jobs\n";	
	$self->runJobs( $jobs, "printLineCounts");
}


=head2

	SUBROUTINE		printBamIdxstats

	PURPOSE

		1. USE bamIdxstats.pl TO PRINT samtools idxstats INFORMATION TO FILE

		2. CREATE ONE FILE FOR EACH INPUT DIRECTORY AND CHROMOSOME

=cut

sub printBamIdxstats {
	my $self	=	shift;

	my $inputdirs	=	$self->get_inputdirs();
	my $outputdir	=	$self->get_outputdir();
	my $window		=	$self->get_window();
	my $referencedir=	$self->get_referencedir();
	my $binlevel	=	$self->get_binlevel();
	my $filename	=	$self->get_filename();
	my $samtools	=	$self->get_samtools();


	#### GET ALL REFERENCE CHROMOSOMES
	my $inputdir = $$inputdirs[0];
	my $chromosomes = $self->getChromosomes($inputdir);

	#### SET EXECUTABLE
	my $executable = "$Bin/../aligners/bamIdxstats.pl";

	#### PRINT ALL RANGE FILES AT ONCE IF FIRST ONE IS MISSING	
	$self->printRangefiles();

	#### REPEAT FOR ALL *bam BIN FILES FOR EACH INPUT DIR AND CHROMOSOME
	my $jobs = [];
	my $counter = 0;
#	foreach my $inputdir ( @$inputdirs )
#	{
#		$counter++;
#		
#		my ($dirname) = $inputdir =~ /([^\/]+)$/;
#	
#		foreach my $chromosome ( @$chromosomes )
#		{
#			
#			#### SET COMMAND
#			my $bamfile = "$outputdir/$chromosome/$filename";
#			$bamfile s/\.bam$/-$dirname.bam/;
#			my $outfile = "$bamfile.idxstats";
#			my $command = qq{$executable \\\n};
#			$command .= qq{ --inputfile $bamfile \\\n};
#			$command .= qq{ --outputfile $outfile\n};
#exit;
#			#### SET LABEL
#			my $label = "printBamIdxstats-$chromosome-$counter";
#	
#			##### SET OUTDIR
#			#my $outdir = "$inputdir/$chromosome/idxstats";
#
#			##### COLLECT JOB
#			my $job = $self->setJob( [ $command ], $label, $outputdir);
#			push @$jobs, $job;
#		}
#	}


	foreach my $inputdir ( @$inputdirs )
	{
		$counter++;

		#### DO ALL CHROMOSOMES AT THE SAME START POINT
		foreach my $chromosome ( @$chromosomes )
		{
			#### SET OUTDIR
			my $outdir = "$outputdir/$chromosome";

			#### CHECK INPUT BAM FILE
			my $bamfile = "$inputdir/$chromosome/$filename";
			print "SNP::unbinSnp    Can't find bamfile: $bamfile\n" and exit if not -f $bamfile;

			#### CREATE CHROMOSOME OUTDIR
			File::Path::mkpath($outdir) if not -d $outdir;
			print "Can't create outdir: $outdir\n" if not -d $outdir;

			#### SET BIN FILES DIRECTORY
			my $bindir = "$inputdir/$chromosome/bins";

			#### CREATE BINNER AND GET CHROMOSOME SIZE
			my $binner = UCSCBin->new({	samtools => $samtools, outputdir => $bindir	});
			my $chromosome_size = $binner->getChromosomeSize($bamfile);

			#### GET ARRAY OF BIN HASHES
			my $bins = $binner->getBins($binlevel, $chromosome_size);

			#### SET BINS BY HIT RANGE	
			my $rangefile = "$inputdir/$chromosome/hit.bam.range";
			$bins = $binner->setBinsByRange($bins, $rangefile);

			#### DO ALL BINS
			foreach my $bin ( @$bins )
			{
				my $commands = [];
				my $binfile = $binner->setBinfile($filename, $bindir, $binlevel, $bin);

				#### SET COMMAND
				my $outfile = "$binfile.idxstats";
				my $command = qq{$executable \\\n};
				$command .= qq{ --inputfile $binfile \\\n};
				$command .= qq{ --outputfile $outfile\n};

				#### SET LABEL
				my $label = "printBamIdxstats-$chromosome-$counter-$bin->{number}";

				##### SET OUTDIR
				#my $outdir = "$inputdir/$chromosome/idxstats";

				##### COLLECT JOB
				my $job = $self->setJob( [ $command ], $label, $outputdir);
				push @$jobs, $job;
			}

		}	####	chromosomes

	}	####	inputdirs


	#### RUN JOBS
	print "\nSNP::printBamIdxstats    Running " , scalar(@$jobs), " jobs\n";	
	$self->runJobs( $jobs, "createTargets");
}





=head2

	SUBROUTINE		createIntervals

	PURPOSE

		1. CREATE INDEL INTERVALS FILE SIMILAR TO THE GATK RealignerTargetCreator

			OUTPUT FILE, FOR INPUT INTO GATK IndelRealigner

		2. CREATE INTERVAL FOR ALL SNPS AND INDEL REGIONS

		3. INCLUDE window NUMBER OF BASES TO LEFT AND RIGHT OF SNP OR INDEL		

	NOTES

		AN EXAMPLE OF 10-COLUMN pileup FORMAT:

			1 		2 				3 		4 		5 		6 		7 		8 		9		
		10
			chr22   41948296        a       C       41      41      39      8       CcCCCcc^GC      #+!)"!++      
	chr22   41953978        G       C       38      38      42      8       Cc.CccCc        !+"!++++      
	chr22   41963553        C       G       72      72      40      24      GgGgggGGGgGgGggGgGGggGGg      !+!+++#!!+!"++$+!&+!!++!

			Column  Definition
		1   Chromosome
		2   Position (1-based)
		3   Reference base at that position
		4   Consensus bases
		5   Consensus quality
		6   SNP quality
		7   Maximum mapping quality
		8   Coverage (# reads aligning over that position)
		9   Read bases
	   10   Base quality values (phred33 scale)

=cut

sub createIntervals {
	my $self	=	shift;
	my $inputfile	=	shift;
	my $outputfile	=	shift;
	my $window		=	shift;


	#### OPEN FILES
	open(FILE, $inputfile) or die "SNP::createIntervals    Can't open inputfile: $inputfile\n";
	open(OUT, ">$outputfile") or die "Can't open outputfile: $outputfile\n";

	my $line = <FILE>;
	my ($chromosome, $currentstart) = $line =~ /^(\S+)\s+(\S+)\s+\S+\s(\S+)/;
	my $currentstop = $currentstart + $window;

	my $counter = 0;
	my $dot = 100000;
	my $merging = 0;
	while ( <FILE> )
	{
		$counter++;
		print "$counter\n" if $counter % $dot == 0;

		my ($position) = $_ =~ /^\S+\s+(\S+)/;

		#### IF OVERLAPS, INCLUDE IT IN THE CURRENT INTERVAL
		if ( ($position - $window) < $currentstop )
		{
			$currentstop = $position + $window;
			$merging = 1;
		}
		#### OTHERWISE, PRINT TO FILE AND START A NEW INTERVAL
		else
		{
			print OUT "$chromosome:$currentstart-$currentstop\n";

			$currentstart = $position - $window;
			$currentstop = $position + $window;
			$merging = 0;
		}

		#last if $counter == 10;
	}

	if ( $merging )
	{
		print OUT "$chromosome:$currentstart-$currentstop\n";
	}

	close(OUT) or die "Can't close outputfile: $outputfile\n";
}


=head2

	SUBROUTINE		createGatkIntervals

	PURPOSE

		1. USE GATK RealignerTargetCreator TO CREATE AN INDEL INTERVALS FILE

			FOR INPUT INTO GATK IndelRealigner

=cut

sub createGatkIntervals {
	my $self	=	shift;

	my $inputdirs	=	$self->get_inputdirs();
	my $outputdir	=	$self->get_outputdir();
	my $referencedir=	$self->get_referencedir();
	my $binlevel	=	$self->get_binlevel();
	my $filename	=	$self->get_filename();
	my $params		=	$self->get_params();


	#### GET REQUIRED VARIABLES
	my $java 		=	$self->get_java();
	my $gatk		=	$self->get_gatk();
	my $samtools	=	$self->get_samtools();

	#### GET ALL REFERENCE CHROMOSOMES
	my $inputdir = $$inputdirs[0];
	my $chromosomes = $self->getChromosomes($inputdir);

	#### CREATE INDEL INTERVALS FILE FOR BIN FILE PER INPUT DIR AND CHROMOSOME
	#### Print intervals for GATK IndelRealigner to target for cleaning
	my $counter = 0;
	my $jobs = [];
	foreach my $inputdir ( @$inputdirs )
	{
		$counter++;

		foreach my $chromosome ( @$chromosomes )
		{	
			#### CREATE CHROMOSOME OUTDIR
			my $bamfile = "$inputdir/$chromosome/$filename";
			my ($dirname) = $inputdir =~ /([^\/]+)$/;
			my $outdir = "$outputdir/$dirname/$chromosome";
			File::Path::mkpath($outdir) if not -d $outdir;
			print "Can't create outdir: $outdir\n" if not -d $outdir;

			#### SET REFERENCE FILE
			my $referencefile = "$referencedir/$chromosome.fa";
			print "Can't find referencefile: $referencefile\n" if not -f $referencefile;

			#### SET BIN FILES DIRECTORY
			my $bindir = "$inputdir/$chromosome/bins";

			#### CREATE BINNER
			my $binner = UCSCBin->new({	samtools => $samtools, outputdir => $bindir	});

			#### GET CHROMOSOME SIZE
			my $chromosome_size = $binner->getChromosomeSize($bamfile);	

			#### GET ARRAY OF BINS
			my $bins = $binner->getBins($binlevel, $chromosome_size);

			#### GET RANGE FOR BINS
			my $rangefile = "$inputdir/$chromosome/hit.bam.range";

			#### SET BINS BY HIT RANGE	
			$bins = $binner->setBinsByRange($bins, $rangefile);

			#### DO EACH BAM BIN FILE
			foreach my $bin ( @$bins )
			{
				#### CHECK INPUT BAM FILE
				my $binfile = $binner->setBinfile($filename, $bindir, $binlevel, $bin);
				print "binfile: $binfile\n";
				print "SNP::createGatkIntervals    Can't find binfile: $binfile\n" and exit if not -f $binfile;

				#### SET OUTPUT FILES
				my ($filename) = $binfile =~ /([^\/]+)$/;
				my $targetsfile = "$outdir/$filename.intervals";
				my $outputfile = "$outdir/$filename";

				######## EXAMPLE:
				####
				#### Emit intervals for the Local Indel Realigner to target for cleaning.
				####
				#### time /usr/local/java/bin/java \
				#### -jar /nethome/syoung/base/apps/gatk/1.0.4705/GenomeAnalysisTK.jar \
				#### -I aln.bam \
				#### -R $REF \
				#### -T RealignerTargetCreator \
				#### -U \
				#### -S SILENT \
				#### -o realignTargets.intervals
				my $command = qq{time $java/bin/java \\\n};
				$command .= qq{-jar $gatk/GenomeAnalysisTK.jar \\\n};
				$command .= qq{-T RealignerTargetCreator \\\n};
				$command .= qq{-I $binfile \\\n};
				$command .= qq{-R $referencefile \\\n};
				$command .= qq{-o $targetsfile \\\n};
				$command .= qq{-U \\\n};
				$command .= qq{-S SILENT \\\n};
				$command .= $params if defined $params;
				$command =~ s/\\+$//;

				print "command: $command\n";

				#### SET LABEL
				my $binnumber = $bin->{number};
				print "\nSNP::createGatkIntervals    bin->{number} not defined\n" if not defined $bin->{number};
				print Dumper $bin and exit if not defined $bin->{number};
				my $label = "realignBam-$chromosome-$counter-$binnumber";

				##### RUN JOBS
				my $job = $self->setJob( [ $command ], $label, $outdir);
				push @$jobs, $job;
			}
		}
	}

	#### RUN JOBS
	print "\nSNP::createGatkIntervals    Running " , scalar(@$jobs), " jobs\n";	
	$self->runJobs( $jobs, "createGatkIntervals");
}





=head2

	SUBROUTINE		createIntervalFiles

	PURPOSE

		1. CREATE INDEL INTERVALS FILES

			Generate an interval file conforming to the GATK RealignerTargetCreator

			output file, suitable for input into the GATK IndelRealigner

		2. CREATE ONE INTERVALS FILE FOR EACH BINNED SNP FILE, PER INPUT DIR,

			PER CHROMOSOME

=cut

sub createIntervalFiles {
	my $self	=	shift;

	my $inputdirs	=	$self->get_inputdirs();
	my $outputdir	=	$self->get_outputdir();
	my $window		=	$self->get_window();
	my $referencedir=	$self->get_referencedir();
	my $binlevel	=	$self->get_binlevel();
	my $filename	=	$self->get_filename();
	my $samtools	=	$self->get_samtools();


	#### GET ALL REFERENCE CHROMOSOMES
	my $inputdir = $$inputdirs[0];
	my $chromosomes = $self->getChromosomes($inputdir);

	#### SET EXECUTABLE
	my $executable = "$Bin/createIntervals.pl";

	### 1. CREATE INDEL INTERVALS FILE
	### Emit intervals for the Local Indel Realigner to target for cleaning
	my $counter = 0;
	my $jobs = [];
	foreach my $inputdir ( @$inputdirs )
	{
		$counter++;

		foreach my $chromosome ( @$chromosomes )
		{	
			#### CREATE CHROMOSOME OUTDIR
			my $bamfile = "$inputdir/$chromosome/$filename";
			my ($dirname) = $inputdir =~ /([^\/]+)$/;
			my $outdir = "$outputdir/$dirname/$chromosome/realign";
			File::Path::mkpath($outdir) if not -d $outdir;
			print "Can't create outdir: $outdir\n" if not -d $outdir;

			#### SET REFERENCE FILE
			my $referencefile = "$referencedir/$chromosome.fa";
			print "Can't find referencefile: $referencefile\n" if not -f $referencefile;

			#### SET BIN FILES DIRECTORY
			my $bindir = "$inputdir/$chromosome/bins";

			#### CREATE BINNER
			my $binner = UCSCBin->new({	samtools => $samtools, outputdir => $bindir	});

			#### GET CHROMOSOME SIZE
			my $chromosome_size = $binner->getChromosomeSize($bamfile);	

			#### GET ARRAY OF BINS
			my $bins = $binner->getBins($binlevel, $chromosome_size);

			#### GET RANGE FOR BINS
			my $rangefile = "$inputdir/$chromosome/hit.bam.range";

			#### SET BINS BY HIT RANGE	
			$bins = $binner->setBinsByRange($bins, $rangefile);

			#### DO EACH BAM BIN FILE
			foreach my $bin ( @$bins )
			{
				#### CHECK INPUT BAM FILE
				my $binfile = $binner->setBinfile($filename, $bindir, $binlevel, $bin);
				print "SNP::createIntervalFiles    Can't find binfile: $binfile\n" and exit if not -f $binfile;

				#### SET OUTPUT FILES
				my ($filename) = $binfile =~ /([^\/]+)$/;
				my $intervalsfile = "$outdir/$filename.intervals";

				#### SET SNP FILE
				my $snpfile = "$outputdir/$chromosome/$filename";
				$snpfile =~ s/\.bam$//;
				$snpfile .= "-$dirname.snp";

				######## EXAMPLE:
				####
				#### Emit intervals for the Local Indel Realigner to target for cleaning.
				####
				my $command = qq{$executable \\\n};
				$command .= qq{ --inputfile $snpfile \\\n};
				$command .= qq{ --outputfile $intervalsfile \\\n};
				$command .= qq{ --window $window \\\n} if defined $window;
				$command =~ s/\\$//;

				#### SET LABEL
				my $binnumber = $bin->{number};
				print "\nSNP::createIntervalFiles    bin->{number} not defined\n" if not defined $bin->{number};
				print Dumper $bin and exit if not defined $bin->{number};
				my $label = "createIntervalFiles-$chromosome-$counter-$binnumber";

				##### RUN JOBS
				my $job = $self->setJob( [ $command ], $label, $outdir);
				push @$jobs, $job;
			}
		}
	}

	#### RUN JOBS
	print "\nSNP::createIntervalFiles    Running " , scalar(@$jobs), " jobs\n";	
	$self->runJobs( $jobs, "createTargets");
}







=head2

	SUBROUTINE		printRangefiles

	PURPOSE

		PRINT ALL RANGE FILES AT ONCE IF FIRST ONE IS MISSING

=cut

sub printRangefiles {
	my $self	=	shift;

	my $inputdirs	=	$self->get_inputdirs();
	my $outputdir	=	$self->get_outputdir();
	my $referencedir=	$self->get_referencedir();
	my $binlevel	=	$self->get_binlevel();
	my $filename	=	$self->get_filename();
	my $params		=	$self->get_params();


	#### GET REQUIRED VARIABLES
	my $java 		=	$self->get_java();
	my $gatk		=	$self->get_gatk();
	my $samtools	=	$self->get_samtools();

	#### GET ALL REFERENCE CHROMOSOMES
	my $inputdir = $$inputdirs[0];
	print "SNP::printRangefiles    inputdir: $inputdir\n";
	my $chromosomes = $self->getChromosomes($inputdir);

	##### PRINT ALL RANGE FILES AT ONCE IF FIRST ONE IS MISSING
	#my $chromosome = $$chromosomes[0];
	#my $bamfile = "$inputdir/$chromosome/$filename";
	#my $rangefile = "$bamfile.range";

	#if ( not -f $rangefile )
	#{

		#### CREATE BINNER
		my $maxjobs = 
		my $binner = UCSCBin->new({
			samtools=> $samtools,
			maxjobs => 	$self->get_maxjobs(),
			queue 	=> 	$self->get_queue(),
			cluster	=>	$self->get_cluster()
		});
		$binner->printRangefiles($inputdirs, $chromosomes);
	#}
}


=head2

	SUBROUTINE		realignBam

	PURPOSE

		REALIGN *bam FILE READS AROUND INDEL OR MISMATCH REGIONS AND

		OUTPUT TO A 'CLEANED' *bam FILE

	NOTES

		THE PROCESS CONSISTS OF TWO STEPS:

			1. CREATE INDEL INTERVALS FILE

				(CARRIED OUT BY THE createIntervals METHOD IN THIS MODULE)

			Emit intervals for the Local Indel Realigner to target for cleaning.

			time /usr/local/java/bin/java \
			-jar /nethome/syoung/base/apps/gatk/1.0.4705/GenomeAnalysisTK.jar \
			-I aln.bam \
			-R $REF \
			-T RealignerTargetCreator \
			-U \
			-S SILENT \
			-o realignTargets.intervals

			2. REALIGN READS IN MERGED INTERVALS (CARRIED OUT BY THIS METHOD)
			Perform local realignment of reads based on misalignments due to the presence of indels. Unlike most mappers, this  walker uses the full alignment context to determine whether an appropriate alternate reference (i.e. indel) exists and updates SAMRecords accordingly.

			time /usr/local/java/bin/java \
			-jar /nethome/syoung/base/apps/gatk/1.0.4705/GenomeAnalysisTK.jar \
			-I aln.bam \
			-R $REF \
			-T IndelRealigner \
			-targetIntervals realignTargets.intervals \
			--out realigned.bam \
			-U \
			-S SILENT


Arguments for IndelRealigner:

 -targetIntervals,--targetIntervals <targetIntervals>          intervals file output from RealignerTargetCreator
 -LOD,--LODThresholdForCleaning <LODThresholdForCleaning>      LOD threshold above which the cleaner will clean
 -entropy,--entropyThreshold <entropyThreshold>                percentage of mismatches at a locus to be considered 
                                                               having high entropy
 -o,--out <out>                                                Output bam
 -compress,--bam_compression <bam_compression>                 Compression level to use for writing BAM files
 --index_output_bam_on_the_fly <index_output_bam_on_the_fly>   Create a BAM index on-the-fly while writing the resulting 
                                                               file.
 -knownsOnly,--useOnlyKnownIndels                              Don't run 'Smith-Waterman' to generate alternate 
                                                               consenses; use only known indels provided as RODs for 
                                                               constructing the alternate references.
 -maxInRam,--maxReadsInRam <maxReadsInRam>                     max reads allowed to be kept in memory at a time by the 
                                                               SAMFileWriter. If too low, the tool may run out of system 
                                                               file descriptors needed to perform sorting; if too high, 
                                                               the tool may run out of memory.
 -maxConsensuses,--maxConsensuses <maxConsensuses>             max alternate consensuses to try (necessary to improve 
                                                               performance in deep coverage)
 -greedy,--maxReadsForConsensuses <maxReadsForConsensuses>     max reads used for finding the alternate consensuses 
                                                               (necessary to improve performance in deep coverage)
 -maxReads,--maxReadsForRealignment <maxReadsForRealignment>   max reads allowed at an interval for realignment; if this 
                                                               value is exceeded, realignment is not attempted and the 
                                                               reads are passed to the output file(s) as-is
 --sortInCoordinateOrderEvenThoughItIsHighlyUnsafe             Should we sort the final bam in coordinate order even 
                                                               though it will be malformed because mate pairs of 
                                                               realigned reads will contain inaccurate information?
 --realignReadsWithBadMates                                    Should we try to realign paired-end reads whose mates map 
                                                               to other chromosomes?
 -noPG,--noPGTag                                               Don't output the usual PG tag in the realigned bam file 
                                                               header. FOR DEBUGGING PURPOSES ONLY. This option is 
                                                               required in order to pass integration tests.
 -noTags,--noOriginalAlignmentTags                             Don't output the original cigar or alignment start tags 
                                                               for each realigned read in the output bam.
 -targetNotSorted,--targetIntervalsAreNotSorted                This tool assumes that the target interval list is 
                                                               sorted; if the list turns out to be unsorted, it will 
                                                               throw an exception.  Use this argument when your interval 
                                                               list is not sorted to instruct the Realigner to first 
                                                               sort it in memory.


=cut

sub realignBam {	
	my $self	=	shift;

	my $inputdirs	=	$self->get_inputdirs();
	my $outputdir	=	$self->get_outputdir();
	my $referencedir=	$self->get_referencedir();
	my $binlevel	=	$self->get_binlevel();
	my $filename	=	$self->get_filename();
	my $params		=	$self->get_params();


	#### GET REQUIRED VARIABLES
	my $java 		=	$self->get_java();
	my $gatk		=	$self->get_gatk();
	my $samtools	=	$self->get_samtools();

	#### GET ALL REFERENCE CHROMOSOMES
	my $inputdir = $$inputdirs[0];
	print "SNP::realignBam    inputdir: $inputdir\n";
	my $chromosomes = $self->getChromosomes($inputdir);

	#### PRINT ALL RANGE FILES AT ONCE IF FIRST ONE IS MISSING
	$self->printRangefiles();

	#### REALIGN READS OVERLAPPING PRECALCULATED INTERVALS
	my $counter = 0;
	my $jobs = [];
	foreach my $inputdir ( @$inputdirs )
	{
		$counter++;

		my ($dirname) = $inputdir =~ /([^\/]+)$/;

		foreach my $chromosome ( @$chromosomes )
		{	
			#### CREATE CHROMOSOME OUTDIR
			my $bamfile = "$inputdir/$chromosome/$filename";
			my ($dirname) = $inputdir =~ /([^\/]+)$/;
			my $outdir = "$outputdir/$dirname/$chromosome";
			File::Path::mkpath($outdir) if not -d $outdir;
			print "Can't create outdir: $outdir\n" if not -d $outdir;

			#### SET REFERENCE FILE
			my $referencefile = "$referencedir/$chromosome.fa";
			print "Can't find referencefile: $referencefile\n" if not -f $referencefile;

			#### SET BIN FILES DIRECTORY
			my $bindir = "$inputdir/$chromosome/bins";

			#### CREATE BINNER
			my $binner = UCSCBin->new({	samtools => $samtools, outputdir => $bindir	});

			#### GET CHROMOSOME SIZE
			my $chromosome_size = $binner->getChromosomeSize($bamfile);	

			#### GET ARRAY OF BINS
			my $bins = $binner->getBins($binlevel, $chromosome_size);

			#### GET RANGE FOR BINS
			my $rangefile = "$inputdir/$chromosome/hit.bam.range";

			#### SET BINS BY HIT RANGE	
			$bins = $binner->setBinsByRange($bins, $rangefile);

			#### DO EACH BAM BIN FILE
			foreach my $bin ( @$bins )
			{
				#### CHECK INPUT BAM FILE
				my $binfile = $binner->setBinfile($filename, $bindir, $binlevel, $bin);
				print "binfile: $binfile\n";
				print "SNP::realignBam    Can't find binfile: $binfile\n" and exit if not -f $binfile;

				#### SET OUTPUT FILES
				my ($filename) = $binfile =~ /([^\/]+)$/;
				my $intervalsfile = "$outputdir/$dirname/$chromosome/realign/$filename.intervals";
				my $outputfile = "$outdir/$filename";

				######## EXAMPLE:
				#### Perform local realignment of reads based on misalignments due to the presence of indels. 
				####
				#### time /usr/local/java/bin/java \
				#### -jar /nethome/syoung/base/apps/gatk/1.0.4705/GenomeAnalysisTK.jar \
				#### -I aln.bam \
				#### -R $REF \
				#### -T IndelRealigner \
				#### -targetIntervals realignTargets.intervals \
				#### --out realigned.bam \
				#### -U \
				#### -S SILENT

				my $command = qq{time $java/bin/java \\\n};
				$command .= qq{-jar $gatk/GenomeAnalysisTK.jar \\\n};
				$command .= qq{-T IndelRealigner \\\n};
				$command .= qq{-I $binfile \\\n};
				$command .= qq{-R $referencefile \\\n};
				$command .= qq{-targetIntervals $intervalsfile \\\n};
				$command .= qq{--out $outputfile \\\n};
				$command .= qq{-U \\\n};
				$command .= qq{-S SILENT \\\n};
				$command .= $params if defined $params;
				$command =~ s/\\+$//;
				print "command: $command\n";

				#### SET LABEL
				my $binnumber = $bin->{number};
				print "\nSNP::realignBam    bin->{number} not defined\n" if not defined $bin->{number};
				print Dumper $bin and exit if not defined $bin->{number};
				my $label = "realignBam-$chromosome-$dirname-$binnumber";

				##### RUN JOBS
				my $job = $self->setJob( [ $command ], $label, $outdir);
				push @$jobs, $job;
			}
		}
	}

	#### RUN JOBS
	print "\nSNP::realignBam    Running " , scalar(@$jobs), " jobs\n";	
	$self->runJobs( $jobs, "realignBam");
}







=head2

	SUBROUTINE		unbinSnp

	PURPOSE

		MERGE SNPS FROM BINNED CUMULATIVE SNP RUN TO CREATE A

		SINGLE SNP FILE

=cut

sub unbinSnp {
	my $self		=	shift;

	my $inputdirs	=	$self->get_inputdirs();
	my $outputdir	=	$self->get_outputdir();
	my $binlevel	=	$self->get_binlevel();
	my $filename	=	$self->get_filename();


	#### FOR EACH CHROMOSOME
	my $inputdir = $$inputdirs[0];
	print "SNP::cumulativeSnp    inputdir: $inputdir\n";
	my $chromosomes = $self->getChromosomes($inputdir);

	#### GET REQUIRED VARIABLES
	my $start 		=	$self->get_start();
	my $samtools	=	$self->get_samtools();

	#### PRINT ALL RANGE FILES AT ONCE IF FIRST ONE IS MISSING	
	$self->printRangefiles();

	#### REPEAT FOR ALL hit.snp FILES
	my $counter = 0;
	my $jobs = [];
	foreach my $inputdir ( @$inputdirs )
	{
		$counter++;

		my ($dirname) = $inputdir =~ /([^\/]+)$/;

		foreach my $chromosome ( @$chromosomes )
		{	
			#### CHECK INPUT BAM FILE
			my $bamfile = "$inputdir/$chromosome/$filename";
			print "SNP::unbinSnp    Can't find bamfile: $bamfile\n" and exit if not -f $bamfile;
			#### CREATE CHROMOSOME OUTDIR
			my $outdir = "$outputdir/$chromosome";
			File::Path::mkpath($outdir) if not -d $outdir;
			print "Can't create outdir: $outdir\n" if not -d $outdir;

			#### SET OUTPUT FILE
			my $outputfile = "$outdir/hit-$dirname.snp";

			#### SET BIN FILES DIRECTORY
			my $bindir = "$inputdir/$chromosome/bins";

			#### CREATE BINNER
			my $binner = UCSCBin->new({	samtools => $samtools, outputdir => $bindir	});

			#### GET CHROMOSOME SIZE
			my $chromosome_size = $binner->getChromosomeSize($bamfile);	

			#### GET ARRAY OF BINS
			my $bins = $binner->getBins($binlevel, $chromosome_size);

			#### SET BINS BY HIT RANGE	
			my $rangefile = "$inputdir/$chromosome/hit.bam.range";
			$bins = $binner->setBinsByRange($bins, $rangefile);

			#### COLLATE JOBS
			my $command = 'cat ';
			foreach my $bin ( @$bins )
			{
				my $binfile = $binner->setBinfile($filename, $bindir, $binlevel, $bin);
				print "binfile: $binfile\n";
				my ($filestub) = $binfile =~ /([^\/]+)$/;
				$filestub =~ s/\.bam$//;
				my $snpfile = "$outdir/$filestub-$dirname.snp";
				print "snpfile: $snpfile\n";
				$command .= " $snpfile ";
			}
			$command .= " > $outputfile";
			print "command: $command\n";

			#### SET LABEL
			my $label = "unbinSnp-$dirname-$chromosome";

			##### RUN JOBS
			my $job = $self->setJob( [ $command ], $label, $outdir);
			push @$jobs, $job;
		}
	}
	#### RUN JOBS
	print "\nSNP::unbinSnp    Running " , scalar(@$jobs), " jobs\n";
	$self->runJobs( $jobs, "unbinSnp");
}




=head2

	SUBROUTINE		mergeBinBam

	PURPOSE

		RUN MERGES IN COHORTS OF BINNED *bam FILE SUB-SECTIONS
		RUN MERGES IN COHORTS OF BINNED *bam FILE SUB-SECTIONS

		FOR Y = 1 TO TOTAL NUMBER OF BINS:

			1. CONCAT hit.binlevelX.numberY.*bam FILE IN N+1 DIRECTORY

				WITH EXISTING cumulative/chrA/hit.binlevelX.numberY-N.bam FILE

			2. PRINT TO cumulative/chrA/hit.binlevelX.numberY-N+1.bam FILE

			3. REPEAT UNTIL ALL DIRECTORIES ARE COMPLETED

=cut

sub mergeBinBam {
	my $self		=	shift;

	#### FILES AND DIRS	
	my $inputdirs 		=	$self->get_inputdirs();
	my $outputdir 		=	$self->get_outputdir();
	my $binlevel		=	$self->get_binlevel();
	my $bindir			=	$self->get_bindir();
	my $filename		=	$self->get_filename();

	#### GET CHROMOSOMES
	my $inputdir = $$inputdirs[0];
	my $chromosomes = $self->getChromosomes($inputdir);


	#### GET REQUIRED VARIABLES
	my $start 	=	$self->get_start();
	my $samtools = $self->get_samtools();
	my $samtools_index = $self->get_samtoolsindex();

	##### SET SLEEP IN BETWEEN JOBS
	my $sleep = $self->get_sleep();
	$sleep = 200 if not defined $sleep;
	#### REPEAT FOR ALL *bam BIN FILES
	#### IF start IS DEFINED, RUN BEGINNING FROM merged.bam-<start> FILE
	my $counter = 0;
	$counter = $start - 1 if defined $start;
	foreach my $inputdir ( @$inputdirs )
	{
		$counter++;

		my ($dirname) = $inputdir =~ /([^\/]+)$/;

		#### COLLATE JOBS
		my $jobs = [];

		#### DO ALL CHROMOSOMES AT THE SAME START POINT
		foreach my $chromosome ( @$chromosomes )
		{
			#### INPUTFILE AND OUTDIR
			my $inputfile = "$inputdir/$chromosome/$filename";
			my $outdir = "$outputdir/$chromosome";


			#### CHECK INPUT FILE		
			print "SNP::mergeBinBam    Can't find inputfile: $inputfile\n" and exit if not -f $inputfile;

			#### CREATE CHROMOSOME OUTDIR
			File::Path::mkpath($outdir) if not -d $outdir;
			print "Can't create outdir: $outdir\n" if not -d $outdir;

			#### SET BIN FILES DIRECTORY
			my $bindir = "$inputdir/$chromosome/bins";

			#### CREATE BINNER AND GET CHROMOSOME SIZE
			my $binner = UCSCBin->new({	samtools => $samtools, outputdir => $bindir	});
			my $chromosome_size = $binner->getChromosomeSize($inputfile);

			#### GET ARRAY OF BIN HASHES
			my $bins = $binner->getBins($binlevel, $chromosome_size);

			#### DO ALL BINS
			foreach my $bin ( @$bins )
			{
				my $commands = [];
				my $binfile = $binner->setBinfile($filename, $bindir, $binlevel, $bin);
				print "SNP::mergeBinBam    Skipping missing binfile (probably off-range): $binfile\n" and next if not -f $binfile;

				my ($filestub) = $binfile =~ /([^\/]+)$/;
				$filestub =~ s/\.bam$//;
				my $previous = $counter - 1;
				my $previousfile = "$outdir/$filestub-$previous.bam";
				my $tempfile = "$outdir/$filestub-$counter.bam";
				my $sortedfile_stub = "$outdir/$filestub-$counter";

				#### JUST COPY IF FIRST ITERATION
				if ( $counter == 1 )
				{
					#### COPY ORIGINAL BAM FILE TO MERGED FILE
					my $copy = "time cp $binfile $tempfile";
					push @$commands, $copy;

					#### SORT MERGED FILE
					my $sort = "time $samtools/samtools sort $tempfile $sortedfile_stub";
					push @$commands, $sort;
				}
				#### OTHERWISE, MERGE WITH PREVIOUS BAM FILE
				else
				{
					my $merge = "time $samtools/samtools merge $tempfile $binfile $previousfile";
					push @$commands, $merge;

					#### SORT MERGED FILE
					my $sort = "time $samtools/samtools sort $tempfile $sortedfile_stub";
					push @$commands, $sort;
				}

				my $main_command = join ";\n", @$commands;

				#### SET LABEL
				my $binnumber = $bin->{number};
				print "\nSNP::mergeBinBam    bin->{number} not defined\n" if not defined $bin->{number};
				print Dumper $bin and exit if not defined $bin->{number};
				my $label = "mergeBinBam-$chromosome-num$binnumber-$counter";

				##### RUN JOBS
				my $job = $self->setJob( [ $main_command ], $label, $outdir);
				push @$jobs, $job;
			}

		}	####	chromosomes	

		#### RUN JOBS
		print "\nSNP::mergeBinBam    Running jobs\n";
		$self->runJobs( $jobs, "cumulativeSnps-$dirname");

	}	#### 	inputdirs
}


=head2

	SUBROUTINE		cumulativeBam

	PURPOSE

		WE ARE PROVIDED A LIST OF DIRECTORIES CONTAINING chr* CHROMOSOME

		SUBDIRECTORIES, EACH CONTAINING A CHROMOSOME *.sam FILE. FOR EACH

		CHROMOSOME, STARTING WITH THE FIRST DIRECTORY AND GOING IN ORDER,

		FOR EACH DIRECTORY:

			1. CONCAT NEW *bam FILE WITH EXISTING merged-N.bam FILE

			2. REPEAT UNTIL ALL *bam FILES MERGED

=cut

sub cumulativeBam {
	my $self			=	shift;

	my $binlevel		=	$self->get_binlevel();
	print "SNP::cumulativeBam    binlevel: $binlevel\n";


	#### STARTING WITH THE FIRST DIRECTORY AND GOING IN ORDER
	####		
	if ( not defined $binlevel )
	{
		$self->mergeBam();
	}
	else
	{
		$self->mergeBinBam();
	}
}


=head2

	SUBROUTINE		cumulativeSnp

	PURPOSE

		WE ARE PROVIDED A LIST OF DIRECTORIES CONTAINING chr* CHROMOSOME

		SUBDIRECTORIES, EACH CONTAINING A CHROMOSOME *.sam FILE. FOR EACH

		CHROMOSOME, STARTING WITH THE FIRST DIRECTORY AND GOING IN ORDER,

		FOR EACH DIRECTORY:

			1. CONCAT NEW *bam FILE WITH EXISTING merged-N.bam FILE

			2. CALL SNPS USING THE merged-N+1.bam FILE

			3. REPEAT UNTIL ALL *bam FILES MERGED

=cut

sub cumulativeSnp {
	my $self			=	shift;

	my $binlevel		=	$self->get_binlevel();
	print "SNP::cumulativeSnp    binlevel: $binlevel\n";


	#### STARTING WITH THE FIRST DIRECTORY AND GOING IN ORDER
	####		
	if ( not defined $binlevel )
	{
		$self->mergeBam();
		$self->bamSnps();
	}
	else
	{
		$self->mergeBinBam();
		$self->binBamToSnp();
	}
}

=head2

	SUBROUTINE		mergeBam

	PURPOSE

		MERGE BAM FILES IN SERIES

=cut

sub mergeBam {
	my $self		=	shift;

	#### FILES AND DIRS	
	my $inputdirs 		=	$self->get_inputdirs();
	my $outputdir 		=	$self->get_outputdir();
	my $filename		=	$self->get_filename();
	my $zipped			=	$self->get_zipped();

	#### GET CHROMOSOMES
	my $inputdir = $$inputdirs[0];
	print "SNP::mergeBam    inputdir: $inputdir\n";
	my $chromosomes = $self->getChromosomes($inputdir);


	#### GET REQUIRED VARIABLES
	my $start 	=	$self->get_start();
	my $samtools = $self->get_samtools();
	my $samtools_index = $self->get_samtoolsindex();

	#### SET SLEEP IN BETWEEN JOBS
	my $sleep = $self->get_sleep();
	$sleep = 1500 if not defined $sleep;


	#### REPEAT FOR ALL INPUT DIRECTORIES	
	#### IF start IS DEFINED, RUN BEGINNING FROM merged.bam-<start> FILE
	my $counter = 0;
	$counter = $start - 1 if defined $start;
	foreach my $inputdir ( @$inputdirs )
	{
		#### COLLATE JOBS
		my $jobs = [];

		$counter++;

		my ($dirname) = $inputdir =~ /([^\/]+)$/;

		foreach my $chromosome ( @$chromosomes )
		{
			#### INPUTFILE AND OUTDIR
			my $inputfile = "$inputdir/$chromosome/$filename";
			my $outdir = "$outputdir/$chromosome";

			#### CHECK INPUT FILE		
			print "SNP::mergeBam    Can't find inputfile: $inputfile\n" and exit if not -f $inputfile;

			#### CREATE CHROMOSOME OUTDIR
			File::Path::mkpath($outdir) if not -d $outdir;
			print "Can't create outdir: $outdir\n" if not -d $outdir;

			#### STORE COMMANDS HERE
			my $commands = [];

			#### UNZIP SAMFILE IF ZIPPED
			if ( $inputfile !~ /\.bam$/ )
			{
				if ( defined $zipped )
				{
					my $zippedfile = $inputfile;
					$inputfile =~ s/\.gz$// if $zipped eq "gzip";
					$inputfile =~ s/\.zip$// if $zipped eq "zip";
					$inputfile =~ s/\.bz2$// if $zipped eq "bz2";

					print "SNP::mergeBam    Can't find zippedfile: $zippedfile\n" and exit if not -f $zippedfile;

					my $unzip = "time zcat $zippedfile > $inputfile" if $zipped eq "gzip" or $zipped eq "zip";
					$unzip = "time bzip2 -dc $zippedfile > $inputfile" if $zipped eq "bz2";
					push @$commands, $unzip;
				}

				#### 1. CONVERT SAM TO BAM
				my $bamfile = $inputfile;
				$bamfile =~ s/\.sam$/.bam/;
				my $convert = "time $samtools/samtools view -bt $samtools_index/$chromosome.fai -o $bamfile $inputfile";
				push @$commands, $convert;
			}

			#### SET FILENAMES
			my $bamfile = $inputfile;
			$bamfile =~ s/\.sam$/.bam/;
			my $previous = $counter - 1;
			my $previousfile = "$outdir/merged-$previous.bam";
			my $tempfile = "$outdir/merged-$dirname.bam.temp";
			my $sortedfile_stub = "$outdir/merged-$dirname";
			my $mergedfile = "$outdir/merged-$dirname.bam";

			#### JUST COPY IF FIRST ITERATION
			if ( $counter == 1 )
			{
				#### COPY ORIGINAL BAM FILE TO MERGED FILE
				my $copy = "time cp $bamfile $outdir/merged-$counter.bam";
				push @$commands, $copy;
			}
			#### OTHERWISE, MERGE WITH PREVIOUS BAM FILE
			else
			{
				my $merge = "time $samtools/samtools merge $tempfile $bamfile $previousfile";
				push @$commands, $merge;
			}

			#### SORT MERGED FILE
			my $sort = "time $samtools/samtools sort $tempfile $sortedfile_stub";
			push @$commands, $sort;

		}	#### chromosomes

		#### RUN JOBS
		$self->runJobs( $jobs, "mergeBam");

	} #### inputdirs

}



=head2

	SUBROUTINE		getChromosomes

	PURPOSE

		RETRIEVE THE LIST OF CHROMOSOME SUBDIRECTORIES

	NOTES

=cut

sub getChromosomes {
	my $self		=	shift;
	my $inputdir	=	shift;

	print "SNP::getChromosomes    inputdir not defined\n" and exit if not defined $inputdir;
	chdir($inputdir);
	my @chromosomes = <chr*>;

	return \@chromosomes;
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

sub new {
 	my $class 		=	shift;
	my $arguments 	=	shift;

	my $self = {};
    bless $self, $class;

	#### INITIALISE THE OBJECT'S ELEMENTS
	$self->initialise($arguments);

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


sub initialise {
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
}


1;


