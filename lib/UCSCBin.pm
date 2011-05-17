package UCSCBin;



=head2

		PACKAGE		UCSCBin

		VERSION		0.01

		PURPOSE

	        BIN ITEMS ACCORDING TO THE UCSC BINNING SCHEMA

		HISTORY

			0.01 BASIC VERSION FOR INITIAL IMPLEMENTATION OF UCSC BINNING
				(Used when chromEnd is less than or equal to 536,870,912 BASES)

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

INPUTDIRS
OUTPUTDIR
REFERENCE
SAMTOOLS
BINLEVEL
EXECUTABLE
SUFFIX
FILENAME

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

#### SIZE THRESHOLD BETWEEN INITIAL AND EXTENDED IMPLEMENTATION
our $MAXSIZE = 536870912;


=head2

	SUBROUTINE		setBinsByRange

	PURPOSE

		MERGE SNPS FROM BINNED CUMULATIVE SNP RUN TO CREATE A

		SINGLE SNP FILE

    INPUTS

        1. SORTED *bam FILE

        2. OUTPUT FILE LOCATION

    OUTPUTS

        1. TAB-SEPARATED VALUES OUTPUT FILE WITH THE FOLLOWING FORMAT:

			#READS	#START	#STOP		<= HEADER LINE
			343040	1500031	51389018	<= DATA LINE

=cut

sub setBinsByRange {
	my $self		=	shift;
	my $bins 		=	shift;
	my $rangefile 	=	shift;


	print "UCSCBin::setBinsByRange    Can't find rangefile: $rangefile\n" and exit if not -f $rangefile;
	open(FILE, $rangefile) or die "Can't open rangefile: $rangefile\n";
	<FILE>;
	my ($readnumber, $hitstart, $hitstop) = <FILE> =~ /^(\S+)\s+(\S+)\s+(\S+)/;	
	print "UCSCBin::setBinsByRange    readnumber not defined\n" and exit if not defined $readnumber;
	print "UCSCBin::setBinsByRange    hitstart not defined\n" and exit if not defined $hitstart;
	print "UCSCBin::setBinsByRange    hitstop not defined\n" and exit if not defined $hitstop;

	#### SET BINS 
	my $counter = 0;
	for ( my $i = 0; $i < @$bins; $i++ )
	{
		my $bin = $$bins[$i];
		$counter++;
		my $binstart = $bin->{start};
		my $binstop = $bin->{stop};

		#### SKIP IF BIN BEFORE START OF HITS
		if ( $binstop < $hitstart or $binstart > $hitstop )
		{
			splice (@$bins, $i, 1);
			$i--;
			next;
		}

		#### ADJUST START/STOP OF BIN TO START/STOP OF HITS
		$bin->{start} = $hitstart if $binstart < $hitstart;
		$bin->{stop} = $hitstop if $binstop > $hitstop;
	}

	return $bins;
}



=head2

	SUBROUTINE		printRangefiles

	PURPOSE

		MERGE SNPS FROM BINNED CUMULATIVE SNP RUN TO CREATE A

		SINGLE SNP FILE

    INPUTS

        1. SORTED *bam FILE

        2. OUTPUT FILE LOCATION

    OUTPUTS

        1. TAB-SEPARATED VALUES OUTPUT FILE WITH THE FOLLOWING FORMAT:

			#READS	#START	#STOP		<= HEADER LINE
			343040	1500031	51389018	<= DATA LINE

=cut

sub printRangefiles {
	my $self		=	shift;
	my $inputdirs	=	shift;
	my $chromosomes =	shift;

	print "UCSCBin::printRangefiles    UCSCBin::printRangefiles(inputdirs, chromosomes)\n";

	#### EXECUTABLE
	my $executable		=	"$Bin/hitRange.pl";

	#### GET HIT RANGES
	my $jobs = [];
	foreach my $inputdir ( @$inputdirs )
	{
		foreach my $chromosome ( @$chromosomes )
		{
			my $bamfile = "$inputdir/$chromosome/hit.bam";

			my $rangefile = "$inputdir/$chromosome/hit.bam.range";
			print "UCSCBin::printRangefiles    rangefile: $rangefile\n";

			my $command = "$executable --inputfile $bamfile --outputfile $rangefile";

			my $label = "binbam-range-$chromosome";
			my $outdir = "$inputdir/$chromosome";

			my $job = $self->setJob( [ $command ], $label, $outdir);
			push @$jobs, $job;
		}
	}
	#### RUN INDEX JOBS
	$self->runJobs($jobs, "binBam-range");
}


=head2

	SUBROUTINE		hitRange

	PURPOSE

		MERGE SNPS FROM BINNED CUMULATIVE SNP RUN TO CREATE A

		SINGLE SNP FILE

    INPUTS

        1. SORTED *bam FILE

        2. OUTPUT FILE LOCATION

    OUTPUTS

        1. TAB-SEPARATED VALUES OUTPUT FILE WITH THE FOLLOWING FORMAT:

			#READS	#START	#STOP		<= HEADER LINE
			343040	1500031	51389018	<= DATA LINE

=cut

sub hitRange {
	my $self		=	shift;
	my $inputfile	=	shift;

	my $samtools	=	$self->get_samtools();
	my $dot			=	$self->get_dot();
	$dot = 1000000 if not defined $dot;

	my $counter = 0;
	my $start = 0;
	my $stop = 0;
	open(FILE, "$samtools/samtools view $inputfile|" ) or die "Can't open inputfile: $inputfile\n";
	my $line = <FILE>;
	($start) = $line =~ /^\S+\s+\S+\s+\S+\s+(\S+)\s+/;
	while ( <FILE> )
	{
		print "UCSCBin::hitRange    $counter\n" if $counter % $dot == 0;
		$counter++;
		last if not defined $_;
		$line = $_;
	}
	($stop) = $line =~ /^\S+\s+\S+\s+\S+\s+(\S+)\s+/;

	return $counter, $start, $stop;
}




=head2

	SUBROUTINE		binBam

    PURPOSE

		SPLIT UP *bam FILES INTO EQUALLY-SIZED CHROMOSOMAL RANGE BINS:

		-	USE UCSC BINNING SYSTEM 
			http://genomewiki.ucsc.edu/index.php/Bin_indexing_system
			Kent & al. The Human Genome Browser at UCSC", doi: 10.1101/gr.229102 . Genome Res. 2002. 12:996-1006

		-   DIVIDE UP *.bam FILE ENTRIES BY CHROMOSOMAL POSITION

		-	WIDTH OF THE BINS IS DEFINED BY THE SPECIFIED BIN LEVEL

		-   DEFAULT: BIN LEVEL 1 (8 x 64MB BINS)

		-   ADD binlevel, binnumber AND totalbins TO BINNED *bam FILE NAME	

	NOTES

		UCSC binning
		http://genomewiki.ucsc.edu/index.php/Bin_indexing_system
		Kent & al. The Human Genome Browser at UCSC", doi: 10.1101/gr.229102 . Genome Res. 2002. 12:996-1006

		Initial implementation		
		Used when chromEnd is less than or equal to 536,870,912 = 229
			bin numbers	bin
		level	#bins	start	end		size
		0		1		0		0		512 Mb
		1		8		1		8		64 Mb
		2		64		9		72		8 Mb
		3		512		73		584		1 Mb
		4		4096	585		4680	128 kb

		Extended implementation
		Used when chromEnd is greater than 536,870,912 = 229 and less than 2,147,483,647 = 231 - 1
			bin numbers	bin
		level	#bins	start	end		size
		0		1		4691	4691	2 Gb
		1		8		4683	4685	512 Mb
		2		64		4698	4721	64 Mb
		3		512		4818	5009	8 Mb
		4		4,096	5778	7313	1 Mb
		5		32,768	13458	25745	128 kb

=cut

sub binBam {
	my $self		=	shift;	

print "UCSCBin::binBam    UCSCBin::binBam()\n";

	#### FILES AND DIRS	
	my $reference 		=	$self->get_reference();
	my $inputdirs 		=	$self->get_inputdirs();
	my $binlevel 		=	$self->get_binlevel();
	my $size 			=	$self->get_size();
	my $outputdir		=	$self->get_outputdir();

	#### GET CHROMOSOMES
	my $indir = $$inputdirs[0];
	my $chromosomes = $self->getChromosomes($indir);

	#### GET SAMTOOLS
	my $samtools = $self->get_samtools();

	#### SORT INPUT BAM FILES
	my $sortjobs = [];
	foreach my $inputdir ( @$inputdirs )
	{
		foreach my $chromosome ( @$chromosomes )
		{
			my $bamfile = "$inputdir/$chromosome/hit.bam";
			print "UCSCBin::binBam    bamfile: $bamfile\n";

			my $filestub = $bamfile;
			$filestub =~ s/\.bam$//;
			my $command = "$samtools/samtools sort $bamfile $filestub";

			my $label = "binbam-sort-$chromosome";
			my $outdir = "$inputdir/$chromosome";

			my $job = $self->setJob( [ $command ], $label, $outdir);
			push @$sortjobs, $job;
		}
	}
	#### RUN SORT JOBS
	$self->runJobs($sortjobs, "binBam-sort");

	#### INDEX INPUT BAM FILES
	my $indexjobs = [];
	foreach my $inputdir ( @$inputdirs )
	{
		foreach my $chromosome ( @$chromosomes )
		{
			my $bamfile = "$inputdir/$chromosome/hit.bam";
			print "UCSCBin::binBam    bamfile: $bamfile\n";

			my $filestub = $bamfile;
			$filestub =~ s/\.bam$//;
			my $command = "$samtools/samtools index $bamfile";

			my $label = "binbam-index-$chromosome";
			my $outdir = "$inputdir/$chromosome";

			my $job = $self->setJob( [ $command ], $label, $outdir);
			push @$indexjobs, $job;
		}
	}

	#### RUN INDEX JOBS
	$self->runJobs($indexjobs, "binBam-index");

	#### GET HIT RANGE FILES
	$self->printRangefiles($inputdirs, $chromosomes);

	#### DO BINS FOR ALL DIRECTORIES AND ALL CHROMOSOMES
	my $binjobs = [];
	foreach my $inputdir ( @$inputdirs )
	{
		foreach my $chromosome ( @$chromosomes )
		{
			my $bamfile = "$inputdir/$chromosome/hit.bam";
			print "UCSCBin::binBam    bamfile: $bamfile\n";

			#### MAKE OUTPUT DIR IF NOT EXISTS
			my $outputdir = "$inputdir/$chromosome/bins";
			print "UCSCBin::binBam    outputdir is a file: $outputdir\n" and exit if -f $outputdir;
			File::Path::mkpath($outputdir) if not -d $outputdir;
			print "UCSCBin::binBam    Can't create output directory: $outputdir\n" and exit if not -d $outputdir;

			my $chromosome_size	=	$self->getChromosomeSize($bamfile);
			print "UCSCBin::binBam    chromosomeSize: $chromosome_size\n";

			#### SET BINS BY HIT RANGE	
			my $bins = $self->getBins($binlevel, $chromosome_size);
			my $rangefile = "$bamfile.range";
			$bins = $self->setBinsByRange($bins, $rangefile);
			print "UCSCBin::binBam    number bins: " , scalar(@$bins), "\n";

			#### SET BIN FILE CREATION JOBS
			my $counter = 0;
			foreach my $bin ( @$bins )
			{
				$counter++;
				#### SET COMMAND		
				my $command = $self->binCommand($bamfile, $outputdir, $reference, $binlevel, $bin);

				#### SET LABEL
				my $label = "binbam-$counter";

				#### SET JOB 
				my $job = $self->setJob([$command], $label, $outputdir);

				push @$binjobs, $job;

			}	#### bins

		}	#### chromosomes

	}	#### inputdirs

	#### RUN BIN JOBS
	$self->runJobs($binjobs, "binBam");

} #	binBam

=head2

	SUBROUTINE		getChromosomes

	PURPOSE

		RETRIEVE THE LIST OF CHROMOSOME SUBDIRECTORIES

	NOTES

=cut

sub getChromosomes {
	my $self		=	shift;
	my $inputdir	=	shift;

	chdir($inputdir);
	my @chromosomes = <chr*>;

	return \@chromosomes;
}



=head2

	SUBROUTINE		binCommand

	PURPOSE

		GENERATE *bam.bai INDEX FILE FOR GIVEN *bam FILE

	NOTES

		http://samtools.sourceforge.net/samtools.shtml

		view	 samtools view [-bchuHS] [-t in.refList] [-o output] [-f reqFlag] [-F skipFlag] [-q minMapQ] [-l library] [-r readGroup] [-R rgFile] <in.bam>|<in.sam> [region1 [...]]

		Extract/print all or sub alignments in SAM or BAM format. If no region is specified, all the alignments will be printed; otherwise only alignments overlapping the specified regions will be output. An alignment may be given multiple times if it is overlapping several regions. A region can be presented, for example, in the following format: ‘chr2’ (the whole chr2), ‘chr2:1000000’ (region starting from 1,000,000bp) or ‘chr2:1,000,000-2,000,000’ (region between 1,000,000 and 2,000,000bp including the end points). The coordinate is 1-based.

		OPTIONS:
		-b	 Output in the BAM format.
		-f INT	 Only output alignments with all bits in INT present in the FLAG field. INT can be in hex in the format of /^0x[0-9A-F]+/ [0]
		-F INT	 Skip alignments with bits present in INT [0]
		-h	 Include the header in the output.
		-H	 Output the header only.
		-l STR	 Only output reads in library STR [null]
		-o FILE	 Output file [stdout]
		-q INT	 Skip alignments with MAPQ smaller than INT [0]
		-r STR	 Only output reads in read group STR [null]
		-R FILE	 Output reads in read groups listed in FILE [null]
		-S	 Input is in SAM. If @SQ header lines are absent, the ‘-t’ option is required.
		-c	 Instead of printing the alignments, only count them and print the total number. All filter options, such as ‘-f’, ‘-F’ and ‘-q’ , are taken into account.
		-t FILE	 This file is TAB-delimited. Each line must contain the reference name and the length of the reference, one line for each distinct reference; additional fields are ignored. This file also defines the order of the reference sequences in sorting. If you run ‘samtools faidx <ref.fa>’, the resultant index file <ref.fa>.fai can be used as this <in.ref_list> file.
		-u	 Output uncompressed BAM. This option saves time spent on compression/decomprssion and is thus preferred when the output is piped to another samtools command.

=cut

sub binCommand {
	my $self		=	shift;
	my $bamfile		=	shift;
	my $outputdir	=	shift;
	my $reference	=	shift;
	my $binlevel	=	shift;
	my $bin			=	shift;

	#### CHECK INPUTS
	print "UCSCBin::binCommand    bamfile not defined\n" and exit if not defined $bamfile;


	#### GET SAMTOOLS	
	my $samtools = $self->get_samtools();
	print "UCSCBin::binCommand    samtools not defined\n" and exit if not defined $samtools;

	#### GET BINFILE
	my $binfile = $self->setBinfile($bamfile, $outputdir, $binlevel, $bin);

	my $region = "$reference:$bin->{start}-$bin->{stop}";
	my $command = "$samtools/samtools view -h -b $bamfile $region > $binfile";
	print "$command\n";

	return $command;
}

=head2

	SUBROUTINE		setBinfile

	PURPOSE

		SET BIN FILE NAME GIVEN BINLEVEL AND BIN HASH

	NOTES

=cut

sub setBinfile {
	my $self		=	shift;
	my $bamfile		=	shift;
	my $outputdir	=	shift;
	my $binlevel	=	shift;				
	my $bin			=	shift;
	my $suffix		=	shift;

	$suffix = "bam" if not defined $suffix;


	my ($filestub) = $bamfile =~ /([^\/]+)$/;
#exit;

	$filestub =~ s/\.bam$//;

	my $binfile = "$outputdir/$filestub.binlevel$binlevel.num$bin->{number}.$suffix";

	return $binfile;
}

=head2

	SUBROUTINE		createIndex

	PURPOSE

		GENERATE *bam.bai INDEX FILE FOR GIVEN *bam FILE

=cut

sub createIndex {
	my $self	=	shift;
	my $bamfile	=	shift;


	#### CHECK INPUTS
	print "UCSCBin::createIndex    bamfile not defined\n" and exit if not defined $bamfile;
	my $samtools = $self->get_samtools();
	print "UCSCBin::createIndex    samtools not defined\n" and exit if not defined $samtools;

	my $indexfile = $bamfile . ".bai";
	return if -f $indexfile;

	my $command = "$samtools/samtools index $bamfile";
	`$command`;

	print "UCSCBin::createIndex    Can't create index file for bamfile: $bamfile\n" and exit if not -f $indexfile;
}


=head2

	SUBROUTINE		getBins

	PURPOSE

		CALL SAMTOOLS idxstats TO GET THE CHROMOSOME SIZE FROM THE BAMFILE

=cut

sub getChromosomeSize {
	my $self	=	shift;
	my $bamfile	=	shift;

	my $samtools = $self->get_samtools();
	print "UCSCBin::getChromosomeSize    samtools not defined\n" and exit if not defined $samtools;
	print "UCSCBin::getChromosomeSize    bamfile not defined\n" and exit if not defined $bamfile;

	#### CREATE INDEX
	$self->createIndex($bamfile);

	my $command = "$samtools/samtools idxstats $bamfile";
	my $idxstats = `$command`;

	#### FORMAT:
	#### chromo	 length	         #reads
	#### chr22   51304566        914819  0
	#### *       0       0       0
	my ($size) = $idxstats =~ /^\S+\s+(\S+)/;

	return $size;
}


=head2

	SUBROUTINE		getBins

	PURPOSE

		RETURN A LIST OF REFERENCE NAMES (N.B.: NOT THE PATHS TO THE FILES)

=cut

sub getBins {
	my $self		=	shift;
	my $binlevel	=	shift;
	my $chromosome_size	=	shift;


	my $bins;	
	if ( $binlevel > 5 )
	{
		$bins = $self->getBinsBySize($binlevel, $chromosome_size);
	}
	else
	{
		$bins = $self->getInitialBins($binlevel, $chromosome_size) if $chromosome_size < $MAXSIZE;
		$bins = $self->getExtendedBins($binlevel, $chromosome_size) if $chromosome_size >= $MAXSIZE;
	}

	return $bins;
}

=head2

	SUBROUTINE		getBinsBySize

	PURPOSE

		RETURN BINS WITH SPECIFIED SIZE

=cut

sub getBinsBySize {
	my $self		=	shift;
	my $binsize		=	shift;
	my $chromosome_size	=	shift;


	#### MINIMUM BIN SIZE
	my $numberbins = int($chromosome_size/$binsize);
	$numberbins++ if $chromosome_size % $binsize > 0;

	#### SET BINS
	my $bins = [];
	my $binstart = 0;
	my $binstop = 0;
	my $number = 0;
	for ( my $i = 0; $i < $numberbins; $i++ )
	{
		$number++;
		$binstart = $binstop + 1;
		$binstop += $binsize;
		$binstop = $chromosome_size if $binstop > $chromosome_size;

		push @$bins, { number => $number, start => $binstart, stop => $binstop };
	}	

	return $bins;
}


=head2

	SUBROUTINE		unbin

	PURPOSE

		CONCATENATE BINNED FILES INTO A SINGLE FILE PER INPUT DIR AND CHROMOSOME

=cut

sub unbin {
	my $self		=	shift;

	my $inputdirs	=	$self->get_inputdirs();
	my $outputdir	=	$self->get_outputdir();
	my $binlevel	=	$self->get_binlevel();
	my $filename	=	$self->get_filename();
	my $suffix		=	$self->get_suffix();


	#### CHECK INPUTS
	print "UCSCBin::unbin    suffix not defined\n" and exit if not defined $suffix;
	print "UCSCBin::unbin    outputdir not defined\n" and exit if not defined $outputdir;
	print "UCSCBin::unbin    inputdirs not defined\n" and exit if not defined $inputdirs;
	print "UCSCBin::unbin    binlevel not defined\n" and exit if not defined $binlevel;
	print "UCSCBin::unbin    filename not defined\n" and exit if not defined $filename;

	#### GET CHROMOSOMES
	my $inputdir = $$inputdirs[0];
	my $chromosomes = $self->getChromosomes($inputdir);

	#### GET REQUIRED VARIABLES
	my $samtools	=	$self->get_samtools();

	#### REPEAT FOR ALL BINNED FILES PER INPUT DIR AND CHROMOSOME
	my $counter = 0;
	my $jobs = [];
	foreach my $inputdir ( @$inputdirs )
	{
		$counter++;

		foreach my $chromosome ( @$chromosomes )
		{	
			#### CHECK INPUT BAM FILE
			my $bamfile = "$inputdir/$chromosome/$filename";
			print "UCSCBin::unbin    Can't find bamfile: $bamfile\n" and exit if not -f $bamfile;
			#### CREATE CHROMOSOME OUTDIR
			my $outdir = "$outputdir/$chromosome";
			File::Path::mkpath($outdir) if not -d $outdir;
			print "Can't create outdir: $outdir\n" if not -d $outdir;

			#### SET OUTPUT FILE
			my $outputfile = "$outdir/hit-$counter.$suffix";

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
			$self->printRangefiles($inputdirs, $chromosomes) if not -f $rangefile;
			$bins = $binner->setBinsByRange($bins, $rangefile);

			#### COLLATE JOBS
			my $command = 'cat ';
			foreach my $bin ( @$bins )
			{
				my $binfile = $binner->setBinfile($filename, $bindir, $binlevel, $bin);
				my ($filestub) = $binfile =~ /([^\/]+)$/;
				$filestub =~ s/\.bam$//;
				my $smallfile = "$outdir/$filestub-$counter.$suffix";
				$command .= " $smallfile ";
			}
			$command .= " > $outputfile";

			#### SET LABEL
			my $label = "unbin-$chromosome-$counter";

			##### RUN JOBS
			my $job = $self->setJob( [ $command ], $label, $outdir);
			push @$jobs, $job;
		}
	}
	#### RUN JOBS
	print "\nUCSCBin::unbin    Running " , scalar(@$jobs), " jobs\n";
	$self->runJobs( $jobs, "unbin");
}



=head2

	SUBROUTINE		getInitialBins

	PURPOSE

		RETURN BINS FOR INITIAL IMPLEMENTATION OF UCSC BINNING SCHEMA

	NOTES

		UCSC binning
		http://genomewiki.ucsc.edu/index.php/Bin_indexing_system
		Kent & al. The Human Genome Browser at UCSC", doi: 10.1101/gr.229102 . Genome Res. 2002. 12:996-1006

		Initial implementation
		Used when chromEnd is less than or equal to 536,870,912 = 2**29
			bin numbers	bin
		level	#bins	start	end		size
		0		1		0		0		512 Mb
		1		8		1		8		64 Mb
		2		64		9		72		8 Mb
		3		512		73		584		1 Mb
		4		4096	585		4680	128 kb

=cut

sub getInitialBins {
	my $self		=	shift;
	my $binlevel	=	shift;
	my $chromosome_size	=	shift;

	#### MINIMUM BIN SIZE
	my $binsize = 2**17;
	my $multiple = (4 - $binlevel);
	$binsize = $binsize << ($multiple * 3);
	my $numberbins = 2**($binlevel * 3);

	#### SET BIN START NUMBER 	
	my $number = 0;
	my $counter = 0;
	while ( $counter < $binlevel )
	{
		$number += 2**($counter * 3);
		$counter++;
	}

	#### SET BINS
	my $bins = [];
	my $binstart = 0;
	my $binstop = 0;
	for ( my $i = 0; $i < $numberbins; $i++ )
	{
		$binstart = $binstop + 1;
		$binstop += $binsize;

		push @$bins, { number => $number, start => $binstart, stop => $binstop };
		$number++;

		last if $binstop >= $chromosome_size;
	}	

	return $bins;
}


=head2

	SUBROUTINE		getExtendedBins

	PURPOSE

		RETURN BINS FOR EXTENDED IMPLEMENTATION OF UCSC BINNING SCHEMA

	NOTES

		UCSC binning
		http://genomewiki.ucsc.edu/index.php/Bin_indexing_system
		Kent & al. The Human Genome Browser at UCSC", doi: 10.1101/gr.229102 . Genome Res. 2002. 12:996-1006

		Extended implementation
		Used when chromEnd is less than or equal to 536,870,912 = 2**29
			bin numbers	bin
		level	#bins	start	end		size
		0		1		4691	4691	2 Gb
		1		8		4683	4685	512 Mb
		2		64		4698	4721	64 Mb
		3		512		4818	5009	8 Mb
		4		4,096	5778	7313	1 Mb
		5		32,768	13458	25745	128 kb

=cut

sub getExtendedBins {
	my $self		=	shift;
	my $binlevel	=	shift;
	#
	#
	##### MINIMUM BIN SIZE
	#my $binsize = 2**17;
	#my $multiple = (4 - $binlevel);
	#$binsize = $binsize << ($multiple * 3);
	#
	#my $numberbins = 2**($binlevel * 3);
	#
	#my $start = 0;
	#my $counter = 0;
	#while ( $counter < $binlevel )
	#{
	#	$start += 2**($counter * 3);
	#	$counter++;
	#}
	#
	#my $end = $start + $numberbins - 1;
	#
	#my $bins = [];
	#my $binstart = 0;
	#my $binstop = 0;
	#for ( my $i = 0; $i < $numberbins; $i++ )
	#{
	#	$binstart = $binstop + 1;
	#	$binstop += $binsize;
	#	
	#	push @$bins, { start => $binstart, stop => $binstop };
	#}	
	#
	#return $bins;

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


#has 'x' => ( isa => 'Int', is => 'ro' );
#has 'y' => ( isa => 'Int', is => 'rw' );

#__PACKAGE__->meta->make_immutable;
#
#sub BUILD
#{
#	my $self	=	shift;
#	my $args	=	shift;
#	
#	#$self->initialise($args);
#	
#	exit;
#	
#}
#

1;

