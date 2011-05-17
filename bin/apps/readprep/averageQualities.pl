#!/usr/bin/perl -w

#### DEBUG
$DEBUG = 1;

#### TIME
my $time = time();

=head2

	APPLICATION     averageQualities

    PURPOSE

        1. RUN averageQuality.pl FOR ALL .fastq/.fastq.gz FILES IN A DIRECTORY

    INPUT

        1. DIRECTORY CONTAINING .fastq OR fastq.gz FILES

    OUTPUT

        1. FOR EACH INPUT FILE:

			1. TAB-SEPARATED *.avqual.tsv FILE

			2. TAB-SEPARATED FILE *.qualstats.tsv FILE

	USAGE

    ./averageQualities.pl  <--inputdir String> <--outputdir String> [--cluster Integer] [--maxjobs Integer] [--limit Integer] [--compress] [--fixed] [--id_length Integer] [--sequence_length Integer] [--compress String] [--dot Integer] [-h]

        --inputdir		       :   /Full/path/to/inputdir 
        --outputdir		       :   /Full/path/to/outputdir 
        --compress             :   Compress with gzip or zip
        --fixed                :   Print FASTA record of fixed length
        --id_length            :   Fixed id length
        --sequence_length      :   Fixed sequence length
        --cluster              :   Run on cluster with 'clusterSubmit.pl'
        --maxjobs                 :   Maximum number of concurrent maxjobs
        --dot                  :   Print counter every 'dot' number of records
        --help                 :   print help info

    EXAMPLES


/nethome/bioinfo/apps/agua/0.5/bin/apps/readprep/averageQualities.pl \
--inputdir /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/SRX000600,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/SRX000601,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/SRX000602,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/SRX000603,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/SRX001539,/scratch/syoung/base/pipeline/solid/NA18507/SRP000239/SRX001540 \
--outputfile /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/samples/100M/averageQualities.txt \
--compress gzip \
--bins -50,-40,-30,-20,-10,0,1,2,3,4,5,6,7,8,9,10,15,20,30,40,50,60,70,80,90,100,200,300 \
--type sanger \
--cluster LSF \
--maxjobs 600 \
--queue small \
--clean \
&> /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/samples/100M/averageQualities.out


=cut

use strict;

#### USE LIBRARY
use FindBin qw($Bin);
use lib "$Bin/../../../lib";
use lib "$Bin/../../../lib/external";

#### FLUSH BUFFER
$| = 1;

#### INTERNAL MODULES
use Monitor;
use ReadPrep;
use SolexaUtil;
use Timer;
use Util;
use Conf::Agua;

#### INITIALISE SolexaUtil OBJECT
my $solexa = SolexaUtil->new();

#### EXTERNAL MODULES
use Term::ANSIColor qw(:constants);
use Getopt::Long;
use IO::Pipe;
use File::Path;
use Data::Dumper;

#### SAVE ARGUMENTS
my @arguments = @ARGV;

#### DEFAULTS
my $binarray = [ 0,1,2,3,4,5,6,7,8,9,10,15,20,30,40,50,60,70,80,90,100,200,300 ];

#### GET OPTIONS
my $inputdir;
my $outputfile;
my $compress;
my $bins;
my $type;
my $clean;
my $min;
my $max;
my $skip;
my $length;
my $dot = 1000000;
my $cluster;
my $maxjobs = 30;
my $sleep = 5;
my $queue;
my $cleanup;
my $help;
if ( not GetOptions (
    'inputdir=s' 	=> \$inputdir,
    'outputfile=s' 	=> \$outputfile,
    'compress=s' 	=> \$compress,
    'bins=s' 		=> \$bins,
    'type=s' 		=> \$type,
    'clean' 		=> \$clean,
    'min=i' 		=> \$min,
    'max=i' 		=> \$max,
    'skip=i' 		=> \$skip,
    'length=i' 		=> \$length,
    'dot=i' 		=> \$dot,
    'cluster=s' 	=> \$cluster,
    'maxjobs=i' 	=> \$maxjobs,
    'queue=s' 		=> \$queue,
    'cleanup' 		=> \$cleanup,
    'sleep=i' 		=> \$sleep,
    'help' 			=> \$help
) )
{ print "Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "inputdir not defined (Use --help for usage)\n" if not defined $inputdir;
die "outputfile not defined (Use --help for usage)\n" if not defined $outputfile;
die "maxjobs not defined\n" if defined $cluster and not defined $maxjobs;
die "type not defined\n" if not defined $type;
die "type must be 'fastq' or 'fasta'\n" if defined $type and $type !~ /^(fastq|fasta)$/;
die "compress must be 'gzip' or 'zip'\n" if defined $compress and $compress !~ /^(gzip|zip)$/;

#### SET EXECUTABLE
my $executable = "$Bin/averageQuality.pl";
print "executable: $executable\n";

#### SET BINARRAY IF NOT DEFAULT
@$binarray = split ",", $bins if defined $bins and $bins;

#### GET DIRECTORIES
my $files = files($inputdir, $compress);
print "No. files: ", scalar(@$files), "\n";
#exit;

#### RUN fastq2fasta.pl ON EVERY FILE
my $commands = [];  #### FOR RUNNING LOCALLY
my $jobs = [];		#### FOR RUNNING ON CLUSTER
my $clusterObject = ReadPrep->new(
	{
		cluster 	=> $cluster,
		queue 		=> $queue,
		cleanup 	=> $cleanup,
		maxjobs     => $maxjobs,
		sleep       => $sleep,
		dot         => $dot
	}
);

#### COLLECT JOBS OR COMMANDS
my $counter = 0;
foreach my $inputfile ( @$files )
{
	$counter++;


	#### SKIP IF OUTPUT FILE EXISTS ALREADY AND clean NOT SPECIFIED
	my $outputfile = $solexa->avqualfile($inputfile, $min, $max, $skip, $length);
	$outputfile .= ".gz" if defined $compress and $compress eq "gzip";
	$outputfile .= ".zip" if defined $compress and $compress eq "zip";
	print "averageQualities.pl    skipping because found outputfile: $outputfile\n"
		and next if -f $outputfile and not defined $clean;
	print "averageQualities.pl    outputfile: $outputfile\n";

	my ($outdir) = $inputfile =~ /^(.+?)\/[^\/]+$/;

	#### SET COMMAND
    my $command = qq{$executable --inputfile $inputfile --type $type --dot $dot};
	$command .= qq{ --compress $compress} if defined $compress;
	$command .= qq{ --bins $bins} if defined $bins;
	$command .= qq{ --min $min} if defined $min;
	$command .= qq{ --max $max} if defined $max;
	$command .= qq{ --skip $skip} if defined $skip;
	$command .= qq{ --length $length} if defined $length;
	$command .= qq{ --clean} if defined $clean;
#exit;
#	
	#### SET LABEL
	my ($filename) = $inputfile =~ /([^\/]+)$/;
	$filename =~ s/\.gz$//;
	$filename =~ s/\.zip$//;
	$filename =~ s/\.fastq$//;
	my $label = "averageQualities-$counter-$filename";

	#### SET JOB
	my $job = $clusterObject->setJob([$command], $label, $outdir);
	push @$jobs, $job;
}

#### RUN JOBS
$clusterObject->runJobs($jobs, "averageQualities");

#### GET COUNTS FOR ALL READS BINNED BY AVERAGE QUALITY
my ($qualbins, $total_average_quality, $totalreads) = qualbins($files, $binarray);

#### GET AVERAGE QUALITIES PER BASE
my ($basequals, $basereads) = basequals($files, $binarray);

#### PRINT TO OUTPUT FILE
open(OUT, ">$outputfile") or die "Can't open outputfile: $outputfile\n";

#### PRINT READ COUNT BY AVERAGE QUALITY STATS
print "total reads\t$totalreads\n";
print "total average quality\t$total_average_quality\n";
print OUT "total reads\t$totalreads\n";
print OUT "total average quality\t$total_average_quality\n";
for ( my $i = 0; $i < @$binarray; $i++ )
{
	print "$$binarray[$i]\t$$qualbins[$i]\n";
	print OUT "$$binarray[$i]\t$$qualbins[$i]\n";
}


#### PRINT BASE QUALITY STATS
print "averageQualities.pl    basereads\t$basereads\n";
print OUT "base reads\t$basereads\n";
for ( my $i = 0; $i < @$basequals; $i++ )
{
	my $readnumber = $i + 1;
	print "$readnumber\t$$basequals[$i]\n";
	print OUT "$readnumber\t$$basequals[$i]\n";
}


close(OUT) or die "Can't close outputfile: $outputfile]\n";




#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "\nRun time: $runtime\n";
print "Completed $0\n";
print Util::datetime(), "\n";
print "****************************************\n\n\n";
exit;

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#									SUBROUTINES
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


=head2

	SUBROUTINE		qualbins

	PURPOSE

		GET COUNTS FOR ALL READS BINNED BY AVERAGE QUALITY

=cut

sub qualbins
{
	my $files	=	shift;
	my $binarray	=	shift;	

	#### COLLATE AVERAGE QUALITY STATS
	my $stats = [];
	#my $counter = 0;
	foreach my $inputfile ( @$files )
	{
		my $statsfile = $inputfile;
		$statsfile =~ s/\.gz$//;
		$statsfile =~ s/\.zip$//;
		$statsfile =~ s/\.fastq$//;
		$statsfile .= ".min$min" if defined $min;
		$statsfile .= ".max$max" if defined $max;
		$statsfile .= ".skip$skip" if defined $skip;
		$statsfile .= ".length$length" if defined $length;
		$statsfile .= ".qualstats";

		print "averageQualities.pl    can't find statsfile: $statsfile\n" and next if not -f $statsfile;

		#### PARSE STATSFILE AND PLACE IN HASH LOADED INTO stats ARRAY
		open(FILE, $statsfile) or die "Can't open statsfile: $statsfile\n";
		$/ = "\n";
		my @lines = <FILE>;
		close(FILE);

		#### 	EXAMPLE *.qualstats FILE FORMAT
		####
		####	head /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/SRX000600/SRR004186_2.qualstats
		####
		####	total reads     6044862
		####	total avg quality       3.62396952431143
		####	bins    -50 -40 -30 -20 -10 0 1 2 3 4 5 6 7 8 9 10 15 20 30 40 50 60 70 80 90 100 200 300
		####	-50     0
		####	-40     0
		####	-30     0
		####	-20     0
		####	-10     0
		####	0       53716
		####	1       368306

		#### EXTRACT total reads, average quality AND bins 
		my ($readcount) = $lines[0] =~ /^total\s+reads\t(\d+)$/;
		my ($avqual) = $lines[1] =~ /^total avg quality\t(\S+)$/;
		my ($bins) = $lines[2] =~ /^bins\t(.+)$/;

		print "averageQualities.pl    skipping because readcount not defined for statsfile: $statsfile\n" and next if not defined $readcount;
		print "averageQualities.pl    skipping because no readcount for statsfile: $statsfile\n" and next if not $readcount;

		#### EXTRACT READ COUNTS FOR AVERAGE QUAL BINS
		my $statsarray = [];
		for ( my $i = 3; $i < $#lines + 1; $i++ )
		{
			my ($count) = $lines[$i] =~ /\t(\d+)$/;
			push @$statsarray, $count;
		}

		my $hash = {};
		$hash->{readcount} = $readcount;
		$hash->{avqual} = $avqual;
		$hash->{bins} = $bins;	
		$hash->{values} = $statsarray;

		push @$stats, $hash;

		#$counter++;
	}

	#### GET GLOBAL TOTAL READS
	my $totalreads = 0;
	foreach my $hash ( @$stats )
	{
		$totalreads += $hash->{readcount} if defined $hash->{readcount};
	}

	#### POPULATE qualbins
	my $qualbins = [];
	foreach my $bin ( @$binarray )
	{
		push @$qualbins, 0;	
	}

	my $total_average_quality = 0;
	foreach my $hash ( @$stats )
	{

		#### ADD AVERAGE QUALITY PROPORTIONALLY TO TOTAL AVERAGE QUALITY
		my $readcount = $hash->{readcount};
		my $average_quality = $hash->{avqual};
		$total_average_quality += ($average_quality *($readcount/$totalreads)) if defined $average_quality;

		my $values = $hash->{values};

		for ( my $i = 0; $i < @$binarray; $i++ )
		{
			$$qualbins[$i] += $$values[$i];
		}
	}

	return $qualbins, $total_average_quality, $totalreads;
}


=head2

	SUBROUTINE		basequals

	PURPOSE

		CALCULATE THE AVERAGE QUALITY PER BASE ACROSS ALL READS

=cut

sub basequals
{
	my $files		=	shift;
	my $binarray	=	shift;

	#### COLLATE AVERAGE QUALITY STATS PER BASE
	my $stats = [];
	my $counter = 0;
	foreach my $inputfile ( @$files )
	{
		my $basequalfile = $inputfile;
		$basequalfile =~ s/\.gz$//;
		$basequalfile =~ s/\.zip$//;
		$basequalfile =~ s/\.fastq$//;
		$basequalfile .= ".min$min" if defined $min;
		$basequalfile .= ".max$max" if defined $max;
		$basequalfile .= ".skip$skip" if defined $skip;
		$basequalfile .= ".length$length" if defined $length;
		$basequalfile .= ".basequal";
		print "averageQualities.pl    can't find basequalfile: $basequalfile\n" and next if not -f $basequalfile;


		####	head /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/SRX000600/SRR004186_1.basequal
		####	total reads     6044862
		####	0       8.13336135713272
		####	1       8.10716158615366
		####	2       8.1007839384919
		####	3       7.80326581483581
		####	4       7.52414397549522
		####	5       7.12680785764836
		####	6       6.80795508648502
		####	7       6.9518604394939
		####	8       6.37849813610302


		#### PARSE STATSFILE AND PLACE IN HASH LOADED INTO stats ARRAY
		open(FILE, $basequalfile) or die "Can't open basequalfile: $basequalfile\n";
		$/ = "\n";
		my @lines = <FILE>;
		close(FILE);

		#### 	EXAMPLE *.qualstats FILE FORMAT
		####
		####	head /scratch/syoung/base/pipeline/solid/NA18507/SRP000239/SRX000600/SRR004186_2.qualstats
		####
		####	total reads     6044862
		####	total avg quality       3.62396952431143
		####	bins    -50 -40 -30 -20 -10 0 1 2 3 4 5 6 7 8 9 10 15 20 30 40 50 60 70 80 90 100 200 300
		####	-50     0
		####	-40     0
		####	-30     0
		####	-20     0
		####	-10     0
		####	0       53716
		####	1       368306

		#### EXTRACT total reads, average quality AND bins 
		my ($readcount) = $lines[0] =~ /^total\s+reads\t(\d+)$/;
		print "averageQualities.pl    skipping because readcount not defined for basequalfile: $basequalfile\n" and next if not defined $readcount;
		print "averageQualities.pl    skipping because no readcount for basequalfile: $basequalfile\n" and next if not $readcount;

		#### EXTRACT READ COUNTS FOR AVERAGE QUAL BINS
		my $basequalarray = [];
		for ( my $i = 1; $i < $#lines + 1; $i++ )
		{
			my ($avqual) = $lines[$i] =~ /\t(\S+)$/;
			push @$basequalarray, $avqual;
		}

		my $hash = {};
		$hash->{readcount} = $readcount;
		$hash->{values} = $basequalarray;

		push @$stats, $hash;

		#last if $counter >= 1;

		#$counter++;
	}


	#### GET GLOBAL TOTAL READS
	my $totalreads = 0;
	foreach my $hash ( @$stats )
	{
		$totalreads += $hash->{readcount} if defined $hash->{readcount};
	}

	#### POPULATE basequals
	my $basequals = [];
	foreach my $bin ( @$binarray )
	{
		push @$basequals, 0;	
	}

	foreach my $hash ( @$stats )
	{
		#### ADD AVERAGE QUALITY PROPORTIONALLY TO TOTAL AVERAGE QUALITY
		my $readcount = $hash->{readcount};
		my $values = $hash->{values};

		for ( my $i = 0; $i < @$values; $i++ )
		{
			$$basequals[$i] += ($$values[$i] *($readcount/$totalreads)) if defined $$values[$i];
		}
	}

	return $basequals, $totalreads;	
}


=head2

	SUBROUTINE		files

	PURPOSE

		GENERATE ARRAY OF HASHES CONTAINING

		INPUTFILE AND OUTPUT FILES

=cut

sub files
{
	my $inputdir		=	shift;
	my $compress		=	shift;

	#### COLLECT ALL THE INPUT FILES
	my $filearray = [];
	my @indirs = split ",", $inputdir;
	for ( my $i = 0; $i < $#indirs + 1; $i++ )
	{
		my $indir = $indirs[$i];
		my $files = Util::files($indir);
		$files = Util::by_suffix($files, "\.fastq.gz") if defined $compress and $compress eq "gzip";
		$files = Util::by_suffix($files, "\.fastq.zip") if defined $compress and $compress eq "zip";
		$files = Util::by_suffix($files, "\.fastq") if not defined $compress;
		die "No files in input directory: $inputdir\n" if not @$files or scalar(@$files) == 0;

		foreach my $file ( @$files )
		{
			push @$filearray, "$indir/$file";
		}
	}

	return $filearray;
}



sub usage
{
	print `perldoc $0`;
    exit;
}
