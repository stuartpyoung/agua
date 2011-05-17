#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();

=head2

    APPLICATION     miniMAQ

    PURPOSE

        WRAPPER SCRIPT FOR RUNNING MAQ ASSEMBLY FOR SMALL INPUT DATA SETS

    INPUT

        1. ASSEMBLY DIRECTORY

        2. FASTA OR FASTQ INPUT FILES

    OUTPUT

        1. MAQ OUTPUT FILES IN ASSEMBLY DIRECTORY

    USAGE

    ./miniMAQ.pl <--inputfile String> <--outputdir String> <--referencefile String> [--convert] [--help]

    --inputfile				:   /full/path/to/input_file.fastq (or '/../read_1.fastq,/../read_2.fastq')
    --outputdir				:   /full/path/to/output_directory
    --referencefile			:   /full/path/to/reference_file
    --clean					:   Clean run (remove old splitfile)
	--convert				:   Convert from Solexa FASTQ to Sanger FASTQ (For sequences generated using GA Pipeline v1.0 and earlier)
    --help                 	:   print help info

    EXAMPLES

/nethome/syoung/base/bin/nextgen/miniMAQ.pl \
--outputdir /mihg/data/NGS/syoung/base/pipeline/run12/lane1 \
--inputfile /mihg/data/NGS/syoung/base/pipeline/run12/lane1/s_1_1_sequence.fastq,/mihg/data/NGS/syoung/base/pipeline/run12/lane1/s_1_2_sequence.fastq \
--referencefile /mihg/data/NGS/syoung/base/pipeline/run12/lane1/CCDS_nucleotide.20080430.fa \
&> /mihg/data/NGS/syoung/base/pipeline/run12/lane1/maq-runlog.txt


/nethome/syoung/base/bin/nextgen/miniMAQ.pl \
--outputdir /mihg/data/NGS/syoung/base/pipeline/run12/lane2 \
--inputfile /mihg/data/NGS/syoung/base/pipeline/run12/lane2/s_2_1_sequence.fastq,/mihg/data/NGS/syoung/base/pipeline/run12/lane2/s_2_2_sequence.fastq \
--referencefile /mihg/data/NGS/syoung/base/pipeline/run12/lane2/CCDS_nucleotide.20080430.fa \
&> /mihg/data/NGS/syoung/base/pipeline/run12/lane2/maq-runlog.txt


/nethome/syoung/base/bin/nextgen/miniMAQ.pl \
--outputdir /mihg/data/NGS/syoung/base/pipeline/run12/lane1 \
--inputfile /mihg/data/NGS/syoung/base/pipeline/run12/lane1/s_1_1_sequence.fastq,/mihg/data/NGS/syoung/base/pipeline/run12/lane1/s_1_2_sequence.fastq \
--referencefile /mihg/data/NGS/syoung/base/pipeline/run12/lane1/CCDS_nucleotide.20080430.fa \
&> /mihg/data/NGS/syoung/base/pipeline/run12/lane1/maq-runlog.txt



/nethome/syoung/base/bin/nextgen/miniMAQ.pl \
--outputdir /nethome/syoung/base/pipeline/comparison/readfiles \
--inputfile /nethome/syoung/base/pipeline/comparison/s_4_sequence.txt \
--referencefile /nethome/syoung/base/pipeline/comparison/mtDNA-AC_000021.fa


/nethome/syoung/base/pipeline/comparison/s_4_sequence.txt


=cut

use strict;

#### EXTERNAL MODULES
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use Getopt::Long;
use FindBin qw($Bin);

#### USE LIBRARY
use lib "$Bin/../../../lib";	

#### INTERNAL MODULES
use Sampler;
use Timer;
use Util;
use Conf::Agua;

##### STORE ARGUMENTS TO PRINT TO FILE LATER
my @arguments = @ARGV;

#### FLUSH BUFFER
$| =1;

#### SET MAQ LOCATION
my $conf = Conf::Agua->new(inputfile=>"$Bin/../../../conf/default.conf");
my $maq = $conf->getKeyValue("agua", 'MAQ');
print "maq: $maq\n";
#exit;
my ($maqbin) = $maq =~ /^(.+)\/[^\/]+$/;

#### DEFAULT MAX FILE LINES (MULTIPLE OF 4 BECAUSE EACH FASTQ RECORD IS 4 LINES)
my $max_lines = 4000000;

#### GET OPTIONS
my $inputfile;
my $readfile;
my $outputdir;
my $convert;
my $clean;
my $referencefile;
my $help;
if ( not GetOptions (
    'inputfile=s' => \$inputfile,
    'outputdir=s' => \$outputdir,
    'referencefile=s' => \$referencefile,
    'convert' => \$convert,
    'clean' => \$clean,
    'help' => \$help
) )
{ print "Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "Input file not defined (Use --help for usage)\n" if not defined $inputfile;
die "Output directory not defined (Use --help for usage)\n" if not defined $outputdir;
die "Reference file not defined (Use --help for usage)\n" if not defined $referencefile;

print "inputfile: $inputfile\n";
print "outputdir: $outputdir\n";
print "referencefile: $referencefile\n";

#### CHECK INPUTS
if ( not defined $inputfile)   {   print "Input file not defined (option i)\n";    usage();    }
if ( not defined $outputdir)   {   print "Output directory not defined (option d)\n";    usage();    }
if ( not defined $referencefile )   {   print "Reference file not defined (option r)\n";    usage();    }

#### MAKE OUTPUT DIR
mkdir($outputdir) if not -d $outputdir;


print "Sleeping for 15 seconds.\n";
sleep(15);
exit;



#### CHECK INPUT FILES
my $inputfiles;
@$inputfiles = split ",", $inputfile;
die "Max two input files (paired end) allowed\n" if scalar(@$inputfiles) > 2;

#### STORE COMMANDS HERE
my $commands = [];

#### CONVERT REFERENCE FILE INTO BINARY FILE
my $referencebinary = $referencefile;
$referencebinary =~ s/\.[^\.]+?$/.bfa/;
push @$commands, "time $maq fasta2bfa $referencefile $referencebinary"
	if not -f $referencebinary;

#### SPLIT INPUT FILE IF > MAX_LINES
my $splitfiles;
my $namesfile = "$outputdir/splitfiles";
print "namesfile: $namesfile\n";
if ( not defined $clean and ( -f $namesfile and not -z $namesfile ) )
{
	print "Getting split files from namesfile: $namesfile ...\n";
    $splitfiles = Sampler::get_splitfiles($namesfile);
}
else
{
	print "Generating split files...\n";
	$splitfiles = Sampler::set_splitfiles(
		{
			'inputfiles' 	=> 	$inputfile,
			'lines'			=> 	$max_lines,
			'splitfile'		=>	$namesfile,
			'outputdir'		=>	$outputdir,
			'suffix'		=> "fastq"
		}
	);	
	print "Names file printed: $namesfile\n";
}

print "inputfiles:\n";
print join "\n", @$inputfiles;
print "\n";
print "splitfiles:\n";
foreach my $splitfile ( @$splitfiles )
{
    print join "\n", @$splitfile;
    print "\n";    
}

#### DO THE MAQ match ALIGNMENT TO PRODUCT THE out.map FILES
for ( my $index = 0; $index < @$splitfiles; $index++ )
{
	my ($outputdir) = $splitfiles->[$index][0] =~ /^(.+?)\/[^\/]+$/;

	my $args;
	$args->{inputfiles} = $$splitfiles[$index];
	$args->{referencefile} = $referencefile;	
	$args->{referencebinary} = $referencebinary;	
	$args->{outputdir} = $outputdir;
	$args->{convert} = $convert;

	my $subcommands = align($args);
	die "No commands returned for splitfiles: $splitfiles->[$index][0], $splitfiles->[$index][1]\n" if not defined $subcommands;
	@$commands = ( @$commands, @$subcommands);
}

#### CHANGE BACK TO OUTPUT DIRECTORY
chdir($outputdir);
push @$commands, "cd $outputdir";


#### GENERATE LIST OF out.map FILES FOR
#### MERGING INTO ONE out.map FILE
my $submaps;
for ( my $index = 0; $index < @$splitfiles; $index++ )
{
	my ($subdir) = $splitfiles->[$index][0] =~ /^(.+?)\/[^\/]+$/;
	my $submap = "$subdir/out.map";
	push @$submaps, $submap;
}


#### CHECK THAT ALL FILES EXIST
my $missedmaps;
for ( my $i = 0; $i < @$submaps; $i++ )
{
	my $submap = $$submaps[$i];
	if ( not -f $submap )
	{
		push @$missedmaps, splice (@$submaps, $i, 1);
	}
}

#### PRINT ANY MISSING out.map FILES
if ( defined $missedmaps and scalar(@$missedmaps) > 0 )
{
	print "Missed maps:\n";
	print join "\n", @$missedmaps;
	print "\n";
}

#### MERGE ALL SUBMAPS
my $outmap = "$outputdir/out.map";
push @$commands, "echo 'out.map: $outmap'";
my $merge_command = qq{time $maq mapmerge $outmap @$submaps}; 
push @$commands, "echo '$merge_command'";
push @$commands, $merge_command;

# 2. Statistics from the alignment
push @$commands, "time $maq mapcheck $referencefile out.map > mapcheck.txt";

#######################
##### DO SNPS
#######################

# 3. Build the mapping assembly
push @$commands, "time $maq assemble consensus.cns $referencebinary out.map 2> assemble.log";

# 4. Extract consensus sequences and qualities
push @$commands, "time $maq cns2fq consensus.cns > cns.fq";

# 5. Extract list of SNPs 
push @$commands, "time $maq cns2snp consensus.cns > cns.snp";


#######################
##### DO INDELS
#######################

#2. rmdup       remove pairs with identical outer coordinates (PE)
push @$commands, "time $maq rmdup out.rmdup out.map";

#3. indelpe     indel calling (PAIRED READS ONLY)
push @$commands, "time $maq indelpe $referencebinary out.rmdup > out.indelpe";

#4. indelsoa    state-of-art homozygous indel detectionll
push @$commands, "time $maq indelsoa $referencebinary out.map > out.indelsoa";

#5. filter indels
push @$commands, "awk '\$5+\$6-\$4 >= 3 \&\& \$4 <= 1' out.indelsoa > out.indelsoa.filter";

#6. SNPfilter    filter SNP predictions
push @$commands, "time $maqbin/scripts/maq.pl SNPfilter -d 1 -s out.indelsoa -F out.indelpe cns.snp &> out.SNPfilter";


######################
#### RUN COMMANDS
######################

foreach my $command ( @$commands )
{
	print "\n";
    print "$command\n";
    print `$command`;
}


#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "\nRun time: $runtime\n";
print "Completed $0\n";
print Timer::datetime(), "\n";
print "****************************************\n\n\n";
exit;

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#									SUBROUTINES
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



=head2

	SUBROUTINE		align

	PURPOSE

		ALIGN SEQUENCE AGAINST REFERENCE

=cut

sub align
{
	my $args		=	shift;

	#### GET VARIABLES
	my $inputfiles		=	$args->{inputfiles};
	my $referencefile	=	$args->{referencefile};
	my $referencebinary	=	$args->{referencebinary};
	my $outputdir		=	$args->{outputdir};
	my $convert			=	$args->{convert};

	foreach my $inputfile ( @$inputfiles )
	{
		Util::checkfile("$inputfile", "Can't find input file: $inputfile\n", \&usage);
		if ( $inputfile !~ /\.(txt|fastq|bfq)$/ )
		{
			print "input file must end in .txt, .fastq or .bfq\n";
			usage();
		}
	}

	#### CHECK REF FILE FORMAT
	if ( $referencefile !~ /(fa|fas|fasta)$/ )
	{
		print "Reference file must end in .fa, .fas or .fasta\n";
		usage();
	}

	#### CHECK FILES AND DIRS
	die "Can't find output dir: $outputdir\n" if not -d $outputdir;
	die "Can't find reference file: $referencefile\n" if not -f $referencefile;

	#### CHANGE TO OUTPUT DIRECTORY
	chdir($outputdir);

	#### KEEP THE COMMANDS HERE
	my $commands;

	########################
	#### PREPARE INPUT FILES
	########################

	# 1. CONVERT SOLEXA QUALITY VALUES TO SANGER FASTQ QUALITY VALUES
	if ( defined $convert )
	{
		push @$commands, "echo 'Converting solexa sequence files to Sanger fastq files...'\n";
		foreach $inputfile ( @$inputfiles )
		{
			my $solexafile = $inputfile . ".solexa";
			push @$commands, "mv $inputfile $solexafile";

			#### JUST IN CASE INPUT FILE IS .txt
			my $fastqfile = $inputfile;
			$fastqfile =~ s/\.txt$/.fastq/;
			push @$commands, "maq sol2sanger $solexafile $fastqfile";

			#### UPDATE NAME OF INPUT FILE
			$inputfile = $fastqfile;
		}
	}


	# 2. CONVERT FASTQ FILE TO BFQ FILE FORMAT
	push @$commands, "echo 'Converting solexa Sanger fastq files to MAQ binary fastq (.bfq) files...'";
	foreach $inputfile ( @$inputfiles )
	{
		my $bfqfile = $inputfile;
		$bfqfile =~ s/\.[^\.]+?$/.bfq/;
		push @$commands, "time $maq fastq2bfq $inputfile $bfqfile";		

		#### UPDATE NAME OF INPUT FILE
		$inputfile = $bfqfile;
	}

	# 3. SET INPUT BINARIES
	my $inputbinaries = join " ", @$inputfiles;

	####################
	#### DO THE ASSEMBLY
	####################

	# 4. Align the reads to the reference
	push @$commands, "time $maq match $outputdir/out.map $referencebinary $inputbinaries";

	return $commands;
}




#=head2
#
#	SUBROUTINE		splitfiles
#	
#	PURPOSE
#	
#		splitfiles SEQUENCE AGAINST REFERENCE
#		
#=cut
#
#sub splitfiles
#{
#	my $inputfiles		=	shift;
#	my $max_lines		=	shift;
#	my $namesfile		=	shift;
#
#
#	open(NAMES, ">$namesfile") or die "Can't open names file: $namesfile\n";
#	
#	my $splitfiles;
#	for ( my $index = 0; $index < @$inputfiles; $index++ )
#	{
#		my $inputfile = $$inputfiles[$index];
#		
#		my $filenumber = 1;
#
#		my ($basedir, $filename) = $inputfile =~ /^(.+?)\/([^\/]+)$/;
#		$filename =~ s/\.([^\.]+?)$//;
#		#print "basedir: $basedir\n";
#		#print "filename: $filename\n";
#
#		#### CREATE SUB DIRECTORY FOR SPLIT FILE
#		my $subdir = "$outputdir/$filenumber";
#		mkdir($subdir) or die "Can't create output directory: $subdir\n" if not -d $subdir;
#
#		#### OPEN OUTPUT FILE 
#		my $outputfile .= "$subdir/$filename.$filenumber.fastq";
#		$splitfiles->[$filenumber - 1][$index] = $outputfile;
#		
#		#print "\$splitfiles->[$filenumber - 1][$index] = $outputfile\n";
#		open(OUTFILE, ">$outputfile") or die "Can't open output file: $outputfile\n";
#			
#		my $counter = 0;
#		open(FILE, $inputfile) or die "Can't open input file: $inputfile\n";
#		while ( <FILE> )
#		{
#			if ( $counter >= $max_lines )
#			{
#				close(OUTFILE);
#				$filenumber++;
#
#				#### CREATE SUB DIRECTORY FOR SPLIT FILE
#				$subdir = "$outputdir/$filenumber";
#				mkdir($subdir) or die "Can't create output directory: $subdir\n" if not -d $subdir;
#
#				$outputfile = "$subdir/$filename.$filenumber.fastq";
#				$splitfiles->[$filenumber - 1][$index] = $outputfile;
#
#
#				#print "\$splitfiles->[$filenumber - 1][$index] = $outputfile\n";
#				open(OUTFILE, ">$outputfile") or die "Can't open output file: $outputfile\n";
#				
#				$counter = 0;
#			}
#
#			$counter++;
#		}
#		close(FILE);
#	}
#	close(NAMES);
#
#	return $splitfiles;
#}


#
#=head2
#
#	SUBROUTINE		lines
#	
#	PURPOSE
#	
#		COUNT THE NUMBER OF LINES IN A FILE
#
#=cut
#
#sub lines
#{
#	my $filename		=	shift;
#	
#	
#	open(FILE, $filename);
#	my $lines = 0;
#	while(<FILE>)
#	{
#		$lines++;
#	}
#	close(FILE);
#	
#	return $lines;
#}
#	


sub usage
{
	print GREEN;
	print `perldoc $0`;
	print RESET;
	exit;
}


