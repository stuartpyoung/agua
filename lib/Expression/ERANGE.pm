package ERANGE;


=head2

		PACKAGE		ERANGE

		PURPOSE

	        WRAPPER SCRIPT FOR RUNNING ERANGE TRANCRIPTOME PREDICTION

=cut

use strict;
use warnings;
use Carp;

require Exporter;
#our @ISA = qw(Exporter Common);
#our @EXPORT_OK = qw();
our $AUTOLOAD;

#### INTERNAL MODULES
use Util;
use Sampler;

#### EXTERNAL MODULES
use Data::Dumper;
use File::Path;
use File::Remove;
#use MD5;

#### SET SLOTS
our @DATA = qw(

LABEL
INPUTFILE
MATEFILE
KNOWNGENE
OUTPUTDIR
CLEAN
SPECIES
FILETYPE

ERANGE
PYTHON
SQLITE
REPMASK

QSTAT
JOBS
CPUS
SLEEP
DOT
QSUB
QUEUE
);
our $DATAHASH;
foreach my $key ( @DATA )	{	$DATAHASH->{lc($key)} = 1;	}



=head2

	SUBROUTINE		runRNAPairedAnalysis

	PURPOSE

		RUN THE ERANGE PAIRED READ ANALYSIS SCRIPT runRNAPairedAnalysis.sh

		NB: USAGE IS AS FOLLOWS -

			runRNAPairedAnalysis.sh genome rdsprefix repeatmaskdb

=cut

sub runRNAPairedAnalysis
{
	my $self			=	shift;


	#### USER INPUTS
	my $label 			=	$self->get_label();
	my $species 		=	$self->get_species();
	my $inputfile 		=	$self->get_inputfile();
	my $matefile 		=	$self->get_matefile();
	my $knowngene 		=	$self->get_knowngene();
	my $outputdir 		=	$self->get_outputdir();
	my $clean 			=	$self->get_clean();
	my $filetype 		=	$self->get_filetype();

	#### EXECUTABLES
	my $python 			=	$self->get_python();
	my $sqlite 			=	$self->get_sqlite();
	my $erange 			=	$self->get_erange();
	my $repmask 		=	$self->get_repmask();

	#### CLUSTER
	my $qstat 			=	$self->get_qstat();
	my $jobs 			=	$self->get_jobs();
	my $sleep 			=	$self->get_sleep();
	my $qsub 			=	$self->get_qsub();
	my $queue 			=	$self->get_queue();

	#### SET PAIRED FLAG
	my $paired = 0;
	$paired = 1 if $matefile;

	#### SET ENVIRONMENT VARIABLES
	#### set PYTHONPATH to point to the parent directory of the Cistematic, e.g.
	$ENV{'PYTHONPATH'} = $erange;
	#### set CISTEMATIC_ROOT to the directory that contains the genome directories
	#### (such as H_sapiens or M_musculus)
	$ENV{'CISTEMATIC_ROOT'} = $erange;
	#### set ERANGEPATH
	$ENV{'ERANGEPATH'} = "$erange/commoncode";

	chdir($outputdir) or die "Can't change to output directory: $outputdir\n";

	my $command = "$erange/commoncode/runRNAPairedAnalysis.sh $species $label $repmask";
	print "ERANGE::runRNAPairedAnalysis    command:\n";
	print "$command\n";

}	#	runRNAPairedAnalysis



=head2

	SUBROUTINE		makeRdsFromEland

	PURPOSE

		CREATE RDS FILE WITH ELAND ALIGNMENT USING makerdsfromeland2.py

		makerdsfromeland2.py USAGE:

			REQUIRED: a label, the eland read file, the output file, the known gene table.

			OPTIONAL: for paired reads, -paired and # of pair.

		NB: ADD THE -index ARGUMENT FOR THE LAST INPUT FILE (OR INPUT + MATE PAIR).

	EXAMPLE

/nethome/apps/bioinfo/python/2.6.1/bin/python /nethome/syoung/base/pipeline/erange2/commoncode/makerdsfromeland2.py s_1_1 /nethome/syoung/base/pipeline/erange2/data/eland/s_1_1_eland_multi.txt /nethome/syoung/base/pipeline/erange2/data/rds/s_1_1.rds -RNA /nethome/syoung/base/pipeline/erange2/data/knownGene/knownGene.hg18 -paired 1 -verbose

=cut

sub makeRdsFromEland
{
	my $self			=	shift;


	#### USER INPUTS
	my $label 			=	$self->get_label();
	my $species 		=	$self->get_species();
	my $inputfile 		=	$self->get_inputfile();
	my $matefile 		=	$self->get_matefile();
	my $knowngene 		=	$self->get_knowngene();
	my $outputdir 		=	$self->get_outputdir();
	my $clean 			=	$self->get_clean();
	my $filetype 		=	$self->get_filetype();

	#### EXECUTABLES
	my $python 			=	$self->get_python();
	my $sqlite 			=	$self->get_sqlite();
	my $erange 			=	$self->get_erange();

	#### CLUSTER
	my $qstat 			=	$self->get_qstat();
	my $jobs 			=	$self->get_jobs();
	my $sleep 			=	$self->get_sleep();
	my $qsub 			=	$self->get_qsub();
	my $queue 			=	$self->get_queue();

	#### SET PAIRED FLAG
	my $paired = 0;
	$paired = 1 if defined $matefile and $matefile !~ /^\s*$/;

	#### CHANGE TO OUTPUT DIR
	chdir($outputdir) or die "Can't change to output directory: $outputdir\n";

	#### SET .rds FILE 		
	my $rdsfile = "$outputdir/$label.rds";

	#### RETURN IF .rds FILE ALREADY EXISTS AND clean FLAG NOT SET
	print "RDS file already exists: $rdsfile\n" and return if -f $rdsfile and not $clean;

	#### REMOVE EXISTING RDS FILE IF clean FLAG IS SET
	if ( $clean and -f $rdsfile )
	{
		my $recursive = 0;
		File::Remove::rm(\$recursive, $rdsfile) or die "Can't remove .rds file: $rdsfile\n";
	}

	#### CREATE RDS FILE
	#### http://woldlab.caltech.edu/rnaseq/README.build-rds
	my $inputfile_command = qq{time $python \\
$erange/commoncode/makerdsfromeland2.py \\
$label $inputfile \\
$rdsfile \\
-RNA $knowngene };

	#### ADD -paired
	$inputfile_command .=  qq{\\\n -paired 1 } if $paired;

	#### ADD -extended
	$inputfile_command .=  qq{\\\n -extended } if $filetype eq "extended";


	print "ERANGE::build_rds    inputfile_command:\n\n$inputfile_command\n\n";
	print `$inputfile_command`;


	#### ADD ITS COMPLEMENTARY MATE FILE IF THIS IS A PAIRED-END ANALYSIS
	if ( $paired )
	{
		my $matefile_command = qq{time $python \\
$erange/commoncode/makerdsfromeland2.py \\
$label $inputfile \\
$rdsfile \\
-RNA $knowngene \\
-paired 2 -append -index };

		#### ADD -extended 
		$matefile_command .=  qq{\\\n -extended } if $filetype eq "extended";
		print "ERANGE::build_rds    matefile_command:\n\n$matefile_command\n\n";
		print `$matefile_command`;

	}
}






=head2

	SUBROUTINE		makeRdsFromBowtie

	PURPOSE

		CREATE RDS FILE WITH ELAND ALIGNMENT USING makerdsfrombowtie2.py

		makerdsfrombowtie2.py USAGE:

			REQUIRED: a label, the eland read file, the output file, the known gene table.

			OPTIONAL: for paired reads, -paired and # of pair.

		NB: ADD THE -index ARGUMENT FOR THE LAST INPUT FILE (OR INPUT + MATE PAIR).

	NOTES

		6. MAPPING READS WITH BOWTIE
		http://woldlab.caltech.edu/rnaseq/README.build-rds

		Bowtie (bowtie-bio.sourceforge.net) is a new read-mapper that 
		is very fast and friendly. Version 0.9.8.1 and higher allow you 
		to control how many multireads are reported. We use the 
		following settings to emulate Eland output:

		$BOWTIEDIR/bowtie zzz -v 2 -k 11 -m 10 -t --best -f s1.32mer.query.txt --unfa s1.unmapped.fa --maxfa s1.repeat.fa s1.zzz.bowtie.txt

		where zzz is the genome prefix that you gave when building the 
		genome. In particular, we ask bowtie to map all multireads up 
		to 11 ("-k") with up to 2 mismatches ("-v" and "--best"), however 
		we will only import all multireads up to 10x multiplicity ("-m").
		Note that bowtie is multithreaded and can use multiple cpu based 
		on the -p flag (e.g. use "-p 4" to use 4 CPUs). Unmapped reads 
		are saved in unmapped.fa for later analysis.

	EXAMPLE

/nethome/apps/bioinfo/python/2.6.1/bin/python /nethome/syoung/base/pipeline/erange2/commoncode/makerdsfrombowtie2.py s_1_1 /nethome/syoung/base/pipeline/erange2/data/eland/s_1_1_eland_multi.txt /nethome/syoung/base/pipeline/erange2/data/rds/s_1_1.rds -RNA /nethome/syoung/base/pipeline/erange2/data/knownGene/knownGene.hg18 -paired 1 -verbose

=cut

sub makeRdsFromBowtie
{
	my $self			=	shift;


	#### USER INPUTS
	my $label 			=	$self->get_label();
	my $species 		=	$self->get_species();
	my $inputfile 		=	$self->get_inputfile();
	my $matefile 		=	$self->get_matefile();
	my $knowngene 		=	$self->get_knowngene();
	my $outputdir 		=	$self->get_outputdir();
	my $clean 			=	$self->get_clean();
	my $filetype 		=	$self->get_filetype();

	#### EXECUTABLES
	my $python 			=	$self->get_python();
	my $sqlite 			=	$self->get_sqlite();
	my $erange 			=	$self->get_erange();

	#### CLUSTER
	my $qstat 			=	$self->get_qstat();
	my $jobs 			=	$self->get_jobs();
	my $sleep 			=	$self->get_sleep();
	my $qsub 			=	$self->get_qsub();
	my $queue 			=	$self->get_queue();

	#### SET PAIRED FLAG
	my $paired = 0;
	$paired = 1 if defined $matefile and $matefile !~ /^\s*$/;

	#### CHANGE TO OUTPUT DIR
	chdir($outputdir) or die "Can't change to output directory: $outputdir\n";

	#### SET .rds FILE 		
	my $rdsfile = "$outputdir/$label.rds";

	#### RETURN IF .rds FILE ALREADY EXISTS AND clean FLAG NOT SET
	print "RDS file already exists: $rdsfile\n" and return if -f $rdsfile and not $clean;

	#### REMOVE EXISTING RDS FILE IF clean FLAG IS SET
	if ( $clean and -f $rdsfile )
	{
		my $recursive = 0;
		File::Remove::rm(\$recursive, $rdsfile) or die "Can't remove .rds file: $rdsfile\n";
	}

	#### CREATE RDS FILE
	#### http://woldlab.caltech.edu/rnaseq/README.build-rds
	my $inputfile_command = qq{time $python \\
$erange/commoncode/makerdsfromeland2.py \\
$label $inputfile \\
$rdsfile \\
-RNA $knowngene };

	#### ADD -paired
	$inputfile_command .=  qq{\\\n -paired 1 } if $paired;

	#### ADD -extended
	$inputfile_command .=  qq{\\\n -extended } if $filetype eq "extended";


	print "ERANGE::build_rds    inputfile_command:\n\n$inputfile_command\n\n";
	print `$inputfile_command`;


	#### ADD ITS COMPLEMENTARY MATE FILE IF THIS IS A PAIRED-END ANALYSIS
	if ( $paired )
	{
		my $matefile_command = qq{time $python \\
$erange/commoncode/makerdsfromeland2.py \\
$label $inputfile \\
$rdsfile \\
-RNA $knowngene \\
-paired 2 -append -index };

		#### ADD -extended 
		$matefile_command .=  qq{\\\n -extended } if $filetype eq "extended";
		print "ERANGE::build_rds    matefile_command:\n\n$matefile_command\n\n";
		print `$matefile_command`;

	}
}



=head2 


http://woldlab.caltech.edu/rnaseq/README.build-rds


This is a description of the sqlite-based read storage 
files and of the scripts designed to import read 
mappings from supported short read mappers.  The code 
should run on any Unix-like system supporting python 2.5 
or better. The code is developed on Linux and MacOS X on 
python 2.5.

This code is made available as open-source, as described 
in the copyright file ERANGE.COPYRIGHT.

1. REQUIREMENTS
2. COMMAND LINE OPTIONS
3. CREATING THE NECESSARY INPUT (RDS) FILES
4. BUILDING EXPANDED GENOMES
5. MAPPING READS WITH ELAND
6. MAPPING READS WITH BOWTIE
7. MAPPING READS WITH BLAT
8. IMPORTING BED FILES
9. MANIPULATING RDS METADATA AND CACHING
10. VISUALIZING THE DATA IN RDS FILES


1. REQUIREMENTS 

See README.chip-seq or README.rna-seq to see the requirements 
for installing and running ERANGE specific to each 
application. 


2. COMMAND LINE OPTIONS

You can find out more about the settings for each script
by typing:

python $ERANGEPATH/<scriptname> 

to see the command line options, where ERANGEPATH is the 
environmental variable set to the path to the directory 
holding the ERANGE scripts. Note that the command line 
options are case sensitive and that they could well 
fail silently. 


3. CREATING THE NECESSARY INPUT (RDS) FILES

Before you can use the rest of the ERANGE scripts to do 
CHiP-seq or RNA-seq analyses, you will need to first 
convert your read mappings to the native ERANGE read 
storage format, which is sqlite-based, and which is 
called RDS (Read DataSet). RDS files consist of four 
tables:
- metadata (tracks required and optional metadata)
- uniqs (stores uniquely mappable-reads)
- multi (stores reads that map equally well to multiple 
locations in the genome)
- splices (stores split reads)

a readDataset python object (in commoncode.py) provides 
the encapsulation of the read database which is accessed 
through specific methods. Since an RDS file is a sqlite3 
database, you can additionally use any of the sqlite-based 
tools to look at the reads in the tables, if you wish to 
do so.

You will need to first map your reads with one of the 
supported read mappers (see next paragraph) against a copy 
of the appropriate genome. For ChIP-seq, it will be your 
genome of interest, whereas for RNA-seq reads should be 
mapped against an expanded genome, which consists of 
chromosomes + splice junctions which depend on the read 
length used. Note that several parts of the code assume 
that your genomic sequences are labelled with the "chr" 
chromosomes prefix. For more information on creating 
expanded genomes, see BUILDING EXPANDED GENOMES.

The currently supported read mappers are:
- Eland (part of the Illumina GA pipeline)
- Bowtie (bowtie-bio.sourceforge.net)
- Blat (from UCSC)

These are described in the sections on MAPPING READS WITH 
ELAND, MAPPING READS WITH BOWTIE, MAPPING READS WITH BLAT. 

For ChIP-seq, you can also import bed files of unique reads 
only using makerdsfrombed.py .

Also see MANIPULATING RDS METADATA AND CACHING to learn about 
some important aspects of working with RDS files.


4. BUILDING EXPANDED GENOMES

For RNA-seq using ELAND or BOWTIE mappings, you will need to build 
an expanded genome consisting of genomic sequences, spike sequences, 
and splice-spanning sequences in order to run ERANGE on your own 
datasets. This expanded genome is specific to the read size used, 
i.e. there will be a different expanded genome for mouse when using 
25bp reads or 32bp reads. For reads longer than 32 bp, we recommend 
using BOWTIE. If your reads are longer than 50bp, consider using 
BLAT instead.

Download the chromosomes from UCSC, as well as the knownGene.txt (or 
equivalent table) and a directory of repeatmask annotations for each 
chromosome (also from UCSC) for your genome of interest.

You will need to build a splice fasta file using the script 
getsplicefa.py, which needs Cistematic, the knownGene table, and a 
paremeter for splice radius, which is 4 bp shorter than the length 
of the reads.

Once you have the splice fasta file, drop it into the same directory 
as well as a fasta file for your spikes. Then use squashGenome 
(part of Eland) or bowtie-build (part of Bowtie), to build the 
expanded genome. Please refer to the documentation for each 
package to run the genome squasher/builder.

You will also build a repeat database using buildrmaskdb.py for use 
in the candidate exon analysis from UCSC repeatmasker annotations.


5. MAPPING READS WITH ELAND

Please refer to the Illumina documentation for the details on 
running squashGenome and Eland. If you do not have access to the 
Illumina pipeline, use bowtie as described in the next section.

For ChIP-seq, you could take the output of the Illumina pipeline, 
e.g. eland_multi.txt or eland_extended.txt and use them as inputs 
for makerdsfromeland2.py .

Once you have run Eland with the --multi option (which we 
colloquially call "eland2") for each RNA-seq lane against the 
expanded genome, combine all of the outputs for one sample into a 
single file e.g. test.comb.eland2

The makerdsfromeland2.py script is used to import the reads 
into RDS:

python makerdsfromeland2.py label infilename outrdsfile [-append] [-RNA ucscGeneModels] 
[propertyName::propertyValue] [-index] [-paired 1 or 2] [-extended] [-verbose] 
[-olddelimiter] [-maxlines num] [-cache numPages]

The first 3 arguments are required:
- label is any label that you wish (a combination flowcell+lane# 
is a good choice) 
- infilename is the output of eland in eland_multi format 
(default) or eland_extended format (with the -extended flag)
- outdbname is the name of the rds file, e.g. test.rds

If the reads are from paired-end runs, enter each eland_multi 
(or extended) file separately with the "-paired 1" or "-paired 2" 
flag, as appropriate.

If entering more than one lane, use -append for all subsequent 
lanes. Upon entering the last lane, use -index to build a read 
index. Refer to MANIPULATING RDS METADATA AND CACHING for 
information on the optional property::value pairs and caching.

For RNA-seq, you must in addition specify the path to knownGene.txt 
using the -RNA flag, e.g.

python $ERANGEPATH/makerdsfromeland2.py myRNAlabel myRNA.eland_multi.txt rnatest.rds -RNA ../mm9/knownGene.txt [more options]


6. MAPPING READS WITH BOWTIE

Bowtie (bowtie-bio.sourceforge.net) is a new read-mapper that 
is very fast and friendly. Version 0.9.8.1 and higher allow you 
to control how many multireads are reported. We use the 
following settings to emulate Eland output:

$BOWTIEDIR/bowtie zzz -v 2 -k 11 -m 10 -t --best -f s1.32mer.query.txt --unfa s1.unmapped.fa --maxfa s1.repeat.fa s1.zzz.bowtie.txt

where zzz is the genome prefix that you gave when building the 
genome. In particular, we ask bowtie to map all multireads up 
to 11 ("-k") with up to 2 mismatches ("-v" and "--best"), however 
we will only import all multireads up to 10x multiplicity ("-m").
Note that bowtie is multithreaded and can use multiple cpu based 
on the -p flag (e.g. use "-p 4" to use 4 CPUs). Unmapped reads 
are saved in unmapped.fa for later analysis.

Once reads are mapped, they can be imported using:

python $ERANGEPATH/makerdsfrombowtie.py testLabel s1.mm9.bowtie.txt bowtietest.rds

The options for the script are:

python makerdsfrombowtie.py label infilename outrdsfile 
[-RNA ucscGeneModels]  [-append]  [-index] [propertyName::propertyValue] 
[-rawreadID] [-verbose] [-cache numPages]

Refer to "MAPPING READS WITH ELAND" for a description of label, 
infilename, outdbname, '-append', '-index', and '-cache'.

****REMEMBER TO USE -index WHEN LOADING THE LAST LANE OF YOUR 
DATASET.****

The script assumes that the read ID are from Illumina, i.e. that 
they have multiple fields separated by ':' and that paired-end 
reads have an additional '/1' or '/2' depending on the end. 
It will by default strip the first part of the readID (up to the 
first ':') and replace it with the label. If you want raw readIDs
because you mapped raw reads that do not have an associated ID or 
an ID that doesn't follow Illumina's conventions, use -rawreadID.

If not using Illumina readIDs, use any identifier of the format

throw_away:uniqueid if unpaired
throw_away:uniqueid/1 and throw_away:uniqueid/2 for paired-ends.

For RNA-seq, you must in addition specify the path to knownGene.txt 
using the -RNA flag, e.g.

python $ERANGEPATH/makerdsfrombowtie.py myRNAlabel myRNA.bowtie.txt rnatest.rds -RNA ../mm9/knownGene.txt [more options]


7. MAPPING READS WITH BLAT

BLAT SUPPORT IN ERANGE IS STILL UNDER DEVELOPMENT AND THE 
SCRIPTS AND SETTINGS BELOW MAY BE OPTIMIZED FURTHER IN 
FUTURE RELEASES OF ERANGE.

Reads longer than 40-50bp can be fruitfully mapped with BLAT 
against the reference genome without needing to provide the 
exon junctions. While BLAT is much slower than BOWTIE, it 
has the great advantage of seeing novel splices (i.e. 
splices not present in knownGene models).

We use the following settings to map 75bp reads with BLAT and 
filter them with pslReps: 

$BLATPATH/blat /tmp/hg18.fa s3_1.query75.txt -out=pslx s3_1.hg18.blat
$BLATPATH/pslReps -minNearTopSize=70 s3_1.hg18.blat s3_1.hg18.blatbetter s3_1.blatpsr

where the binaries are in $BLATPATH anywhere on your system.

Once the reads have been filtered, the makerdsfromblat.py 
script is used to import the mapped reads (in the example 
above s3_1.hg18.blatbetter) into RDS:

python makerdsfromblat.py label infilename outrdsfile [-append] [-index] [propertyName::propertyValue] 
[-rawreadID] [-forceRNA]  [-flag] [-strict minSpliceLen] [-spliceonly] [-verbose] [-cache numPages]

If you are using BLAT for RNA-seq, please be sure to use
-forceRNA in order to import spliced reads and consider 
using -strict to require a minimum length of bases on 
each side of the splice. 

You can combine BOWTIE and BLAT by mapping reads with BOWTIE 
first, and then using BLAT to map the unmapped reads. In 
that case, you may want to only load the spliced reads 
using the -spliceonly flag. To track those reads in the RDS 
file, use -flag ; you can then retrieve those reads using 
the options "-flag blat -flagLike" with the makebedfromrds.py 
script.


8. IMPORTING BED FILES

If you do not have the raw read data, you can import unique 
reads only using the script makerdsfrombed.py . Note that 
this is not particularly useful for RNA-seq since you will 
have neither the multireads nor the spliced reads.

The command line options are similar to those for other 
scripts described in part 5-7:

python makerdsfrombed.py label bedfile outrdsfile [-append] [-index] [propertyName::propertyValue] [-cache numPages]


9. MANIPULATING RDS METADATA AND CACHING

One of the advantages of RDS over bed, is the possibility of 
attaching arbitrary sets of annotations with the data, which 
are then carried along. Both the makerds* scripts and 
rdsmetadata.py allows you to both enter key::value 
combinations. Entering a key multiple times will cause the 
same instance to be recorded multiple times, which is 
appropriate in some settings (e.g. to enter flowcell info). 
In addition rdsmetadata.py allows you to inspect various 
attributes of your RDS files such as # of reads and size 
of the default cache size.

Sqlite files have a certain amount of RAM set aside as cache 
for lookups, indexes, etc.... where the amount is measured in 
1.5kb pages. Each RDS instance come with a default of 100000 
pages (150MB) of cache, which is needlessly small in most 
situations. Whenever appropriate, try using more cache (e.g. 
750000 pages on a 2GB RAM machine, much more if more RAM is 
available) for a significant speed increase in indexing and 
lookups. You can change the default value for each RDS file 
by using the -defaultcache option of rdsmetadata.py.

Note that sqlite can be very slow over NFS. Wherever 
possible, copy your RDS file locally before running an I/O 
intensive script.


10. VISUALIZING THE DATA IN RDS FILES

You can output bed-files of the raw reads using 
makebedfromrds.py. A more practical way to look at the data 
might be to ouput it as a WIG file using makewiggle.py . 

Note that UCSC has a hard limit on the size of their WIG files 
and you will likely need to break the WIGs on a per-chromosome 
basis for mammalian genomes.

RELEASE HISTORY

version 3.01   February 2009 - bug fixes
version 3.0    January  2009 - added logging to buildrdsfrom*
version 3.0rc1 December 2008 - added blat support

=cut


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

	#### INITIALISE THE OBJECT'S ELEMENTS
	$self->initialise($arguments);

	$self->{_sqlite} = "sqlite3" if not defined $self->{_sqlite};

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

    #### VALIDATE USER-PROVIDED ARGUMENTS
	($arguments) = $self->validate_arguments($arguments, $DATAHASH);	

	#### LOAD THE USER-PROVIDED ARGUMENTS
	foreach my $key ( keys %$arguments )
	{		
		#### LOAD THE KEY-VALUE PAIR
		$self->value($key, $arguments->{$key});
	}
}


=head2

	SUBROUTINE		value

	PURPOSE

		SET A PARAMETER OF THE self OBJECT TO A GIVEN value

    INPUT

        1. parameter TO BE SET

		2. value TO BE SET TO

    OUTPUT

        1. THE SET parameter INSIDE THE BioSVG OBJECT

=cut
sub value
{
    my $self		=	shift;
	my $parameter	=	shift;
	my $value		=	shift;

	$parameter = lc($parameter);

    if ( defined $value)
	{	
		$self->{"_$parameter"} = $value;
	}
}

=head2

	SUBROUTINE		validate_arguments

	PURPOSE

		VALIDATE USER-INPUT ARGUMENTS BASED ON

		THE HARD-CODED LIST OF VALID ARGUMENTS

		IN THE data ARRAY
=cut

sub validate_arguments
{
	my $self		=	shift;
	my $arguments	=	shift;
	my $DATAHASH	=	shift;

	my $hash;
	foreach my $argument ( keys %$arguments )
	{
		if ( $self->is_valid($argument, $DATAHASH) )
		{
			$hash->{$argument} = $arguments->{$argument};
		}
		else
		{
			warn "'$argument' is not a known parameter\n";
		}
	}

	return $hash;
}


=head2

	SUBROUTINE		is_valid

	PURPOSE

		VERIFY THAT AN ARGUMENT IS AMONGST THE LIST OF

		ELEMENTS IN THE GLOBAL '$DATAHASH' HASH REF

=cut

sub is_valid
{
	my $self		=	shift;
	my $argument	=	shift;
	my $DATAHASH	=	shift;

	#### REMOVE LEADING UNDERLINE, IF PRESENT
	$argument =~ s/^_//;

	#### CHECK IF ARGUMENT FOUND IN '$DATAHASH'
	if ( exists $DATAHASH->{lc($argument)} )
	{
		return 1;
	}

	return 0;
}





=head2

	SUBROUTINE		AUTOLOAD

	PURPOSE

		AUTOMATICALLY DO 'set_' OR 'get_' FUNCTIONS IF THE

		SUBROUTINES ARE NOT DEFINED.

=cut

sub AUTOLOAD {
    my ($self, $newvalue) = @_;


    my ($operation, $attribute) = ($AUTOLOAD =~ /(get|set)(_\w+)$/);

    # Is this a legal method name?
    unless($operation && $attribute) {
        croak "Method name $AUTOLOAD is not in the recognized form (get|set)_attribute\n";
    }
    unless( exists $self->{$attribute} )
	{
        #croak "No such attribute '$attribute' exists in the class ", ref($self);
		return undef;
    }

    # Turn off strict references to enable "magic" AUTOLOAD speedup
    no strict 'refs';

    # AUTOLOAD accessors
    if($operation eq 'get') {

		return $self->{$attribute};

        # define subroutine
        *{$AUTOLOAD} = sub { shift->{$attribute} };

    # AUTOLOAD mutators
    }elsif($operation eq 'set') {

        # define subroutine
		print "ERANGE::AUTOLOAD    Setting attribute: $attribute\n";
        *{$AUTOLOAD} = sub { shift->{$attribute} = shift; };

        # set the new attribute value
        $self->{$attribute} = $newvalue;
    }

    # Turn strict references back on
    use strict 'refs';

    # return the attribute value
    return $self->{$attribute};
}


# When an object is no longer being used, this will be automatically called
# and will adjust the count of existing objects
sub DESTROY {
    my($self) = @_;

	#if ( defined $self->{_databasehandle} )
	#{
	#	my $dbh =  $self->{_databasehandle};
	#	$dbh->disconnect();
	#}

#    my($self) = @_;
}




1;

