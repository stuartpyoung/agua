#!/usr/bin/perl -w

$DEBUG = 1;

=head2

    TEST        pileupToSnp

    PURPOSE

		1. CHECK IF SNP IN GENE (WITHIN RIGHT OR LEFT MARGINS)

			IGNORE 'N' POSITIONS

			CHECK IF POSITION IN ccdsGene


		2. IF IN GENE, CALCULATE POSITION OF SNP

		**        UPSTREAM (MARGIN)
			5' UTR
			NON-CODING
			EXONIC
			INTRONIC
			3' UTR
		**        DOWNSTREAM (MARGIN)


		3. IF CODING, CALCULATE EFFECT OF SNP 

			MISSENSE
			SYNONYMOUS
			STOP
			START


	NOTES

		INPUT FORMAT
		------------

		emacs /scratch/syoung/base/pipeline/SRA/NA18507/maq/maq1/chr22/out.filter.bkp

			...
			chr22   16017040        N       A       36      0       49      3       a-2nna-2nna     8_A
			chr22   16017040        *       -NN/-NN 0       0       49      3       -NN     *       2       1       0           0       0
			chr22   16038081        N       T       15      0       51      3       CT-1NT  629
			chr22   16038081        *       -N/-N   0       0       51      3       -N      *       1       2       0           0       0
			chr22   16051270        t       G       36      36      60      3       Ggg     IAI
			chr22   16054071        *       +G/+G   53      159     44      5       +G      *       1       4       0           0       0
			chr22   16061884        t       G       39      39      31      4       Gg^!G^!G        15FF
			chr22   16062193        C       T       36      36      60      3       TTT     I0C
			chr22   16062195        G       T       36      36      60      3       T$TT    I20
			chr22   16066324        c       G       36      36      51      3       G$g^8g  I6I


		OUTPUT FORMAT
		-------------

		EXTENDED PILEUP CONTAINING SNP ANNOTATION

			PILEUP COLUMNS

			1. chromosome                      ('chromosome')
			2. 1-based coordinate
			3. reference base
			4. consensus base                   ('variant')
			5. consensus quality                Phred-scaled probability that the consensus is wrong
			6. SNP quality                      Phred-scaled probability that the consensus is identical to the reference. (For SNP calling, SNP quality is of more importance.)
			7. max. mapping quality of reads covering sites
			8. number of reads covering site    ('depth')
			9. read bases
			10. phred33 base qualities

			ADDITIONAL COLUMNS

			11. variantfrequency                (calculated from readbases in column 9)
			12. ccds
			13. ccdstype                        (5-utr|non-coding|exonic|intronic|3-utr)
			14. ccdsstrand                      (+|-)
			15. ccdsnumber                      (CCDS1.1)
			16. ccdsstart
			17. referencecodon
			18. variantcodon
			19. referenceaa
			20. variantaa
			21. effect                          (missense|synonymous|stop|start)
			22. snp                             ('rs_***')
			22. snptype                         ('dbsnp')
			23. snpscore


		NB: ^* CHARACTERS IN READ BASES COLUMNS, E.G.: 

			chr22   16066324        c       G       36      36      51      3       G$g^8g  I6I

			'^'     start of a read segment which is a contiguous subsequence on the read separated by `N/S/H' CIGAR operations. The ASCII of the character following `^' minus 33 gives the mapping quality. 

			'$'     end of a read segment.

		Start and end markers of a read are largely inspired by Phil Green's CALF format. These markers make it possible to reconstruct the read sequences from pileup.

	EXAMPLES

./pileupToSnp.pl \
--inputtype pileup \
--inputfile /scratch/syoung/base/pipeline/SRA/NA18507/maq/maq1/chr22/out.filter \
--outputfile /scratch/syoung/base/pipeline/SRA/NA18507/maq/maq1/chr22/out.filter.snp


=cut

use strict;

#### USE FINDBIN
use FindBin qw($Bin);

#### USE LIBS
use lib "$Bin/../../../lib";
use lib "$Bin/../../../lib/external";
use lib "$Bin/../../../lib/external/lib64/perl5/site_perl/5.8.8/x86_64-linux-thread-multi";

#### INTERNAL MODULES
use Filter::SNP;
use DBaseFactory;
use Feature;

#### EXTERNAL MODULES
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use Getopt::Long;
use File::Path;
use File::Copy;

##### STORE ARGUMENTS TO PRINT TO FILE LATER
my @arguments = @ARGV;

#### GET OPTIONS
my $inputfile;
my $inputtype;
my $outputfile;
my $dbdir;
my $dbsnp = "snp130";
my $help;
if ( not GetOptions (
    'inputfile=s'  	=> \$inputfile,
    'inputtype=s'   => \$inputtype,
    'outputfile=s'  => \$outputfile,
    'dbdir=s'   	=> \$dbdir,
    'dbsnp=s'   	=> \$dbsnp,
    'help'          => \$help
) )
{ print "Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "inputfile not defined (option --help for usage)\n" if not defined $inputfile;
die "outputfile not defined (Use --help for usage)\n" if not defined $outputfile;
die "inputtype not defined (Use --help for usage)\n" if not defined $inputtype;

#### DEBUG
print "pileupToSnp.pl    inputfile: $inputfile\n";
print "pileupToSnp.pl    outputfile: $outputfile\n";
print "pileupToSnp.pl    inputtype: $inputtype\n";

#### INITIALISE DATABASE HANDLE TO THESE TABLES:
#### ccdsGene - exonStarts, exonStops
#### ccdsSeq - sequence
#### snp130 - SNP positions
my $dbtype = "SQLite";
my $testdir = "$Bin/t/04-Filter-SNP";
my $dbfile = "$testdir/dbfile/filtersnp.dbl";
my $database = "filtersnp";
my $dbobject = 	DBaseFactory->new( $dbtype,
	{
		'DBFILE'	=>	$dbfile,
		'DATABASE'	=>	$database
	}
) or print "Can't create database: $database: $!\n" and exit;

#### SET DBDIR CONTAINING <DBSNP>-chr*.dbl DB FILES
$dbdir = "$testdir/dbfile" if not defined $dbdir;

#### INITIALISE Filter::SNP OBJECT 
my $filterSNP = Filter::SNP->new( { 'DBOBJECT' => $dbobject } );
$filterSNP->annotate(
    {
        inputfile   =>  $inputfile,
        outputfile  =>  $outputfile,
        inputtype   =>  $inputtype,
        dbobject    =>  $dbobject,
		dbdir		=>	$dbdir,
		dbsnp		=>	$dbsnp
    }
);


################################################################################
################################################################################
########################           SUBROUTINES          ########################
################################################################################
################################################################################


