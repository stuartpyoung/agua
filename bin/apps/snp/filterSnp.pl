#!/usr/bin/perl -w

#### DEBUG
$DEBUG = 1;

#### TIME
my $time = time();

=head2

    APPLICATION     filterSNP

    PURPOSE

        1. FILTER A LIST OF SNPS BASED ON USER-DEFINED CRITERIA

        2. OUTPUT THE NUMBER OF SNPS PASSING EACH SUCCESSIVE FILTER

        3. GENERATE A SET OF DATABASE TABLES, ONE FOR EACH FILTER LEVEL

        4. FILTERS:

			1. species
			2. chromosome
			3. variant
			4. depth
			5. sense
			6. exonic
			7. dbsnp

    INPUT

        1. INPUT FILE

        2. FASTA OR FASTQ INPUT FILES

    OUTPUT

        1. FILTERED OUTPUT FILES:

			JSON 

    USAGE

./filterSNP.pl <--inputfile String> <--filetype String> <--depth String> <--species String> <--lines Integer> [--sense] [--dbsnp] [--help]    

    --inputfile				:   Full path to input SNP file
	--filetype				: 	Type of SNP file '454', 'casava', 'maq' (DEFAULT: 454)
    --species				:   'human' or 'mouse'
    --chromosome			:   Results for one chromosome, e.g., 'chr1'. Whole genome if omitted.
    --variant				:   Variant frequency as a fraction (e.g., use '0.47' for 47%)
    --depth					:   Minimum read coverage depth (DEFAULT: 10)
    --exonic				:   Exonic (exons only) or intronic (introns only).
    --dbsnp					:   'only' or 'non'. Return only hits in dbSNP or not in dbSNP.
	--sense					:   'sense' or 'missense': Return only synonymous or missense SNPs
    --help                 	:   print help info





COMPLETE JSON ARGUMENTS:

{"username":"admin","sessionId":"9999999999.9999.999","project":"Project1","workflow":"Workflow1","report":"Report1","mode":"filterReport","class":"Report::SNP","editor":false,"reportCombo":"Report1","workflowCombo":"Workflow1","projectCombo":"Project1","chromosomeCombo":"chr1","speciesCombo":"Human","senseCombo":"Missense","exonicCombo":"Exonic","dbsnpCombo":"Only dbsnp","chromosomeCheckbox":true,"senseCheckbox":false,"exonicCheckbox":false,"dbsnpCheckbox":false,"depthCheckbox":true,"variantCheckbox":true,"fileInput":"Project1/Workflow1/454HCDiffs-headers-SNP.txt","variantInput":"47%","depthSpinner":"10.00","totalResult":"-","chromosomeResult":"-","variantResult":"-","depthResult":"-","senseResult":"-","exonicResult":"-","dbsnpResult":"-"}

REDUCED TO COMMAND LINE ARGUMENTS:

--inputfile Project1/Workflow1/454HCDiffs-headers-SNP.txt \
--chromosomeCombo chr1 \
--speciesCombo Human \
--senseCombo Missense \
--exonicCombo Exonic \
--dbsnpCombo Only dbsnp \
--chromosomeCheckbox true \
--senseCheckbox 	false \
--exonicCheckbox 	false \
--dbsnpCheckbox 	false \
--depthCheckbox 	true \
--variantCheckbox 	true \
--variantInput 47% \
--depthSpinner 10.00 \
--editor \:false

WHICH WILL BE PARSED FROM THE FOLLOWING:



    EXAMPLES

perl /nethome/bioinfo/apps/agua/0.4/bin/apps/filterSNP.pl \
--inputfile /nethome/syoung/.agua/Project1/Workflow1/454HCDiffs-headers-SNP.txt \
--species human \
--chromosome chr1 \
--sense missense \
--exonic exonic \
--variant 0.47 \
--depth 10


	NOTES

		PRINT JUST THE TABLE HEADERS AND THEN AN ARRAY OF
		LINES TO REDUCE TRANSPORT:
		   outputResult : {
			   headers: [ 'column1', 'column2', ... ],
			   data: [
						   [ 1,2,3,... ],
						   [ ... ],
						   ...
			   ]
		   }

		AT THE CLIENT END, PARSE THIS INTO THE FOLLOWING
		FORMAT:
		PRINT OUTPUT DATA TO JSON FILE FOR THIS REPORT
		sense output_rows INTO AN ARRAY CONTAINING ARRAYS OF HASHES

		   {
			   identifier :"id",
			   label : "id",
			   items : [ array of hashes ]
		   }


		454HCDiffs.txt FORMAT:

			name    chromosome      ccdsstart       ccdsstop        referencenucleotide     variantnucleotide       depth   variantfrequency        chromosomestart chromosomestop  sense   referencecodon    variantcodon    referenceaa     variantaa       strand  snp     score   strand
			CCDS3.1 chr1    770     770     C       T       3       100%    881093  881093  missense        GCC     GTC     Alanine Valine  -                       -
			CCDS3.1 chr1    780     780     G       A       3       100%    879243  879243  synonymous      CTG     CTA     Leucine Leucine -                       -
			CCDS3.1 chr1    790     790     C       G       3       100%    879233  879233  missense        CTG     GTG     Leucine Valine  -                       -

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
use Report::SNP;
use Timer;
use Util;
use Conf::Agua;

##### STORE ARGUMENTS TO PRINT TO FILE LATER
my @arguments = @ARGV;

#### FLUSH BUFFER
$| =1;

#### SET filterSNP LOCATION
my $conf = Conf::Agua->new(inputfile=>"$Bin/../../../conf/default.conf");



############# 	CHECK THIS!!!!

my $filterSNP = $conf->getKeyValue("agua", 'filterSNP');


my $qstat = $conf->getKeyValue("cluster", 'QSTAT');
my $qsub = $conf->getKeyValue("cluster", 'QSUB'); #### /usr/local/bin/qsub
my ($filterSNPbin) = $filterSNP =~ /^(.+)\/[^\/]+$/;

#### DEFAULT MAX FILE LINES (MULTIPLE OF 4 BECAUSE EACH FASTQ RECORD IS 4 LINES)
my $maxlines = 4000000;

##### GET OPTIONS
my $inputfile;
my $species;
my $chromosome;
my $variant;
my $depth;
my $sense;
my $exonic;
my $dbsnp;

my $help;
if ( not GetOptions (
    'inputfile=s' => \$inputfile,
    'species=s' => \$species,
    'chromosome=s' => \$chromosome,
    'depth=s' => \$depth,
    'sense=s' => \$sense,
    'exonic=s' => \$exonic,
    'dbsnp=s' => \$dbsnp,
    'help' => \$help
) )
{ print "filterSNP.pl    Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
print "input file 1 not defined (Use --help for usage)\n" if not defined $inputfile;
print "depth directory not defined (Use --help for usage)\n" if not defined $depth;
print "species not defined (Use --help for usage)\n" if not defined $species;

print "filterSNP.pl    inputfile: $inputfile\n";
print "filterSNP.pl    depth: $depth\n";

#### CHECK INPUTS
print "filterSNP.pl    input file not defined (option --inputfile)\n" if not defined $inputfile;
print "filterSNP.pl    Output directory not defined (option d)\n" if not defined $depth;
print "filterSNP.pl    Reference file not defined (option r)\n" if not defined $species;

#{"username":"admin","sessionId":"9999999999.9999.999","project":"Project1","workflow":"Workflow1","report":"Report1","mode":"filterReport","class":"Report::SNP","editor":false,"reportCombo":"Report1","workflowCombo":"Workflow1","projectCombo":"Project1","chromosomeCombo":"chr1","speciesCombo":"Human","senseCombo":"Missense","exonicCombo":"Exonic","dbsnpCombo":"Only dbSNP","chromosomeCheckbox":true,"senseCheckbox":false,"exonicCheckbox":false,"dbsnpCheckbox":false,"depthCheckbox":true,"variantCheckbox":true,"fileInput":"Project1/Workflow1/454HCDiffs-headers-SNP.txt","variantInput":"47%","depthSpinner":"10.00","totalResult":"-","chromosomeResult":"-","variantResult":"-","depthResult":"-","senseResult":"-","exonicResult":"-","dbsnpResult":"-"}

my $mappings = {
	inputfile => "fileInput",
	chromosome => "chromosomeCombo",
	species => "speciesCombo",
	sense => "senseCombo",
	exonic => "exonicCombo",
	dbsnp => "dbSNPCombo Only",
	chromosome => "chromosomeCheckbox",
	variant => "variantInput 47%",
	depth => "depthSpinner 10.00",
	editor => "editor"
};

#### CONVERT ARGUMENTS TO HASH
my $arguments = {};
mapInput($arguments, $mappings, $sense, "sense");
mapInput($arguments, $mappings, $exonic, "exonic");
mapInput($arguments, $mappings, $dbsnp, "dbsnp");
mapInput($arguments, $mappings, $depth, "depth");
mapInput($arguments, $mappings, $variant, "variant");


$arguments->{exonicCheckbox} = 1 if defined $arguments->{exonic};
$arguments->{dbsnpCheckbox} = 1 if defined $arguments->{dbsnp};
$arguments->{depthCheckbox} = 1 if defined $arguments->{depth};
$arguments->{variantCheckbox} = 1 if defined $arguments->{variant};


my $report = Report::SNP->new(
	{
		inputfile 	=> $inputfile,
		depth 		=> $depth,
		dbsnp 		=> $dbsnp,
		splitfile 	=> $splitfile,
		qstat 		=> $qstat,
		jobs 		=> $jobs,
		sleep 		=> $sleep,
		sense 		=> $sense,
		dot 		=> $dot,
		qsub 		=> $qsub,
		queue 		=> $queue,
		filterSNP 	=> $filterSNP,
		filterSNPbin => $filterSNPbin
	}
);

print $runfilterSNP->snpcall();


#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "filterSNP.pl    Run time: $runtime\n";
print "filterSNP.pl    Completed $0\n";
print "filterSNP.pl    ";
print Timer::datetime(), "\n";
print "filterSNP.pl    ****************************************\n\n\n";
exit;

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#									SUBROUTINES
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

=head2

	SUBROUTINE 		mapInput

	PURPOSE

		MAP THE COMMAND LINE INPUTS TO THE filterReport JSON INPUTS

=cut

sub mapInput
{
	my $arguments		=	shift;
	my $mappings		=	shift;
	my $value			=	shift;
	my $name			=	shift;

	my $checkbox = $name . "Checkbox";
	$arguments->{$mappings->{$name}} = $value if defined $value;
	$arguments->{$checkbox} = 1 if defined $value;
	$arguments->{$checkbox} = 0 if not defined $value;	
}


sub usage
{
	print GREEN;
	print `perldoc $0`;
	print RESET;
	exit;
}


