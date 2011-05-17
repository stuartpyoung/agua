#!/usr/bin/perl -w

#### DEBUG
$DEBUG = 1;

#### TIME
my $time = time();

=head2

    APPLICATION     filterFile

    PURPOSE

        1. FILTER A TSV FILE BASED ON VALUES IN COLUMNS

        2. PRINT THE FILTERED LINES TO FILE

    INPUT

        1. INPUT FILE

        2. FASTA OR FASTQ INPUT FILES

    OUTPUT

        1. FILTERED OUTPUT FILES:

			JSON 

    USAGE

./filterFile.pl <--inputfile String> <--columns String> <--values String> <--types String> [--help]    

    --inputfile				:   Full path to input SNP file
	--columns				: 	Type of SNP file '454', 'casava', 'maq' (DEFAULT: 454)
    --values				:   'human' or 'mouse'
    --types					:   Results for one values, e.g., 'chr1'. Whole genome if omitted.
    --help                 	:   print help info

	EXAMPLES

perl /nethome/bioinfo/apps/agua/0.4/bin/apps/filterFile.pl \
--inputfile /scratch/syoung/base/pipeline/SRA/NA18507/maq/maq1/chr22/out.filter \
--outputfile /scratch/syoung/base/pipeline/SRA/NA18507/maq/maq1/chr22/out.filter.fil \
--columns 1,8 \
--values chr22,9 \
--types "eq,>"


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


=cut

use strict;

#### EXTERNAL MODULES
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use Getopt::Long;
use FindBin qw($Bin);

#### USE LIBRARY
use lib "$Bin/../../../lib";	
use lib "$Bin/../../../lib/external";	

#### INTERNAL MODULES
use Filter;
use Timer;
use Util;
use Conf::Agua;

##### STORE ARGUMENTS TO PRINT TO FILE LATER
my @arguments = @ARGV;

##### GET OPTIONS
my $inputfile;
my $outputfile;
my $columns;
my $values;
my $types;
my $help;
if ( not GetOptions (
    'inputfile=s' => \$inputfile,
    'outputfile=s' => \$outputfile,
    'columns=s' => \$columns,
    'values=s' => \$values,
    'types=s' => \$types,
    'help' => \$help
) )
{ print "filterFile.pl    Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
print "inputfile not defined (Use --help for usage)\n" if not defined $inputfile;
print "columns not defined (Use --help for usage)\n" if not defined $columns;
print "values not defined (Use --help for usage)\n" if not defined $values;
print "types not defined (Use --help for usage)\n" if not defined $types;

#### DEBUG
print "filterFile.pl    inputfile: $inputfile\n";
print "filterFile.pl    columns: $columns\n";
print "filterFile.pl    values: $values\n";
print "filterFile.pl    types: $types\n";

#### GET ARRAYS
my @cols = split ",", $columns;
my @vals = split ",", $values;
my @typs = split ",", $types;

#### INSTANTIATE FILTER OBJECT AND RUN FILTER
my $filter = Filter->new();
$filter->filterColumns(
	{
		inputfile 	=> $inputfile,
		outputfile 	=> $outputfile,
		columns 	=> \@cols,
		values 		=> \@vals,
		types 		=> \@typs
	}
);


#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "filterFile.pl    Run time: $runtime\n";
print "filterFile.pl    Completed $0\n";
print "filterFile.pl    ";
print Timer::datetime(), "\n";
print "filterFile.pl    ****************************************\n\n\n";
exit;

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#									SUBROUTINES
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

sub usage
{
	print GREEN;
	print `perldoc $0`;
	print RESET;
	exit;
}


