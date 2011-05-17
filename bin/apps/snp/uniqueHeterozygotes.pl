#!/usr/bin/perl -w

#### DEBUG
$DEBUG = 1;

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     uniqueHeterozygotes

=head2  PURPOSE

    PRINT THE HETEROZYGOTE SNPS FOR EACH INDIVIDUAL 

    THAT WHILE ALL OTHER MEMBERS OF THE GROUP HAVE 

    HOMOZYGOTE OR UNKNOWN GENOTYPES AT THAT POSITION

=head2  INPUTS

    1. TSV INPUTFILE CONTAINING HAPMAP ENTRIES

    2. OUTPUTFILE

    3. LIST OF IDS


=head2  OUTPUTS

    1. TSV OUTPUTFILE CONTAINING ONLY THOSE SNPS WHERE THE INDIVIDUAL

        HAS THE ONLY HETEROZYGOTE GENOTYPE AMONGST ALL THE GROUP

        MEMBERS.

        THE OUTPUT FILE HAS THESE FIELDS:

        snp_id  chromosome  position    strand  person_id genotype
        rs2075511	chr16	15725642	+	NA12156	AC
        rs16967494	chr16	15728364	+	NA12156	CT
        rs1050113	chr16	15746535	+	NA12156	AG
        rs2272554	chr16	15757705	+	NA12156	AG
        rs4781689	chr16	15772973	+	NA12878	AG
        rs1050111	chr16	15824698	+	NA12878	AG

=head2    USAGE

    ./uniqueHeterozygotes.pl <--inputfile String> <--outputfile String> <--ids String>
		[--help]

        --inputfile     File containing Hapmap entries downloaded from Hapmap site
                        http://hapmap.ncbi.nlm.nih.gov/biomart/martview/bc7cd0b892fcaab14ad78d7568588a40
        --outputfile    File containing Hapmap entries downloaded from 
        --ids           Comma-separated list of ids of all samples in the group


=head2	NOTES

        CAUTION: THE HAPMAP SAMPLE GROUPS

        THE INDIVIDUAL SNPS SEEM TO BE ALWAYS PREDICTED USING FIXED GROUPS

        OF INDIVIDUALS. IF AN INDIVIDUAL DOES NOT BELONG TO THE SAME GROUP

        AS THE OTHERS FOR SNP SURVEYING PURPOSES, THIS WILL RESULT IN ALL

        OF ITS HETEROZYGOUS SNPS PASSING THE UNIQUE HETEROZYGOTE TEST BECAUSE

        THEY ARE UNIQUE WITH RESPECT TO THE OTHER INDIVIDUALS OF THE GROUP

        (BECAUSE THE OTHER MEMBERS OF THE GROUP ARE ABSENT FROM ANY SNPS CALLED

        FOR THE INDIVIDUAL)

=head2  EXAMPLES

RUN IT WITH YOUR OWN DATA:

/nethome/bioinfo/apps/agua/0.5/bin/apps/uniqueHeterozygotes.pl \
--inputfile /nethome/uozomaro/hapmap/unique-heterozygotes-INPUT.txt \
--outputfile /nethome/uozomaro/hapmap/unique-heterozygotes-OUTPUT.txt \
--ids NA12156,NA12878,NA12878,NA18507,NA18517,NA18555,NA18956,NA19129,NA19240


OR USE THE TEST DATA:

/nethome/bioinfo/apps/agua/0.5/bin/apps/uniqueHeterozygotes.pl \
--inputfile /nethome/bioinfo/apps/agua/0.5/bin/apps/t/Hapmap/hapmap-ceph-uzoezi-all.txt \
--outputfile /mihg/users/uozomaro/hapmap/unique-heterozygotes-OUTPUT.txt \
--ids NA12156,NA12878,NA12878,NA18507,NA18517,NA18555,NA18956,NA19129,NA19240

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
use Hapmap;

#### GET OPTIONS
my $outputfile;
my $inputfile;
my $ids;
my $help;
print "uniqueHeterozygotes.pl    Use option --help for usage instructions.\n" and exit if not GetOptions (	
	#### JBROWSE
    'inputfile=s'	=> \$inputfile,
    'outputfile=s'	=> \$outputfile,
    'ids=s' 	    => \$ids,
    'help' 			=> \$help
);

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "outputfile not defined (Use --help for usage)\n" if not defined $outputfile;
die "inputfile not defined (Use --help for usage)\n" if not defined $inputfile;
die "Can't find inputfile: $inputfile\n" if not -f $inputfile;
die "ids not defined (Use --help for usage)\n" if not defined $ids;

#### DEBUG

#### INSTANTIATE Hapmap OBJECT
my $hapmap = Hapmap->new(
	{
	}
);

print "uniqueHeterozygotes.pl    Doing uniqueHeterozygotes()\n";
$hapmap->uniqueHeterozygotes($inputfile, $outputfile, $ids);
print "uniqueHeterozygotes.pl    Completed uniqueHeterozygotes()\n";


#### PRINT RUN TIME
print "uniqueHeterozygotes.pl    ****************************************\n\n\n";
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

