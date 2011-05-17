#!/usr/bin/perl -w

$DEBUG = 1;

=head2

    TEST        replaceString

    PURPOSE

		1. REPLACE A STRING IN A FILE AND PRINT TO A NEW FILE

	NOTES

		USE cat AND sed

	EXAMPLES

/nethome/bioinfo/apps/agua/0.4/bin/apps/replaceString.pl \
--inputfile /nethome/syoung/base/pipeline/query/snp130_template.sql \
--outputfile /nethome/syoung/base/pipeline/query/snp130_chr1.sql,\
/nethome/syoung/base/pipeline/query/snp130_chr1.sql,\
/nethome/syoung/base/pipeline/query/snp130_chr2.sql,\
/nethome/syoung/base/pipeline/query/snp130_chr3.sql,\
/nethome/syoung/base/pipeline/query/snp130_chr4.sql,\
/nethome/syoung/base/pipeline/query/snp130_chr5.sql,\
/nethome/syoung/base/pipeline/query/snp130_chr6.sql,\
/nethome/syoung/base/pipeline/query/snp130_chr7.sql,\
/nethome/syoung/base/pipeline/query/snp130_chr8.sql,\
/nethome/syoung/base/pipeline/query/snp130_chr9.sql,\
/nethome/syoung/base/pipeline/query/snp130_chr10.sql,\
/nethome/syoung/base/pipeline/query/snp130_chr11.sql,\
/nethome/syoung/base/pipeline/query/snp130_chr12.sql,\
/nethome/syoung/base/pipeline/query/snp130_chr13.sql,\
/nethome/syoung/base/pipeline/query/snp130_chr14.sql,\
/nethome/syoung/base/pipeline/query/snp130_chr15.sql,\
/nethome/syoung/base/pipeline/query/snp130_chr16.sql,\
/nethome/syoung/base/pipeline/query/snp130_chr17.sql,\
/nethome/syoung/base/pipeline/query/snp130_chr18.sql,\
/nethome/syoung/base/pipeline/query/snp130_chr19.sql,\
/nethome/syoung/base/pipeline/query/snp130_chr20.sql,\
/nethome/syoung/base/pipeline/query/snp130_chr21.sql,\
/nethome/syoung/base/pipeline/query/snp130_chr22.sql,\
/nethome/syoung/base/pipeline/query/snp130_chrX.sql,\
/nethome/syoung/base/pipeline/query/snp130_chrY.sql \
--target chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \
--query chr1



/nethome/bioinfo/apps/agua/0.4/bin/apps/utils/replaceString.pl \
--global \
--inputfile /nethome/syoung/base/pipeline/dbsnp/command-snp130-template.txt \
--outputfile /nethome/syoung/base/pipeline/dbsnp/command-snp130-chr2.txt \
--target chr2 \
--query chr1


=cut

use strict;

#### USE FINDBIN
use FindBin qw($Bin);

#### USE LIBS
use lib "$Bin/../../../lib";
use lib "$Bin/../../../lib/external";

#### INTERNAL MODULES
use FileTools;

#### EXTERNAL MODULES
use Term::ANSIColor qw(:constants);
use Getopt::Long;

#### GET OPTIONS
my $inputfile;
my $outputfile;
my $target;
my $query;
my $global;
my $help;
if ( not GetOptions (
    'inputfile=s'  	=> \$inputfile,
    'outputfile=s'  => \$outputfile,
    'query=s'   	=> \$query,
    'target=s'   	=> \$target,
    'global'   		=> \$global,
    'help'          => \$help
) )
{ print "Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "inputfile not defined (option --help for usage)\n" if not defined $inputfile;
die "outputfile not defined (Use --help for usage)\n" if not defined $outputfile;

#### DEBUG

my @targs = split ",", $target;
my @outfiles = split ",", $outputfile;
print "No. targets (", $#targs + 1 ,") does not equal number of outputfiles (", $#outfiles + 1 , ")\n" and exit if $#targs != $#outfiles;

#### INSTANTIATE FileTools OBJECT
my $filetool = FileTools->new();

for ( my $i = 0; $i < $#outfiles + 1; $i++ )
{
	$filetool->replaceString(
		{
			inputfile	=>	$inputfile,
			outputfile	=>	$outfiles[$i],
			query		=>	$query,
			target		=>	$targs[$i],
			global		=>	$global
		}
	);


	my $command = "cat $inputfile | sed 's/$query/$targs[$i]/' > $outfiles[$i]";
	$command = "cat $inputfile | sed 's/$query/$targs[$i]/g' > $outfiles[$i]" if defined $global;

	print `$command`;
}

print "Completed $0\n";
exit;
