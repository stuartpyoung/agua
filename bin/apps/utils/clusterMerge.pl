#!/usr/bin/perl -w

$DEBUG = 1;

=head2

    TEST        clusterMerge

    PURPOSE

		CARRY OUT A DISTRIBUTED MERGE ON A CLUSTER OF MULTIPLE INPUT FILES

		INTO A SINGLE OUTPUT FILE:

			1. ASSUMES FILE CAN BE MERGED IN ANY ORDER

			2. AT EACH STEP, FILES ARE MERGED IN PAIRS

			3. USE DIFFERENT MERGE FUNCTIONS FOR DIFFERENT FILE TYPES:

					map (MAQ MERGE)

					sam, txt, cat (CONCAT)

					bam (SAMTOOLS MERGE)

	EXAMPLES

/nethome/bioinfo/apps/agua/0.4/bin/apps/clusterMerge.pl \
--inputdir /scratch/syoung/base/pipeline/SRA/NA18507/maq/maq1/chr22 \
--pattern "*/out.sam" \
--outputfile /scratch/syoung/base/pipeline/SRA/NA18507/maq/maq1/chr22/out.sam \
--type sam


/nethome/bioinfo/apps/agua/0.4/bin/apps/clusterMerge.pl \
--inputdir /scratch/syoung/base/pipeline/SRA/NA18507/maq/maq1/chr22 \
--pattern "*/out.map" \
--outputfile /scratch/syoung/base/pipeline/SRA/NA18507/maq/maq1/chr22/out.sam \
--type maq




=cut

use strict;

#### USE FINDBIN
use FindBin qw($Bin);

#### USE LIBS
use lib "$Bin/../../../lib";
use lib "$Bin/../../../lib/external";
use lib "$Bin/../../../lib/external/lib64/perl5/site_perl/5.8.8/x86_64-linux-thread-multi";

#### INTERNAL MODULES
use Cluster;

#### EXTERNAL MODULES
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use Getopt::Long;

##### STORE ARGUMENTS TO PRINT TO FILE LATER
my @arguments = @ARGV;

#### SET MAQ LOCATION
my $conf = Conf::Agua->new(inputfile=>"$Bin/../../../conf/default.conf");
my $maq = $conf->getKeyValue("agua", 'MAQ');

#### GET OPTIONS
my $inputdir;
my $pattern;
my $type;
my $outputfile;
my $help;
if ( not GetOptions (
    'inputdir=s'  	=> \$inputdir,
    'pattern=s'   	=> \$pattern,
    'type=s'   		=> \$type,
    'outputfile=s'  => \$outputfile,
    'help'          => \$help
) )
{ print "Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "inputdir not defined (option --help for usage)\n" if not defined $inputdir;
die "outputfile not defined (Use --help for usage)\n" if not defined $outputfile;
die "type not defined (Use --help for usage)\n" if not defined $type;

#### DEBUG
print "clusterMerge.pl    inputdir: $inputdir\n";
print "clusterMerge.pl    outputfile: $outputfile\n";
print "clusterMerge.pl    type: $type\n";

#### INITIALISE Filter::SNP OBJECT 
my $cluster = Cluster->new();


#my $infiles = $cluster->listFiles($inputdir, $pattern);

my $infiles;
@$infiles = <$inputdir/$pattern>;
print "No input files found in inputdir: $inputdir\n... and with pattern: $pattern\n\n" and exit if not defined $infiles or not @$infiles;

my $mergesub;

#### CONCAT FILES SUBROUTINE
$mergesub = sub {
	my $firstfile	=	shift;
	my $secondfile	=	shift;

	my $command;
	$command .= "cp $firstfile $firstfile.temp;" if $firstfile !~ /temp$/;
	$firstfile = "$firstfile.temp" if $firstfile !~ /temp$/;
	$command = "cat $secondfile >> $firstfile;";
	$command .= "rm -fr $secondfile;" if $secondfile =~ /temp$/;

	return $command;
} if $type eq "sam" or $type eq "cat" or $type eq "txt";


#### maq mapmerge SUBROUTINE
print "clusterMerge.pl    maq not defined and type: $type\n" if $type eq "map" and not defined $maq;
$mergesub = sub {
	my $firstfile	=	shift;
	my $secondfile	=	shift;

	my $tempfile = "$firstfile.temp";	
	return "$maq/maq mapmerge $tempfile $firstfile $secondfile;";

} if $type eq "map";

print "clusterMerge.pl    mergesub not defined. Exiting\n" if not defined $mergesub;

#### DO MERGE
$cluster->merge(
    {
        inputfiles   =>  $infiles,
        outputfile  =>  $outputfile,
		mergesub	=>	$mergesub
    }
);



################################################################################
################################################################################
########################           SUBROUTINES          ########################
################################################################################
################################################################################


