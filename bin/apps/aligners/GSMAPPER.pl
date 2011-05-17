#!/usr/bin/perl -w

#### DEBUG
$DEBUG = 1;

#### TIME
my $time = time();

=head2

    APPLICATION     GSMAPPER

    PURPOSE

            THIS APPLICATION WILL ALIGN SAMPLES FROM DIFFERENT INDIVIDUALS

            (ONE .fna FILE PER INDIVIDUAL) AGAINST A SINGLE REFERENCE FILE 

            AND PUT THE ALIGNMENT OUTPUTS FOR THE DIFFERENT INDIVIDUALS INTO 

            SEPARATE OUTPUT DIRECTORES. FOR EACH INPUT .fna FILE IT WILL:

                1. ALIGN 454 READS AGAINST A REFERENCE SEQUENCE

                2. GENERATE MULTIPLE GS MAPPER OUTPUT FILES

            ** NB **: IT DOES NOT POOL MULTIPLE .fna FILES FROM ONE INDIVIDUAL 

            AND PUT THEM INTO A SINGLE OUTPUT DIRECTORY

    INPUT

        1. FASTA .fna FILES AND QUALITY .qual FILES

        2. FASTA REFERENCE FILE

        3. OUTPUT DIRECTORY

        4. (optional) CLUSTER QUEUE OPTIONS

        5. (optional) ADDITIONAL GS MAPPER OPTIONS

    OUTPUT

        1. GS MAPPER GENERATES A 'mapping' FOLDER CONTAINING THE ALIGNMENT 

            OUTPUT FILES, E.G., 454HCDiffs.txt AND 454AllDiffs.txt

        2. ONE MAPPING FOR EACH .fna FILE IN SEPARATE OUTPUT DIRECTORIES

    USAGE

    ./GSMAPPER.pl  <--inputfiles String> <--outputdirs String> <--referencefile String> [--queue String] [--options String] [-h]

        --inputfiles        :   COMMA-SEPARATED LIST OF /full/path/to/input1.fna FILES
        --outputdirs        :   COMMA-SEPARATED LIST OF /full/path/to/output_directories
        --referencefile        :   /full/path/to/reference.fasta FILE
        --queue                :    CLUSTER QUEUE OPTIONS
        --options            :    ADDITIONAL GS MAPPER OPTIONS
                                (DEFAULT OPTIONS: -ace -pairt -nrm )
        --help                 :   print help info

    EXAMPLES

cd /nethome/syoung/base/bin/exome

perl GSMAPPER.pl \
--inputfiles /nethome/syoung/base/pipeline/nimblegen-run2/SID10032/ccds/SID10032_454Reads.fna,/nethome/syoung/base/pipeline/nimblegen-run2/SID10033/ccds/SID10033_454Reads.fna,/nethome/syoung/base/pipeline/nimblegen-run2/SID10034/ccds/SID10034_454Reads.fna,/nethome/syoung/base/pipeline/nimblegen-run2/SID10035/ccds/SID10035_454Reads.fna,/nethome/syoung/base/pipeline/nimblegen-run2/SID10036/ccds/SID10036_454Reads.fna,/nethome/syoung/base/pipeline/nimblegen-run2/SID10037/ccds/SID10037_454Reads.fna,/nethome/syoung/base/pipeline/nimblegen-run2/SID10039/ccds/SID10039_454Reads.fna,/nethome/syoung/base/pipeline/nimblegen-run2/SID10082/ccds/SID10082_454Reads.fna \
--referencefile /nethome/syoung/base/pipeline/nimblegen-gsmapper/CCDS_nucleotide.20080430.fa \
--outputdirs /nethome/syoung/base/pipeline/nimblegen-run2/SID10032/ccds,/nethome/syoung/base/pipeline/nimblegen-run2/SID10033/ccds,/nethome/syoung/base/pipeline/nimblegen-run2/SID10034/ccds,/nethome/syoung/base/pipeline/nimblegen-run2/SID10035/ccds,/nethome/syoung/base/pipeline/nimblegen-run2/SID10036/ccds,/nethome/syoung/base/pipeline/nimblegen-run2/SID10037/ccds,/nethome/syoung/base/pipeline/nimblegen-run2/SID10039/ccds,/nethome/syoung/base/pipeline/nimblegen-run2/SID10082/ccds 


=cut

use strict;

#### USE LIBRARY
use FindBin qw($Bin);
use lib "$Bin/../../../lib";

#### INTERNAL MODULES
use Timer;
use Util;
use Conf::Agua;

#### EXTERNAL MODULES
use Term::ANSIColor qw(:constants);
use Getopt::Long;
use Data::Dumper;

#### GET CONFIGURATION PARAMETERS
my $application = $0;
my ($path) = $application =~ /^(.+?)\/([^\/]+)$/;
if ( not defined $path )
{
    $path = "";
}
else
{
    $path .= "/";
}
my $conf = Conf::Agua->new(inputfile=>"$Bin/../../../conf/default.conf", 'is_file');
my $bin = $conf->getKeyValue("agua", 'BIN');
my $runMapping = $conf->getKeyValue("agua", 'RUNMAPPING');
#exit;

#### SAVE ARGUMENTS, MINUS 'imagedir' ARGUMENT AND VALUE
my @arguments = @ARGV;

#### GET OPTIONS
my $options = '';	
my $inputfiles;	
my $outputdirs;
my $referencefile;
my $queue = '';
my $help;
GetOptions (
	'options=s' => \$options,
    'inputfiles=s' => \$inputfiles,
    'outputdirs=s' => \$outputdirs,
    'referencefile=s' => \$referencefile,
    'queue=s' => \$queue,
    'help' => \$help
)
or usage();

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
if ( not defined $inputfiles)   {   die "Input file not defined (option --inputfiles)\n";	}
if ( not defined $outputdirs)   {   die "Output dir not defined (option --outputdirs)\n";	}
if ( not defined $referencefile)   {   die "Reference file not defined (option --referencefile)\n";	}

#### CHECK REFERENCE FILE
if ( not -f $referencefile )    {   die "Can't find reference file: $referencefile";    }
if ( -z $referencefile )    {   die "Reference file is empty: $referencefile";    }

#### GET INPUT FILES
my @infiles = split ",", $inputfiles;

#### GET OUTPUT FILES
my @outdirs = split ",", $outputdirs;

if ( $#outdirs != $#infiles )
{
	print "The number of output directories (" , $#outdirs + 1, ") does not equal the number of input files (", $#infiles + 1, ")\n";
	exit;
}

#### DO EACH INPUT FILE
for ( my $i = 0; $i < $#infiles + 1; $i++ )
{
    my $inputfile = $infiles[$i];
	my $outputdir = $outdirs[$i];

	#### CREATE OUTPUT DIR
	if ( not -d $outputdir )
	{
		if ( -f $outputdir )
		{
			die "Output dir is a file: $outputdir";
		}
		`mkdir -p $outputdir`;
	}

    #### CHECK INPUT FILE
    if ( not -f $inputfile )    {    die "Could not find input file: $inputfile\n";   }

    ### 2. CREATE SHELL SCRIPT 
    my $scriptfile = $inputfile;
    $scriptfile =~ s/\.fna$//;    
    $scriptfile .= "-runMapping.sh";

    my $shellscript = qq{#!/bin/sh

HOST=`hostname`
echo \$HOST

echo "Doing runMapping.pl..."
echo "time $runMapping -ace -o $outputdir -pairt -nrm -ref $referencefile -read $inputfile $options";
cd $outputdir;
echo `pwd`;

time $runMapping -ace -o $outputdir \\
-pairt -nrm -ref $referencefile \\
-read $inputfile $options

echo `/usr/local/bin/qstat -f \$PBS_ID`
};



    open(SHFILE, ">$scriptfile") or die "Can't open script file: $scriptfile\n";
    print SHFILE $shellscript;
    close(SHFILE);
    #exit;

    ### 3. RUN QSUB COMMAND
    my $qsub_command = "qsub $queue $scriptfile";
	print "Qsub command:\n";
}


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

sub usage
{
	print `perldoc $0`;
}



