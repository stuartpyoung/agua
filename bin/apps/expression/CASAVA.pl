#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();

=head2

    APPLICATION     CASAVA

    VERSION         0.01

    PURPOSE

        WRAPPER SCRIPT FOR RUNNING CASAVA:

            1. RUN ONE LANE AT A TIME (MATE PAIRED)

            2. OUTPUT TO A USER-SPECIFIED DIRECTORY

    INPUT

        1. LIST OF 'INPUTFILES' (s_*_1_export.txt)

        2. LIST OF 'MATEFILES' (s_*_2_export.txt)

        3. OUTPUT DIRECTORY CONTAINING THE SEQUENCES

    OUTPUT

        1. SNP CALLS WRITTEN TO cN.snp.txt FILE FOR

            EACH CHROMOSOME, EACH WITHIN THE SPECIFIC

            CHROMOSOME DIRECTORY:

                outputdir/build/Parsed*/cN

            WHERE N = 1,2,3,..,X,Y

    NOTES

        1. IN THE CASE OF MULTIPLE INPUT FILES, THE READS IN ALL

            FILES MUST BE THE SAME LENGTH BECAUSE A SINGLE s_*_pair.xml

            FILE WILL BE USED FOR THE READS MERGED INTO ONE FILE

    USAGE

./CASAVA.pl <--inputfiles String> <--matesfiles String> <--sequencedir String>
            <--rundir String> <--outputdir String> <--referencedir String>
            <--lane Integer> [--clean]  [--queue String]
             [--jobs]  [--cpus]  [--help]

    --inputfiles    :   Comma-separated *export.txt file names (e.g., 's_1_1_export.txt,s_2_1_export.txt')
    --matefiles	    :   Comma-separated mate *export.txt file names
    --sequencedir   :	Location of directory containing *export.txt files
    --rundir	    :	Location of base run directory
    --outputdir	    :   Create this directory and write output files to it
    --referencedir  :   Location of squashed genome reference files
    --lane	    :   Lane number N to extract s_N_pairs.xml file from run directory
    --clean	    :   Overwrite existing intermediate (merged) input file
    --queue	    :   Cluster queue options
    --jobs	    :   Max. number of concurrent cluster jobs (DEFAULT: 30)
    --cpus	    :   Max. number of cpus per job (DEFAULT: 8)
    --help          :   print help info

    EXAMPLES

perl /nethome/bioinfo/apps/agua/0.4/bin/apps/CASAVA.pl \
--inputfiles s_1_1_export.txt \
--matefiles s_1_2_export.txt \
--sequencedir /p/NGS/syoung/base/pipeline/run12/1.4.0/090508_HWI-EAS185_0004_3031UAAXX_Duan4_Jia3/Data/C1-104_Firecrest1.4.0_12-07-2009_syoung/Bustard1.4.0_12-07-2009_syoung/GERALD_12-07-2009_syoung \
--rundir /p/NGS/syoung/base/pipeline/run12/1.4.0/090508_HWI-EAS185_0004_3031UAAXX_Duan4_Jia3 \
--referencedir /nethome/syoung/base/pipeline/human-genome/illumina \
--outputdir /nethome/syoung/base/pipeline/casava/sample2 \
--lane 1

=cut

use strict;

#### EXTERNAL MODULES
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use Getopt::Long;
use FindBin qw($Bin);
use File::Path;
use File::Copy;

#### USE LIBRARY
use lib "$Bin/../../../lib";	
use lib "$Bin/../../../lib/external";	

#### INTERNAL MODULES
use CASAVA;
use FileTools;
use Timer;
use Util;
use Conf::Agua;

##### STORE ARGUMENTS TO PRINT TO FILE LATER
my @arguments = @ARGV;

#### FLUSH BUFFER
$| =1;

#### GET CONF 
my $conf = Conf::Agua->new(inputfile=>"$Bin/../../../conf/default.conf");
my $casava = $conf->getKeyValue("applications", 'CASAVA');
my $qstat = $conf->getKeyValue("cluster", 'QSTAT');
my $qsub = $conf->getKeyValue("cluster", 'QSUB'); #### /usr/local/bin/qsub
print "CASAVA.pl    casava: $casava\n";

#### GET OPTIONS
my $inputfiles;
my $matefiles;
my $sequencedir;
my $rundir;
my $outputdir;
my $referencedir;
my $lane;
my $clean;
my $parallel;

#### CLUSTER OPTIONS
my $jobs = 30;
my $cpus = 8;
my $sleep = 5;
my $queue = "-q gsmall";
my $dot = 1;

my $help;
if ( not GetOptions (
    'inputfiles=s'  => \$inputfiles,
    'matefiles=s'   => \$matefiles,
    'sequencedir=s' => \$sequencedir,
    'rundir=s'      => \$rundir,
    'outputdir=s'   => \$outputdir,
    'referencedir=s' => \$referencedir,
    'lane=i'        => \$lane,
    'clean'        => \$clean,
    'parallel'        => \$parallel,
    'queue=s'       => \$queue,
    'jobs=s'        => \$jobs,
    'cpus=s'        => \$cpus,
    'help'          => \$help
) )
{ print "Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "Inputfiles not defined (option --help for usage)\n" if not defined $inputfiles;
die "Matefiles not defined (option --help for usage)\n" if not defined $matefiles;
die "Output directory not defined (Use --help for usage)\n" if not defined $outputdir;
die "Reference directory not defined (Use --help for usage)\n" if not defined $referencedir;

#### DEBUG
print "inputfiles: $inputfiles\n";
print "outputdir: $outputdir\n";
print "referencedir: $referencedir\n";

#### CONVERT INPUTFILES AND MATEFILES INTO ARRAYS
my (@ins, @mates);
@ins = split ",", $inputfiles;
@mates = split ",", $matefiles if defined $matefiles;
@mates = () if not defined $matefiles;
print "CASAVA.pl    no. files: ", $#ins + 1, "\n";
print "CASAVA.pl    ins: @ins\n";
print "CASAVA.pl    mates: @mates\n";

#### CHECK ALL INPUT FILES ARE *export.txt FILES
foreach my $infile ( @ins )
{   
    print "Input file is not a *_export.txt file: $infile\n" and exit if $infile !~ /export\.txt$/;
}
foreach my $matefile ( @mates )
{   
    print "Mate file is not a *_export.txt file: $matefile\n" and exit if $matefile !~ /export\.txt$/;
}

#### CHECK LENGTH OF ARRAYS MATCHES AND ORDER OF FILES 
print "CASAVA.pl    Checking paired files match up.\n";
my $filetool = FileTools->new();
my $mates_paired = $filetool->matesPaired(\@ins, \@mates);
print "CASAVA.pl    Paired files match\n" if $mates_paired;
print "CASAVA.pl    Paired files do not match. Exiting\n" if not $mates_paired;

#### SET MERGED FILE NAMES
my $suffix = "_export.txt";
my $inputfile = "$outputdir/s_$lane" . "_1" . $suffix;
my $matefile = "$outputdir/s_$lane" . "_2" . $suffix;

#### MERGEFILES
my @infilepaths;
my @matefilepaths;
for my $i (0..$#ins)  {    $infilepaths[$i] = "$sequencedir/$ins[$i]"; }
for my $i (0..$#mates)  {  $matefilepaths[$i] = "$sequencedir/$mates[$i]"; }
$filetool->mergeFiles($inputfile, \@infilepaths) unless (-f $inputfile and not $clean); 
$filetool->mergeFiles($matefile, \@matefilepaths) unless ((-f $matefile and not $clean) or not defined $matefiles); 

#### MAKE OUTPUT DIR IF NOT EXISTS
print "CASAVA.pl    Output directory is a file: $outputdir\n" and exit if -f $outputdir;
File::Path::mkpath($outputdir) if not -d $outputdir;
print "CASAVA.pl    Can't create output directory: $outputdir\n" if not -d $outputdir;

#### SET SOURCE LANE FOR s_N_pair.xml FILE AS FIRST INPUT FILE'S LANE
my $firstfile = $ins[0];
print "CASAVA.pl    firstfile: $firstfile\n";
my ($firstlane) = $firstfile =~ /^s_(\d+)/;
print "CASAVA.pl    firstlane: $firstlane\n";

#### INSTANTIATE CASAVA
my $casavaObject = CASAVA->new(
	{
		#### INPUTS (FROM USER)
		inputfile   => $inputfile,
		matefile    => $matefile,
		sequencedir => $sequencedir,
		rundir      => $rundir,
		outputdir   => $outputdir,
		referencedir => $referencedir,
		lane        => $lane,
		firstlane   => $firstlane,

		#### EXECUTABLES (FROM CONF)
		casava => $casava,

		#### CLUSTER (PRESETS AND USER)
		qstat => $qstat,
		jobs => $jobs,
		cpus => $cpus,
		sleep => $sleep,
		dot => $dot,
		qsub => $qsub,
		queue => $queue
	}
);


#### CREATE RDS FILE CONTAINING ALL INPUT SEQUENCES
print "CASAVA.pl    Doing copyFiles()\n";
$casavaObject->copyFiles();

#### RUN RNA PAIRED TRANSCRIPTOME ANALYSIS
print "CASAVA.pl    Doing run()\n";
$casavaObject->run();

#### COLLECT SNP FILES
$casavaObject->collectSNPs();

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "CASAVA.pl    Run time: $runtime\n";
print "CASAVA.pl    Completed $0\n";
print "CASAVA.pl    ";
print Timer::datetime(), "\n";
print "CASAVA.pl    ****************************************\n\n\n";
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


