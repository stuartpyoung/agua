#!/usr/bin/perl -w

#### DEBUG
$DEBUG = 1;

#### TIME
my $time = time();

=head2

    APPLICATION     ERANGE

    PURPOSE

        WRAPPER SCRIPT FOR RUNNING ERANGE TRANSCRIPTOME PREDICTION

    INPUT

        1. LABEL OF EXPERIMENT

        2. NAME OF SPECIES

        3. ELAND 'multi' OR 'extended' OUTPUT FILES

		4. OUTPUT DIRECTORY (WILL BE CREATED)

		5. 'knowngene.txt' FILE (LIST OF GENES DOWNLOADED FROM UCSC)        

    OUTPUT

        1. ERANGE OUTPUT FILES IN ASSEMBLY DIRECTORY

    USAGE

    ./ERANGE.pl <--label String> <--species String> <--inputfiles String> [--matefiles String] <--outputdir String> <--knowngene String> [--clean] [--help]

    --label					:   Name of experiment, e.g., 'sample1'
    --species				:   Name of species, e.g., 'human' or 'mouse'
    --inputfiles			:   Location of 'multi.txt' or 'extended.txt' ELAND output file
    --matefiles				:   Location of mate pair 'multi' or 'extended' ELAND output file
ead_2.fastq')
    --outputdir				:   Full path to output directory (to be created if not present)
    --knowngene				:   Location of UCSC knowngene.txt file
    --outputdir				:   Full path to output directory (to be created if not present)
    --clean					:   Set this flag to delete label.rds file if already exists

    --queue					:   Cluster queue options
    --jobs					:   Max. number of concurrent cluster jobs
    --cpus					:   Max. number of cpus per job
    --help                 	:   print help info

    EXAMPLES


cd /nethome/syoung/base/pipeline/erange2/data/rds/paired

perl /nethome/bioinfo/apps/agua/0.4/bin/apps/ERANGE.pl \
--label sample2 \
--species human \
--outputdir /nethome/syoung/base/pipeline/erange2/data/rds/paired \
--inputfiles /nethome/syoung/base/pipeline/erange2/data/reads/s_1_1_eland_extended.txt \
--matefiles /nethome/syoung/base/pipeline/erange2/data/reads/s_1_2_eland_extended.txt \
--knowngene /nethome/syoung/base/pipeline/erange2/data/knownGene/knownGene.hg18




perl /nethome/bioinfo/apps/agua/0.4/bin/apps/ERANGE.pl \
--label sample2 \
--species human \
--outputdir /nethome/syoung/base/pipeline/erange2/data/rds/paired \
--inputfiles /nethome/syoung/base/pipeline/erange2/data/reads/s_1_1_eland_extended.txt,/nethome/syoung/base/pipeline/erange2/data/reads/s_2_1_eland_extended.txt \
--matefiles /nethome/syoung/base/pipeline/erange2/data/reads/s_1_2_eland_extended.txt,/nethome/syoung/base/pipeline/erange2/data/reads/s_2_2_eland_extended.txt \
--knowngene /nethome/syoung/base/pipeline/erange2/data/knownGene/knownGene.hg18



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
use ERANGE;
use FileTools;
use Timer;
use Util;
use Conf::Agua;

#### EXTERNAL MODULES

##### STORE ARGUMENTS TO PRINT TO FILE LATER
my @arguments = @ARGV;

#### FLUSH BUFFER
$| =1;

#### GET CONF
my $conf = Conf::Agua->new(inputfile=>"$Bin/../../../conf/default.conf");
my $python = $conf->getKeyValue("applications", 'PYTHON');
my $erange = $conf->getKeyValue("agua", 'ERANGE');
my $sqlite = $conf->getKeyValue("applications", 'SQLITE');
my $repmask = $conf->getKeyValue("applications", 'REPMASK');
my $qstat = $conf->getKeyValue("cluster", 'QSTAT');
my $qsub = $conf->getKeyValue("cluster", 'QSUB'); #### /usr/local/bin/qsub
print "ERANGE.pl    erange: $erange\n";

#### GET OPTIONS
my $label;
my $species;
my $inputfiles;
my $matefiles = '';
my $knowngene;
my $outputdir;
my $clean;

#### CLUSTER OPTIONS
my $cpus = 1;
my $jobs = 30;
my $sleep = 5;
my $queue = "-q gsmall";
my $dot = 1;

my $help;
if ( not GetOptions (

    'label=s' 		=> \$label,
    'species=s' 	=> \$species,
    'inputfiles=s' 	=> \$inputfiles,
    'matefiles=s' 	=> \$matefiles,
    'knowngene=s' 	=> \$knowngene,
    'outputdir=s' 	=> \$outputdir,
    'clean' 		=> \$clean,

    'queue=s' 		=> \$queue,
    'cpus=s' 		=> \$cpus,
    'jobs=s' 		=> \$jobs,
    'help' 			=> \$help
) )
{ print "ERANGE.pl    Use option --help for usage instructions.\n";  exit;    };

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "Input files not defined (Use --help for usage)\n" if not defined $inputfiles;
die "Output directory not defined (Use --help for usage)\n" if not defined $outputdir;
die "Reference file not defined (Use --help for usage)\n" if not defined $knowngene;
print "ERANGE.pl    inputfiles: $inputfiles\n";
print "ERANGE.pl    matefiles: $matefiles\n" if defined $matefiles;
print "ERANGE.pl    outputdir: $outputdir\n";
print "ERANGE.pl    knowngene: $knowngene\n";

#### CHECK INPUTS
if ( not defined $inputfiles)   {   print "ERANGE.pl    Input file 1 not defined (option i)\n";    usage();    }
if ( not defined $outputdir)   {   print "ERANGE.pl    Output directory not defined (option d)\n";    usage();    }
if ( not defined $knowngene )   {   print "ERANGE.pl    Reference file not defined (option r)\n";    usage();    }

#### MAKE OUTPUT DIR IF NOT EXISTS
print "ERANGE.pl    Output directory is a file: $outputdir\n" and exit if -f $outputdir;
File::Path::mkpath($outputdir) if not -d $outputdir;
print "ERANGE.pl    Can't create output directory: $outputdir\n" if not -d $outputdir;

#### CONVERT INPUTFILES AND MATEFILES INTO ARRAYS
my @ins = split ",", $inputfiles;
my @mates;
@mates = split ",", $matefiles if defined $matefiles;
@mates = () if not defined $matefiles;
print "ERANGE.pl    ins: @ins\n";
print "ERANGE.pl    mates: @mates\n";

#### CHECK LENGTH OF ARRAYS MATCHES AND ORDER OF FILES 
my $filetool = FileTools->new();
my $mates_paired = $filetool->matesPaired(\@ins, \@mates);
print "ERANGE.pl    mates_paired: $mates_paired\n";

#### MAKE OUTPUT DIR
mkdir($outputdir) if not -d $outputdir;

#### CHECK INPUT FILE TYPE TO DETERMINE IF WE NEED TO RUN ELAND
#### OR BOWTIE TO ALIGN THE SEQUENCES
my $filetype = 'sequence';
$filetype = "extended" if $ins[0] =~ /extended\.txt$/;
$filetype = "multi" if $ins[0] =~ /multi\.txt$/;

#### SET MERGED FILE SUFFIX
my $suffix = '_sequence.txt';
$suffix = '_eland_multi.txt' if $filetype eq "multi";
$suffix = "_eland_extended.txt" if $filetype eq "extended";

#### MERGE INPUT FILES IF THEY ARE ELAND OUTPUT FILES
my $inputfile = "$outputdir/$label" . "_1" . $suffix;
my $matefile = "$outputdir/$label" . "_2" . $suffix;
print "ERANGE.pl    no. files: ", $#ins + 1, "\n";
$filetool->mergeFiles($inputfile, \@ins) unless (-f $inputfile and not $clean); 
$filetool->mergeFiles($matefile, \@mates) unless ((-f $matefile and not $clean) or not defined $matefiles); 

#### INSTANTIATE ERANGE
my $erangeObject = ERANGE->new(
	{
		#### INPUTS (FROM USER)
		label => $label,
		species => $species,
		inputfile => $inputfile,
		matefile => $matefile,
		knowngene => "$knowngene",
		outputdir => $outputdir,
		clean => $clean,
		filetype => $filetype,

		#### EXECUTABLES (FROM CONF)
		python => $python,
		erange => $erange,
		sqlite => $sqlite,
		repmask => $repmask,

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
print "ERANGE.pl    RUNNING makeRdsFromEland()\n";
$erangeObject->makeRdsFromEland();

#### RUN RNA PAIRED TRANSCRIPTOME ANALYSIS
$erangeObject->runRNAPairedAnalysis();


#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "ERANGE.pl    Run time: $runtime\n";
print "ERANGE.pl    Completed $0\n";
print "ERANGE.pl    ";
print Timer::datetime(), "\n";
print "ERANGE.pl    ****************************************\n\n\n";
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






__END__

Step 8: Build rds format files from eland output

It is necessary to create the output file BEFORE generating the output data. So you might use sqlite3 to create the output file and then ask erange to make the rds file.

You will need to build an rds format file for each read. So for read s_1_1_eland_multi.txt you might use the following:

mkdir rds
sqlite3 s_1_1.rds "create table dummy (integer int)"

/nethome/bioinfo/apps/python/2.6.1/bin/python /nethome/syoung/base/pipeline/erange2/commoncode/makerdsfromeland2.py s_1_1 /nethome/syoung/base/pipeline/erange2/data/eland/s_1_1_eland_multi.txt /nethome/syoung/base/pipeline/erange2/data/rds/s_1_1.rds -RNA /nethome/syoung/base/pipeline/erange2/data/knownGene/knownGene.hg18 -paired 1 -verbose


The output will be put in s_1_1.rds. The arguments to makerdsfromeland2.py are a label, the eland read file, the output file, the known gene table, and, for paired reads, -paired and # of pair.
For the last read you need to run the above code with the argument -index.



Step 9: Build a repeatmask database using buildrmaskdb.py

Make sure that your input file (downloaded from UCSC) has a name starting with rmsk, otherwise this
won't work. So you might need to rename the file you downloaded.


/nethome/apps/bioinfo/python/2.6.1/bin/python /nethome/syoung/base/pipeline/erange2/commoncode/buildrmaskdb.py /nethome/syoung/base/pipeline/erange2/data/repMask/RepMask3.2.7.hg18 /nethome/syoung/base/pipeline/erange2/data/repMask/rmask.db


Step 10: Run paired read analysis; runRNAPairedAnalysis.sh


Form of function call: runRNAPairedAnalysis.sh genome rdsprefix repeatmaskdb



# preliminary: set PYTHONPATH to point to the parent directory of the Cistematic, e.g.
export PYTHONPATH=/nethome/syoung/base/pipeline/erange2
# preliminary: set CISTEMATIC_ROOT to the directory that contains the genome directories (such as H_sapiens or M_musculus), e.g.
export CISTEMATIC_ROOT=/nethome/syoung/base/pipeline/erange2
# preliminary: set ERANGEPATH, e.g. 
export ERANGEPATH=/nethome/syoung/base/pipeline/erange2/commoncode

echo $PYTHONPATH
echo $CISTEMATIC_ROOT
echo $ERANGEPATH



1. MOVE TO rds FILES DIR AND RUN

screen -S erange

cd /nethome/syoung/base/pipeline/erange2/data/rds

#/nethome/apps/bioinfo/python/2.6.1/bin/python /nethome/syoung/base/pipeline/erange2/commoncode/runRNAPairedAnalysis.sh human s_1_1 /nethome/syoung/base/pipeline/erange2/data/repMask/rmask.db

/nethome/syoung/base/pipeline/erange2/commoncode/runRNAPairedAnalysis.sh human s_1_1_1000 /nethome/syoung/base/pipeline/erange2/data/repMask/rmask.db




2. DO IT WITH A SHELL SCRIPT AND RENAMED INPUT FILE

cd /nethome/syoung/base/pipeline/erange2/data/rds
cp s_1_1_1000.rds s_1_1_1000-shell.rds

emacs s_1_1_1000.sh

#!/bin/bash

#PBS -N test-er

#PBS -j oe

# SET ENVIRONMENT VARIABLES
# preliminary: set PYTHONPATH to point to the parent directory of the Cistematic, e.g.
export PYTHONPATH=/nethome/syoung/base/pipeline/erange2
# preliminary: set CISTEMATIC_ROOT to the directory that contains the genome directories (such as H_sapiens or M_musculus), e.g.
export CISTEMATIC_ROOT=/nethome/syoung/base/pipeline/erange2
# preliminary: set ERANGEPATH, e.g. 
export ERANGEPATH=/nethome/syoung/base/pipeline/erange2/commoncode

echo $PYTHONPATH
echo $CISTEMATIC_ROOT
echo $ERANGEPATH

cd /nethome/syoung/base/pipeline/erange2/data/rds

/nethome/syoung/base/pipeline/erange2/commoncode/runRNAPairedAnalysis.sh human s_1_1_1000-shell /nethome/syoung/base/pipeline/erange2/data/repMask/rmask.db

########  END OF SHELL SCRIPT


qsub s_1_1_1000.sh

    227486

    Mon Feb  8 01:37:31 EST 2010
    qstat
    Job id                    Name             User            Time Use S Queue
    ------------------------- ---------------- --------------- -------- - -----
    227486.kronos             test-er          syoung                 0 R default        


FINISHED IN 5 MINS

checkjob 227486

AName: test-er
State: Completed 
Complete Time:  Mon Feb  8 01:41:11
  Completion Code: 0
Creds:  user:syoung  group:bioinfo  account:bioinfo  class:default
WallTime:   00:04:09 of 4:00:00
SubmitTime: Mon Feb  8 01:37:02
  (Time Queued  Total: 00:05:02  Eligible: 00:00:00)

Total Requested Tasks: 1

Req[0]  TaskCount: 1  Partition: base  
Memory >= 0  Disk >= 0  Swap >= 0
NodeCount:  1

Allocated Nodes:
[n05:1]

IWD:            /home/syoung/base/pipeline/erange2/data/rds
Executable:     /opt/moab/spool/moab.job.6hEiLj

Execution Partition:  base
StartPriority:  0


MANUAL AND CLUSTER OUTPUT FILES ARE IDENTICAL

[syoung@u01 rds]$ diff  s_1_1_1000.uniqs.count  s_1_1_1000-shell.uniqs.count




3. RUN WHOLE LANE ON CLUSTER



cd /nethome/syoung/base/pipeline/erange2/data/rds
emacs s_1_1.sh

#!/bin/bash

#PBS -N test-er

#PBS -j oe

# SET ENVIRONMENT VARIABLES
# preliminary: set PYTHONPATH to point to the parent directory of the Cistematic, e.g.
export PYTHONPATH=/nethome/syoung/base/pipeline/erange2
# preliminary: set CISTEMATIC_ROOT to the directory that contains the genome directories (such as H_sapiens or M_musculus), e.g.
export CISTEMATIC_ROOT=/nethome/syoung/base/pipeline/erange2
# preliminary: set ERANGEPATH, e.g. 
export ERANGEPATH=/nethome/syoung/base/pipeline/erange2/commoncode

echo $PYTHONPATH
echo $CISTEMATIC_ROOT
echo $ERANGEPATH

cd /nethome/syoung/base/pipeline/erange2/data/rds

/nethome/syoung/base/pipeline/erange2/commoncode/runRNAPairedAnalysis.sh human s_1_1 /nethome/syoung/base/pipeline/erange2/data/repMask/rmask.db

########  END OF SHELL SCRIPT


chmod 755 *.sh
qsub s_1_1.sh; sleep 5; date; qstat


qsub s_1_1.sh

    227495

sleep 20; date; qstat
    Mon Feb  8 01:45:12 EST 2010

    Job id                    Name             User            Time Use S Queue
    ------------------------- ---------------- --------------- -------- - -----
    227495.kronos             test-er          syoung                 0 R default 




