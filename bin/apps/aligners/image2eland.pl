#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();

=head2

	APPLICATION     image2eland

    PURPOSE

        **** DUMMY FILE TO TEST WORKFLOW ****


        WRAPPER TO RUN goat_pipeline.py TO CARRY OUT THE FOLLOWING STEPS:

            1. PROCESS THE IMAGE FILES (Firecrest)

            2. DO THE BASE CALLING (Bustard)

            3. ALIGN AGAINST A REFERENCE GENOME (Eland)

    INPUT

        1. IMAGE FILE DIRECTORIES INSIDE Images DIRECTORY

                Images --------+ 
                                L001
                                L002
                                L003
                                ...
                                L008----+
                                        C1.1  
                                        C2.1  
                                        C3.1 
                                        ...
                                        C36.1 ----+
                                                    s_1_1_a.tif
                                                    s_1_1_c.tif
                                                    s_1_1_g.tif
                                                    s_1_1_t.tif
                                                    ...
                                                    s_1_300_t.tif 

        2. LOCATION OF PARENT DIRECTORY

        3. NUMBER OF CPUs TO USE FOR THE DATA PROCESSING

        4. OUTPUT DIRECTORY

    OUTPUT

        1. GERALD CONFIGURATION FILE, WITH FORMAT:

            ANALYSIS eland_pair
            ELAND_GENOME /store/home/jhoffman/myGenome
            USE_BASES nY26

        2 chmod 660 .params TO MAKE FILE WRITEABLE BY ANYONE REDOING BASE CALLS

            chmod 660 /store/data/pipeline_in/080814_HWI-EAS185_0001_SeqCapture_Barcoding_RNA_JH/Data/.params

        3. Eland OUTPUT *sorted.txt FILES IN GERALD.pl DIRECTORY

        4. COPY *sorted.txt FILES TO OUTPUT DIRECTORY

    NOTES

        IF NOT SPECIFIED BY THE USER, THE APPLICATION CREATES A GERALD

        CONFIGURATION FILE NAMED geraldfile.txt IN THE OUTPUT DIRECTORY:

        ANALYSIS        eland_extended
        ELAND_GENOME    /home/syoung/base/pipeline/run2-human-mitochondria/data
        READ_LENGTH     30
        USE_BASES       nY*

        WHERE THE USER SPECIFIED FILE CAN CONTAIN ANY OF THE FOLLOWING:

            QUALITY_FORMAT      Two options: 'numeric' or 'symbolic'
            ANALYSIS            Type of alignment to be performed (default (phageAlign),
                                none ('8:ANALYSIS none' will ignore lane 8), sequence,
                                sequence_pair, eland, eland_extended, eland_pair)
            ELAND_GENOME        Name of the directory containing a compressed (“squashed”)
                                reference for Eland 
            READ_LENGTH         Read length to use for alignment (has to be less than or
                                equal to the number of sequencing cycles).
            USE_BASES nY[YY…]   Specifies which of the bases to align to the reference sequence.
                                “n” means the first base will not be used (maybe part of the
                                primer). Number of Y’s has to be equal to the read length.
            CONTAM_FILE         Switch contaminant filtering on
            CONTAM_DIR          Specifies location of CONTAM_FILE

    USAGE

image2eland.pl --type single 
    <--rundir> 
    <--outputdir> 
    <--referencefile>
    <--readlength>
    <--tiles>
    [--geraldfile]
    [--cpus]
    [--help]

        --rundir        /FULL/PATH/TO/SOLEXA/RUN/DIRECTORY        
        --outputdir     /FULL/PATH/TO/DESTINATION/DIRECTORY (COPY SELECTED RUN FILES TO HERE)
        --referencefile /FULL/PATH/TO/REFERENCE.fasta FILE
        --readlength    USE SELECTED NUMBER OF READ BASES
        --tiles         PROCESS IMAGES FOR THESE TILES ONLY
        --geraldfile    /FULL/PATH/TO/GERALD_CONFIG.txt FILE
        --cpus          NUMBER OF CPUs TO USE
        --help          PRINT OUT THIS HELP INFORMATION

    EXAMPLES


TEST PAIRED END


image2eland.pl --type paired \
    --rundir /store/data/pipeline_in/workflow1,/store/data/pipeline_in/workflow2 \
    --outputdir /mihg/users/syoung/base/pipeline/workflow2/eland \
    --referencefile /store/home/syoung/base/pipeline/human-mtdna/human-mtDNA-AC_000021.fasta \
    --readlength 26 \
    --tiles s_6_015

    Run time: 00:22:32
    Completed /home/syoung/base/bin/nextgen/image2eland.pl
    1:00AM, 11 October 2008
    ****************************************

    [syoung@solexa01 nextgen]$ cd /mihg/users/syoung/base/pipeline/workflow2/eland
    [syoung@solexa01 eland]$ ll
    total 21330
    -rw-rw-rw-+ 1 syoung users     132 Oct 11  2008 geraldfile.txt
    drwxrwxrwx  2 syoung users      12 Oct 11 01:05 intfiles
    -rw-rw-rw-+ 1 syoung users   29891 Oct 11 01:04 make.error
    -rw-rw-rw-+ 1 syoung users     816 Oct 11  2008 makefile.error
    -rw-rw-rw-+ 1 syoung users       0 Oct 11  2008 makefile.out
    -rw-rw-rw-+ 1 syoung users  766788 Oct 11 01:04 make.out
    -rw-rw-rw-+ 1 syoung users 9662662 Oct 11 01:05 s_6_1_sequence.txt
    -rw-rw-rw-+ 1 syoung users  692691 Oct 11 01:05 s_6_1_sorted.txt
    -rw-rw-rw-+ 1 syoung users 9662662 Oct 11 01:05 s_6_2_sequence.txt
    -rw-rw-rw-+ 1 syoung users  729083 Oct 11 01:05 s_6_2_sorted.txt
    drwxrwxrwx  2 syoung users      12 Oct 11 01:05 seqfiles


TEST SINGLE READ

image2eland.pl --type single \
    --rundir /store/data/pipeline_in/080805_HWI-EAS185_0006_20GVYAAXX_Jia_cDNA2_mtDNA_Total_Small_Medium_read2_JH \
    --geraldfile=/home/syoung/base/pipeline/run2lane6-test/eland/phix-config.txt \
    --outputdir /home/syoung/base/pipeline/run2lane6-test/eland \
    --referencefile /store/home/syoung/base/pipeline/phix/phiFasta.fa \
    --readlength 26 \
    --tiles s_5_015

image2eland.pl --type single \
    --rundir /store/data/pipeline_in/workflow1 \
    --geraldfile /mihg/users/syoung/base/pipeline/workflow1/eland/geraldfile.txt \
    --outputdir /mihg/users/syoung/base/pipeline/workflow1/eland \
    --referencefile /store/home/syoung/base/pipeline/human-mtdna/human-mtDNA-AC_000021.fasta \
    --readlength 26 \
    --tiles s_6_015


    Run time: 00:12:34
    Completed /home/syoung/base/bin/nextgen/image2eland.pl
    0:22AM, 11 October 2008
    ****************************************


image2eland.pl --type single \
    --rundir /store/data/pipeline_in/workflow1 \
    --geraldfile=/home/syoung/base/pipeline/run2lane6-test/eland/geraldfile.txt \
    --outputdir /home/syoung/base/pipeline/workflow1/eland \
    --referencefile /store/home/syoung/base/pipeline/phix/phiFasta.fa \
    --readlength 26 \
    --tiles s_6_015

image2eland.pl --type single     --rundir /store/data/pipeline_in/workflow1     --outputdir /home/syoung/base/pipeline/workflow1/eland     --referencefile /store/data/pipeline_in/workflow1/human-mtdna/human-mtDNA-AC_000021.fasta     --readlength 26     --tiles s_6_015

    Run time: 00:12:23
    Completed /home/syoung/base/bin/nextgen/image2eland.pl
    5:47PM, 10 October 2008
    ****************************************


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
my $conf = Conf::Agua->new(inputfile=>"$Bin/../../../conf/default.conf", 'is_file');

my $rootdir = $conf->getKeyValue("agua", 'ROOTDIR');
my $goat_pipeline = $conf->getKeyValue("agua", 'GOAT_PIPELINE');
my $squash = $conf->getKeyValue("agua", 'SQUASH');

#### DEFAULT VARIABLES
my $DEFAULT_CPUs = 8;

#### ADD TO PATH
#`export PATH=/usr/local/Pipeline/Gerald:/usr/local/Pipeline/Goat:\$PATH`;

#### GET OPTIONS
my $rundir;
my $outputdir;	
my $referencefile;
my $type;
my $readlength;
my $tiles;
my $geraldfile;
my $help;
my $cpus;
GetOptions (
    'rundir=s' => \$rundir,
    'outputdir=s' => \$outputdir,
    'referencefile=s' => \$referencefile,
    'type=s' => \$type,
    'readlength=s' => \$readlength,
    'tiles=s' => \$tiles,
    'geraldfile=s' => \$geraldfile,
    'cpus=s' => \$cpus,
    'help' => \$help
)
or usage();
if ( defined $help )    {   usage();    }


##### TEST DELAY
#my $delay = 4;
#sleep($delay);



#### CHECK REQUIRED VARIABLES ARE DEFINED
Util::check_defined($outputdir, "\n*** Option --outputdir not defined\n", \&exit) or die;
Util::check_defined($referencefile, "\n*** Option --referencefile not defined\n") or die;
Util::check_defined($type, "\n*** Option --type not defined\n") or die;
Util::check_defined($rundir, "\n*** Option --rundir not defined\n") or die;
Util::check_defined($readlength, "\n*** Option --readlength not defined\n") or die;
Util::check_defined($tiles, "\n*** Option --tiles not defined\n") or die;

#### CHECK TYPE
print "Checking type\n";
if ( $type !~ /^(paired|single)$/ )
{
    die "Option --type ('$type') should be either 'paired' or 'single'\n";
}

#### CHECK Images DIRECTORY
my $rundir2 = '';
if ( $type =~ /^paired$/ )
{
    ($rundir, $rundir2) = $rundir =~ /^([^,]+),([^,]+)$/;
    print "Rundir 1: $rundir\n";
    print "Rundir 2: $rundir2\n";
}

print "Checking rundir\n";
if ( not Util::checkdir($rundir, "Run directory not found\n") )
{
    print "$rundir\n\n";
    die "Quitting...\n";
}

#### CHECK OUTPUT DIRECTORY
print "Checking outputdir\n";
if ( not Util::checkdir($outputdir, "Creating output directory...\n") )
{
    my $create_directory = "mkdir -p $outputdir";
    print "$create_directory\n";
    print `$create_directory`;

    if ( not Util::checkdir($outputdir, "Could not created output directory: $@\n" ) )
    {
        print "Could not cr eate output directory:\n\n$outputdir\n\n";
        print "Quitting...\n";
    }
}

#### CHECK REFERENCE FILE
print "Checking referencefile\n";
if ( not Util::checkfile($referencefile, "Reference file not found or empty"))
{
    print "\n$referencefile\n\n";
    die "Quitting...";
}

#### CHECK READLENGTH
print "Checking readlength\n";
if ( $readlength !~ /^\d+$/ )
{
    die "Option --readlength ('$readlength') is not numeric\n";
}

#### SET CPUs TO DEFAULT IF NOT DEFINED
if ( not defined $cpus )
{
    $cpus = $DEFAULT_CPUs;
}

print "Output dir: $outputdir\n";
print "Reference file: $referencefile\n";
print "Type: $type\n";
print "Images dir: $rundir\n";
print "Read length: $readlength\n";
print "Tiles: $tiles\n";
print "CPUs: $cpus\n";

#### GET ELAND GENOME DIRECTORY
my ($eland_genome) = $referencefile =~ /^(.+)\/[^\/]+$/;
print "Genome dir: $eland_genome\n";

#### IF NOT DEFINED, SET DEFAULT GERALD CONFIG FILE
if (not defined $geraldfile )
{
    $geraldfile = "$outputdir/geraldfile.txt";

    #### CREATE GERALD CONFIG FILE 
    #my $analysis = "sequence";
    my $analysis = "eland_extended";
    if ( $type =~ /^paired$/)
    {
        #$analysis = "sequence_pair";
        $analysis = "eland_pair";
    }

    my $gerald_config = qq{
ANALYSIS        $analysis
ELAND_GENOME    $eland_genome
READ_LENGTH     $readlength
USE_BASES       nY*
    };
    Util::print_file($geraldfile, $gerald_config);
}
print "Gerald configfile: $geraldfile\n";


##### SQUASH REFERENCE GENOME
#my $squashfile1 = "$referencefile.2bpb";
#my $squashfile2 = "$referencefile.vld";
#if ( not -f $squashfile1 or -z $squashfile1
#    or not -f $squashfile2 or -z $squashfile2 )
#{
#
#    my $squash_command = "$squash $eland_genome $referencefile";
#    my $squash_result = `$squash_command`;
#    if ( $squash_result =~ /Error/)
#    {
#        usage();
#    }
#}    
#
##### MOVE TO THE RUN DIRECTORY
#chdir($rundir) or die "Can't move to run directory:\n\n$rundir\n\n";
#
###### CREATE THE Makefile
#my $makefile_command = qq{time $goat_pipeline \\
#    --GERALD=$geraldfile \\
#    $rundir $rundir2 \\
#    --tiles=$tiles \\
#    --make};
#
#
##### REDIRECT OUTPUT TO FILE
#use IO::Handle;
#my $makefile_out = "$outputdir/makefile.out";
#my $makefile_error = "$outputdir/makefile.error";
#
##### SAVE OLD STDOUT AND STDERR
#my ($makefile_oldout, $makefile_olderr);
#open $makefile_oldout, ">&STDOUT" or die "Can't open old STDOUT\n";
#open $makefile_olderr, ">&STDERR" or die "Can't open old STDERR\n";
#
##### REDIRECT
#open OUTPUT, '>', "$makefile_out" or die $!;
#open ERROR,  '>', "$makefile_error"  or die $!;
#STDOUT->fdopen( \*OUTPUT, 'w' ) or die $!;
#STDERR->fdopen( \*ERROR,  'w' ) or die $!;
#
#$| = 1;
#
##### RESTORE OLD STDOUT AND STDERR, IF DEFINED
#if ( defined $makefile_oldout )
#{
#    open STDOUT, ">&", $makefile_oldout;
#}
#if ( defined $makefile_olderr )
#{
#    open STDERR, ">&", $makefile_olderr;
#}
#
#
#my $makefile_results = Util::contents($makefile_out);
#
#
#my $makefile;
#if ( $makefile_results =~ /(\S+\/Makefile)/ms )
#{
#    $makefile = $1;
#}
#else
#{
#}
#
##### EXTRACT Firecrest AND Bustard DIRECTORIES FROM MAKEFILE PATH
#my ($firecrest_directory) = $makefile =~ /^(.+Firecrest[^\/]+)/;
#my ($bustard_directory) = $makefile =~ /^(.+Bustard[^\/]+)/;
#my ($gerald_directory) = $makefile =~ /^(.+GERALD[^\/]+)/;
#
##### MOVE TO *** Firecrest *** DIRECTORY AND RUN make
#chdir($firecrest_directory) or die "Can't chdir to data directory:\n\n$firecrest_directory\n\n";
#
##### PREPARE make COMMAND
#my $make_command = "make recursive -j$cpus";
#
##### REDIRECT OUTPUT TO FILE
#my $make_out = "$outputdir/make.out";
#my $make_error = "$outputdir/make.error";
#
##### SAVE OLD STDOUT AND STDERR
#my ($make_oldout, $make_olderr);
#open $make_oldout, ">&STDOUT" or die "Can't open old STDOUT\n";
#open $make_olderr, ">&STDERR" or die "Can't open old STDERR\n";
#
##### REDIRECT
#open OUTPUT, '>', "$make_out" or die $!;
#open ERROR,  '>', "$make_error"  or die $!;
#STDOUT->fdopen( \*OUTPUT, 'w' ) or die $!;
#STDERR->fdopen( \*ERROR,  'w' ) or die $!;
#
#$| = 1;
#
##### RESTORE OLD STDOUT AND STDERR, IF DEFINED
#if ( defined $make_oldout )
#{
#    open STDOUT, ">&", $make_oldout;
#}
#if ( defined $make_olderr )
#{
#    open STDERR, ">&", $make_olderr;
#}
#
#
##### COPY FILES TO DIRECTORIES
##### 1. COPY *.int.txt.gz FILES
#my $int_directory = "$outputdir/intfiles";
#copy_files($int_directory, $firecrest_directory, "\*int\.txt\.gz");
#
##### 2. COPY .*seq.txt FILES
#my $seq_directory = "$outputdir/seqfiles";
#copy_files($seq_directory, $bustard_directory, "\*seq\.txt");
#
##### 3. COPY *sequence.txt FILES
#copy_files($outputdir, $gerald_directory, "\*sequence\.txt");
#
##### 4. COPY *sorted.txt FILES
#copy_files($outputdir, $gerald_directory, "\*sorted\.txt");

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

=head2

    SUBROUTINE      copy_files

    PURPOSE

        1. MAKE DESTINATION DIRECTORY IF NOT FOUND
        2. COPY FILES FROM SOURCE TO DESTINATION
        3. REPORT PROGRESS

=cut

sub copy_files
{
    my $destination_dir     =   shift;
    my $source_dir          =   shift;
    my $filename            =   shift;

    if ( not -d $destination_dir )
    {
        `mkdir $destination_dir`;
    }

    if ( not Util::checkdir($destination_dir) )
    {
        die "Could not create destination dir:\n\n$destination_dir\n\n";
    }

    print "Copying '$filename' files...\n";
    my $copy_command = "cp $source_dir/$filename $destination_dir";
    print "Copy command:\n$copy_command\n\n";
    `$copy_command`;
    print "Finished copying files\n\n";

    return 1;
}



sub usage
{
	print GREEN <<"EOF";

	APPLICATION     image2eland

    PURPOSE

        WRAPPER TO RUN goat_pipeline.py TO CARRY OUT THE FOLLOWING STEPS:

            1. PROCESS THE IMAGE FILES (Firecrest)

            2. DO THE BASE CALLING (Bustard)

            3. ALIGN AGAINST A REFERENCE GENOME (Eland)

    INPUT

        1. IMAGE FILE DIRECTORIES INSIDE Images DIRECTORY

                Images --------+ 
                                L001
                                L002
                                L003
                                ...
                                L008----+
                                        C1.1  
                                        C2.1  
                                        C3.1 
                                        ...
                                        C36.1 ----+
                                                    s_1_1_a.tif
                                                    s_1_1_c.tif
                                                    s_1_1_g.tif
                                                    s_1_1_t.tif
                                                    ...
                                                    s_1_300_t.tif 

        2. LOCATION OF PARENT DIRECTORY

        3. NUMBER OF CPUs TO USE FOR THE DATA PROCESSING

        4. OUTPUT DIRECTORY

    OUTPUT

        1. GERALD CONFIGURATION FILE, WITH FORMAT:

            ANALYSIS eland_pair
            ELAND_GENOME /store/home/jhoffman/myGenome
            USE_BASES nY26

        2 chmod 660 .params TO MAKE FILE WRITEABLE BY ANYONE REDOING BASE CALLS

            chmod 660 /store/data/pipeline_in/080814_HWI-EAS185_0001_SeqCapture_Barcoding_RNA_JH/Data/.params

        3. Eland OUTPUT *sorted.txt FILES IN GERALD.pl DIRECTORY

        4. COPY *sorted.txt FILES TO OUTPUT DIRECTORY

    NOTES

        IF NOT SPECIFIED BY THE USER, THE APPLICATION CREATES A GERALD

        CONFIGURATION FILE NAMED geraldfile.txt IN THE OUTPUT DIRECTORY:

        ANALYSIS        eland_extended
        ELAND_GENOME    /home/syoung/base/pipeline/run2-human-mitochondria/data
        READ_LENGTH     30
        USE_BASES       nY*

        WHERE THE USER SPECIFIED FILE CAN CONTAIN ANY OF THE FOLLOWING:

            ANALYSIS            Type of alignment to be performed (default (phageAlign),
                                none ('8:ANALYSIS none' will ignore lane 8), sequence,
                                sequence_pair, eland, eland_extended, eland_pair)
            ELAND_GENOME        Name of the directory containing a compressed (“squashed”)
                                reference for Eland 
            READ_LENGTH         Read length to use for alignment (has to be less than or
                                equal to the number of sequencing cycles).
            USE_BASES nY[YY…]   Specifies which of the bases to align to the reference sequence.
                                “n” means the first base will not be used (maybe part of the
                                primer). Number of Y’s has to be equal to the read length.
            CONTAM_FILE         Switch contaminant filtering on
            CONTAM_DIR          Specifies location of CONTAM_FILE

    USAGE

        image2eland.pl --type single 
            <--rundir> 
            <--outputdir> 
            <--referencefile>
            <--readlength>
            <--tiles>
            [--geraldfile]
            [--cpus]
            [--help]

        --rundir        /FULL/PATH/TO/SOLEXA/RUN/DIRECTORY        
        --outputdir     /FULL/PATH/TO/DESTINATION/DIRECTORY (COPY SELECTED RUN FILES TO HERE)
        --referencefile /FULL/PATH/TO/REFERENCE.fasta FILE
        --readlength    USE SELECTED NUMBER OF READ BASES
        --tiles         PROCESS IMAGES FOR THESE TILES ONLY
        --geraldfile    /FULL/PATH/TO/GERALD_CONFIG.txt FILE
        --cpus          NUMBER OF CPUs TO USE
        --help          PRINT OUT THIS HELP INFORMATION

    EXAMPLES

image2eland.pl --type single \
    --rundir /store/data/pipeline_in/080805_HWI-EAS185_0006_20GVYAAXX_Jia_cDNA2_mtDNA_Total_Small_Medium_read2_JH \
    --geraldfile=/home/syoung/base/pipeline/run2lane6-test/eland/geraldfile.txt \
    --outputdir /home/syoung/base/pipeline/runlin2lane6-test/eland \
    --referencefile /store/home/syoung/base/pipeline/phix/phiFasta.fa \
    --readlength 26 \
    --tiles s_5_015

image2eland.pl --type single \
    --rundir /store/data/pipeline_in/080805_HWI-EAS185_0006_20GVYAAXX_Jia_cDNA2_mtDNA_Total_Small_Medium_read2_JH \
    --geraldfile=/home/syoung/base/pipeline/run2lane6-test/eland/phix-config.txt \
    --outputdir /home/syoung/base/pipeline/run2lane6-test/eland \
    --referencefile /store/home/syoung/base/pipeline/phix/phiFasta.fa \
    --readlength 26 \
    --tiles s_5_015



=cut

EOF

	print RESET;

	exit;
}



