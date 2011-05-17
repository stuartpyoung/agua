#!/usr/bin/perl -w

#### DEBUG
$DEBUG = 1;

#### TIME
my $time = time();


my $delay = 30;
print "Sleeping $delay seconds...\n";
sleep($delay);




=head2

	APPLICATION     eland2ace


        **** DUMMY EXECUTABLE TO TEST WORKFLOW ****


    PURPOSE

        1. Convert Eland output *_sorted.txt files to .ace format contig files

        .2 If the reference sequence file is specified, the contig consensus sequence will be identical to the reference sequence for the matched regions. If not specified, the contig consensus is the high quality value read for each position in the multiple alignment of reads.

    INPUT

        1. *sorted.txt Eland OUTPUT FILE

        2. REFERENCE GENOME FASTA FILE

    OUTPUT

        1. INTERMEDIATE '*.single.txt' FILE

        2. PRINT REFERENCE GENOME FASTA TO A SEQUENCE-ONLY FILE FOR FAST SEQUENCE

            RETRIEVAL USING seek

        3. CREATE acefiles DIRECTORY INSIDE INPUT FILE PARENT DIRECTORY

        4. AN .ace FILE FOR EACH 'CONTIG' (DEFINED AS EACH AREA WHERE THERE

            ARE OVERLAPPING READS) IN NUMBERED FOLDERS INSIDE acefiles DIRECTORY


    NOTES

        ELAND FORMAT

            s_L_export.txt AND s_L_sorted.txt

            Database friendly export format
            – tab delimited fields:
                0 – Machine (as parsed from run folder name.
                1 – Run Number (as parsed from run folder name).
                2 – Lane.
                3 – Tile.
                4 – X Coordinate of cluster.
                5 – Y Coordinate of cluster.
                6 – Index String (blank for a non-indexed run).
                7 – Read Number ("1" or "2" for paired read, blank for a single read).
                8 – Read.
                9 – Quality String - in symbolic ASCII format (ASCII chracter code =
                quality value + 64) by default, set QUALITY_FORMAT --numeric in
                the GERALD config file to get numeric values instead.
                10 – Match Chromosome - name of chromosome match was to OR code
                indicating why no match was done.
                11 – Match Contig (blank if no match found) - gives contig name if there is
                a match and the match chromosome is split into contigs.
                12 – Match Position (always with respect to forward strand, numbering
                starts at 1).
                13 – Match Strand ("F" for forward or "R" for reverse, blank if no match).
                14 – Match Descriptor - concise description of alignment. A numeral
                denotes a run of matching bases, a letter denotes substituation of a
                nucleotide, so e.g. for a 35 base read, "35" denotes an exact match
                and "32C2" denotes substitution of a "C" at the 33rd position.
                15 – Single Read Alignment Score - alignment score of single read match
                (if a paired read, gives alignment score of read if it were to be treated
                as a single read).
                16 – Paired Read Alignment Score - alignment score of read pair
                (alignment score of a paired read and its partner, taken as a pair.
                Blank for a single read run).
                17 – Partner Chromosome - not blank only if read is paired and its partner
                aligns to another chromosome, in which case it gives the name of the
                chromosome.
                18 – Partner Contig - not blank only if read is paired and its partner aligns
                to another chromosome and that partner is split into contigs.
                19 – Partner Offset - if a paired read's partner hits to the same
                chromosome (as it will in the vast majority of cases) and contig (if the
                chromosome is split into contigs) then this number added to Match
                Position gives the alignment position of the read's partner.
                20 – Partner Strand - which strand did the paired read's partner hit to("F"
                for forward or "R" for reverse, blank if no match).
                21 – Filtering. Did the read pass quality filtering? "Y" for yes, "N" for no.

        EXAMPLE

            INPUT FILES:

                head s_6_1_sorted.txt
                HWI-EAS185_6_20GVYAAXX_Jia_cDNA2_mtDNA_Total_Small_Medium_read2_JH              6       161     501     257             1       TGGCTTAATATGCTTGGCACGTTCGT      TTTTTTTTTTTTTTTTTTTTTTTTTT        phiFasta.fa             241     F       26      64      0                       0       N
                HWI-EAS185_6_20GVYAAXX_Jia_cDNA2_mtDNA_Total_Small_Medium_read2_JH              6       28      873     812             1       ATTGCGATAAACGGTCACATTAAATT      TTTTTTTTTTTTTTTTTTTTTTTTTT        phiFasta.fa             2223    R       26      64      0                       -126    F
                HWI-EAS185_6_20GVYAAXX_Jia_cDNA2_mtDNA_Total_Small_Medium_read2_JH              6       102     635     260             1       ATTGCTTTTGATGCCGACCCTAAATT      TTTTTTTTTTTTTTTTTTTTTTTTTT        phiFasta.fa             2629    F       26      64      0                       0       N
                HWI-EAS185_6_20GVYAAXX_Jia_cDNA2_mtDNA_Total_Small_Medium_read2_JH              6       179     767     363             1       AGACTTGCCACCAAGTCCAACCAAAT      TTTTTTTTTTTTTTTTTTTTTTTTTT        phiFasta.fa             3232    R       26      64      0                       0       N


            OUTPUT FILE:

                head /home/syoung/base/pipeline/run2-mtdna-eland/assembly/s_6_1_sorted.paired.txt

                HWI-EAS185_6_20GVYAAXX_Jia_cDNA2_mtDNA_Total_Small_Medium_read2_JH              6       38      97      55              1       TGAGTGGTTAATAGGGTGATAGACCTGTGATCCAT     ZZZZZZZZZZZZZZZZZZZZZZZZZZZVYZZZZZV     human-mtDNA-AC_000021.fasta             -2      R       32NNN   77      0                       16441   F
                HWI-EAS185_6_20GVYAAXX_Jia_cDNA2_mtDNA_Total_Small_Medium_read2_JH              6       38      97      55              2       CTACTCTCCTCGCTCCGGGCCCATCACACTTGGGG     ZZZZZZZZZZZZZZZXVPZZZZWCRHZCKPNWKWV     human-mtDNA-AC_000021.fasta             16439   F       24A10   1       0                       -16441  R
                HWI-EAS185_6_20GVYAAXX_Jia_cDNA2_mtDNA_Total_Small_Medium_read2_JH              6       40      693     765             1       TGAGTGGTTAATAGGGTGATAGACCTGTGATCAAT     TZZZTZZZZZZZZYYTYTZYYOLYZZTYYGPZRQP     human-mtDNA-AC_000021.fasta             -2      R       32NNN   38      0                       16457   F
                HWI-EAS185_6_20GVYAAXX_Jia_cDNA2_mtDNA_Total_Small_Medium_read2_JH              6       40      693     765             2       GGGCCCATAACACTTGGGGGTAGCTAAAGTGAACT     ZZZZZZZZZSZZZSZZZSZZWOZKVJORRJSRINS     human-mtDNA-AC_000021.fasta             16455   F       35      35      0                       -16457  R
                HWI-EAS185_6_20GVYAAXX_Jia_cDNA2_mtDNA_Total_Small_Medium_read2_JH              6       65      586     329             1       TGAGTGGTTAATAGGGTGATAGACCTGTGATCCAT     ZZZZZZZZZZZZZZZZZZZZZZZZZZYSSZZXZXZ     human-mtDNA-AC_000021.fasta             -2      R       32NNN   66      0                       16449   F
                HWI-EAS185_6_20GVYAAXX_Jia_cDNA2_mtDNA_Total_Small_Medium_read2_JH              6       65      586     329             2       CTCGCTCCGGGCCCCTAACACTTGGGGGTAGCTAA     ZZZZZZZZZZZZZZEWPTTPTPRPINTJCCEEHPV     human-mtDNA-AC_000021.fasta             16447   F       14A20   6       0                       -16449  R
                HWI-EAS185_6_20GVYAAXX_Jia_cDNA2_mtDNA_Total_Small_Medium_read2_JH              6       118     815     10              1       GTGAGTGGTTAATAGGGTGATAGACCTGTGATCCA     ZZZZZZZZZZZZZZZZYYZZZZZZ\ZZZZZXXLSP     human-mtDNA-AC_000021.fasta             -1      R       33NN    77      0                       16443   F
                HWI-EAS185_6_20GVYAAXX_Jia_cDNA2_mtDNA_Total_Small_Medium_read2_JH              6       118     815     10              2       CTCTCCTCGCTCCGGGCCCATAACAATCGGGGGTA     WZZZZZZZZZZZZZZZZVGZWZZZZNNEZPWWVDV     human-mtDNA-AC_000021.fasta             16442   F       25C1T7  3       0                       -16443  R
                HWI-EAS185_6_20GVYAAXX_Jia_cDNA2_mtDNA_Total_Small_Medium_read2_JH              6       118     816     8               1       GTGAGTGGTTAATAGGGTGATAGACCTGTGATCCA     ZZZZZZZZZZZZZZZZYYZZZZZZUZZZZZXXLSP     human-mtDNA-AC_000021.fasta             -1      R       33NN    73      0                       16443   F



    USAGE

    ./eland2ace.pl <--inputfile String> <--referencefile String> [-h]

    -i inputfile            :   /full/path/to/eland_export.txt file (comma-separated if multiple files)
    -r referencefile        :   /full/path/to/reference_genome.fasta file (must end in .fa, .fas or .fasta)
    -h help                 :   print help info

    EXAMPLES

eland2ace.pl --inputfile /home/syoung/base/pipeline/run2lane6-test/eland/s_6_1_sorted.txt --referencefile /home/syoung/base/pipeline/run2-human-mitochondria/data/human-mtDNA-AC_000021.fasta


/store/home/syoung/base/pipeline/run2-mtdna-eland/assembly/s_6_1_sorted.txt

head /store/home/syoung/base/pipeline/run2-mtdna-eland/assembly/s_6_1_sorted.txt

HWI-EAS185_6_20GVYAAXX_Jia_cDNA2_mtDNA_Total_Small_Medium_read2_JH              6       15      132     181             1       TGAGTGGTTAATAGGGTGATAGACCTGTGATCCAT     ZZZZZZZZZZZZZZZZZZZZZZZZZZYZYZZZZRZ       human-mtDNA-AC_000021.fasta             -2      R       32NNN   80      0                       0       N
HWI-EAS185_6_20GVYAAXX_Jia_cDNA2_mtDNA_Total_Small_Medium_read2_JH              6       18      48      86              1       TGAGTGGTTAATAGGGTGATAGACCTGTGATCCAT     ZZZZZZZZZZZZZYZZZZZZZZZZZZZZZYZZZZX       human-mtDNA-AC_000021.fasta             -2      R       32NNN   80      0                       0       N


head /home/syoung/base/pipeline/human1-eland/assembly/human1-cdna_embl-eland.txt

>HWI-EAS185_1_JIA_cDNA_JH:3:1:129:540   TGGAGATGGCGACTAGTGGACAACAGAACAATAT      U0      1       0       0       Homo_sapiens.NCBI36.49.cdna.known.fas/ENST00000272348   225     F..
>HWI-EAS185_1_JIA_cDNA_JH:3:1:113:576   GGGTACCTCTATAGTGCTAGAGTGCCCTAACGAT      U1      0       1       0       Homo_sapiens.NCBI36.49.cdna.known.fas/ENST00000313486   1620    R..       16G



    USAGE

    ./eland2ace.pl <-i inputfile> [-r referencefile] [-h]

    -i inputfile            :   /full/path/to/eland_sorted.txt file (comma-separated if multiple files)
    -r referencefile        :   /full/path/to/reference_genome.fasta file (must end in .fa, .fas or .fasta). 
    -h help                 :   print help info

    EXAMPLES

./eland2ace.pl -i ~/base/pipeline/run2-mtdna-eland/assembly/s_6_1_sorted.txt,~/base/pipeline/run2-mtdna-eland/assembly/s_7_1_sorted.txt,~/base/pipeline/run2-mtdna-eland/assembly/s_8_1_sorted.txt -r /home/syoung/base/pipeline/run2-human-mitochondria/data/human-mtDNA-AC_000021.fasta 

=cut

use strict;

#### USE LIBRARY
use FindBin qw($Bin);
use lib "$Bin/../../../lib";


#### INTERNAL MODULES
#use Alignment::Ace;
#use SolexaUtil;
use Timer;
use Util;
use Conf::Agua;

#### EXTERNAL MODULES
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use Tie::File;
#use Tie::File::Hashify;

#### INITIALISE SolexaUtil OBJECT
#my $solexa = SolexaUtil->new();

##### STORE ARGUMENTS TO PRINT TO FILE LATER
my @arguments = @ARGV;

#### GET OPTIONS
my $inputfile;
my $referencefile;
my $help;
GetOptions (
    'inputfile=s' => \$inputfile,
    'referencefile=s' => \$referencefile,
    'help' => \$help
)
or usage();
if ( defined $help )    {   usage();    }

print "Reference file: $referencefile\n" if defined $referencefile;

#### FLUSH BUFFER
$| =1;


my $delay = 4;
print "Sleeping $delay seconds...\n";
sleep($delay);


#
##### CHECK REFERENCE FILE
#if ( defined $referencefile )
#{
#    if ( $referencefile !~ /(fa|fas|fasta)$/ )
#    {
#        usage();
#    }
#    if ( not -f $referencefile )
#    {
#        usage();
#    }    
#}
#
#
##### CHECK INPUT FILES
#if ( not defined $inputfile )
#{
#    usage();
#}
#my @inputfiles = split ",", $inputfile;
#foreach my $infile ( @inputfiles )
#{
#    if ( not -f $infile )
#    {
#        die "Could not find input file: $infile\n";
#        usage();
#    }
#}
#
##### PRINT REFERENCE GENOME FASTA TO A SEQUENCE-ONLY FILE TO BE USED FOR
##### FAST SEQUENCE RETRIEVAL USING seek
##### CHECK REFERENCE FILE
#my $reference_length;
#my $sequence_filehandle;
#if ( defined $referencefile )
#{
#
#    my ($sequencefile) = $referencefile =~ /^(.+)\.(fa|fas|fasta)$/;
#    $sequencefile .= ".sequence.fasta";
#    
#    if ( not -f $sequencefile )
#    {
#        open(SEQUENCE, ">$sequencefile") or die "Can't open reference 'sequence-only' file: $sequencefile\n";
#        open(REFERENCE, $referencefile) or die "Can't open reference file: $referencefile\n";
#        while ( <REFERENCE> )
#        {
#            if ( $_ =~ /^>/ ) {   next;   }
#            $_ =~ s/\s+//g;
#        }
#    }
#    else
#    {
#    }
#    
#    #### GET REFERENCE LENGTH
#    my ($dev, $ino, $mode, $nlink, $uid, $gid, $rdev, $size, $atime, $mtime, $ctime, $blksize, $blocks) = stat($sequencefile);
#    $reference_length = $size - 1;
#    
#    #### GET SEQUENCE FILEHANDLE
#    open($sequence_filehandle, $sequencefile) or die "Can't open sequence 'sequence-only' file: $sequencefile\n";
#}
#
#
##### DO .ace FILE USING PAIRED FILE FOR EACH
#foreach my $infile( @inputfiles )
#{
#    my $outfile = $inputfile;
#    $outfile =~ s/\.txt$//;
#    $outfile .= ".ace";
#
#
#    #### COUNT CONTIG NUMBER AND TOTAL BASE LENGTH
#    my $number_contigs = 0;
#    my $total_bases = 0;
#
#    #### INITIALISE Alignment::Ace CONTIG OBJECT
#    my $contig = Alignment::Ace->new();
#    
#    #### OPEN OUTPUT FILE
#    open(OUTFILE, ">$outfile") or die "Can't open output file: $outfile\n";
#
#    #### SET COUNTERS
#    my $counter = 0;
#    my $contig_counter = 1;
#    my $total_reads = 0;
#    my $contig_start;
#    my $current_position;
#
#    #### GO THROUGH INPUT FILE PUTTING ALL CONTIGUOUS READS INTO DISTINCT CONTIGS
#    open(INFILE, $infile) or die "Can't open intermediate output file: $infile\n";
#    while ( <INFILE> )
#    {            
#        if ( $_ =~ /^#/ )
#        {
#            next;
#        }
#
#        if ( $counter % 1000 == 0 )  {   print "$counter\n"; }
#
#        #### REMOVE TRAILING ENDLINE
#        $_ =~ s/\s$//;
#        
#        #### GET READ START/STOP 
#        my @elements = split "\t", $_;
#        my $read_start = $elements[12];
#        my $read_sequence = $elements[8];
#        my $read_strand = $elements[13];
#        
#        #### SET CURRENT POSITION IF THIS IS THE FIRST READ
#        if ( $counter == 0 )
#        {
#            $contig_start = $read_start;
#            $current_position = $read_start + length($read_sequence);
#        }
#
#        
#        #### ADD *_sorted.txt LINE TO CONTIG
#        $contig->add_sorted($_, $sequence_filehandle, $reference_length);
#
#        
#        #### END THIS CONTIG AND START A NEW ONE IF THIS READ DOES NOT OVERLAP
#        #### THE CURRENT MOST DOWNSTREAM READ STOP POSITION
#        if ( $counter !=  0 and $current_position < $read_start - 1 )
#        {
#            #print Dumper $contig;
#
#            #### SET CONTIG NAME WITH POSITIONS
#            $contig->{contigname} = "Contig$contig_counter\[$contig_start..$current_position\]";
#
#            #### PRINT CURRENT CONTIG TO OUTPUT FILE\
#
#            #### ADD CONTIG LENGTH TO TOTAL SEQUENCE
#            #### NB: MUST BE AFTER $contig->as_text();
#            $total_reads += $contig->{numberreads};
#
#            #### START NEW CONTIG
#            $contig = undef;
#            $contig_counter++;
#
#            #### CHECK THAT A NEW LINE OF DATA IS AVAILABLE
#            $_ = <INFILE>;               
#            if ( defined $_ )
#            {
#                $contig = Alignment::Ace->new();
#                #$contig->{contigname} = "Contig$contig_counter";
#                
#                #### ADD *_sorted.txt LINE TO NEW CONTIG
#                $contig->add_sorted($_, $sequence_filehandle, $reference_length);
#
#                #### GET READ START/STOP 
#                my @elements = split "\t", $_;
#                my $read_start = $elements[12];
#                my $read_sequence = $elements[8];
#                my $read_strand = $elements[13];
#
#                #### SET CONTIG START
#                $contig_start = $read_start;
#            }
#
#            $counter++;
#        }
#
#        #### SET CURRENT POSITION            
#        $current_position = $read_start + length($read_sequence);
#
#        $counter++;
#        
#    }   #   while <INFILE>
#
#    #### IF THE CONTIG IS DEFINED THEN PRINT IT OUT
#    if ( defined $contig )
#    {
#        #### SET CONTIG NAME WITH POSITIONS
#        $contig->{contigname} = "Contig$contig_counter\[$contig_start..$current_position\]";
#
#        #### PRINT CURRENT CONTIG TO OUTPUT FILE\
#
#        #### ADD CONTIG LENGTH TO TOTAL SEQUENCE
#        #### NB: MUST BE AFTER $contig->as_text();
#        $total_reads += $contig->{numberreads};
#
#        #### CLOSE FILES        
#        close(INFILE);
#        close(OUTFILE);
#        
#        
#        #### MOVE OUTPUT FILE
#        my $tempfile = "outfile.tmp";
#        `mv $outfile $tempfile`;
#        
#        #### ADD .ace FILE HEADER 'AS ....' TO TOP OF FILE
#        open(OUTFILE, ">$outfile") or die "Can't open output file: $outfile\n";
#        open(TEMPFILE, "$tempfile") or die "Can't open temp file: $tempfile\n";
#        while ( <TEMPFILE> )
#        {
#        }
#
#    }    
#}
#

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
	print GREEN <<"EOF";

	APPLICATION     eland2ace

    PURPOSE

        CONVERT ELAND PAIRED READ *_sorted.txt FILES TO .ace FORMAT

    INPUT

        1. COMMA-SEPARATED INPUT FILE NAMES (ONLY s_N_1_sorted.txt FILES. THE

            APPLICATION WILL SEARCH FOR THE s_N_2_sorted.txt FILE IN THE SAME FOLDER)

    OUTPUT

        1. AN .ace FILE FOR EACH 'CONTIG' DEFINED AS EACH AREA WHERE THERE

            ARE OVERLAPPING READS

    USAGE

    ./eland2ace.pl <-i inputfile> [-h]

    -i inputfile            :   /full/path/to/eland_sorted.txt file (comma-separated if multiple files)
    -m max_distance         :   max distance between paired reads (reads allocated as 'paired' if
                                inter-read distance < max distance, and 'distant' if otherwise)
    -o output_type          :   print .ace files for these comma-separated read types (paired|distant|unpaired)
    -r referencefile        :   /full/path/to/reference_genome.fasta file (must end in .fa, .fas or .fasta)
    -h help                 :   print help info

    EXAMPLES

    ./eland2ace.pl -i ~/base/pipeline/run2-mtdna-eland/assembly/s_6_1_sorted.txt,~/base/pipeline/run2-mtdna-eland/assembly/s_7_1_sorted.txt,~/base/pipeline/run2-mtdna-eland/assembly/s_8_1_sorted.txt -o paired -r /home/syoung/base/pipeline/run2-human-mitochondria/data/human-mtDNA-AC_000021.fasta -m 400 -c 16500


eland2ace.pl -i ~/base/pipeline/run2-mtdna-eland/assembly/s_6_1_sorted.txt -o paired -m 400 -c 16500


EOF

	print RESET;

	exit;
}


