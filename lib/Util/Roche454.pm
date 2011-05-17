package Roche454;

#### DEBUG
$DEBUG = 1;

=head2

	PACKAGE		Roche454

    VERSION:        0.01

    PURPOSE

        CALCULATE COVERAGE AND RELATED STATISTICS FOR ROCHE 454 SEQUENCES 

        BASED ON 'GS MAPPER' REFERENCE ASSEMBLER OUTPUT

            1. PER SEQUENCE STATISTICS

                - Average coverage of CDS targets (% OF SEQUENCE COVERED)
                - Maximal coverage, minimal coverage (% OF SEQUENCE COVERED)
                - Average coverage per chromosome
                - Graph with average coverage per each gene

                    EXTRACT FROM FILE:

                    454RefStatus.txt

                    Reference       Num Unique      Pct of All      Pct of  Pct Coverage
                    Accession       Matching Reads  Unique Matches  All Reads       of Reference    Description
                    CCDS6193.1|Hs36.3|chr8  3071    0.3%    0.1%    100.00%
                    CCDS8731.1|Hs36.3|chr12 2993    0.2%    0.1%    20.56%
                    CCDS12128.1|Hs36.3|chr19        2724    0.2%    0.1%    90.37%
                    CCDS8755.1|Hs36.3|chr12 2300    0.2%    0.1%    98.33%
                    CCDS35212.1|Hs36.3|chrX 2077    0.2%    0.1%    84.33%
                    CCDS34869.1|Hs36.3|chr8 1752    0.1%    0.1%    10.22%
                    CCDS9868.1|Hs36.3|chr14 1518    0.1%    0.1%    97.56%
                    CCDS14785.1|Hs36.3|chrXY        1364    0.1%    0.1%    7.96%


                - Number of aligned and unaligned reads (on-target/ off-target)

                    EXTRACT FROM FILE:

                    454MappingQC.xls

                    Num. Reads:     2405723
                    Num. Bases:     712282210

                    >>> Mapped Reads:   1200713 49.91%  1475521 61.33%
                    >>> Mapped Bases:   173799488       24.40%  267745758       37.59%
                    Inf. Read Error:        3.50%   14.6    6089688 173799488
                    Exp. Read Error:        0.81%   20.9    1399662

                    Last 100 Base IRE:      3.58%   14.5    1919243 53609629
                    Last 50 Base IRE:       4.44%   13.5    1111174 25033961
                    ...

                    TO GENERATE FILE

                    stats-sequence-coverage.tsv


            2. PER BASE STATISTICS

                - Average coverage of CDS targets (FOLD PER BASE IN TOTAL SEQUENCE)
                - Maximal coverage, minimal coverage (FOLD PER BASE IN TOTAL SEQUENCE)
                - What did we miss of the targeted region

                    454AlignmentInfo.tsv <==
                    Position        Reference       Consensus       Quality Score   Depth   Signal  StdDeviation
                    >CCDS2.2|Hs36.3|chr1    1
                    1       A       A       64      3       1.00    0.00
                    2       T       T       64      3       1.00    0.00
                    3       G       G       64      3       1.00    0.00
                    4       T       T       64      3       1.00    0.00
                    5       C       C       64      3       2.00    0.00
                    6       C       C       64      3       2.00    0.00
                    7       A       A       64      3       2.00    0.00
                    8       A       A       64      3       2.00    0.00

                GENERATE FILES

                    <reference_file>.length
                        CCDS2.2 chr1    2046

                TO GENERATE

                    stats-base-coverage.tsv

                        CCDS2.2 chr1  2046  1..72  3.4
                        CCDS2.2 chr1  2046  158..240   5.4


                USE stats-base-coverage.tsv TO EXTRACT:

                - How much sequence below certain cut-off (e.g. 5x)
                - Average BASE coverage per chromosome

                - Graph with average coverage per each gene

                TO GENERATE

                    stats-base-coverage-summary.tsv

                - What did we miss of the targeted region

                USE FILES

                    stats-base-coverage.tsv
                    <reference_file>

                TO GENERATE .FASTA FILE

                    stats-unmapped-regions.fasta


            - Homopolymericity of unmapped regions and mapped regions

                1. count polymers as 1 unique nucleotide

                2. all unique nucleotides / all nucleotide


                EXTRACT homopolymericity DIFFERENCE FOR EACH READ FROM FILES:

                    stats-unmapped-regions.fasta

                    454AllContigs.fna

                        >contig00001  CCDS2.2|Hs36.3|chr1, 158..305  length=147   numreads=2
                        ccagcgatggtgacagcgacgggagtggccccacctgtgggcggcggccaggcttgaagc
                        aggaggatggtccgcacatccgtatcatgaagagaaGGTccacacccactgggacgtgaa
                        catctctttccgagaggcgtcctgcag
                        >contig00002  CCDS2.2|Hs36.3|chr1, 518..767  length=247   numreads=4
                        agaggcgctgctgctgccgcgggagctggggcccagcatggccccggaggaccattaccg
                        ccggcttgtgtcagcactgagcgaggccagcacctttgaggaccctcagcgcctctacca
                        cctgggcctcccagCCACGGTGAGGACCCACCCTGGCATGAtccccctcatcaccTCCCc
                        ACCcAGATCTCCTGAGGGTCCGGCAGGAGGTGGCGGCTGCAGCTCTGAGGGGCCCCAGTG
                        GCCTGGA

                GENERATE FILE:

                    stats-homopolymer-diffs.txt


    INPUT

        1. REFERENCE FASTA FILE (TO CALCULATE REFERENCE SEQUENCE lengthhash)

        2. 454RefStatus.txt GS MAPPER OUTPUT FILE

            ==> 454RefStatus.txt <==
            Reference       Num Unique      Pct of All      Pct of  Pct Coverage
            Accession       Matching Reads  Unique Matches  All Reads       of Reference    Description
            CCDS6193.1|Hs36.3|chr8  3071    0.3%    0.1%    100.00%
            CCDS8731.1|Hs36.3|chr12 2993    0.2%    0.1%    20.56%
            CCDS12128.1|Hs36.3|chr19        2724    0.2%    0.1%    90.37%
            CCDS8755.1|Hs36.3|chr12 2300    0.2%    0.1%    98.33%
            CCDS35212.1|Hs36.3|chrX 2077    0.2%    0.1%    84.33%
            CCDS34869.1|Hs36.3|chr8 1752    0.1%    0.1%    10.22%
            CCDS9868.1|Hs36.3|chr14 1518    0.1%    0.1%    97.56%
            CCDS14785.1|Hs36.3|chrXY        1364    0.1%    0.1%    7.96%


            ==> 454MappingQC.xls <==
            Num. Reads:     2405723
            Num. Bases:     712282210

            Mapped Reads:   1200713 49.91%  1475521 61.33%
            Mapped Bases:   173799488       24.40%  267745758       37.59%
            Inf. Read Error:        3.50%   14.6    6089688 173799488
            Exp. Read Error:        0.81%   20.9    1399662

            Last 100 Base IRE:      3.58%   14.5    1919243 53609629
            Last 50 Base IRE:       4.44%   13.5    1111174 25033961


            6.4.2.9 454MappingQC.xls
            This file contains a number of detailed metrics regarding the mapping results; its format
            and structure are intended for reading by MS Excel or a similar spreadsheet program, so
            that the metrics can be easily visualized using Excel’s charting/graphing tools. The file is in
            tab-delimited format, and contains six main sections (Note: The section titles given here
            are not displayed in the file.)
            ( Summary Statistics
            a. Num. Reads – the number of input reads used in the mapping computation
            b. Num. Bases – the number of bases in the input reads
            c. Mapped Reads – the number and percentage of reads that uniquely mapped to
            the reference, followed by the number and percentage of reads that uniquely or
            multiply mapped
            d. Mapped Bases – the number and percentage of bases that uniquely mapped to
            the reference, followed by the number and percentage of reads that uniquely or
            multiply mapped
            e. Inf. Read Error – the “inferred read error” percentage and quality score (calculated
            as the number of read alignment differences over the number of mapped bases),
            along with the counts of the number of read alignment differences and mapped
            bases
            f. Exp. Read Error – the expected read error computed from the input read quality
            scores, given as a percentage, quality score and expected number of alignment
            differences. This is computed by summing the expected number of errors for each
            quality score value (i.e. number of bases with a quality score times the accuracy
            rate of that quality score).
            g. Last 100 Base IRE – the “inferred read error” numbers, using only the last (3’)
            100 bases of each read
            h. Last 50 Base IRE – the “inferred read error” numbers, using only the last (3’)
            50 bases of each read
            i. Last 20 Base IRE – the “inferred read error” numbers, using only the last (3’)
            20 bases of each read
            j. Genome Size – the number of bases in the reference
            k. Num. Large Contigs – the number of large contigs reported in the
            454LargeContigs.fna file
            l. Num. Large Contig Bases – the number of bases in the large contigs
            m. Avg. Depth – the average alignment depth (i.e. how many reads aligned to each
            position of the reference)
            n. Avg. Map Length – the average length of the alignment of a read (the read’s “map
            length”)
            ) Read Error Histogram
            a. This section breaks down the number of reads based on their read status
            and/or number of alignment differences, and displays percentages and counts per
            category.
            b. Reads that did not map to the reference (“Unmapped”), reads where only part
            of the sequence mapped to the reference (“Partial”) and reads that mapped to
            multiple locations in the reference (“Multiple”) are displayed as one category each.
            c. Reads that mapped fully to the reference (meaning every base of the read occurred
            in the alignment) are then divided by the number of alignment differences found,
            and percentages and number of reads having 0, 1, 2, …, 9 or 10+ (10 or more)
            errors are shown.


            ==> 454AlignmentInfo.tsv <==
            Position        Reference       Consensus       Quality Score   Depth   Signal  StdDeviation
            >CCDS2.2|Hs36.3|chr1    1
            1       A       A       64      3       1.00    0.00
            2       T       T       64      3       1.00    0.00
            3       G       G       64      3       1.00    0.00
            4       T       T       64      3       1.00    0.00
            5       C       C       64      3       2.00    0.00
            6       C       C       64      3       2.00    0.00
            7       A       A       64      3       2.00    0.00
            8       A       A       64      3       2.00    0.00


            6.4.2.10 454AlignmentInfo.tsv

            The 454AlignmentInfo.tsv file contains position-by-position summary information
            about the alignment of the reference sequence and the consensus generated from the
            aligned reads, listed one nucleotide per line (in a tab-delimited format). For a sample
            454AlignmentInfo.tsv file, see section 13.4.11.

            The columns of each line contain the following information:
            1. Position – the position in the reference sequence (or the position of the last reference
            base if the alignment has a gap for the reference)
            2. Reference – the alignment character for the reference, either a reference base or a
            gap (‘-’)
            3. Consensus – the alignment character of the consensus of the aligned reads, either a
            nucleotide or a gap (‘-’)
            4. Quality Score – the quality score of the consensus base, or the quality score giving the
            probability of an undercall at that location, in the case of a gap for the consensus
            5. Depth – the number of reads that align at that position in the alignment
            6. Signal – the average signalof the read flowgrams, for the aligned flows that
            correspond to that position in the alignment
            7. StdDevation – the standard deviation of the read flowgram signals at the
            corresponding flows

            This file only contains regions of the reference sequence where the aligned depth is greater
            than zero and, prior to each region of lines, a header line beginning with a ‘>’ displays the
            accession number of the sequence in the reference and starting position of the next region
            with a non-zero aligned depth.

    OUTPUT

        1. FILE CONTAINING A LIST OF SNPS THAT PASS FILTER CRITERIA

            stats-sequence-coverage.tsv


    USAGE

    ./capturestats.pl  <--mapdir String> <--reference String>  [-h]

        --mapdir               :   /full/path/to/mapping_directory
        --reference            :   /full/path/to/chromosome_positions.txt file
        --help                 :   print help info

    EXAMPLES

ON SOLEXA

./capturestats.pl --mapdir /home/syoung/base/pipeline/nimblegen-gsmapper/P_2009_01_09_03_26_14_runMapping/mapping --reference /home/syoung/base/pipeline/nimblegen-gsmapper/CCDS_nucleotide.20080430.fa


=cut 

use strict;
use warnings;
use Carp;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our $AUTOLOAD;

#### INTERNAL MODULES
use Util;

#### EXTERNAL MODULES
use POSIX;
use Data::Dumper;

#### DEFAULT PARAMETERS
our @DATA = qw(
	MAPDIR
    REFERENCE
    SEQSTATSFILE
    BASESTATSFILE
    BASESUMMARYFILE
    MAPPEDFILE
    UNMAPPEDFILE
    MAPPEDHOMOPOLYFILE
    UNMAPPEDHOMOPOLYFILE
    LENGTHHASH
);
our $DATAHASH;
foreach my $key ( @DATA )	{	$DATAHASH->{lc($key)} = 1;	}

our $SEQSTATSFILE = "stats-sequence.tsv";
our $BASESTATSFILE = "stats-base.tsv";
our $BASESUMMARYFILE = "stats-base-summary.tsv";
our $MAPPEDFILE = "stats-mapped.fasta";
our $UNMAPPEDFILE = "stats-unmapped.fasta";
our $MAPPEDHOMOPOLYFILE = "stats-mapped-homopoly.tsv";
our $UNMAPPEDHOMOPOLYFILE = "stats-unmapped-homopoly.tsv";



=head2

	SUBROUTINE		new

	PURPOSE

		CREATE A NEW self OBJECT

=cut

sub new
{
    my $class 		=	shift;
	my $arguments 	=	shift;


	my $self = {};
    bless $self, $class;

	#### INITIALISE THE OBJECT'S ELEMENTS
	$self->initialise($arguments);


    return $self;
}


=head2

	SUBROUTINE		initialise

	PURPOSE

		INITIALISE THE self OBJECT

=cut

sub initialise
{
    my $self		=	shift;
	my $arguments	=	shift;


    #### LOAD DEFAULTS
    $self->{_seqstatsfile} = $SEQSTATSFILE;
    $self->{_basestatsfile} = $BASESTATSFILE;
    $self->{_basesummaryfile} = $BASESUMMARYFILE;
    $self->{_mappedfile} = $MAPPEDFILE;
    $self->{_unmappedfile} = $UNMAPPEDFILE;
    $self->{_mappedhomopolyfile} = $MAPPEDHOMOPOLYFILE;
    $self->{_unmappedhomopolyfile} = $UNMAPPEDHOMOPOLYFILE;


    #### VALIDATE USER-PROVIDED ARGUMENTS
	($arguments) = $self->validate_arguments($arguments);	

    #### PROCESS USER-PROVIDED ARGUMENTS
	foreach my $key ( keys %$arguments )
	{
		$self->value($key, $arguments->{$key});
	}

    #### CHECK REQUIRED ARGUMENTS
    if ( not defined $self->{_mapdir} )
    {
        die "Mapping directory not defined. Exiting...\n";
    }
    if ( not defined $self->{_reference} )
    {
        die "Reference file not defined. Exiting...\n";
        exit;
    }
}


=head2

    SUBROUTINE      sequence_stats

    PURPOSE

        #### 1. PER SEQUENCE STATISTICS
        ####
        ####    - Average coverage of CDS targets (% OF SEQUENCE COVERED)
        ####    - Maximal coverage, minimal coverage (% OF SEQUENCE COVERED)
        ####    - Average coverage per chromosome
        ####    - Graph with average coverage per each gene
        ####
        ####        EXTRACT FROM FILE:
        ####
        ####        454RefStatus.txt
        ####    
        ####        Reference       Num Unique      Pct of All      Pct of  Pct Coverage
        ####        Accession       Matching Reads  Unique Matches  All Reads       of Reference    Description
        ####        CCDS6193.1|Hs36.3|chr8  3071    0.3%    0.1%    100.00%
        ####        CCDS8731.1|Hs36.3|chr12 2993    0.2%    0.1%    20.56%
        ####        CCDS12128.1|Hs36.3|chr19        2724    0.2%    0.1%    90.37%
        ####        CCDS8755.1|Hs36.3|chr12 2300    0.2%    0.1%    98.33%
        ####        CCDS35212.1|Hs36.3|chrX 2077    0.2%    0.1%    84.33%
        ####        CCDS34869.1|Hs36.3|chr8 1752    0.1%    0.1%    10.22%
        ####        CCDS9868.1|Hs36.3|chr14 1518    0.1%    0.1%    97.56%
        ####        CCDS14785.1|Hs36.3|chrXY        1364    0.1%    0.1%    7.96%
        ####
        ####
        ####    - Number of aligned and unaligned reads (on-target/ off-target)
        ####
        ####        EXTRACT FROM FILE:
        ####
        ####        454MappingQC.xls
        ####
        ####        Num. Reads:     2405723
        ####        Num. Bases:     712282210
        ####        
        ####        >>> Mapped Reads:   1200713 49.91%  1475521 61.33%
        ####        >>> Mapped Bases:   173799488       24.40%  267745758       37.59%
        ####        Inf. Read Error:        3.50%   14.6    6089688 173799488
        ####        Exp. Read Error:        0.81%   20.9    1399662
        ####        
        ####        Last 100 Base IRE:      3.58%   14.5    1919243 53609629
        ####        Last 50 Base IRE:       4.44%   13.5    1111174 25033961
        ####        ...
        ####
        ####        TO GENERATE FILE
        ####        
        ####        stats-sequence-coverage.tsv
        ####        
        ####

=cut

sub sequence_stats
{
    my $self        =   shift;


    my $sequence_statsfile = $self->{_seqstatsfile};
    my $mapdir = $self->{_mapdir};
    my $reference = $self->{_reference};

    #### OPEN OUTPUT FILE
    open(SEQSTATS, ">$mapdir/$sequence_statsfile") or die "Can't open sequence coverage file: $sequence_statsfile\n";

    ############################
    #### GET AVERAGE COVERAGE
    ############################

    #### SET VARIABLES
    my $average_coverage;           #### (% OF TARGET SEQUENCE COVERED)
    my $max_coverage = 0;           #### (% OF TARGET SEQUENCE COVERED)
    my $min_coverage = 100;     #### (% OF TARGET SEQUENCE COVERED)
    my $chromosome_coverage;
    my @coverage_bins;              #### ARRAY OF BINS (10%, 20%, 30%, ..., 100%)

    #### INITIALISE COVERAGE BINS
    for my $i ( 0..9 )  {   $coverage_bins[$i] = 0; }

    #### OPEN 454RefStatus.txt FILE
    open(REFSTATS, "$mapdir/454RefStatus.txt") or die "Can't open '454RefStatus.txt' file in mapping directory: $mapdir\n";
    <REFSTATS>;  #### REMOVE COLUMN HEADER LINES
    <REFSTATS>;
    my $total = 0;
    while ( <REFSTATS> )
    {
        ####   FORMAT    
        ####        CCDS6193.1|Hs36.3|chr8  3071    0.3%    0.1%    100.00%
        ####        CCDS8731.1|Hs36.3|chr12 2993    0.2%    0.1%    20.56%

        my ($chromosome, $coverage) = $_ =~ /^\S+\|(\S+)\s+.+\s+(\S+)%\s*$/;

        #### DO CHROMOSOME COVERAGE
        if ( not defined $chromosome )
        {
            print "Chromosome not defined: $_\n";
            exit;
        }
        $chromosome_coverage->{$chromosome}->{total} += $coverage;
        $chromosome_coverage->{$chromosome}->{count}++;

        #### DO MIN/MAX
        if ( $coverage < $min_coverage )    {   $min_coverage = $coverage; }
        if ( $coverage > $max_coverage )    {   $max_coverage = $coverage; }

        #### GET INDEX FOR COVERAGE BINS
        my $bin_index = int($coverage/10);
        if ( $bin_index == 10 ) {   $bin_index = 9; }
        $coverage_bins[$bin_index] += $coverage;

        $total++;
        last if $total == 100;
    }

    #### PRINT TO stats-sequence-coverage FILE
    for ( my $i = 0; $i < $#coverage_bins + 1; $i++ )
    {
        print SEQSTATS ($i + 1) * 10, "% coverage: " , $coverage_bins[$i] / $total, "\n";
    }
    print SEQSTATS "Max coverage: $max_coverage\n";
    print SEQSTATS "Min coverage: $min_coverage\n";
    print SEQSTATS "Average chromosome coverage:\n";

    my @keys = keys %$chromosome_coverage;
    @keys = sort {
        my ($aa) = $a =~ /(\d+)$/;
        $aa = 999 if not defined $aa;
        my ($bb) = $b =~ /(\d+)$/;
        $bb = 999 if not defined $bb;
        $aa <=> $bb
    } @keys;
    foreach my $key ( @keys )
    {
        print SEQSTATS "$key: ";
        print SEQSTATS sprintf  "%.3f", $chromosome_coverage->{$key}->{total} / $chromosome_coverage->{$key}->{count};
        print SEQSTATS "\n";
    }

    ###############################
    #### GET NUMBER ON/OFF TARGET
    ###############################

    #### OPEN 454RefStatus.txt FILE
    #### FORMAT:
    ####        Num. Reads:     2405723
    ####        Num. Bases:     712282210
    ####        
    ####        >>> Mapped Reads:   1200713 49.91%  1475521 61.33%
    ####        >>> Mapped Bases:   173799488       24.40%  267745758       37.59%

    print "Opening file:\n\n$mapdir/454MappingQC.xls\n\n";
    open(MAPPINGQC, "$mapdir/454MappingQC.xls") or die "Can't open '454MappingQC.txt' file in mapping directory: $mapdir\n";
    $/ = "END OF FILE";
    my $content = <MAPPINGQC>;
    #### VARIABLES
    my $status;
    ($status->{uniquemappedreads}, $status->{uniquemappedreadspct}, $status->{mappedreads}, $status->{mappedreadspct})
            = $content =~ /Mapped Reads:\s+(\S+)\s+(\S+)%\s+(\S+)\s+(\S+)%\s*$/ms;
    ($status->{uniquemappedbases}, $status->{uniquemappedbasespct}, $status->{mappedbases}, $status->{mappedbasespct})
            = $content =~ /Mapped Bases:\s+(\S+)\s+(\S+)%\s+(\S+)\s+(\S+)%\s*$/ms;
    print SEQSTATS "Unique mapped reads: $status->{uniquemappedreads} ($status->{uniquemappedreadspct})\n";
    print SEQSTATS "Mapped reads: $status->{mappedreads} ($status->{mappedreadspct})\n";
    print SEQSTATS "Unique mapped bases: $status->{uniquemappedbases} ($status->{uniquemappedbasespct})\n";
    print SEQSTATS "Mapped bases: $status->{mappedbases} ($status->{mappedbasespct})\n";

    #### CLOSE SEQUENCE STATS FILE
    close(SEQSTATS);

    ##### REPORT FILE PRINTED

    return "$mapdir/$sequence_statsfile";
}




=head2

    SUBROUTINE      base_stats

    PURPOSE

        #### 2. PER BASE STATISTICS
        ####
        ####    - Average coverage of CDS targets (FOLD PER BASE IN TOTAL SEQUENCE)
        ####    - Maximal coverage, minimal coverage (FOLD PER BASE IN TOTAL SEQUENCE)
        ####    - What did we miss of the targeted region?
        ####    - How much sequence below certain cut-off (e.g. 5x)?
        ####    - Average BASE coverage per chromosome
        ####    - Graph with average coverage per each gene
        ####    
        ####    
        ####    
        ####    - Average coverage of CDS targets (FOLD PER BASE IN TOTAL SEQUENCE)
        ####    - Maximal coverage, minimal coverage (FOLD PER BASE IN TOTAL SEQUENCE)
        ####        454AlignmentInfo.tsv <==
        ####        Position        Reference       Consensus       Quality Score   Depth   Signal  StdDeviation
        ####        >CCDS2.2|Hs36.3|chr1    1
        ####        1       A       A       64      3       1.00    0.00
        ####        2       T       T       64      3       1.00    0.00
        ####        3       G       G       64      3       1.00    0.00
        ####        4       T       T       64      3       1.00    0.00
        ####        5       C       C       64      3       2.00    0.00
        ####        6       C       C       64      3       2.00    0.00
        ####        7       A       A       64      3       2.00    0.00
        ####        8       A       A       64      3       2.00    0.00
        ####
        ####    GENERATE FILES
        ####    
        ####        <reference_file>.length
        ####            CCDS2.2 chr1    2046
        ####    
        ####    TO GENERATE
        ####    
        ####        stats-base-coverage.tsv
        ####        
        ####            CCDS2.2 chr1  2046  1..72  3.4
        ####            CCDS2.2 chr1  2046  158..240   5.4
        ####            
        ####    
        ####    USE stats-base-coverage.tsv TO EXTRACT:
        ####
        ####    - How much sequence below certain cut-off (e.g. 5x)?
        ####    - Average BASE coverage per chromosome
        ####    - Graph with average coverage per each gene
        ####    
        ####    TO GENERATE
        ####    
        ####        stats-base-coverage-summary.tsv
        ####    
        ####    - What did we miss of the targeted region
        ####        
        ####    USE FILES
        ####    
        ####        stats-base-coverage.tsv
        ####        <reference_file>
        ####        
        ####    TO GENERATE .FASTA FILE
        ####    
        ####        stats-unmapped-regions.fasta
        ####
        ####
        ####- Homopolymericity of unmapped regions and mapped regions
        ####
        ####    1. count polymers as 1 unique nucleotide
        ####    
        ####    2. all unique nucleotides / all nucleotide
        ####
        ####
        ####    EXTRACT homopolymericity DIFFERENCE FOR EACH READ FROM FILES:
        ####    
        ####        stats-unmapped-regions.fasta
        ####    
        ####        454AllContigs.fna
        ####        
        ####            >contig00001  CCDS2.2|Hs36.3|chr1, 158..305  length=147   numreads=2
        ####            ccagcgatggtgacagcgacgggagtggccccacctgtgggcggcggccaggcttgaagc
        ####            aggaggatggtccgcacatccgtatcatgaagagaaGGTccacacccactgggacgtgaa
        ####            catctctttccgagaggcgtcctgcag
        ####            >contig00002  CCDS2.2|Hs36.3|chr1, 518..767  length=247   numreads=4
        ####            agaggcgctgctgctgccgcgggagctggggcccagcatggccccggaggaccattaccg
        ####            ccggcttgtgtcagcactgagcgaggccagcacctttgaggaccctcagcgcctctacca
        ####            cctgggcctcccagCCACGGTGAGGACCCACCCTGGCATGAtccccctcatcaccTCCCc
        ####            ACCcAGATCTCCTGAGGGTCCGGCAGGAGGTGGCGGCTGCAGCTCTGAGGGGCCCCAGTG
        ####            GCCTGGA
        ####
        ####    GENERATE FILE:
        ####    
        ####        stats-homopolymer-diffs.txt
        ####        

=cut

sub base_stats
{
    my $self        =   shift;

    print "++++ Roche454::base_stats\n";

    my $mapdir = $self->{_mapdir};
    my $basestatsfile = $self->{_basestatsfile};

    if ( not defined $self->{_lengthhash} )
    {
        my $lengthfile = $self->generate_lengthfile();
        if ( not defined $lengthfile )
        {
            die "Roche454::base_stats. Reference sequence length file not defined. Exiting...\n";
        }
        $self->load_lengthfile($lengthfile);
    }

    my $lengthhash = $self->{_lengthhash};
    if ( not defined $self->{_lengthhash} )
    {
        die "Length hash not defined. Exiting...\n";
    }


    #### OPEN BASE COVERAGE OUTPUT FILE
    ####
    ####        stats-base-coverage.tsv
    ####
    ####        id    chromsome length  start..stop avg_cvg bin1 bin2 bin3 ... bin10
    ####        CCDS2.2 chr1  2046  1..72  3.4 
    ####        CCDS2.2 chr1  2046  158..240   5.4
    open(BASESTATS, ">$mapdir/$basestatsfile") or die "Can't open base coverage file: $mapdir/$basestatsfile\n";

    #### OPEN ALIGNMENT INFO INPUT FILE
    #### FORMAT:
    ####        Position        Reference       Consensus       Quality Score   Depth   Signal  StdDeviation
    ####        >CCDS2.2|Hs36.3|chr1    1
    ####        1       A       A       64      3       1.00    0.00
    ####        2       T       T       64      3       1.00    0.00
    ####        3       G       G       64      3       1.00    0.00
    ####        4       T       T       64      3       1.00    0.00
    ####        5       C       C       64      3       2.00    0.00
    ####        6       C       C       64      3       2.00    0.00
    ####        7       A       A       64      3       2.00    0.00
    ####        8       A       A       64      3       2.00    0.00

    print "Opening file:\n\n$mapdir/454AlignmentInfo.tsv\n\n";
    open(ALIGNINFO, "$mapdir/454AlignmentInfo.tsv") or die "Can't open '454AlignmentInfo.tsv' file in mapping directory: $mapdir\n";
    $/ = ">";
    <ALIGNINFO>; #### REMOVE THE TOP COLUMN HEADER LINE
    my $bins = [1, 3, 5, 10, 15, 20, 25, 50, 100];

    #### FLUSH BUFFER
    $| = 1;

    print "Calculating base statistics...\n";
    while ( <ALIGNINFO> )
    {
        next if $_ =~ /^\s*$/;

        #### REMOVE FASTA HEADER
        my ($header) = $_ =~ /^([^\n]+)/ms;
        chop($_);
        $_ =~ s/^([^\n]+)//ms;
        $_ =~ s/\s+$//;
        my ($id, $chromosome) = $header =~ /^([^\|]+)\|[^\|]+\|(\S+)/ms;

        #### SUM COVERAGES
        my $binner = BinData->new( { 'BINS' => $bins } );
        my $total = 0;
        my $count = 0;
        my $start;
        my $stop;
        my @lines = split "\n", $_;
        for ( my $i = 0; $i < $#lines + 1; $i++ )
        {
            next if $lines[$i] =~ /^\s*$/;

            my ($position, $base_coverage) = $lines[$i] =~ /^(\S+)\s+\S+\s+\S+\s+\S+\s+(\S+)/;
            if ( not defined $start )  {   $start = $position; }
            if ( $i == $#lines )  {   $stop = $position; }
            $binner->add($base_coverage);
            $total += $base_coverage;
            $count++;        
        }
        my $average = $total/$count;
        my $length = $lengthhash->{$id}->{length};

        print BASESTATS "$id\t$chromosome\t$length\t$start..$stop\t$average\t";
        my $databins = $binner->get_data_bins();
        foreach my $bindata ( @$databins )
        {
            print BASESTATS "$bindata\t";
        }
        print BASESTATS "\n";

    }
    close(ALIGNINFO);

    #### CLOSE AND REPORT BASE COVERAGE STATS FILE
    close(BASESTATS);

    return "$mapdir/$basestatsfile";
}


=head2

    SUBROUTINE      generate_lengthfile

    PURPOSE

        ###    GENERATE REFERENCE.length FILE
        ###    
        ###        <reference_file>.length
        ###            CCDS2.2 chr1    2046

=cut

sub generate_lengthfile
{
    my $self        =   shift;

    my $reference = $self->{_reference};
    return if not defined $reference;

    my $lengthfile = $reference . ".length";
    open(LENGTHFILE, ">$lengthfile") or die "Can't open reference length file: $lengthfile\n";

    #### OPEN REFERENCE FILE AND CALCULATE lengthhash
    open(REFERENCE, $reference) or die "Can't open reference file: $reference\n";
    $/ = ">";
    <REFERENCE>;
    while ( <REFERENCE> )
    {
        chop($_);
        my $fasta = $_;
        my ($header) = $_ =~ /^([^\n]+)/ms;
        $_ =~ s/^([^\n]+)\n//ms;
        $_ =~ s/\s+//g;
        my $length = length($_);
        my ($id, $chromosome) = $header =~ /^([^\|]+)\|[^\|]+\|(\S+)/ms;

#exit;

        if ( not defined $id or not defined $chromosome or not defined $_ )
        {
            print "Problem getting id, chromosome or sequence of fasta: $fasta\n";
            exit;
        }
        print LENGTHFILE "$id\t$chromosome\t$length\n";
    }

    close(LENGTHFILE);

    return $lengthfile;
}


=head2

    SUBROUTINE      generate_lengthfile

    PURPOSE

        ####     LOAD REFERENCE.length FILE INTO MEMORY
        ####    
        ####        <reference_file>.length
        ####            CCDS2.2 chr1    2046

=cut



sub load_lengthfile
{
    my $self        =   shift;
    my $reference_lengthfile      =   shift;

    $/ = "\n";
    open(LENGTHFILE, $reference_lengthfile) or die "Can't open reference length file: $reference_lengthfile\n";
    my $lengthhash;
    while ( <LENGTHFILE> )
    {
        $_ =~ /^(\S+)\s+(\S+)\s+(\S+)/;
        $lengthhash->{$1}->{chromosome} = $2;
        $lengthhash->{$1}->{length} = $3;
    }

    $self->{_lengthhash} = $lengthhash;
}


=head2

    SUBROUTINE      base_summary

    PURPOSE

        ####    USE BASE STATS FILE TO EXTRACT:
        ####
        ####    - How much sequence below certain cut-off (e.g. 5x)
        ####    - Average BASE coverage per chromosome
        ####
        ####    - Graph with average coverage per each gene
        ####    
        ####    TO GENERATE
        ####    
        ####        stats-base-summary.tsv

=cut

sub base_summary
{
    my $self        =   shift;

    my $mapdir = $self->{_mapdir};    
    my $basestatsfile = $self->{_basestatsfile};    

    #### OPEN BASE STATS FILE
    open(BASESTATS, "$mapdir/$basestatsfile") or die "Can't open base coverage file: $mapdir/$basestatsfile\n";    

    #### SET BASE SUMMARY FILE AND OPEN
    my $base_summaryfile = $self->{_basesummaryfile};
    open(BASESUMMARY, ">$mapdir/$base_summaryfile") or die "Can't open base coverage file: $mapdir/$base_summaryfile\n";

        #### FORMAT:
    ####    CCDS2.2 chr1    2046    1..72   3       0       0       72      0       0       0       0       0       0       0
    ####    CCDS2.2 chr1    2046    158..305        1.04054054054054        0       142     6       0       0       0       0       0       0       0
    ####    CCDS2.2 chr1    2046    518..767        1.41832669322709        0       148     103     0       0       0       0       0       0       0
    ####    CCDS2.2 chr1    2046    1796..2046      3.25396825396825        0       0       147     105     0       0       0       0       0       0
    ####    CCDS3.1 chr1    2250    92..420 8.0326409495549 0       14      28      109     84      102     0       0       0       0
    ####    CCDS3.1 chr1    2250    455..2166       8.07884288145207        0       90      143     274     675     567     14      0       0       0
    ####    CCDS4.1 chr1    1836    81..603 2.23282442748092        0       253     136     126     9       0       0       0       0       0
    ####    CCDS4.1 chr1    1836    786..1534       4.62714097496706        0       90      163     218     286     2       0       0       0       0
    ####    CCDS4.1 chr1    1836    1684..1836      2       0       0       160     0       0       0       0       0       0       0
    ####    CCDS6.1 chr1    498     3..498  6.0281124497992 0       0       157     7       334     0       0       0       0       0
    ####    CCDS7.1 chr1    765     202..765        2.89575971731449        0       0       510     47      9       0       0       0       0       0
    ####    CCDS8.1 chr1    1215    288..409        1       0       123     0       0       0       0       0       0       0       0
    ####    CCDS8.1 chr1    1215    549..978        3.6183908045977 0       155     86      33      161     0       0       0       0       0
    ####    CCDS8.1 chr1    1215    1038..1215      5.76404494382022        0       0       0       89      89      0       0       0       0       0
    ####    CCDS11.1        chr1    834     1..88   1       0       89      0       0       0       0       0       0       0       0
    ####    CCDS11.1        chr1    834     264..370        1       0       110     0       0       0       0       0       0       0       0
    ####    CCDS11.1        chr1    834     475..635        4       0       0       0       165     0       0       0       0       0       0


    my $summary;
    my $previous_id;
    my $current_id;
    my $first_id = "true";
    #my $print_counter = 0;

    $/ = "\n";
    while ( <BASESTATS> )
    {
        #### STORE FIRST ID AS 'previous_id', THEN COMPARE WITH SUBSEQUENT
        #### LINES TO MAKE SURE WE ARE STILL READING THE SAME TARGET
        #### FORMAT:
        ####    CCDS2.2 chr1    2046    1..72   3       0       0       72      0       0       0       0       0       0       0
        ####    CCDS2.2 chr1    2046    158..305        1.04054054054054        0       142     6       0       0       0       0       0       0       0
        $_ =~ s/\s+$//;

        my @elements = split "\t", $_;
        my $id = shift @elements;
        my $chromosome = shift @elements;
        my $length = shift @elements;
        my $startstop = shift @elements;
        my $base_coverage = shift @elements;

        if ( not defined $id or not defined $chromosome or not defined $length or not defined $startstop or not defined $base_coverage )
        {
            print "Problem with line: $_\n";
            exit;
        }

        #### GET TARGET LENGTH
        my ($start, $stop) = $startstop =~ /^(\d+)\.\.(\d+)$/;
        my $targetlength = $stop - $start + 1;

        #### IF FIRST LINE, SET PREVIOUS ID = CURRENT ID
        if ( $first_id =~ /^true$/ )
        {
            $first_id = "false";
            ($previous_id) = $id;
            if ( not defined $previous_id )
            {
                print "ID not defined in line: $_\n";
                exit;
            }
            $current_id = $previous_id;
        }
        else
        {
            ($current_id) = $id;
        }

        #### CHECK IF WE ARE DOING A NEW GENE
        if ( $current_id !~ /^$previous_id$/ )
        {
            #### PRINT CUMULATIVE STATS TO FILE IF DEFINED
            if ( defined $summary )
            {
                print BASESUMMARY $summary->{id} , "\t";
                print BASESUMMARY $summary->{chromosome} , "\t";
                print BASESUMMARY $summary->{coverage} , "\t";
                print BASESUMMARY $summary->{length} , "\t";
                print BASESUMMARY $summary->{coveredlength} , "\t";
                my $pct_covered = sprintf "%.3f", $summary->{coveredlength} / $summary->{length};
                print BASESUMMARY $pct_covered , "\t";
                my $bins = join "\t", @{$summary->{bins}};
                print BASESUMMARY $bins, "\t";
                print BASESUMMARY "\n";

                #### SET CUMULATIVE SUMMARY TO UNDEF
                $summary = undef;

                #last if $print_counter == 2;
            }

            $previous_id = $current_id;
        }


        if ( not defined $summary->{id} )
        {
            #### ADD THIS ENTRY TO CUMULATIVE SUMMARY
            $summary->{id} = $id;
            $summary->{chromosome} = $chromosome;
            $summary->{length} = $length;
            $summary->{coveredlength} = 0;
            $summary->{coverage} = 0;
            $summary->{bins} = [];
        }

        #### ADD THIS ENTRY TO CUMULATIVE SUMMARY
        $summary->{coverage} = $summary->{coverage} * ( $summary->{coveredlength} / ($summary->{coveredlength} + $targetlength) )
                                + $base_coverage * ( $targetlength / ($summary->{coveredlength} + $targetlength) );

        my $bins = $summary->{bins};
        for ( my $i = 0; $i < $#elements + 1; $i++ )
        {
            last if $elements[$i] =~ /^\s*$/;

            my $bin = $$bins[$i];
            if ( not defined $bin )
            {
                $bin = 0;
            }
            $$bins[$i] = $bin * ( $summary->{coveredlength} / ($summary->{coveredlength} + $targetlength) )
                    + $elements[$i] * ( $targetlength / ($summary->{coveredlength} + $targetlength) );
        }
        $summary->{coveredlength} += $targetlength;


    }
    close(BASESTATS);

    #### CLOSE AND REPORT BASE COVERAGE STATS FILE
    close(BASESUMMARY);

    return "$mapdir/$base_summaryfile";
}



=head2

    SUBROUTINE      mappedfile

    PURPOSE

        ####    
        ####    - What did we miss of the targeted region?
        ####        
        ####    USE FILES
        ####    
        ####        stats-base-coverage.tsv
        ####        <reference_file>
        ####        
        ####    TO GENERATE .FASTA FILE
        ####    
        ####        stats-unmapped-regions.fasta
        ####
        ####
        ####- Homopolymericity of unmapped regions and mapped regions
        ####
        ####    1. count polymers as 1 unique nucleotide
        ####    
        ####    2. all unique nucleotides / all nucleotide
        ####
        ####

=cut

sub mappedfile
{
    my $self        =   shift;

    my $mapdir = $self->{_mapdir};    
    my $basestatsfile = $self->{_basestatsfile};    


    #### GET INDEX OF REFERENCE FILE FASTA RECORD POSITIONS
    if ( not defined $self->{_indexer} )
    {
        $self->index_reference();
        #$self->load_index();
    }
    if ( not defined $self->{_indexer} )
    {
        die "Can't index reference or load index file. Exiting...\n";
    }
    my $indexer = $self->{_indexer};


    ### OPEN UNMAPPED REGIONS OUTPUT FILE
    my $unmappedfile = "stats-unmapped.fasta";
    open(UNMAPPED, ">$mapdir/$unmappedfile") or die "Can't open unmapped regions file: $mapdir/$unmappedfile\n";

    ### OPEN MAPPED REGIONS OUTPUT FILE
    my $mappedfile = "stats-mapped.fasta";
    open(MAPPED, ">$mapdir/$mappedfile") or die "Can't open mapped regions file: $mapdir/$mappedfile\n";

    ### OPEN INPUTFILES
    $/ = "\n";
    open(BASESTATS, "$mapdir/$basestatsfile") or die "Can't open base coverage file: $mapdir/$basestatsfile\n";
    my $previous_id;
    my $current_id;
    my $mapped;
    my $first_id = "true";

    my $counter = 0;
    while ( <BASESTATS> )
    {
        next if $_ =~ /^\s*$/;

        #### FORMAT:
        ####    CCDS2.2 chr1    2046    1..72   3       0       0       72      0       0       0       0       0       0       0
        ####    CCDS2.2 chr1    2046    158..305        1.04054054054054        0       142     6       0       0       0       0       0       0       0
        ####    CCDS2.2 chr1    2046    518..767        1.41832669322709        0       148     103     0       0       0       0       0       0       0
        ####    CCDS2.2 chr1    2046    1796..2046      3.25396825396825        0       0       147     105     0       0       0       0       0       0


        my @elements = split "\t", $_;
        my $id = shift @elements;
        my $chromosome = shift @elements;
        my $length = shift @elements;
        my $startstop = shift @elements;
        #my $base_coverage = shift @elements;

        #### IF FIRST LINE, SET PREVIOUS ID = CURRENT ID
        if ( $first_id =~ /^true$/ )
        {
            $first_id = "false";
            ($previous_id) = $id;
            if ( not defined $previous_id )
            {
                print "ID not defined in line: $_\n";
                exit;
            }
            $current_id = $previous_id;
        }
        else
        {
            ($current_id) = $id;
        }

        #### CHECK IF WE ARE DOING A NEW GENE
        if ( $current_id !~ /^$previous_id$/ )
        {
            $counter++;    

            #### PRINT CUMULATIVE STATS TO FILE IF DEFINED
            if ( defined $mapped )
            {
                my $id = $mapped->{id};
                my $hits = $mapped->{hits};
                my $length = $mapped->{length};            
                my $chromosome = $mapped->{chromosome};

                my $record = $indexer->record($id);
                $record =~ s/^[^\n]+\n//;
                $record =~ s/\s+//g;

                foreach my $startstop ( @$hits )
                {
                    #### GET TARGET LENGTH
                    my ($start, $stop) = $startstop =~ /^(\d+)\.\.(\d+)$/;
                    my $targetlength = $stop - $start + 1;
                    my $target_sequence = substr($record, $start - 1, $targetlength);

                    my $header = ">$id $chromosome $start..$stop length=$length\n";
                    print MAPPED $header;
                    print MAPPED $target_sequence, "\n";
                }

                my $unmapped = $self->unmapped($hits, $length);
                foreach my $startstop ( @$unmapped )
                {
                    #### GET TARGET LENGTH
                    my ($start, $stop) = $startstop =~ /^(\d+)\.\.(\d+)$/;
                    my $targetlength = $stop - $start + 1;
                    my $target_sequence = substr($record, $start - 1, $targetlength);

                    my $header = ">$id $chromosome $start..$stop length=$length\n";
                    print UNMAPPED $header;
                    print UNMAPPED $target_sequence, "\n";
                }

                #### SET CUMULATIVE SUMMARY TO UNDEF
                $mapped = undef;            
            }
            #last if $counter == 2;

            $previous_id = $current_id;
        }


        #### INITIALISE mapped WITH THE FIRST LINE FOR THIS ID    
        if ( not defined $mapped ->{id} )
        {
            #### ADD THIS ENTRY TO CUMULATIVE SUMMARY
            $mapped->{id} = $id;
            $mapped->{chromosome} = $chromosome;
            $mapped->{hits} = [];
            $mapped->{length} = $length;
        }

        push @{$mapped->{hits}}, $startstop;


    }
    close(BASESTATS);

    #### CLOSE AND REPORT OUTPUT FILES
    close(MAPPED);
    close(UNMAPPED);

    return "$mapdir/$mappedfile", "$mapdir/$unmappedfile";    
}



=head2

    SUBROUTINE      index_reference

    PURPOSE

        BUILD AN INDEX OF THE REFERENCE FASTA FILE RECORDS

        AND LOAD THE INDEXED RECORDS INTO THE INDEXER OBJECT

=cut


sub index_reference
{
    my $self        =   shift;

    my $reference = $self->{_reference};

    #### GENERATE REFERENCE FILE INDEX
    my $reference_indexfile = "$reference.index";
    my $indexer = Indexer->new(
        {
            'INPUTFILE' => $reference,
            'ID_REGEX' => "^>*([^\\|]+)",
            'RECORD_TYPE' => 'FASTA'
        }
    );
    $indexer->build_index();

    $self->{_indexer} = $indexer;
}



=head2

    SUBROUTINE

    PURPOSE

        ####    EXTRACT homopolymericity DIFFERENCE FOR EACH READ FROM FILES:
        ####    
        ####        stats-unmapped-regions.fasta
        ####    
        ####        454AllContigs.fna
        ####        
        ####            >contig00001  CCDS2.2|Hs36.3|chr1, 158..305  length=147   numreads=2
        ####            ccagcgatggtgacagcgacgggagtggccccacctgtgggcggcggccaggcttgaagc
        ####            aggaggatggtccgcacatccgtatcatgaagagaaGGTccacacccactgggacgtgaa
        ####            catctctttccgagaggcgtcctgcag
        ####            >contig00002  CCDS2.2|Hs36.3|chr1, 518..767  length=247   numreads=4
        ####            agaggcgctgctgctgccgcgggagctggggcccagcatggccccggaggaccattaccg
        ####            ccggcttgtgtcagcactgagcgaggccagcacctttgaggaccctcagcgcctctacca
        ####            cctgggcctcccagCCACGGTGAGGACCCACCCTGGCATGAtccccctcatcaccTCCCc
        ####            ACCcAGATCTCCTGAGGGTCCGGCAGGAGGTGGCGGCTGCAGCTCTGAGGGGCCCCAGTG
        ####            GCCTGGA
        ####
        ####    GENERATE MAPPED HOMOPOLY AND UNMAPPED HOMOPOLY FILES
        ####    
        ####        stats-mapped-homopoly.tsv
        ####        stats-unmapped-homopoly.tsv
        ####        

=cut

sub homopolies
{
    my $self        =   shift;

    my @types = ( "mapped", "unmapped" );
    my @outputfiles;
    foreach my $type ( @types )
    {
        push @outputfiles, $self->homopoly($type);
    }

    return @outputfiles;
}


sub homopoly
{
    my $self        =   shift;
    my $type        =   shift;

    my $mapdir = $self->{_mapdir};    
    my $basestatsfile = $self->{_basestatsfile};    

    my $homopolyfile;
    my $fastafile;
    if ( $type =~ /^mapped$/i )
    {
        $homopolyfile = $self->{_mappedhomopolyfile};
        $fastafile = $self->{_mappedfile};
    }
    else
    {
        $homopolyfile = $self->{_unmappedhomopolyfile};
        $fastafile = $self->{_unmappedfile};
    }

    #my $homopolyfile = "stats-homopolymericity-mapped.tsv";
    open(HOMOPOLY, ">$mapdir/$homopolyfile") or die "Can't open homopolymericity file: $homopolyfile\n";

    #### DO MAPPED REGIONS OUTPUT FILE
    open(FASTA, "$mapdir/$fastafile") or die "Can't open mapped regions file: $mapdir/$fastafile\n";

    print "Doing type '$type' homopolymericity file...\n";
    my $average_mapped_homopolymericity = 0;
    my $average_g = 0;
    my $average_t = 0;
    my $average_a = 0;
    my $average_c = 0;
    my $total_mapped = 0;
    $/ = ">";
    <FASTA>;
    while ( <FASTA> )
    {
        next if $_ =~ /^\s*$/;
        $total_mapped++;

        chop($_);
        $_ =~ s/\s+$//;
        my ($header, $sequence) = $_ =~ /^([^\n]+)\n(.+)$/ms;
        my ($id, $chromosome, $startstop) = $header =~ /^(\S+)\s+(\S+)\s+(\S+)\s+/;
        $sequence =~ s/\s+//g;
        my $length = length($sequence);

        #my ($counter) = $sequence =~ s/^([A|[G|[T|[C)//g;
        my $counter = 0;
        my $g = 0;
        my $c = 0;
        my $a = 0;
        my $t = 0;
        while ($sequence =~ /G/ig) { $g++; }
        while ($sequence =~ /C/ig) { $c++; }
        while ($sequence =~ /A/ig) { $a++; }
        while ($sequence =~ /T/ig) { $t++; }
        $g = $g / $length;
        $c = $c / $length;
        $a = $a / $length;
        $t = $t / $length;
        $average_g += $g;
        $average_c += $c;
        $average_a += $a;
        $average_t += $t;

        while ( $sequence =~ s/^A+//g
               or $sequence =~ s/^G+//g
               or $sequence =~ s/^T+//g
               or $sequence =~ s/^C+//g) { $counter++ }
        my $homopolymericity = sprintf "%.3f", $counter / $length;
        $average_mapped_homopolymericity += $homopolymericity;    
        print HOMOPOLY "$id\t$chromosome\t$startstop\t$length\t$homopolymericity\n";
    }

    $average_g = sprintf "%.3f", $average_g / $total_mapped;
    $average_c = sprintf "%.3f", $average_c / $total_mapped;
    $average_a = sprintf "%.3f", $average_a / $total_mapped;
    $average_t = sprintf "%.3f", $average_t / $total_mapped;
    my $gc = $average_g + $average_c;
    my $at = $average_a + $average_t;
    $average_mapped_homopolymericity = sprintf "%.3f", $average_mapped_homopolymericity / $total_mapped;
    print HOMOPOLY "Average_mapped_homopolymericity: $average_mapped_homopolymericity\n";
    print HOMOPOLY "Average g: $average_g\n";
    print HOMOPOLY "Average c: $average_c\n";
    print HOMOPOLY "Average a: $average_a\n";
    print HOMOPOLY "Average t: $average_t\n";
    print HOMOPOLY "Average gc: $gc\n";
    print HOMOPOLY "Average at: $at\n";

    #### CLOSE AND REPORT OUTPUT FILES
    close(FASTA);
    close(HOMOPOLY);
    print "Type '$type' homopolymericity file printed:\n\n$mapdir/$homopolyfile\n\n";

    return "$mapdir/$homopolyfile";
}




#sub unmapped_homopoly
#{
#    my $self        =   shift;
#
#    my $mapdir = $self->{_mapdir};
#    my $unmappedfile = $self->{_unmappedfile};
#    my $unmapped_homopolyfile = $self->{_unmappedhomopolyfile};
#    
#    #my $unmapped_homopolyfile = "stats-homopolymericity-unmapped.tsv";
#    open(UNMAPPEDHOMOPOLY, ">$mapdir/$unmapped_homopolyfile") or die "Can't open homopolymericity file: $unmapped_homopolyfile\n";
#    
#    #### DO UNMAPPED REGIONS OUTPUT FILE
#    open(UNMAPPED, "$mapdir/$unmappedfile") or die "Can't open unmapped regions file: $mapdir/$unmappedfile\n";
#    
#    my $average_unmapped_homopolymericity = 0;
#    my $average_g = 0;
#    my $average_t = 0;
#    my $average_a = 0;
#    my $average_c = 0;
#    my $total_unmapped = 0;
#    $/ = ">";
#    <UNMAPPED>;
#    while ( <UNMAPPED> )
#    {
#        next if $_ =~ /^\s*$/;
#        $total_unmapped++;
#        
#        chop($_);
#        $_ =~ s/\s+$//;
#        my ($header, $sequence) = $_ =~ /^([^\n]+)\n(.+)$/ms;
#        my ($id, $chromosome, $startstop) = $header =~ /^(\S+)\s+(\S+)\s+(\S+)\s+/;
#        $sequence =~ s/\s+//g;
#        my $length = length($sequence);
#        #print "header: $header\n";
#        #print "sequence: $sequence\n";
#        
#        #my ($counter) = $sequence =~ s/^([A|[G|[T|[C)//g;
#        my $counter = 0;
#        my $g = 0;
#        my $c = 0;
#        my $a = 0;
#        my $t = 0;
#        while ($sequence =~ /G/ig) { $g++; }
#        while ($sequence =~ /C/ig) { $c++; }
#        while ($sequence =~ /A/ig) { $a++; }
#        while ($sequence =~ /T/ig) { $t++; }
#        $g = $g / $length;
#        $c = $c / $length;
#        $a = $a / $length;
#        $t = $t / $length;
#        $average_g += $g;
#        $average_c += $c;
#        $average_a += $a;
#        $average_t += $t;
#        
#        while ( $sequence =~ s/^A+//g
#               or $sequence =~ s/^G+//g
#               or $sequence =~ s/^T+//g
#               or $sequence =~ s/^C+//g) { $counter++ }
#        #print "counter: $counter\n";
#        my $homopolymericity = sprintf "%.3f", $counter / $length;
#        #print "homopolymericity: $homopolymericity\n";
#        $average_unmapped_homopolymericity += $homopolymericity;    
#    }
#    
#    $average_g = sprintf "%.3f", $average_g / $total_unmapped;
#    $average_c = sprintf "%.3f", $average_c / $total_unmapped;
#    $average_a = sprintf "%.3f", $average_a / $total_unmapped;
#    $average_t = sprintf "%.3f", $average_t / $total_unmapped;
#    my $gc = $average_g + $average_c;
#    my $at = $average_a + $average_t;
#    $average_unmapped_homopolymericity = sprintf "%.3f", $average_unmapped_homopolymericity / $total_unmapped;
#    
#    #### CLOSE AND REPORT OUTPUT FILES
#    close(UNMAPPED);
#    close(UNMAPPEDHOMOPOLY);
#
#}

=head2

    SUBROUTINE      unmapped

    PURPOSE

        GENERATE START/STOPS OF UNMAPPED REGIONS GIVEN TARGET HITS AND

        TOTAL LENGTH OF SEQUENCE

=cut

sub unmapped
{
    my $self        =   shift;
    my $hits        =   shift;
    my $length      =   shift;

    if ( not defined $hits or not @$hits )   {   return; }
    if ( not defined $length or not $length =~ /^\d+$/ )
    {
        return;
    }
    my $unmapped; 
    my ($firststart) = $$hits[0] =~ /^(\d+)/;
    if ( $firststart > 1 )
    {
        push @$unmapped, 1 . ".." . ($firststart - 1);
    }

    for ( my $i = 1; $i <@$hits; $i++ )
    {
        my ($previousstop) = $$hits[$i - 1] =~ /^\d+\.\.(\d+)$/;
        my ($currentstart) = $$hits[$i] =~ /^(\d+)\.\.\d+$/;

        push @$unmapped, ($previousstop + 1) . ".." . ($currentstart - 1);
    }

    my ($laststop) = $$hits[scalar(@$hits) - 1] =~ /^\d+\.\.(\d+)$/;
    if ( $laststop < $length )
    {
        push @$unmapped, ($laststop + 1) . ".." . $length;
    }

    return $unmapped;
}



=head2

	SUBROUTINE		value

	PURPOSE

		SET A PARAMETER OF THE self OBJECT TO A GIVEN value

    INPUT

        1. parameter TO BE SET

		2. value TO BE SET TO

    OUTPUT

        1. THE SET parameter INSIDE THE self OBJECT

=cut

sub value
{
    my $self		=	shift;
	my $parameter	=	shift;
	my $value		=	shift;

	$parameter = lc($parameter);

    if ( not defined $value)	{	return;	}
	$self->{"_$parameter"} = $value;
}

=head2

	SUBROUTINE		validate_arguments

	PURPOSE

		VALIDATE USER-INPUT ARGUMENTS BASED ON

		THE HARD-CODED LIST OF VALID ARGUMENTS

		IN THE data ARRAY
=cut

sub validate_arguments
{
	my $self		=	shift;
	my $arguments	=	shift;

	my $hash;
	foreach my $argument ( keys %$arguments )
	{
		if ( $self->is_valid($argument) )
		{
			$hash->{$argument} = $arguments->{$argument};
		}
		else
		{
			warn "'$argument' is not a known parameter\n";
		}
	}

	return $hash;
}


=head2

	SUBROUTINE		is_valid

	PURPOSE

		VERIFY THAT AN ARGUMENT IS AMONGST THE LIST OF

		ELEMENTS IN THE GLOBAL '$DATAHASH' HASH REF

=cut

sub is_valid
{
	my $self		=	shift;
	my $argument	=	shift;

	#### REMOVE LEADING UNDERLINE, IF PRESENT
	$argument =~ s/^_//;

	#### CHECK IF ARGUMENT FOUND IN '$DATAHASH'
	if ( exists $DATAHASH->{lc($argument)} )
	{
		return 1;
	}

	return 0;
}


=head2

	SUBROUTINE		AUTOLOAD

	PURPOSE

		AUTOMATICALLY DO 'set_' OR 'get_' FUNCTIONS IF THE

		SUBROUTINES ARE NOT DEFINED.

=cut

sub AUTOLOAD {
    my ($self, $newvalue) = @_;


    my ($operation, $attribute) = ($AUTOLOAD =~ /(get|set)(_\w+)$/);

    # Is this a legal method name?
    unless($operation && $attribute) {
        croak "Method name $AUTOLOAD is not in the recognized form (get|set)_attribute\n";
    }
    unless( exists $self->{$attribute} or $self->is_valid($attribute) )
	{
        #croak "No such attribute '$attribute' exists in the class ", ref($self);
		return;
    }

    # Turn off strict references to enable "magic" AUTOLOAD speedup
    no strict 'refs';

    # AUTOLOAD accessors
    if($operation eq 'get') {
        # define subroutine
        *{$AUTOLOAD} = sub { shift->{$attribute} };

    # AUTOLOAD mutators
    }elsif($operation eq 'set') {
        # define subroutine4

        *{$AUTOLOAD} = sub { shift->{$attribute} = shift; };

        # set the new attribute value
        $self->{$attribute} = $newvalue;
    }

    # Turn strict references back on
    use strict 'refs';

    # return the attribute value
    return $self->{$attribute};
}



# When an object is no longer being used, this will be automatically called
# and will adjust the count of existing objects
sub DESTROY {
    my($self) = @_;
}



1;




