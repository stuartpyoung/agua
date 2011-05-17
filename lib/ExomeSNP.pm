package ExomeSNP;

#### DEBUG

=head2

	PACKAGE		ExomeSNP

    VERSION:        0.01

    PURPOSE

        IDENTIFY SNPs THAT BELONG TO dbSNP AND DETERMINE

        WHETHER THEY ARE SYNONYMOUS OR NON-SYNONYMOUS

    INPUT

        1. MYSQL TABLES:

            ccdsGene    -   ALL CCDS GENES

                +--------------+------------------------------------+------+-----+---------+-------+
                | Field        | Type                               | Null | Key | Default | Extra |
                +--------------+------------------------------------+------+-----+---------+-------+
                | bin          | smallint(5) unsigned               | NO   |     | 0       |       | 
                | name         | varchar(255)                       | NO   | MUL |         |       | 
                | chrom        | varchar(255)                       | NO   | MUL |         |       | 
                | strand       | char(1)                            | NO   |     |         |       | 
                | txStart      | int(10) unsigned                   | NO   |     | 0       |       | 
                | txEnd        | int(10) unsigned                   | NO   |     | 0       |       | 
                | cdsStart     | int(10) unsigned                   | NO   |     | 0       |       | 
                | cdsEnd       | int(10) unsigned                   | NO   |     | 0       |       | 
                | exonCount    | int(10) unsigned                   | NO   |     | 0       |       | 
                | exonStarts   | longblob                           | NO   |     |         |       | 
                | exonEnds     | longblob                           | NO   |     |         |       | 
                | score        | int(11)                            | YES  |     | NULL    |       | 
                | name2        | varchar(255)                       | NO   | MUL |         |       | 
                | cdsStartStat | enum('none','unk','incmpl','cmpl') | NO   |     | none    |       | 
                | cdsEndStat   | enum('none','unk','incmpl','cmpl') | NO   |     | none    |       | 
                | exonFrames   | longblob                           | NO   |     |         |       | 
                +--------------+------------------------------------+------+-----+---------+-------+

                SELECT * FROM ccdsGene LIMIT 1\G

                        bin: 592
                        name: CCDS30551.1
                       chrom: chr1
                      strand: +
                     txStart: 945415
                       txEnd: 980224
                    cdsStart: 945415
                      cdsEnd: 980224
                   exonCount: 36
                  exonStarts: 945415,947443,960519,965907,966415,966720,967198,968481,968780,969065,969351,969576,970403,970601,970975,971206,971402,971639,972062,972569,972815,973018,973254,974109,974478,974808,975145,975475,975669,975968,976495,976695,976970,978995,979690,980066,
                    exonEnds: 945616,947705,960567,966123,966640,966945,967405,968700,968975,969266,969500,969682,970520,970766,971119,971331,971508,971978,972200,972697,972930,973138,973608,974302,974694,975038,975280,975572,975834,976080,976612,976888,977058,979220,979794,980224,
                       score: 0
                       name2: 
                 cdsStartStat: cmpl
                  cdsEndStat: cmpl
                  exonFrames: 0,0,1,1,1,1,1,1,1,1,1,0,1,1,1,1,0,1,1,1,0,1,1,1,2,2,1,1,2,2,0,0,1,2,2,1,

                select chrom, min(txStart), max(txStart) from ccdsGene group by chrom;
                +-------+--------------+--------------+
                | chrom | min(txStart) | max(txStart) |
                +-------+--------------+--------------+
                | chr1  |        58953 |    247178159 | 
                | chr10 |        82996 |    135217783 | 
                | chr11 |       183099 |    133757041 | 
                | chr12 |       117699 |    132231151 | 
                | chr13 |     18646002 |    114107419 | 
                | chr14 |     18447593 |    105066216 | 
                | chr15 |     18848597 |    100279867 | 
                | chr16 |        37429 |     88652259 | 
                | chr17 |        63642 |     78630980 | 
                | chr18 |       148698 |     76018644 | 
                | chr19 |        61678 |     63765250 | 
                | chr2  |        31607 |    242460581 | 
                | chr20 |        16350 |     62361762 | 
                | chr21 |      9928775 |     46881291 | 
                | chr22 |     14828823 |     49584070 | 
                | chr3  |       336459 |    199171489 | 
                | chr4  |       321772 |    191183709 | 
                | chr5  |       193422 |    180726893 | 
                | chr6  |       237539 |    170728571 | 
                | chr7  |       506606 |    158516067 | 
                | chr8  |       106085 |    146248681 | 
                | chr9  |       106799 |    139730906 | 
                | chrX  |       140854 |    154880618 | 
                | chrY  |       140854 |     57739818 | 
                +-------+--------------+--------------+


            ccdsSNP     -   INTERSECTION OF dbSNP 129 AND CCDS GENES

                +-----------------+------------------+------+-----+---------+-------+
                | Field           | Type             | Null | Key | Default | Extra |
                +-----------------+------------------+------+-----+---------+-------+
                | chromosome      | varchar(255)     | YES  |     | NULL    |       | 
                | chromosomeStart | int(10) unsigned | YES  |     | NULL    |       | 
                | chromosomeEnd   | int(10) unsigned | YES  |     | NULL    |       | 
                | snp             | varchar(20)      | NO   | PRI |         |       | 
                | score           | int(10)          | YES  |     | NULL    |       | 
                | strand          | char(1)          | YES  |     | NULL    |       | 
                +-----------------+------------------+------+-----+---------+-------+

                select chromosome, min(chromosomeStart), max(chromosomeStart) from ccdsSNP group by chromosome;
               +------------+----------------------+----------------------+
               | chromosome | min(chromosomeStart) | max(chromosomeStart) |
               +------------+----------------------+----------------------+
               | chr1       |                58996 |            247178652 | 
               | chr10      |                83040 |            135222828 | 
               | chr11      |               183111 |            133759013 | 
               | chr12      |               118160 |            132243019 | 
               | chr13      |             18646014 |            114109500 | 
               | chr14      |             18447613 |            105067111 | 
               | chr15      |             18848691 |            100280447 | 
               | chr16      |                37540 |             88652371 | 
               | chr17      |                97057 |             78636327 | 
               | chr18      |               156818 |             76019389 | 
               | chr19      |               232404 |             63774450 | 
               | chr2       |                35894 |            242463731 | 
               | chr20      |                24770 |             62367108 | 
               | chr21      |              9928785 |             46907829 | 
               | chr22      |             14829713 |             49584412 | 
               | chr3       |               336507 |            199235925 | 
               | chr4       |               327612 |            191185332 | 
               | chr5       |               193531 |            180620045 | 
               | chr6       |               256937 |            170735424 | 
               | chr7       |               506655 |            158627934 | 
               | chr8       |               106245 |            146250282 | 
               | chr9       |               106831 |            139848768 | 
               | chrX       |               140859 |            154893028 | 
               | chrY       |               539533 |             57752001 | 
               +------------+----------------------+----------------------+


            snp129      -   ALL dbSNP 129

                -----------------------------------------------+------+-----+---------+-------+
                | bin        | smallint(5) unsigned            | NO   |     | 0       |       | 
                | chrom      | varchar(31)                     | NO   | MUL |         |       | 
                | chromStart | int(10) unsigned                | NO   |     | 0       |       | 
                | chromEnd   | int(10) unsigned                | NO   |     | 0       |       | 
                | name       | varchar(15)                     | NO   | MUL |         |       | 
                | score      | smallint(5) unsigned            | NO   |     | 0       |       | 
                | strand     | enum('+','-')                   | YES  |     | NULL    |       | 
                | refNCBI    | blob| NO   |     |         |       | 
                | refUCSC    | blob| NO   |     |         |       | 
                | observed   | varchar(255)| NO   |     |         |       | 
                | molType    | enum('genomic','cDNA')| YES  |     | NULL    |       | 
                | class      | enum('unknown','single','in-del','het','microsatellite','named','mixed','mnp','insertion','deletion')| NO   |     | unknown |       | 
                | valid      | set('unknown','by-cluster','by-frequency','by-submitter','by-2hit-2allele','by-hapmap')| NO   |     | unknown |       | 
                | avHet      | float| NO   |     | 0       |       | 
                | avHetSE    | float| NO   |     | 0       |       | 
                | func       | set('unknown','coding-synon','intron','cds-reference','near-gene-3','near-gene-5','nonsense','missense','frameshift','untranslated-3','untranslated-5','splice-3','splice-5') | NO   |     | unknown |       | 
                | locType    | enum('range','exact','between','rangeInsertion','rangeSubstitution','rangeDeletion')| YES  |     | NULL    |       | 
                | weight     | int(10) unsigned | NO   |     | 0       |       | 
                -----------------------------------------------+------+-----+---------+-------+            

            codingExons -   ALL CCDS GENES IN THE NIMBLEGEN CAPTURE ARRAY

                +----------------------------------+-------------+------+-----+---------+-------+
                | Field                            | Type        | Null | Key | Default | Extra |
                +----------------------------------+-------------+------+-----+---------+-------+
                | CCDS_ID                          | varchar(20) | NO   | PRI |         |       | 
                | GENE_SYMBOL                      | varchar(20) | YES  |     | NULL    |       | 
                | DESCRIPTION                      | text        | YES  |     | NULL    |       | 
                | REFSEQ                           | varchar(20) | YES  |     | NULL    |       | 
                | UCSC_GENE_ID                     | varchar(20) | YES  |     | NULL    |       | 
                | ENSEMBL                          | varchar(20) | YES  |     | NULL    |       | 
                | CHROMOSOME                       | varchar(10) | YES  |     | NULL    |       | 
                | STRAND                           | varchar(1)  | YES  |     | NULL    |       | 
                | CDS_START                        | int(12)     | YES  |     | NULL    |       | 
                | CDS_END                          | int(12)     | YES  |     | NULL    |       | 
                | EXON_COUNT                       | int(6)      | YES  |     | NULL    |       | 
                | ARRAY_COVERAGE                   | varchar(6)  | YES  |     | NULL    |       | 
                | ARRAY_COVERAGE_W_100BP_EXTENSION | varchar(6)  | YES  |     | NULL    |       | 
                +----------------------------------+-------------+------+-----+---------+-------+

                select * from codingExons order by CDS_START limit 10;
                +-------------+-------------+---------------------------------------------------------------+--------------+--------------+-----------------+------------+--------+-----------+---------+------------+----------------+----------------------------------+
                | CCDS_ID     | GENE_SYMBOL | DESCRIPTION                                                   | REFSEQ       | UCSC_GENE_ID | ENSEMBL         | CHROMOSOME | STRAND | CDS_START | CDS_END | EXON_COUNT | ARRAY_COVERAGE | ARRAY_COVERAGE_W_100BP_EXTENSION |
                +-------------+-------------+---------------------------------------------------------------+--------------+--------------+-----------------+------------+--------+-----------+---------+------------+----------------+----------------------------------+
                | CCDS12989.2 | DEFB125     | "defensin, beta 125"                                          | NM_153325    | uc002wcw.1   | ENST00000382410 | chr20      | +      |     16350 |                               |       | 100%
                | CCDS42645.1 | FAM110C     | "family with sequence similarity 110, member C"               | NM_001077710 | uc002qvt.1   | ENST00000327669 | chr2       | -      |     31607 |                                 |     | 0%
                | CCDS10395.1 | POLR3K      | "polymerase (RNA) III (DNA directed) polypeptide K, 12.3 kDa" | NM_016310    | uc002cfi.1   | ENST00000293860 | chr16      | -      |     37429 |                               |       | 100%
                | CCDS10396.1 | C16orf33    | chromosome 16 open reading frame 33                           | NM_024571    | uc002cfj.2   | ENST00000293861 | chr16      | +      |     43989 |                               |       | 100%
                | CCDS32344.1 | RHBDF1      | rhomboid 5 homolog 1 (Drosophila)                             | NM_022450    | uc002cfl.2   | ENST00000262316 | chr16      | -      |     48338 |                               |       | 100%
                | CCDS30547.1 | OR4F5       | "olfactory receptor, family 4, subfamily F, member 5"         | NM_001005484 | uc001aal.1   | ENST00000326183 | chr1       | +      |     58953 |                               |       | 100%
                | CCDS32854.1 | OR4F17      | "olfactory receptor, family 4, subfamily F, member 17"        | NM_001005240 | uc002loc.1   | ENST00000318050 | chr19      | +      |     61678 |                               |       | 100%
                | CCDS10994.1 | RPH3AL      | rabphilin 3A-like (without C2 domains)                        | NM_006987    | uc002frd.1   | ENST00000323434 | chr17      | -      |     63642 |  1                            |       | 100%
                | CCDS32345.1 | MPG         | N-methylpurine-DNA glycosylase                                | NM_001015054 | uc002cfm.1   | ENST00000397817 | chr16      | +      |     68308 |                               |       | 100%
                | CCDS32346.1 | MPG         | N-methylpurine-DNA glycosylase                                | NM_002434    | uc002cfn.1   | ENST00000219431 | chr16      | +      |     69291 |                               |       | 100%
                +-------------+-------------+---------------------------------------------------------------+--------------+--------------+-----------------+------------+--------+-----------+---------+------------+----------------+----------------------------------+


        2. NIMBLEGEN SNP PIPELINE OUTPUT FILES, E.G.:

            >Reference      Start   End     Ref     Var     Total   Var
            >Accno           Pos    Pos     Nuc     Nuc     Depth   Freq
            >CCDS3.1|Hs36.3|chr1    778     783     CTGGTG  GTGCTAT 5       60%
            >CCDS3.1|Hs36.3|chr1    1016    1016    C       T       4       100%
            >CCDS3.1|Hs36.3|chr1    1078    1078    T       C       4       100%
            >CCDS3.1|Hs36.3|chr1    1084    1084    T       C       4       100%
            >CCDS3.1|Hs36.3|chr1    1102    1102    C       T       4       100%
            >CCDS3.1|Hs36.3|chr1    1166    1166    G       A       5       80%
            >CCDS3.1|Hs36.3|chr1    1182    1182    T       C       5       100%
            >CCDS3.1|Hs36.3|chr1    1359    1359    G       A       6       100%

        3. FILTER CRITERIA

                % FREQUENCY
                QUALITY
                HETEROZYGOTE/HOMOZYGOTE
                SYNONYMOUS/NON-SYNONYMOUS

    OUTPUT

        1. FILE CONTAINING A LIST OF SNPS THAT PASS FILTER CRITERIA

    USAGE

    ./captureSNP.pl  <--inputfile String> <--positionfile String> <--outputfile String> [-h]

        --inputfile            :   /full/path/to/input_SNP.txt file
        --positionfile         :   /full/path/to/chromosome_positions.txt file
        --outputfile           :   /full/path/to/output_SNP.txt file
        --help                 :   print help info

    EXAMPLES

ON SOLEXA

./captureSNPs.pl --inputfile /home/syoung/base/pipeline/nimblegen-gsmapper/P_2009_01_09_03_26_14_runMapping/mapping/snps/454HCDiffs-headers.txt --outputfile /home/syoung/base/pipeline/nimblegen-gsmapper/P_2009_01_09_03_26_14_runMapping/mapping/snps/snps.out --positionfile /home/syoung/base/pipeline/human-genome/chromosome_positions.txt

ON KRONOS

./captureSNPs.pl --inputfile /nethome/syoung/base/pipeline/nimblegen-gsmapper/P_2009_01_09_03_26_14_runMapping/mapping/snps/454HCDiffs-headers.txt --outputfile /nethome/syoung/base/pipeline/nimblegen-gsmapper/P_2009_01_09_03_26_14_runMapping/mapping/snps/snps.out --positionfile /nethome/syoung/base/pipeline/human-genome/chromosome_positions.txt


 select * from ccdsSNP limit 1\G
*************************** 1. row ***************************
     chromosome: chr1
chromosomeStart: 58996
  chromosomeEnd: 58997
            snp: rs1638318
          score: 0
         strand: +
1 row in set (0.01 sec)

rs1638318 is A/G

TEST FILE

emacs /nethome/syoung/base/pipeline/nimblegen-gsmapper/P_2009_01_09_03_26_14_runMapping/mapping/snps/test.txt

>Reference      Start   End     Ref     Var     Total   Var
>Accno           Pos    Pos     Nuc     Nuc     Depth   Freq
>CCDS3.1|Hs36.3|chr1    58996     598997     A  G 5       60%

./captureSNPs.pl --inputfile /nethome/syoung/base/pipeline/nimblegen-gsmapper/P_2009_01_09_03_26_14_runMapping/mapping/snps/test.txt --outputfile /nethome/syoung/base/pipeline/nimblegen-gsmapper/P_2009_01_09_03_26_14_runMapping/mapping/snps/snps.out --positionfile /nethome/syoung/base/pipeline/human-genome/chromosome_positions.txt

=cut 

use strict;
use warnings;
use Carp;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our $AUTOLOAD;

use strict;

#### EXTERNAL MODULES
use POSIX;
use Data::Dumper;

#### DEFAULT PARAMETERS
our @DATA = qw(
	DBOBJECT
);
our $DATAHASH;
foreach my $key ( @DATA )	{	$DATAHASH->{lc($key)} = 1;	}


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
    #$self->{_seqstatsfile} = $SEQSTATSFILE;

    #### VALIDATE USER-PROVIDED ARGUMENTS
	($arguments) = $self->validate_arguments($arguments);	

    #### PROCESS USER-PROVIDED ARGUMENTS
	foreach my $key ( keys %$arguments )
	{
		$self->value($key, $arguments->{$key});
	}    
}



=head2

    SUBROUTINE      dbSNP_entry

    PURPOSE

        1. FIND OUT IF THE SNP IS ALREADY ENTERED IN dbSNP

        2. IF SO, RETRIEVE RELATED INFO ABOUT THE dbSNP ENTRY

=cut

sub linehash
{
    my $self        =   shift;
    my $line        =   shift;

    return if $line =~ /^#/ or $line =~ /^\s*$/;

    my $hash;
    my @elements = split " ", $line;
    ($hash->{name}, $hash->{chromosome}) = $elements[0] =~ /^>([^\|]+).+?([^\|]+)$/;
    $hash->{ccdsstart} = $elements[1];
    $hash->{ccdsstop} = $elements[2];
    $hash->{referencenucleotide} = $elements[3];
    $hash->{variantnucleotide} = $elements[4];
    $hash->{depth} = $elements[5];
    ($hash->{variantfrequency}) = $elements[6] =~ /^(\S+)/;

    #### RESET LINEHASH
    $self->{_linehash} = undef;
    $self->{_linehash} = $hash;

    return $hash;
}

=head2

    SUBROUTINE      chromosome_startstop

    PURPOSE

        GET CHROMOSOME START/STOP OF SNP BASED ON THE SNP POSITION

        IN GENE, GENE START/STOP AND EXON START/STOPS

=cut

sub chromosome_startstop
{
    my $self        =   shift;

    my $dbobject = $self->{_dbobject};
    my $linehash = $self->{_linehash};


    #### GET SNP CCDS START
    my $snp_ccds_start = $linehash->{ccdsstart};

    #### GET txStart FROM ccdsGene TABLE
    my $query = "SELECT txStart, txEnd, strand, exonStarts, exonEnds FROM ccdsGene WHERE name='$linehash->{name}'";
    my $ccds = $dbobject->queryhash($query);
    return if ( not defined $ccds or not $ccds );
    my $transcription_start = $ccds->{txStart};
    my $transcription_stop = $ccds->{txEnd};

    my $strand = $ccds->{strand};

    #if ( $strand !~ /^(\-|\+)$/ )
    #{
    #    die "Strand not +/-: $strand for CCDS $linehash->{name}\n";
    #}    

    #my $snp_gene_start = $linehash->{};


    #### SET CHROMOSOME-RELATIVE START/STOP (REM: ZERO-INDEXED SO SUBTRACT 1)

    #### GET CHROMOSOME-SPECIFIC START SITE OF SNP WITHIN CCDS
    my $snp_chromosome_start = my $snp_chromosome_stop = $self->snp_startstop($ccds, $snp_ccds_start);

    $linehash->{chromosomestart} = $snp_chromosome_start;
    $linehash->{chromosomestop} = $snp_chromosome_stop;

    $self->{_linehash} = $linehash;

    return $snp_chromosome_start;
}



=head2

    SUBROUTINE      snp_startstop

    PURPOSE

        CALCULATE THE CHROMOSOME POSITION OF A SNP WITHIN A GENE

        GIVEN THE SNP'S GENE-SPECIFIC START/STOP AND THE GENE'S

        EXON START/STOPS

=cut

sub snp_startstop
{
    my $self    =   shift;
    my $ccds    =   shift;  #### ENTRY IN ccdsGene TABLE
    my $snp_ccds_start    =   shift;


    #### VERIFY SNP CCDS POSITION
    my @exon_starts = split /,/, $ccds->{exonStarts};
    my @exon_stops = split /,/, $ccds->{exonEnds};

    #### CHECK IF EQUAL NUMBERS OF STARTS AND STOPS
    if ( $#exon_starts != $#exon_stops )
    {
        die "Unequal number of exon_starts (", $#exon_starts + 1, ") and exon_stops (" , $#exon_stops + 1, ")\n";
    }

    #### GET LENGTHS OF EXONS
    my $lengths;
    for ( my $i = 0; $i < $#exon_starts + 1; $i++ )
    {
        #### MAKE SURE THE STOP COMES AFTER THE START
        if ( ($exon_stops[$i] - $exon_starts[$i]) < 0 )
        {
            die "exon_stop ($exon_stops[$i]) is before exon_start ($exon_starts[$i])\n"
        }

        push @$lengths, $exon_stops[$i] - $exon_starts[$i];
    }

    my $total_length = $self->sum($lengths);


    #### GET CHROMOSOME POSITION OF SNP
    my $strand = $ccds->{strand};

    #### DO POSITIVE STRAND
    if ( $strand =~ /^\+$/ )
    {
        #### 1. GO THROUGH THE EXONS UNTIL THE EXON'S CHROMOSOME START IS GREATER THAN
        ####    THE SNP'S POSITION ON THE CHROMOSOME

        my $cumulative_length = 0;
        my $exon_chromosome_start = $exon_starts[0];
        my $exon_counter = 0;
        while ( $cumulative_length <= $snp_ccds_start
               and $exon_counter < $#exon_starts + 1 )
        {

            #### ADD EXON LENGTH TO CUMULATIVE LENGTH
            $cumulative_length += $$lengths[$exon_counter];
            #### SET NEW EXON CHROMOSOME START
            $exon_chromosome_start = $exon_starts[$exon_counter];

            #### INCREMENT EXON COUNTER
            $exon_counter++;
        }        

        #### DECREMENT 1 FOR PREVIOUS EXON (ZERO INDEXED)
        $exon_counter--;

        #### 2. SUBTRACT THE PREVIOUS EXON'S CHROMOSOME POSITION FROM THE SNP'S CHROMOSOME
        ####    POSITION 
        ####    (BACKTRACK ONE TO GET THE PREVIOUS EXON)
        my $previous_exon_start = $exon_starts[$exon_counter];

        #### CALCULATE DIFFERENCE TO ADD TO START OF EXON
        #### TO GET CHROMOSOME POSITION OF SNP
        my $difference = $snp_ccds_start - ($cumulative_length - $$lengths[$exon_counter]);


        my $snp_chromosome_start = $exon_starts[$exon_counter] + $difference;

        return $snp_chromosome_start;
    }

    #### DO NEGATIVE STRAND
    elsif ( $strand =~ /^\-$/ )
    {
        #### SWAP THE EXON START AND STOP ARRAYS AND REVERSE THEM
        my @temp = @exon_starts;
        @exon_starts = reverse @exon_stops;
        @exon_stops = reverse @temp;
        @$lengths = reverse @$lengths;

        #### 

        my $cumulative_length = 0;
        my $exon_counter = 0;
        while ( $cumulative_length <= $snp_ccds_start
               and $exon_counter < $#exon_starts + 1 )
        {

            #### ADD EXON LENGTH TO CUMULATIVE LENGTH
            $cumulative_length += $$lengths[$exon_counter];

            #### INCREMENT EXON COUNTER
            $exon_counter++;
        }

        ####    NB: DECREMENT exon_counter BY 1 FOR PREVIOUS EXON (ZERO INDEXED)
        $exon_counter--;

        $cumulative_length = $cumulative_length - $$lengths[$exon_counter];


        #### 2. ADD THE CUMULATIVE LENGTH OF THE PRECEDING EXONS PLUS THE
        ####    DIFFERENCE BETWEEN THE EXON STOP AND THE SNP CHROMOSOME POSITION
        my $difference = $snp_ccds_start - $cumulative_length;

        my $snp_chromosome_start = $exon_stops[$exon_counter - 1] - $difference;

        return $snp_chromosome_start;
    }
}



=head2

    SUBROUTINE      is_SNP

    PURPOSE


=cut

sub is_SNP
{
    my $self        =   shift;


    my $linehash = $self->{_linehash};

    #### ONLY ONE BASE BETWEEN START AND STOP
    if ( abs($linehash->{ccdsstart} - $linehash->{ccdsstop}) > 1 )
    {
        return 0;
    }

    #### LENGTH OF REFERENCE AND VARIANT MUST BE 1
    return 0 if length($linehash->{referencenucleotide}) > 1;
    return 0 if length($linehash->{variantnucleotide}) > 1;

    #### NEITHER NUCLEOTIDE IS '-'
    return 0 if ($linehash->{referencenucleotide}) =~ /^\-$/;
    return 0 if ($linehash->{variantnucleotide}) =~ /^\-$/;

    return 1;
}


=head2

    SUBROUTINE      dbSNP_entry

    PURPOSE

        1. FIND OUT IF THE SNP IS ALREADY ENTERED IN dbSNP

        2. IF SO, RETRIEVE RELATED INFO ABOUT THE dbSNP ENTRY

=cut

sub dbsnp
{
    my $self        =   shift;


    my $dbobject = $self->{_dbobject};
    return if not defined $dbobject;

    my $linehash = $self->get_linehash();

    #### GET CHROMOSOME START
    my $chromosome_start = $self->chromosome_startstop();
    return if ( not defined $linehash->{chromosomestart} );

    #### GET CHROMOSOME AND CHANGE chrXY TO chrY
    my $chromosome = $linehash->{chromosome};

    #### chrXY Used by illumina for diploid snps on XY
    #### http://www.obiba.org/genobyte/apidocs/org/obiba/genobyte/model/Chromosome.html
    #### SET TO chrY TO MATCH WITH ENTRIES IN ccdsGene AND ccdsSNP TABLES (chrX, chrY)
    if ( $chromosome =~ /^chrXY$/ )
    {
        $chromosome = "chrY";
    }


    #### TEST: SET MARGIN WITHIN WHICH TO FIND dbSNP AROUND OUR PREDICTED SNP
    #### MARGIN  = 0: dbSNP MUST FALL ON SAME BASE AS PREDICTED SNP
    my $margin = 5;
    my $upper = $chromosome_start + $margin;
    my $lower = $chromosome_start - $margin;

    #### CHECK IF SNP BELONGS IN dbSNP
    my $query = qq{SELECT * FROM ccdsSNP
    WHERE chromosomeStart <= $upper
    AND chromosomeStart >= $lower
    AND chromosome = '$chromosome'};

    my $snps = $dbobject->queryhash($query);
    if ( defined $snps)
    {
    }
    else
    {
    }

    return $snps;    
}



=head2

    SUBROUTINE      sum

    PURPOSE

    1. GET FRAME OF SNP AND HENCE WHETHER sum OR NON-sum

    2. RETURN

=cut

sub sum
{
    my $self            =   shift;
    my $array           =   shift;

    my $sum = 0;
    foreach my $value ( @$array )
    {
        $sum += $value;
    }

    return $sum;
}



=head2

    SUBROUTINE      synonymous

    PURPOSE

    1. GET FRAME OF SNP AND HENCE WHETHER SYNONYMOUS OR NON-SYNONYMOUS

    2. RETURN

=cut

sub synonymous
{
    my $self            =   shift;    

    my $dbobject = $self->{_dbobject};
    return if not defined $dbobject;

    my $linehash        =   $self->get_linehash();
    return if not defined $linehash;


    my $chromosome              =   $linehash->{chromosome};
    my $snp_chromosome_start    =   $linehash->{chromosomestart};
    my $variant_nucleotide      =   $linehash->{variantnucleotide};
    my $snp_ccds_start          =   $linehash->{ccdsstart};
    my $name                    =   $linehash->{name};

    #### RETURN IF NO SNP CHROMOSOME START
    return if not defined $snp_chromosome_start;

    #### ADD 1 TO CCDS START
    #$snp_ccds_start++;


    #### CONFIRM IDENTITY OF SPANNING CCDS
    #### GET ALL CCDS THAT SPAN THIS SNP
    my $query = qq{SELECT * FROM ccdsGene
    WHERE txStart <= $snp_chromosome_start
    AND txEnd >= $snp_chromosome_start
    AND chrom = '$chromosome'
    AND name = '$name'};

    #### GET RESULT
    my $ccds = $dbobject->queryhash($query);

    #### RETURN IF NOT DEFINED
    return if not defined $ccds;


    #### GET CCDS SEQUENCE
    $query = qq{SELECT sequence FROM ccdsSeq WHERE id='$name'};
    my $sequence = $dbobject->query($query);

    #### GET CCDS STRAND
    $query = qq{SELECT strand FROM ccdsGene WHERE name='$name'};
    my $strand = $dbobject->query($query);
    if ( $strand !~ /^(\-|\+)$/ )
    {
        die "Strand not +/-: $strand for CCDS $name\n";
    }    

    #### GET REFERENCE BASE
    my $reference_nucleotide = substr($sequence, $snp_ccds_start, 1);

    #### SET SNP FRAME    
    my $snp_frame = $snp_ccds_start % 3;
    if ( $snp_frame == 0 )  {   $snp_frame = 3; }

    #### SET CODON CCDS START
    my $codon_ccds_start = $snp_ccds_start - $snp_frame;

    #### GET REFERENCE AND VARIENT CODONS
    my $variant_codon = my $reference_codon = substr($sequence, $codon_ccds_start, 3);
    substr($variant_codon, $snp_frame - 1, 1, $variant_nucleotide);

    my $reference_aa = $self->codon2aa($reference_codon, "long");

    my $variant_aa = $self->codon2aa($variant_codon, "long");
    if ( not defined $variant_aa )
    {
        print "NO VARIANT AA FOR variant codon: $variant_codon\n";
        return;
    }


    my $sense = "synonymous";
    if ( $reference_aa ne $variant_aa )
    {
        $sense = "missense";
    }
    $linehash->{sense} = $sense;
    $linehash->{referencecodon} = $reference_codon;
    $linehash->{variantcodon} = $variant_codon;
    $linehash->{referenceaa} = $reference_aa;
    $linehash->{variantaa} = $variant_aa;
    $linehash->{strand} = $strand;

    return $linehash;
}




=head2

	SUBROUTINE		codon2aa

	PURPOSE

		TRANSLATE FROM A CODON TRIPLET TO AN AMINO ACID

    INPUT

        1. THREE-NUCLEOTIDE CODON

        2. OUTPUT TYPE (threeletter, oneletter, long)

    OUTPUT

        1. AMINO ACID NAME, 3-LETTER OR 1-LETTER SYMBOL

=cut

sub codon2aa
{
    my $self		=	shift;
	my $codon   	=	shift;
    my $type        =   shift;

    return if not defined $codon;
    return '' if not $codon;
    return if not $type =~ /^(threeletter|oneletter|long)$/;

    if ( not defined $self->{_codons} )
    {
        my $hash;
        while ( <DATA> )
        {
            $_ =~ /^(\S+)\s+(\S+)\s+(\S+)\s+(\w+)/;
            $hash->{$1}->{threeletter} = $2;
            $hash->{$1}->{oneletter} = $3;
            $hash->{$1}->{long} = $4;
        }
        $self->{_codons} = $hash;
    }
    my $codons = $self->{_codons};

    return $codons->{uc($codon)}->{$type};
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

__DATA__
TTT Phe  F Phenylalanine
TTC Phe F Phenylalanine
TTA Leu L Leucine
TTG Leu L Leucine
TCT Ser S Serine
TCC Ser S Serine
TCA Ser S Serine
TCG Ser S Serine
TAT Tyr Y Tyrosine
TAC Tyr Y Tyrosine
TAA Ochre X Stop
TAG Amber X Stop
TGT Cys C Cysteine
TGC Cys C Cysteine
TGA Opal X Stop
TGG Trp W Tryptophan
CTT Leu L Leucine
CTC Leu L Leucine
CTA Leu L Leucine
CTG Leu L Leucine
CCT Pro P Proline
CCC Pro P Proline
CCA Pro P Proline
CCG Pro P Proline
CAT His H Histidine
CAC His H Histidine
CAA Gln Q Glutamine
CAG Gln Q Glutamine
CGT Arg R Arginine
CGC Arg R Arginine
CGA Arg R Arginine
CGG Arg R Arginine
ATT Ile I Isoleucine
ATC Ile I Isoleucine
ATA Ile I Isoleucine
ATG Met M Methionine, Start
ACT Thr T Threonine
ACC Thr T Threonine
ACA Thr T Threonine
ACG Thr T Threonine
AAT Asn N Asparagine
AAC Asn N Asparagine
AAA Lys K Lysine
AAG Lys K Lysine
AGT Ser S Serine
AGC Ser S Serine
AGA Arg R Arginine
AGG Arg R Arginine
GTT Val V Valine
GTC Val V Valine
GTA Val V Valine
GTG Val V Valine
GCT Ala A Alanine
GCC Ala A Alanine
GCA Ala A Alanine
GCG Ala A Alanine
GAT Asp D Aspartic acid
GAC Asp D Aspartic acid
GAA Glu E Glutamic acid
GAG Glu E Glutamic acid
GGT Gly G Glycine
GGC Gly G Glycine
GGA Gly G Glycine
GGG Gly G Glycine


