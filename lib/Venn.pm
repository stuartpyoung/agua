package Venn;


=head2

	PACKAGE		Venn

	PURPOSE

		PERFORM BASIC FILE FILTERING TASKS, SUCH AS:

			1. match: FIND MATCHING AND NON-MATCHING TRANSCRIPTS IN TWO

				.GTF FILES

			2. filterTranscripts: SELECT ONLY TRANSCRIPT FEATURES FROM A .GTF FILE

=cut

use strict;
use warnings;
use Carp;

require Exporter;
our @ISA = qw(Exporter Cluster);
our @EXPORT_OK = qw();
our $AUTOLOAD;

#### INTERNAL MODULES
use Cluster;
use Util;
use Filter;

#### EXTERNAL MODULES
use File::Remove;
use File::Copy;
use File::Path;

{
	use Data::Dumper;	
}

#### FLUSH BUFFER
#$| = 1;

#### SET SLOTS
our @DATA = qw(
	TYPE
	MODE
	JSON
	ARGS
);
our $DATAHASH;
foreach my $key ( @DATA )	{	$DATAHASH->{lc($key)} = 1;	}


=head2

	SUBROUTINE		collate

	PURPOSE

		COLLATE THE RESULTS OF A SERIES OF VENN COMPARISON

		IN TERMS OF LINE COUNTS PER OUTPUT FILE ('*-AND-*'

		AND '*-NOT-*' FILES) AND PRINT TO THE TABLE TO A

		FILE, E.G.:

		(COLUMN HEADINGS ADDED AFTERWARDS)

		Sample ID	Sample 	Different to	Shared with 	Sample 33  		
		(Billion 	SNPs 	Sample 33 (A)	Sample 33 (B)	ONLY (C)	A + B	B + C
		Reads)									
		0.1			24284	20720			3563	1		6904		24284	20467
		0.3			18250	14918			3331			17136		18250	20467
		0.4			21206	16801			4404			16063		21206	20467
		0.5			21206	16801			4404			16063		21206	20467
		0.6			22896	17261			5634			14833		22896	20467
		...
		3.2			20619	843				19775			692			20619	20467
		3.3			20468	0				20467			0			20468	20467


=cut

sub collate {
	my $self		=	shift;

	return $self->collateBins() if defined $self->get_binlevel(); 

	return $self->collateSimple();
}

sub collateSimple {
	my $self		=	shift;
	print "Venn::collateSimple    Venn::collateSimple()\n";

	my $outputfile	=	$self->get_outputfile();
	my $replicates	= 	$self->get_replicates(); 
	my $querylabel	=	$self->get_querylabel();
	my $targetlabel	=	$self->get_targetlabel();
	my $suffix		=	$self->get_suffix();

	my ($outputdir) = $outputfile =~ /^(.+?)\/([^\/]+)$/;
	print "Venn::collateSimple    outputdir: $outputdir\n";
	print "Venn::collateSimple    replicates: $replicates\n";
	print "Venn::collateSimple    querylabel: $querylabel\n";
	print "Venn::collateSimple    targetlabel: $targetlabel\n";

	my $repliques = $self->getReps($replicates);	
	print "Venn::collateSimple    repliques: @$repliques\n";

	my $outputs = [];
	my $querylabels = [];
	my $targetlabels = [];
	foreach my $replique ( @$repliques )
	{
		my $target = $targetlabel;
		$target =~ s/%REPLICATE%/$replique/ig;
		my $query = $querylabel;
		$query =~ s/%REPLICATE%/$replique/ig;

		my ($queryonly, $targetonly, $both) = $self->setOutputFiles($outputdir, $query, $target, $suffix);
		my $hash;
		$hash->{queryonly}->{filename} = $queryonly;
		$hash->{targetonly}->{filename} = $targetonly;
		$hash->{both}->{filename} = $both;
		push @$outputs, $hash;
	}

	foreach my $output ( @$outputs )
	{
		$output->{targetonly}->{lines} = $self->lineCount($output->{targetonly}->{filename});
		$output->{queryonly}->{lines} = $self->lineCount($output->{queryonly}->{filename});
		$output->{both}->{lines} = $self->lineCount($output->{both}->{filename});
	}
#exit;

	open(OUT, ">$outputfile") or die "Can't open outputfile: $outputfile\n";
	print OUT "$querylabel\t$targetlabel\t$querylabel ONLY (A)\tshared(B)\t$targetlabel ONLY(C)\tA + B\tB + C\n";
	foreach my $output ( @$outputs )
	{
		####	FORMAT:
		####	
		####	Sample ID	Sample 	Different to	Shared with 	Sample 33  		
		####	(Billion 	SNPs 	Sample 33 (A)	Sample 33 (B)	ONLY (C)	A + B	B + C
		####	Reads)									
		####	0.1			24284	20720			3563	1		6904		24284	20467
		####	0.3			18250	14918			3331			17136		18250	20467
		####	0.4			21206	16801			4404			16063		21206	20467
		####	0.5			21206	16801			4404			16063		21206	20467
		####	0.6			22896	17261			5634			14833		22896	20467

		my ($queryfile) = $output->{queryonly}->{filename} =~ /^.+?\/([^\/]+)-NOT-.+$/;
		my ($targetfile) = $output->{targetonly}->{filename} =~ /^.+?\/([^\/]+)-NOT-.+$/;
		my $a = $output->{targetonly}->{lines};
		my $b = $output->{both}->{lines};
		my $c = $output->{queryonly}->{lines};
		my $a_plus_b = $a + $b;
		my $b_plus_c = $b + $c;
		my $line = "$queryfile\t$targetfile\t$a\t$b\t$c\t$a_plus_b\t$b_plus_c";
		print OUT "$line\n";
	}
	close(OUT);
	print "Venn::collateSimple    outputfile printed:\n\n$outputfile\n\n";
}

=head2

	SUBROUTINE		collateBins

	PURPOSE

		COLLATE THE RESULTS OF A SERIES OF VENN COMPARISON SPLIT INTO

		BINS BY AGLOMERATING ALL OF THE TOTALS FOR ALL OF THE BINS 

=cut

sub collateBins {
	my $self		=	shift;

	print "Venn::collateBins    Venn::collateBins()\n";

	my $outputfile	=	$self->get_outputfile();
	my $replicates	= 	$self->get_replicates(); 
	my $querylabel	=	$self->get_querylabel();
	my $targetlabel	=	$self->get_targetlabel();
	my $querydir	=	$self->get_querydir();
	my $targetdir	=	$self->get_targetdir();
	my $queryindex	=	$self->get_queryindex();
	my $targetindex	=	$self->get_targetindex();
	my $suffix		=	$self->get_suffix();
	my $bamfile 	=	$self->get_bamfile();
	my $binlevel	=	$self->get_binlevel();
	my $samtools	=	$self->get_samtools();

	my ($outputdir) = $outputfile =~ /^(.+?)\/([^\/]+)$/;
	print "Venn::collateBins    outputdir: $outputdir\n";
	print "Venn::collateBins    replicates: $replicates\n";
	print "Venn::collateBins    querylabel: $querylabel\n";
	print "Venn::collateBins    targetlabel: $targetlabel\n";
	print "Venn::collateBins    suffix: $suffix\n";

	#### GET REPLIQUES
	my $repliques = $self->getReps($replicates);	
	print "Venn::collateBins    repliques: @$repliques\n";

	#### GET ARRAY OF BINS
	my $binner = UCSCBin->new({	samtools	=> $samtools	});
	my $chromosome_size = $binner->getChromosomeSize($bamfile);	
	my $bins = $binner->getBins($binlevel, $chromosome_size);
	print "Venn::collateBins    No. bins: ", scalar(@$bins), "\n";

	my $outputs = [];
	my $querylabels = [];
	my $targetlabels = [];
	foreach my $replique ( @$repliques )
	{
		my $target = $targetlabel;
		$target =~ s/%REPLICATE%/$replique/ig;
		my $query = $querylabel;
		$query =~ s/%REPLICATE%/$replique/ig;

		my $targetnumber = $targetindex;
		$targetnumber =~ s/%REPLICATE%/$replique/ig;
		my $querynumber = $queryindex;
		$querynumber =~ s/%REPLICATE%/$replique/ig;

		my $queryonly_count 	=	0;
		my $targetonly_count 	=	0;
		my $both_count			=	0;

		foreach my $bin ( @$bins )
		{
			#### SET FILE NAMES
			my ($filename) = $bamfile =~ /([^\/]+)$/;
			my $binfile = $binner->setBinfile($filename, $outputdir, $binlevel, $bin);
			my ($filestub) = $binfile =~ /([^\/]+)$/;
			$filestub =~ s/\.bam//;
			my $prefix = $filestub . ".";		

			#### IN THE CASE WHERE ONE OF THE FILES IS MISSING, ADD
			#### THE LINE COUNT OF THE OTHER FILE TO ITS 'only' TOTAL
			my $queryfile = "$querydir/$filestub-$querynumber.$suffix";
			my $targetfile = "$targetdir/$filestub-$targetnumber.$suffix";

			if ( not -f $queryfile )
			{
				next if not -f $targetfile;

				$targetonly_count += $self->lineCount($targetfile) and next;
			}
			if ( not -f $targetfile )
			{
				next if not -f $queryfile;

				$queryonly_count += $self->lineCount($queryfile) and next;
			}


			my ($queryonly, $targetonly, $both) = $self->setOutputFiles($outputdir, $query, $target, $suffix, $prefix);
			$queryonly_count 	+= $self->lineCount($queryonly);
			$targetonly_count 	+= $self->lineCount($targetonly);
			$both_count 		+= $self->lineCount($both);

		}


		#### PUSH TOTALS FOR THIS REPLIQUE ONTO outputs
		my $hash;
		$hash->{querynumber} = $querynumber;
		$hash->{targetnumber} = $targetnumber;
		$hash->{queryonly} = $queryonly_count;
		$hash->{targetonly} = $targetonly_count;
		$hash->{both} = $both_count;
		push @$outputs, $hash;
	}

#exit;

	open(OUT, ">$outputfile") or die "Can't open outputfile: $outputfile\n";
	print OUT "query\ttarget\ttarget_only(A)\tshared(B)\tquery_only(C)\tA + B\tB + C\n";
	foreach my $output ( @$outputs )
	{
		####	FORMAT:
		####	
		####	Sample ID	Sample 	Different to	Shared with 	Sample 33  		
		####	(Billion 	SNPs 	Sample 33 (A)	Sample 33 (B)	ONLY (C)	A + B	B + C
		####	Reads)									
		####	0.1			24284	20720			3563	1		6904		24284	20467
		####	0.3			18250	14918			3331			17136		18250	20467
		####	0.4			21206	16801			4404			16063		21206	20467
		####	0.5			21206	16801			4404			16063		21206	20467
		####	0.6			22896	17261			5634			14833		22896	20467

		my ($querynumber) = $output->{querynumber};
		my ($targetnumber) = $output->{targetnumber};
		my $a = $output->{targetonly};
		my $b = $output->{both};
		my $c = $output->{queryonly};
		my $a_plus_b = $a + $b;
		my $b_plus_c = $b + $c;
		my $line = "$querynumber\t$targetnumber\t$a\t$b\t$c\t$a_plus_b\t$b_plus_c";
		#my $line = "$a\t$b\t$c\t$a_plus_b\t$b_plus_c";
		print "$line\n";
		print OUT "$line\n";
	}
	close(OUT);
	print "Venn::collateBins    outputfile printed:\n\n$outputfile\n\n";	
}


sub replicateLabels {
	my $self		=	shift;
	my $replicates 	=	shift;
	my $filestub	=	shift;

	#### RUN replicates TIMES WITH values VALUES
	my $labels = [];
	for ( my $counter = 0; $counter < @$replicates; $counter++ )
	{
		my $replicate = $$replicates[$counter];
		my $label = $filestub;
		$label =~ s/%REPLICATE%/$replicate/ig;
		push @$labels, $label;
	}

	return $labels;
}

=head2

	SUBROUTINE		lineCount

	PURPOSE

		RETURN THE NUMBER OF LINES IN A FILE

=cut

sub lineCount {
	my $self	=	shift;
	my $filename = shift;

	# SET LINE-END TO '\n'
	$/ = "\n";

	my $counter = 0;
	open(FILE, $filename) or return 0;
	while ( my $line = <FILE> )
	{
		$counter++;
	}
	close(FILE);

	return $counter;
}	



=head2

	SUBROUTINE		getReps

	PURPOSE

		RETURN replicates STRING PARSED INTO AN ARRAY

=cut

sub getReps {
	my $self		=	shift;
	my $replicates =	shift;

	my $replicates_array;
	if ( $replicates =~ /^(\d+)\-(\d+)$/ )
	{
		my $start = $1;
		my $end = $2;
		for my $time ( $start .. $end )
		{
			push @$replicates_array, $time;
		}
	}
	else
	{
		@$replicates_array = split ",", $replicates;
	}

	return $replicates_array;
}


=head2

	SUBROUTINE		gtfSummary	

	PURPOSE

		RETURN THE FOLLOWING SUMMARY OF A GTF-FORMAT LINE:

			chromosome <tab> start <tab> stop <tab> id

=cut

sub gtfSummary {
	my $self		=	shift;
	my $line		=	shift;


	my ($chromosome, $start, $stop, $id) =
		$line =~ /([^\t]+)\t[^\t]+\t[^\t]+\t([^\t]+)\t([^\t]+)\t[^\t]+\t[^\t]+\t[^\t]\tgene_id\s+"([^"]+)"/;


	return "$chromosome\t$start\t$stop\t$id";
}


=head2

	SUBROUTINE		getAnnotation

	PURPOSE

		1. EXTRACT ANNOTATION FOR TRANSCRIPTS IN GTF FILE


	INPUTS

		1. ARRAY OF THE LINES OF A GTF-FORMAT FILE

			CONTAINING TRANSCRIPTS WITH knownGene IDS

		2. ARRAY OF THE LINES OF A TSV-FORMAT FILE

			CONTAINING UCSC kgXref KNOWN GENES EXTERNAL

			REFERENCES TABLE DATA (TO BE READ INTO MEMORY)

	OUTPUTS

		1. AN ARRAY OF EXTENDED GTF-FORMAT LINES WITH

			ADDITIONAL XREF INFORMATION FIELDS 

	NOTES

		INPUT FILE CONTENTS:

		head /scratch/syoung/base/pipeline/bixby/run1/tophat/analysis1/chrY/transcripts-only.gtf		

		chrY    Cufflinks       transcript      6622    31738   1000    -       .       gene_id "uc009uyu.1"; transcript_id "uc009uyu.1"; RPKM "53.9340766766"; frac "1.000000"; conf_lo "0.000000"; conf_hi "0.000000"; cov "0.607097";
		chrY    Cufflinks       transcript      61650   133852  1000    -       .       gene_id "uc009uyv.2"; transcript_id "uc009uyv.2"; RPKM "36.8056246792"; frac "1.000000"; conf_lo "0.000000"; conf_hi "0.000000"; cov "0.414295";



		XREF TABLE FILE CONTENTS:

		head /nethome/bioinfo/data/annotation/ucsc/knownGene/kgXref-table

		#kgID	mRNA	spID	spDisplayID	geneSymbol	refseq	protAcc	description
		uc007afh.1	NM_008866	P97823	LYPA1_MOUSE	Lypla1	NM_008866	NP_032892	lysophospholipase 1
		uc007afg.1	AK050549	P97823-2	P97823-2	Lypla1	NM_008866	NP_032892	lysophospholipase 1
		uc007afi.1	NM_011541	P10711	TCEA1_MOUSE	Tcea1	NM_011541	NP_035671	transcription elongation factor A (SII) 1

=cut

sub getAnnotation {
	my $self		=	shift;
	my $lines		=	shift;
	my $annohash	=	shift;

	print "Venn::getAnnotation    Venn::getAnnotation(lines)\n";
	print "Venn::getAnnotation    lines: $lines\n";
	print "Venn::getAnnotation    annohash: $annohash\n";

	foreach my $line ( @$lines )
	{
		my ($id) = $line =~ /(\S+)$/;
		print "id: $id\n";
		my $annotation = $annohash->{$id};
		print "annotation: $annotation\n";
		$line .= "\t$annotation" if defined $annotation;
	}

	return $lines;	
}



=head2

	SUBROUTINE		filterTranscripts

	PURPOSE

		1. FILTER 'transcript' LINES FROM .GTF FILE

=cut

sub filterTranscripts {
	my $self		=	shift;
	my $args		=	shift;

	my $inputfile 	=	$args->{inputfile};
	my $outputfile 	=	$args->{outputfile};
	if ( not defined $outputfile )
	{
		$outputfile = $inputfile;
		$outputfile =~ s/^(.+)(\.gtf)/$1-only$2/g;
	}




	#### FILTER ALL LINES
	open(FILE, $inputfile) or die "Can't open inputfile: $inputfile\n";
	my $lines = [];
	while ( <FILE> )
	{
		next if $_ =~ /^\s*$/;

		my @elements = split "\t", $_;
		push @$lines, $_ if $elements[2] eq "transcript"
	}
	close(FILE) or die "Can't close inputfile: $inputfile\n";

	#### SORT FILE BY START POSITION
	@$lines = sort by_start @$lines;

	#### PRINT TO OUTPUT FILE 
	open(OUT, ">$outputfile") or die "Can't open outputfile: $outputfile\n";
	for my $line ( @$lines )	{	print OUT $line;	}
	close(OUT) or die "Can't close outputfile: $outputfile\n";


	return $outputfile;
}

#### SORT BY START POSITION
sub by_start {
	my ($aa) = $a =~ /^\S+\s+\S+\s+transcript\s+(\d+)/;
	my ($bb) = $b =~ /^\S+\s+\S+\s+transcript\s+(\d+)/;

	$aa <=> $bb;	
}

=head2

	SUBROUTINE

	PURPOSE

		SET DEFAULT OUTPUT FILES

=cut

sub setOutputFiles {
	my $self		=	shift;
	my $outputdir	=	shift;
	my $querylabel	=	shift;
	my $targetlabel	=	shift;
	my $suffix		=	shift;
	my $prefix	=	shift;
	$prefix = '' if not defined $prefix;

	#### SET OUTPUT FILES
	my $queryonly = "$outputdir/$prefix$querylabel-NOT-$targetlabel.$suffix";
	my $targetonly = "$outputdir/$prefix$targetlabel-NOT-$querylabel.$suffix";
	my $both = "$outputdir/$prefix$querylabel-AND-$targetlabel.$suffix";

	return ($queryonly, $targetonly, $both);
}

=head2

	SUBROUTINE		match

	PURPOSE

		1. LOAD QUERY FILE INTO MEMORY AND SCROLL THROUGH

			TARGET FILE LOOKING FOR MATCHING AND NON-MATCHING

			TRANSCRIPTS 

	INPUTS

		1. QUERY AND TARGET .gtf FORMAT FILES WITH 'transcript' IN 3RD FIELD

			OF ALL LINES INDICATING THAT THE LINE IS A TRANSCRIPT

	OUTPUTS

		1. RETURN AN ARRAY OF MATCHING TRANSCRIPT LINES

		2. RETURN AN ARRAY OF THE NON-MATCHING TRANSCRIPT LINES

=cut



sub match {
	my $self		=	shift;
	my $queryfile	=	shift;
	my $targetfile	=	shift;





	#### OPEN TARGET FILE
	open(TARGET, $targetfile) or die "Can't open target file: $targetfile\n";
	my @queries = <TARGET>;
	close(TARGET);

	my $matched = [];
	my $unmatched = [];

	#### OPEN QUERY FILE
	open(QUERY, $queryfile) or die "Can't open target file: $targetfile\n";
	while ( <QUERY> )
	{
		next if $_ =~ /^\s*$/ or /^#/;
		my ($queryposition) = $_ =~ /^\S+\s\S+\s+\S+\s+(\S+)/;
		my $match = 0;
		foreach my $target ( @queries )
		{
			my ($targetposition) = $target =~ /^\S+\s\S+\s+\S+\s+(\S+)/;

			if ( $targetposition == $queryposition )
			{
				push @$matched, $target;
				$match = 1;
				last;
			}
		}

		push @$unmatched, $_ if not $match;
	}
	close(QUERY);

	return $matched, $unmatched;
}







################################################################################
##################			HOUSEKEEPING SUBROUTINES			################
################################################################################

=head2

	SUBROUTINE		initialise

	PURPOSE

		INITIALISE THE self OBJECT WITH USER-INPUT ARGUMENT VALUES

=cut

sub initialise {
    my $self		=	shift;
	my $arguments	=	shift;

    #### VALIDATE USER-PROVIDED ARGUMENTS
	($arguments) = $self->validate_arguments($arguments, $DATAHASH);	

	#### LOAD THE USER-PROVIDED ARGUMENTS
	foreach my $key ( keys %$arguments )
	{		
		#### LOAD THE KEY-VALUE PAIR
		$self->value($key, $arguments->{$key});
	}
}


=head2

	SUBROUTINE		new

	PURPOSE

		CREATE THE NEW self OBJECT AND INITIALISE IT, FIRST WITH DEFAULT 

		ARGUMENTS, THEN WITH PROVIDED ARGUMENTS

        THE PASSWORD FOR THE DEFAULT ROOT USER CAN BE SET ON THE MYSQL

        COMMAND LINE:

        use myEST;
        insert into users values('admin', 'mypassword', now());

=cut

sub new {
 	my $class 		=	shift;
	my $arguments 	=	shift;

	my $self = {};
    bless $self, $class;

	#### INITIALISE THE OBJECT'S ELEMENTS
	$self->initialise($arguments);

    return $self;
}


=head2

	SUBROUTINE		value

	PURPOSE

		SET A PARAMETER OF THE self OBJECT TO A GIVEN value

    INPUT

        1. parameter TO BE SET

		2. value TO BE SET TO

    OUTPUT

        1. THE SET parameter INSIDE THE BioSVG OBJECT

=cut
sub value {
    my $self		=	shift;
	my $parameter	=	shift;
	my $value		=	shift;

	$parameter = lc($parameter);

    if ( defined $value)
	{	
		$self->{"_$parameter"} = $value;
	}
}

=head2

	SUBROUTINE		validate_arguments

	PURPOSE

		VALIDATE USER-INPUT ARGUMENTS BASED ON

		THE HARD-CODED LIST OF VALID ARGUMENTS

		IN THE data ARRAY
=cut

sub validate_arguments {
	my $self		=	shift;
	my $arguments	=	shift;
	my $DATAHASH	=	shift;

	my $hash;
	foreach my $argument ( keys %$arguments )
	{
		if ( $self->is_valid($argument, $DATAHASH) )
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

sub is_valid {
	my $self		=	shift;
	my $argument	=	shift;
	my $DATAHASH	=	shift;

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

    return if not defined $newvalue or not $newvalue;

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

	#if ( defined $self->{_databasehandle} )
	#{
	#	my $dbh =  $self->{_databasehandle};
	#	$dbh->disconnect();
	#}

#    my($self) = @_;
}

1;
