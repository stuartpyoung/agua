package Venn::Snp;


=head2

	PACKAGE		Venn::Snp

	PURPOSE

		PERFORM BASIC FILE FILTERING TASKS, SUCH AS:

			1. match: FIND MATCHING AND NON-MATCHING TRANSCRIPTS IN TWO

				.GTF FILES

			2. filterTranscripts: SELECT ONLY TRANSCRIPT FEATURES FROM A .GTF FILE

=cut

use strict;
use warnings;
use Carp;

#### USE LIB FOR INHERITANCE
use FindBin qw($Bin);
use lib "$Bin/..";

require Exporter;
our @ISA = qw(Venn);
use Venn;
our $AUTOLOAD;

#### INTERNAL MODULES
use Cluster;
use UCSCBin;
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

OUTPUTFILE

QUERYFILE
TARGETFILE
QUERYLABEL
TARGETLABEL
QUERYINDEX
TARGETINDEX
OUTPUTDIR

REPLICATES
SAMTOOLS
INPUTDIRS
OUTPUTDIR

QUERYDIR
TARGETDIR
QUERYLABEL
TARGETLABEL
WINDOW
INPUTFILE
INPUTFILES
OUTPUTFILE
FILENAME

BINLEVEL
BAMFILE
SUFFIX
QUERYINDEX
TARGETINDEX


PREFIX
SUFFIX


COMMAND
CLUSTER
QUEUE
WALLTIME
CPUS
QSTAT
QSUB
MAXJOBS
SLEEP
CLEANUP
VERBOSE
TEMPDIR
DOT


);
our $DATAHASH;
foreach my $key ( @DATA )	{	$DATAHASH->{lc($key)} = 1;	}




=head2

	SUBROUTINE		snpToVenn

	PURPOSE

		1. RUN snpToSav.pl FOR ALL BIN FILES PER INPUT DIR AND CHROMOSOME:

            - ANNOTATE WHETHER OR NOT SNP IS FOUND IN dbSNP

            - ANNOTATE EFFECT OF SNP IF FOUND IN CCDS GENE

	NOTES

	/nethome/bioinfo/apps/agua/0.5/bin/apps/venn/snpVenn.pl \
	--queryfile /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland/cumulative3/chr22/hit-2.sav \
	--targetfile /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/bowtie/cumulative3/chr22/hit-2.sav \
	--querylabel eland-2 \
	--targetlabel bowtie-2 \
	--queryindex 1 \
	--targetindex 33 \
	--outputdir /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/venn/eland-bowtie

		-rw-rw-rw-  1 syoung root  129M Jan 14 08:18 hit.binlevel500000.num33-10.bam
		-rw-rw-rw-  1 syoung root  1.2M Jan 17 05:41 hit.binlevel500000.num33-10.sav
		-rw-rw-rw-  1 syoung root  858K Jan 19 23:59 hit.binlevel500000.num33-10.snp

	NB: REQUIRES THAT THE *bam FILE THAT ORIGINATED THE *snp BIN FILES

		IS NAMED 'hit.bam' AND IS AVAILABLE AT THE SPECIFIED LOCATION

		IN ORDER TO LOOKUP THE CHROMOSOME SIZE FOR BIN GENERATION

=cut

sub snpToVenn {
	my $self			=	shift;
	print "Venn::Snp::run    START snpToVenn()       ", Timer::current_datetime(), "\n";

	my $outputdir 		=	$self->get_outputdir();
	my $bamfile 		=	$self->get_bamfile();
	my $querydir 		=	$self->get_querydir();
	my $targetdir 		=	$self->get_targetdir();
	my $suffix			=	$self->get_suffix();
	my $queryindex		=	$self->get_queryindex();
	my $targetindex		=	$self->get_targetindex();
	my $binlevel		=	$self->get_binlevel();
	my $querylabel		=	$self->get_querylabel();
	my $targetlabel		=	$self->get_targetlabel();
	my $samtools		=	$self->get_samtools();

	#### GET CLUSTER ARGS TO PASS THROUGH TO snpToSav.pl
	my $cluster			=	$self->get_cluster();
	my $queue			=	$self->get_queue();
	my $maxjobs			=	$self->get_maxjobs();
	my $walltime		=	$self->get_walltime();

	require Venn;
	my $venn = Venn->new();


	#### SET EXECUTABLE TO CONVERT PILEUP FORMAT TO ANNOTATED SNP FORMAT
	my $executable = "$Bin/../venn/snpToVenn.pl";

	#### CHECK INPUTS

	print "Venn::Snp::snpToVenn    outputdir: $outputdir\n" and exit if not defined $outputdir;
	print "Venn::Snp::snpToVenn    querydir: $querydir\n" and exit if not defined $querydir;
	print "Venn::Snp::snpToVenn    targetdir: $targetdir\n" and exit if not defined $targetdir;
	print "Venn::Snp::snpToVenn    querylabel: $querylabel\n" and exit if not defined $querylabel;
	print "Venn::Snp::snpToVenn    targetlabel: $targetlabel\n" and exit if not defined $targetlabel;
	print "Venn::Snp::snpToVenn    suffix: $suffix\n" and exit if not defined $suffix;
	print "Venn::Snp::snpToVenn    queryindex: $queryindex\n" and exit if not defined $queryindex;
	print "Venn::Snp::snpToVenn    targetindex: $targetindex\n" and exit if not defined $targetindex;



	File::Path::mkpath($outputdir) if not -d $outputdir;
	print "Can't create outputdir: $outputdir\n" if not -d $outputdir;

	#### GET ARRAY OF BINS
	my $binner = UCSCBin->new({	samtools	=> $samtools	});
	my $chromosome_size = $binner->getChromosomeSize($bamfile);	
	my $bins = $binner->getBins($binlevel, $chromosome_size);
#exit;

	#### DO EACH BAM BIN FILE
	my $jobs = [];
	foreach my $bin ( @$bins )
	{
		#### SET FILE NAMES
		my ($filename) = $bamfile =~ /([^\/]+)$/;
		my $binfile = $binner->setBinfile($filename, $querydir, $binlevel, $bin);
		my ($filestub) = $binfile =~ /([^\/]+)$/;
		$filestub =~ s/\.bam//;
		my $prefix = $filestub . ".";		

		#### GET BIN NUMBER
		my $binnumber = $bin->{number};
		print "Venn::Snp::snpToVenn    binnumber not defined: \n" if not defined $binnumber;
		print Dumper $bin and exit if not defined $binnumber;				

		my $queryfile = "$querydir/$filestub-$queryindex.$suffix";
		my $targetfile = "$targetdir/$filestub-$targetindex.$suffix";

		my $stdout = "$outputdir/stdout/$filestub-$queryindex-$targetindex.stdout.txt";
		print "Venn::Snp::snpToVenn    Skipping missing queryfile: $queryfile\n" and next if not -f $queryfile;
		print "Venn::Snp::snpToVenn    Skipping missing targetfile: $targetfile\n" and next if not -f $targetfile;

		my $command = qq{/usr/bin/perl $executable \\\n};
		$command .= qq{ --queryfile $queryfile \\\n};
		$command .= qq{ --targetfile $targetfile \\\n};
		$command .= qq{ --querylabel $querylabel \\\n};
		$command .= qq{ --targetlabel $targetlabel \\\n};
		$command .= qq{ --suffix $suffix \\\n};
		$command .= qq{ --prefix $prefix \\\n};
		$command .= qq{ --outputdir $outputdir \\\n};
		$command .= qq{ --stdout $stdout \n};

		#### SET LABEL
		print Dumper $bin and exit if not defined $bin->{number};
		my $label = "snpToVenn-$binnumber";

		my $checkfile = $self->setOutputFiles($outputdir, $querylabel, $targetlabel, $suffix, "$filestub.");

		#### SET UNIQUE JOB LABEL TO BE USED FOR SCRIPTFILE NAME
		my $joblabel = "$label-$filestub-$queryindex-$targetindex";

		##### RUN JOBS
		my $job = $self->setJob( [ $command ], $joblabel, $outputdir);
		$job->{checkfile} = $checkfile;
		push @$jobs, $job;
	} #### bins

	#### RUN COMMANDS IN PARALLEL
	$self->runJobs( $jobs, "snpToVenn-$queryindex-$targetindex" );	

	print "Venn::Snp::run    END snpToVenn()       ", Timer::current_datetime(), "\n";
}

=head2

	SUBROUTINE		compare

	PURPOSE

		CALCULATE THE VENN INTERSECTION OF TWO SETS OF SNPS

=cut
sub compare {
	my $self		=	shift;
	my $queryfile	=	$self->get_queryfile();
	my $targetfile	=	$self->get_targetfile();
	my $querylabel	=	$self->get_querylabel();
	my $targetlabel	=	$self->get_targetlabel();
	my $outputdir	=	$self->get_outputdir();
	my $suffix		=	$self->get_suffix();
	my $prefix		=	$self->get_prefix();


	#### SET LABELS
	if ( not defined $querylabel )
	{
		($querylabel)	=	$queryfile =~ /([^\/]+)$/ if not defined $querylabel;
		$querylabel =~ s/\.[^\.]+//;
	}
	if ( not defined $targetlabel )
	{
		($targetlabel)	=	$targetfile =~ /([^\/]+)$/;
		$targetlabel =~ s/\.[^\.]+//;
	}

	#### SET OUTPUT FILES
	my $queryonly;
	my $targetonly;
	my $both;
	($queryonly, $targetonly, $both) = $self->setOutputFiles($outputdir, $querylabel, $targetlabel, $suffix, $prefix);

	#### DO QUERY AND BOTH
	my $querycount = 0;
	my $bothcount = 0;

	my ($query_only_lines, $target_only_lines, $query_both_lines, $target_both_lines);

	##### DO QUERY
	($query_both_lines, $query_only_lines) = $self->match($queryfile, $targetfile);

	#### DO TARGET
	($target_both_lines, $target_only_lines) = $self->match($targetfile, $queryfile);

	#### THESE TWO SHOULD BE IN AGREEMENT
	@$target_both_lines = sort @$target_both_lines;
	@$query_both_lines = sort @$query_both_lines;	
	#

	print "Venn::Snp::compare    target both != query both\n" if @$target_both_lines != @$query_both_lines;

	#### PRINT OUTPUT FILES
	open(QUERYONLY, ">$queryonly") or die "Can't open queryonly file: $queryonly\n";
	foreach my $line ( @$query_only_lines )	{	print QUERYONLY $line; }
	close(QUERYONLY);
	open(TARGETONLY, ">$targetonly") or die "Can't open targetonly file: $targetonly\n";
	foreach my $line ( @$target_only_lines )	{	print TARGETONLY $line; }
	close(TARGETONLY);
	open(BOTH, ">$both") or die "Can't open both file: $both\n";
	foreach my $line ( @$query_both_lines )	{	print BOTH $line; }
	close(BOTH);

	#### PRINT INFO
	print "\n\nVenn::Snp::compare    output files printed:\n\n";
	print "Venn::Snp::compare    $queryonly\n";
	print "Venn::Snp::compare    $targetonly\n";
	print "Venn::Snp::compare    $both\n";
	print "\n\n";
}



#### SORT BY START POSITION
sub by_start {
	my ($aa) = $a =~ /^\S+\s+\S+\s+transcript\s+(\d+)/;
	my ($bb) = $b =~ /^\S+\s+\S+\s+transcript\s+(\d+)/;

	$aa <=> $bb;	
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
	my @targets = [];
	if ( $targetfile =~ /\.gz$/ )
	{
		my $pipe_command = "zcat $targetfile |";
		open(TARGET, $pipe_command);
		@targets = <TARGET>;
		close(TARGET);
	}
	else
	{
		open(TARGET, $targetfile) or die "Venn::Snp::matchCan't open target file: $targetfile\n";
		@targets = <TARGET>;
		close(TARGET);
	}

	my $matched = [];
	my $unmatched = [];

	#### OPEN QUERY FILE
	open(QUERY, $queryfile) or die "Venn::Snp::matchCan't open target file: $targetfile\n";
	my $counter = 0;
	my $dot = 1000;
	while ( <QUERY> )
	{
		next if $_ =~ /^\s*$/ or /^#/;

		$counter++;
		print "$counter\n" if $counter % $dot == 0;

#next if $counter < 1475;

		my ($queryposition) = $_ =~ /^\S+\s+(\S+)/;
		if ( not defined $queryposition )
		{
			next;
		}
		my $match = 0;
		for ( my $i = 0; $i < @targets; $i++ )
		{

#next if $i < 1475;


			my $target = $targets[$i];
			my ($targetposition) = $target =~ /^\S+\s+(\S+)/;
			if ( not defined $targetposition )
			{
				next;
			}
			$targetposition =~ s/#.+$//;


			if ( $targetposition eq $queryposition )
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
