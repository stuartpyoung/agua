package CASAVA;


=head2

		PACKAGE		CASAVA

		PURPOSE

	        WRAPPER SCRIPT FOR RUNNING CASAVA SNP PREDICTION

=cut

use strict;
use warnings;
use Carp;

require Exporter;
#our @ISA = qw(Exporter Common);
#our @EXPORT_OK = qw();
our $AUTOLOAD;

#### INTERNAL MODULES
use Util;
use Sampler;

#### EXTERNAL MODULES
use Data::Dumper;
use File::Path;
use File::Remove;
#use MD5;

our $VERSION = 0.01;

#### SET SLOTS
our @DATA = qw(

INPUTFILE
MATEFILE
SEQUENCEDIR
RUNDIR
OUTPUTDIR
REFERENCEDIR
LANE
FIRSTLANE

CASAVA

QSTAT
JOBS
CPUS
SLEEP
DOT
QSUB
QUEUE
);
our $DATAHASH;
foreach my $key ( @DATA )	{	$DATAHASH->{lc($key)} = 1;	}



=head2

	SUBROUTINE		copyFiles

	PURPOSE

		COPY THE INFO FILES REQUIRED BY CASAVA FROM

		THE RUN DIRECTORY TO THE OUTPUT DIRECTORY

=cut

sub copyFiles
{
	my $self			=	shift;


	#### USER INPUTS
	my $inputfile 		=	$self->get_inputfile();
	my $matefile 		=	$self->get_matefile();
	my $sequencedir 	=	$self->get_sequencedir();
	my $rundir 			=	$self->get_rundir();
	my $outputdir 		=	$self->get_outputdir();
	my $referencedir 	=	$self->get_referencedir();
	my $lane 			=	$self->get_lane();
	my $firstlane 		=	$self->get_firstlane();

	##### EXECUTABLES
	#my $casava 			=	$self->get_casava();
	#
	##### CLUSTER
	#my $qstat 			=	$self->get_qstat();
	#my $jobs 			=	$self->get_jobs();
	#my $sleep 			=	$self->get_sleep();
	#my $qsub 			=	$self->get_qsub();
	#my $queue 			=	$self->get_queue();

	#### SET PAIRED FLAG
	my $paired = 0;
	$paired = 1 if defined $matefile and $matefile !~ /^\s*$/;

	#### CHANGE TO OUTPUT DIR
	chdir($outputdir) or die "Can't change to output directory: $outputdir\n";

	#### GET PAIR FILE FOR FIRST INPUT FILE
	my $pair_xml = "s_" . $firstlane . "_pair.xml";
	File::Copy::copy("$sequencedir/$pair_xml", "$outputdir/$pair_xml") or die "Can't copy pair_xml file: $sequencedir/$pair_xml\n";

	#### GET Summary.htm
	File::Copy::copy("$sequencedir/Summary.htm", "$outputdir/Summary.htm") or die "Can't copy pair_xml file: $sequencedir/Summary.htm\n";

	#### GET config.xml
	File::Copy::copy("$sequencedir/config.xml", "$outputdir/config.xml") or die "Can't copy pair_xml file: $sequencedir/config.xml\n";

}	#	copyFiles



=head2

	SUBROUTINE		shellScript

	PURPOSE

		CREATE .sh FILE

=cut


sub shellScript
{
	my $self			=	shift;


	#### USER INPUTS
	my $inputfile 		=	$self->get_inputfile();
	my $matefile 		=	$self->get_matefile();
	my $rundir 			=	$self->get_rundir();
	my $outputdir 		=	$self->get_outputdir();
	my $referencedir 	=	$self->get_referencedir();
	my $lane 			=	$self->get_lane();
	my $firstlane 		=	$self->get_firstlane();

	#### EXECUTABLES
	my $casava 			=	$self->get_casava();

	#### CLUSTER
	my $qstat 			=	$self->get_qstat();
	my $jobs 			=	$self->get_jobs();
	my $cpus 			=	$self->get_cpus();
	my $sleep 			=	$self->get_sleep();
	my $qsub 			=	$self->get_qsub();
	my $queue 			=	$self->get_queue();

	#### CHANGE TO OUTPUT DIR
	chdir($outputdir) or die "Can't change to output directory: $outputdir\n";

	#### SET SCRIPT FILE AND REMOVE IF PRESENT
	my $shellscript = "$outputdir/casava-$lane.sh";
	print "CASAVA::shellScript    shellscript: $shellscript\n";
	my $recursive = 0;
	File::Remove::remove(\$recursive, $shellscript) if -f $shellscript;

	#### REMOVE BUILD DIRECTORY IF IT EXISTS
	my $builddir = "$outputdir/casava-$lane";
	$recursive = 1;
	File::Remove::remove(\$recursive, $builddir) if -d $builddir;

	my $run_command = $self->runCommand();

my $script = qq{
#!/bin/sh

#PBS -N casava-$lane
#PBS -l ncpus=$cpus
#PBS -j oe

cd $outputdir

$run_command

#### END OF SCRIPT
};
	open(SHFILE, ">$shellscript") or die "Can't open script file: $shellscript\n";
	print SHFILE $script;
	close(SHFILE);

	print "CASAVA::shellScript    script:\n";
	print "$script\n";
	print "\n";

	return $shellscript;

}	#	shellScript



=head2

	SUBROUTINE		run

	PURPOSE

		CREATE .sh FILE

=cut

sub run
{
	my $self			=	shift;


	my $run_command = $self->runCommand();


	#### RUN
	#`$run_command`;
}



=head2

	SUBROUTINE		runCommand

	PURPOSE

		CREATE .sh FILE

=cut


sub runCommand
{
	my $self			=	shift;


	#### USER INPUTS
	my $inputfile 		=	$self->get_inputfile();
	my $matefile 		=	$self->get_matefile();
	my $rundir 			=	$self->get_rundir();
	my $outputdir 		=	$self->get_outputdir();
	my $referencedir 	=	$self->get_referencedir();
	my $lane 			=	$self->get_lane();
	my $firstlane 		=	$self->get_firstlane();

	#### EXECUTABLES
	my $casava 			=	$self->get_casava();

	##### CLUSTER
	#my $qstat 			=	$self->get_qstat();
	#my $jobs 			=	$self->get_jobs();
	#my $cpus 			=	$self->get_cpus();
	#my $sleep 			=	$self->get_sleep();
	#my $qsub 			=	$self->get_qsub();
	#my $queue 			=	$self->get_queue();

	#### GET RUN NAME
	my ($run) = $rundir =~ /([^\/]+)$/;

	my $command = qq{$casava/run.pl \\
--exportDir=$outputdir \\
--genomeSize=$casava/features/human_genome_size.xml \\
--lanes=$lane \\
--projectDir=$outputdir/casava-$lane \\
--runId=$run \\
--refSequences=$referencedir \\
--workflowAuto \\
--snpCovCutoff=-1 
};
	return $command;

}	#	runCommand



=head2

	SUBROUTINE		collectSNPs

	PURPOSE

		COLLECT SNPS FILES GENERATED BY RUN INTO

		A 'SNPs' DIRECTORY IN THE OUTPUT DIRECTORY

=cut

sub collectSNPs
{
	my $self			=	shift;


	#### USER INPUTS
	#my $inputfile 		=	$self->get_inputfile();
	#my $matefile 		=	$self->get_matefile();
	#my $rundir 			=	$self->get_rundir();
	my $outputdir 		=	$self->get_outputdir();
	#my $referencedir 	=	$self->get_referencedir();
	#my $lane 			=	$self->get_lane();
	#my $firstlane 		=	$self->get_firstlane();

	##### EXECUTABLES
	#my $casava 			=	$self->get_casava();
	#
	##### CLUSTER
	#my $qstat 			=	$self->get_qstat();
	#my $jobs 			=	$self->get_jobs();
	#my $cpus 			=	$self->get_cpus();
	#my $sleep 			=	$self->get_sleep();
	#my $qsub 			=	$self->get_qsub();
	#my $queue 			=	$self->get_queue();
	#


	#### CHANGE TO OUTPUT DIR
	chdir($outputdir) or die "Can't change to output directory: $outputdir\n";

	#### CREATE OUTPUT SUBDIR FOR SNPs
	my $snpdir = "$outputdir/snp-text";
	print "CASAVA::collectSNPs    SNP directory already exists as a file: $snpdir. Exiting\n" and exit if -f $snpdir;
	File::Path::mkpath($snpdir) if not -d $snpdir;
	print "CASAVA::collectSNPs    Can't create snpdir: $snpdir\n" and exit if not -d $snpdir;

	#### GET 'TXT' SNP FILES IN 'Parsed*' DIRECTORY
	my @textfiles = <Parsed*/c*/c*.snp.txt>;

	#### COLLECT SNPs
	foreach my $textfile ( @textfiles )
	{
		print "CASAVA::collectSNPs    Copying file: $outputdir/$textfile to snp directory: $snpdir";
		File::Copy::copy("$outputdir/$textfile", $snpdir) or die "Can't copy file: $outputdir/$textfile to snp directory: $snpdir";
	}

	#### GET 'GFF' SNP FILES IN 'Parsed*' DIRECTORY
	my @gfffiles = <Parsed*/c*/c*.snp.txt>;

	#### COLLECT SNPs
	foreach my $gfffile ( @gfffiles )
	{
		print "CASAVA::collectSNPs    Copying file: $outputdir/$gfffile to snp directory: $snpdir";
		File::Copy::copy("$outputdir/$gfffile", $snpdir) or die "Can't copy file: $outputdir/$gfffile to snp directory: $snpdir";
	}

}	#	collectSNPs


################################################################################
##################			HOUSEKEEPING SUBROUTINES			################
################################################################################


=head2

	SUBROUTINE		new

	PURPOSE

		CREATE THE NEW self OBJECT AND INITIALISE IT, FIRST WITH DEFAULT 

		ARGUMENTS, THEN WITH PROVIDED ARGUMENTS

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

		INITIALISE THE self OBJECT:

			1. LOAD THE DATABASE, USER AND PASSWORD FROM THE ARGUMENTS

			2. FILL OUT %VARIABLES% IN XML AND LOAD XML

			3. LOAD THE ARGUMENTS

=cut


sub initialise
{
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

	SUBROUTINE		value

	PURPOSE

		SET A PARAMETER OF THE self OBJECT TO A GIVEN value

    INPUT

        1. parameter TO BE SET

		2. value TO BE SET TO

    OUTPUT

        1. THE SET parameter INSIDE THE BioSVG OBJECT

=cut
sub value
{
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

sub validate_arguments
{
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

sub is_valid
{
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
    unless( exists $self->{$attribute} )
	{
        #croak "No such attribute '$attribute' exists in the class ", ref($self);
		return undef;
    }

    # Turn off strict references to enable "magic" AUTOLOAD speedup
    no strict 'refs';

    # AUTOLOAD accessors
    if($operation eq 'get') {

		return $self->{$attribute};

        # define subroutine
        *{$AUTOLOAD} = sub { shift->{$attribute} };

    # AUTOLOAD mutators
    }elsif($operation eq 'set') {

        # define subroutine
		print "CASAVA::AUTOLOAD    Setting attribute: $attribute\n";
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

