#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();

=head2

    APPLICATION     elandHits.pl

    VERSION         0.01

    PURPOSE

        1. COUNT THE NUMBER AND PERCENTAGE OF READS ALIGNED

			BY ELAND

    INPUT

		1. INPUT FASTQ-FORMAT READ FILE

        2. ELAND OUTPUT *export.txt FILE

    OUTPUT

        1. OUTPUT THESE STATISTICS:

			NUMBER OF READS

			NUMBER OF HITS

			% HITS

    USAGE

    ./elandHits.pl <--readsfile String> <--exportfile String> [--help]

		--readsfile		:	Location of FASTQ format reads file
		--exportfile	:   Location of ELAND output *export file
		--stdout   		:   Write STDOUT to this file
		--help          :   print help info

    EXAMPLES


/nethome/bioinfo/apps/agua/0.5/bin/apps/elandHits.pl \
--readsfile /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/rerun/eland/SRX001539/10M/simpleheader/srx001539-1.reads_1.sequence.txt \
--exportfile /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/rerun/eland/SRX001539/manual/1/reanalysis_1_export.txt



=cut

use strict;

#### EXTERNAL MODULES
use Data::Dumper;
use Getopt::Long;
use File::Path;
use FindBin qw($Bin);
use Term::ANSIColor qw(:constants);

#### USE LIBRARY
use lib "$Bin/../../../lib";	
use lib "$Bin/../../../lib/external";	

#### INTERNAL MODULES
use FileTools;


#### GET OPTIONS
use Getopt::Long;
my $readsfile;
my $exportfile;
my $stdout;
my $zipped;
my $dot = 1000000;
my $help;
GetOptions (
    'readsfile=s' 	=> \$readsfile,
    'exportfile=s' 	=> \$exportfile,
    'stdout=s' 		=> \$stdout,
    'zipped=s' 		=> \$zipped,
    'dot=i' 		=> \$dot,
    'help' 			=> \$help             
) or die "No options specified. Use --help for usage\n";
usage() if defined $help;

#### CHECK INPUTS
die "readsfile is not defined (use --help for options)\n" if not defined $readsfile;
die "exportfile is not defined (use --help for options)\n" if not defined $exportfile;
die "Can't find readsfile: $readsfile\n" if not -f $readsfile;
die "Can't find exportfile: $exportfile\n" if not -f $exportfile;
die "readsfile is empty: $readsfile\n" if -z $readsfile;
die "exportfile is empty: $exportfile\n" if -z $exportfile;

#### PRINT TO STDOUT IF DEFINED stdout
if ( defined $stdout )
{
	print "BOWTIE.pl    Printing STDOUT to stdoutfile:\n\n$stdout\n\n";

	my ($stdout_path) = $stdout =~ /^(.+)\/[^\/]+$/;
	File::Path::mkpath($stdout_path) if not -d $stdout_path;
	open(STDOUT, ">$stdout") or die "Can't open STDOUT file: $stdout\n" if defined $stdout;
}

#### GET READS COUNT
my $reads = FileTools::records($readsfile, "\n", $zipped);
$reads = $reads / 4;
my $export_lines = FileTools::records($exportfile, "\n", $zipped);
print "\n";
print "input reads: $reads\n";
print "export lines: $export_lines\n";
print "input reads and export lines don't match\n" if $reads != $export_lines;


my $misses  = `cut -f 11 < $exportfile | grep -c NM`;
$misses =~ s/\s*$//;
print "misses: $misses\n";
my $quality_rejects = `cut -f 11 < $exportfile | grep -c QC`;
$quality_rejects =~ s/\s*$//;
print "quality_rejects: $quality_rejects\n";
my $hits = $export_lines - $misses - $quality_rejects;
print "hits: $hits\n";
my $pct_hits = $hits/$reads;
print "% hits: $pct_hits\n";


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#									SUBROUTINES
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#=head2
#
#	SUBROUTINE		loop
#	
#	PURPOSE
#	  
#        1. REPEATEDLY EXECUTE AN APPLICATION, CHANGING THE VALUE OF
#		
#			A PARTICULAR PARAMETER EVERY TIME
#		
#		2. REPEAT FOR A SPECIFIED NUMBER OF REPLICATES
#
#    INPUT
#
#        1. readsfile AND ITS ARGUMENTS (WITH STRING '%VALUE%' 
#		
#			IN THE PLACE OF THE ACTUAL PARAMETER VALUE)
#		
#		2. PARAMETER TO BE CHANGED
#		
#		3. COMMA-SEPARATED LIST OF VALUES FOR THE PARAMETER
#		
#		4. COMMA-SEPARATED LIST OF REPLICATES
#        
#    OUTPUT
#    
#        1. OUTPUTS OF EACH RUN OF THE readsfile USING A
#		
#			DIFFERENT VALUE FOR THE PARAMETER EACH TIME
#
#=cut
#
#sub loop
#{
#	my $self		=	shift;
#	my $readsfile	=	shift;
#	my $arguments	=	shift;
#	my $parameter	=	shift;
#	my $values		=	shift;
#	my $replicates	=	shift;
#	
#	#### RUN replicates TIMES WITH values VALUES
#	for ( my $counter = 0; $counter < @$replicates; $counter++ )
#	{
#		my $replicate = $$replicates[$counter];
#		#print "replicate: $replicate\n";
#		
#		for ( my $i = 0; $i < @$values; $i++ )
#		{
#			my $instance_args;
#			@$instance_args = @$arguments;
#			
#			my $value = $$values[$i];
#	
#			#### SUBSTITUTE replicate FOR ONE OR MORE '%REPLICATE%' STRINGS IN ALL ARGUMENTS
#			$instance_args = fill_in($instance_args, "%REPLICATE%", $replicate);
#	
#			#### SUBSTITUTE value FOR ONE OR MORE '%VALUE%' STRINGS IN ALL ARGUMENTS
#			$instance_args = fill_in($instance_args, "%VALUE%", $value);
#	
#			#### SUBSTITUTE parameter FOR ONE OR MORE '%PARAMETER%' STRINGS IN ALL ARGUMENTS
#			$instance_args = fill_in($instance_args, "%PARAMETER%", $parameter);
#	
#			my $command = "$readsfile @$instance_args";
#			
#			`$command`;		
#		}
#	}
#}
#


#=head2
#
#	SUBROUTINE		set_parameter
#	
#	PURPOSE
#	
#		SET THE VALUE OF A PARAMETER IN arguments
#
#=cut
#
#sub set_parameter
#{	
#	my $arguments		=	shift;
#	my $parameter			=	shift;
#	my $value			=	shift;
#
#	
#	for ( my $i = 0; $i < @$arguments; $i++ )
#	{
#		if ( "--$parameter" eq $$arguments[$i] )
#		{
#			$$arguments[$i + 1] = $value;
#			return $arguments;
#		}	
#	}
#	
#	return $arguments;
#}
#
#
#
#=head2
#
#	SUBROUTINE		fill_in
#	
#	PURPOSE
#	
#		SUBSTITUTE counter FOR ONE OR MORE '%REPLICATE%' STRINGS IN ALL ARGUMENTS
#
#=cut
#
#sub fill_in
#{	
#	my $arguments		=	shift;
#	my $pattern			=	shift;
#	my $value			=	shift;
#	
#
#	
#	foreach my $argument ( @$arguments )
#	{
#		$argument =~ s/$pattern/$value/ig;
#	}
#
#	return $arguments;
#}
#
#
#=head2
#
#	SUBROUTINE		get_argument
#	
#	PURPOSE
#	
#		EXTRACT AN ARGUMENT FROM THE ARGUMENT ARRAY AND RETURN IT 
#
#=cut
#
#sub get_argument
#{	
#	my $arguments		=	shift;
#	my $name			=	shift;
#
#
#	
#	my $argument;
#	for ( my $i = 0; $i < @$arguments; $i++ )
#	{
#		if ( "--$name" eq $$arguments[$i] )
#		{
#			$argument = $$arguments[$i + 1];
#			splice @$arguments, $i, 2;
#			return $argument;
#		}
#	}
#	
#	return;
#}
#


sub usage
{
	print GREEN;
	print `perldoc $0`;
	print RESET;
	exit;
}


