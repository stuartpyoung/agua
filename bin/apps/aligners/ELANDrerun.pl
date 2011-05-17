#!/usr/bin/perl -w

#### DEBUG
$DEBUG = 1;

#### TIME
my $time = time();
my $duration = 0;
my $current_time = $time;

=head2

    APPLICATION     XELAND

    VERSION         0.01

    PURPOSE

        CHECK OUTPUT FILES FROM XELAND ASSEMBLY

	EXAMPLES

/nethome/syoung/0.5/bin/apps/aligners/ELANDrerun.pl \
--replicates 1-33 \
--inputtype fastq \
--label eland \
--cluster LSF \
--referencedir /nethome/bioinfo/data/sequence/chromosomes/human/hg19/eland/chr22 \
--outputdir /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland \
--min 0.5 \
--max 2


mv /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland/2/chr22/42/out.sam /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/sampled/200bp/chr22/eland/2/chr22/42/out.sam.bkp
=cut


use strict;

#### USE LIBRARY
use Scalar::Util qw(weaken);
use FindBin qw($Bin);
use lib "$Bin/../../../lib/external";
use lib "$Bin/../../../lib";
BEGIN {
    unshift @INC, "/nethome/syoung/base/pipeline/moose/tmp/lib64/perl5/site_perl/5.8.8/x86_64-linux-thread-multi";
    unshift @INC, "/nethome/syoung/0.5/lib/external/perl5-32/site_perl/5.8.8";
    unshift @INC, "/nethome/syoung/0.5/lib/external/perl5-64/site_perl/5.8.8/x86_64-linux-thread-multi";
    unshift @INC, "/nethome/syoung/0.5/lib/external/perl5-32/5.8.8";    
}

#### EXTERNAL MODULES
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use Getopt::Long;
use FindBin qw($Bin);
use File::Path;


#### INTERNAL MODULES
use XELAND;
use Timer;
use Util;
use Conf::Agua;

##### STORE ARGUMENTS FOR PRINTING TO USAGE FILE LATER
my @arguments = @ARGV;
unshift @arguments, $0;

#### FLUSH BUFFER
$| =1;

#### SET XELAND LOCATION
my $conf = Conf::Agua->new(inputfile=>"$Bin/../../../conf/default.conf");
my $casava = $conf->getKeyValue("applications", 'CASAVA');

#### GET OPTIONS
my $select;
my $min;
my $max;
my $rundir;
my $outputdir;
my $replicates;
my $referencedir;
my $label;
my $stdout;
my $inputtype;
my $cluster;
my $verbose;

my $help;
if ( not GetOptions (

	#### GENERAL
    'select=s' 		=> \$select,
    'min=s' 		=> \$min,
    'max=s' 		=> \$max,
    'outputdir=s' 	=> \$outputdir,
    'replicates=s' 	=> \$replicates,
    'referencedir=s' => \$referencedir,
    'label=s' 		=> \$label,
    'stdout=s' 		=> \$stdout,
    'inputtype=s' 	=> \$inputtype,
    'cluster=s' 	=> \$cluster,
    'verbose' 		=> \$verbose,
    'help' 			=> \$help
) )
{ print "Use option --help for usage instructions.\n";  exit;    };

#### MAKE OUTPUT DIR IF NOT EXISTS
print "XELANDrerun.pl    Output directory is a file: $outputdir\n" and exit if -f $outputdir;
File::Path::mkpath($outputdir) if not -d $outputdir;
print "Could not create output directory: $outputdir\n" if not -d $outputdir;

#### PRINT TO STDOUT IF DEFINED stdout
if ( defined $stdout )
{
	my ($stdout_path) = $stdout =~ /^(.+)\/[^\/]+$/;
	File::Path::mkpath($stdout_path) if not -d $stdout_path;
	open(STDOUT, ">$stdout") or die "Can't open STDOUT file: $stdout\n" if defined $stdout;
}

#### PRINT HELP
if ( defined $help )	{	usage();	}

#### CHECK INPUTS
die "min not defined (Use --help for usage)\n" if not defined $min;
die "max not defined (Use --help for usage)\n" if not defined $max;
die "outputdir not defined (Use --help for usage)\n" if not defined $outputdir;
die "replicates not defined (Use --help for usage)\n" if not defined $replicates;
die "referencedir not defined (Use --help for usage)\n" if not defined $referencedir;
die "label not defined (Use --help for usage)\n" if not defined $label;

#### DEBUG
print "outputdir: $outputdir\n";
print "referencedir: $referencedir\n";


my $eland = XELAND->new(
	{
		min 		=> $min,
		max 		=> $max,
		inputtype 	=> $inputtype,
		outputdir 	=> $outputdir,
		replicates 	=> $replicates,
		referencedir=> $referencedir,
		conf 		=> $conf,
		casava		=> $casava,
		label 		=> $label,
		verbose 	=> $verbose,
		cluster 	=> $cluster
	}
);

my ($completed, $sublabels, $missingfiles, $dubiousfiles);
($completed, $label, $sublabels, $missingfiles, $dubiousfiles) = $eland->rerun($select);	

#### SEND JOB COMPLETION SIGNAL
print "\n------------------------------------------------------------\n";
print "---[completed $label: $completed $sublabels]---";
if ( scalar(@$missingfiles) > 0 )
{
	print "\n";
	print scalar(@$missingfiles);
	print " missing file:\n" if scalar(@$missingfiles) == 1;
	print " missing files:\n" if scalar(@$missingfiles) != 1;
	print @$missingfiles if scalar(@$missingfiles) > 0;	
}
if ( scalar(@$dubiousfiles) > 0 )
{
	print "\n";
	print scalar(@$dubiousfiles);
	print " dubious file:\n" if scalar(@$dubiousfiles) == 1;
	print " dubious files:\n" if scalar(@$dubiousfiles) != 1;
	print "lowerbound\tupperbound\tsizemultiple\taveragesize\tfilesize\tlocation\n";
	print @$dubiousfiles if scalar(@$dubiousfiles) > 0;	
}
print "\n------------------------------------------------------------\n";

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "\nRun time: $runtime\n";
print "Completed $0 @arguments\n";
print Timer::datetime(), "\n";
print "****************************************\n\n\n";

#### PRINT TO STDOUT IF DEFINED stdout
close(STDOUT) or die "Can't close STDOUT file\n";
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

