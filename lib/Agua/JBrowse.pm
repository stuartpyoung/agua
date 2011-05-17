use MooseX::Declare;

=head2

	PACKAGE		JBrowse

	PURPOSE

		THIS MODULE ENABLES THE FOLLOWING USE SCENARIOS:


		1. USER GENERATES JBROWSE JSON FILES IN OWN FOLDER

		    1.1 CREATE USER-SPECIFIC VIEW FOLDER

			1.2 GENERATE JBrowse FEATURES IN PARALLEL ON A CLUSTER

		        plugins/view/jbrowse/users/__username__/__project__/__view__/data

			1.3 COPY BY ln -s ALL STATIC FEATURE TRACK SUBFOLDERS

				(I.E., AT THE chr*/FeatureDir LEVEL) TO THE USER'S VIEW FOLDER.

				NB: __NOT__ AT A HIGHER LEVEL BECAUSE WE WANT TO BE ABLE

				TO ADD/REMOVE DYNAMIC TRACKS IN THE USER'S VIEW FOLDER WITHOUT

				AFFECTING THE PARENT DIRECTORY OF THE STATIC TRACKS.

				THE STATIC TRACKS ARE MERELY INDIVIDUALLY 'BORROWED' AS IS.


		2. USER ADDS/REMOVES TRACKS TO/FROM VIEW

			2.1 ADD/REMOVE FEATURE trackData.json INFORMATION TO/FROM data/trackInfo.js

			2.2 RERUN generate-names.pl AFTER EACH ADD/REMOVE


=cut

use strict;
use warnings;
use Carp;

#### EXTERNAL MODULES
use FindBin qw($Bin);

#### USE LIB
use lib "$Bin/..";

#### USES ROLES
use Agua::Cluster::Checker;
use Agua::Cluster::Cleanup;
use Agua::Cluster::Jobs;
use Agua::Cluster::Loop;
use Agua::Cluster::Usage;
use Agua::Cluster::Util;

class Agua::JBrowse with (Agua::Cluster::Checker,
	Agua::Cluster::Cleanup,
	Agua::Cluster::Jobs,
	Agua::Cluster::Loop,
	Agua::Cluster::Usage,
	Agua::Cluster::Util)
{


#### INTERNAL MODULES
use Sampler;

use Data::Dumper;
use File::Path;
use File::Copy;
use JSON 2;
use IO::File;
use Fcntl ":flock";
use POSIX qw(ceil floor);

# FLAGS
has 'compress'	=> ( isa => 'Bool|Undef', is => 'rw', default => 0 );
# INTS
has 'chunksize'	=> ( isa => 'Int|Undef', is => 'rw', default => 0 );
has 'sortmem'	=> ( isa => 'Int|Undef', is => 'rw', default => 0 );
# STRINGS
has 'configfile'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'inputdir'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'outputdir'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'filename'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'filetype'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'label'			=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'key'			=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'refseqfile'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'trackfile'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'jbrowse'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'species'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'build'			=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'username'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'htmlroot'		=> ( isa => 'Str|Undef', is => 'rw', default => '' );
has 'chromofile'	=> ( isa => 'Str|Undef', is => 'rw', default => '' );
# OBJECTS
has 'config'		=> ( isa => 'HashRef', is => 'rw');
has 'jsonparser'	=> ( isa => 'JSON', is => 'rw');

	####/////}

method BUILD ($hash) {

	$self->setJsonParser();
}

=head2

	SUBROUTINE		gtfToGff

	PURPOSE

		CONVERT GTF TO GFF FILES

	INPUTS

		1. REFSEQS FILE CONTAINING AN ENTRY FOR EACH REFERENCE CHROMOSOME

=cut

method gtfToGff ($args) {


	#### DO SINGLE FILE IF SPECIFIED	
	$self->_gtfToGff($args)
		and return if not defined $args->{inputdir};

	#### DO ALL FILES IN INPUT DIRECTORY
	opendir(DIR, $args->{inputdir}) or die "Can't open inputdir directory: $args->{inputdir}\n";
	my @infiles = readdir(DIR);
	close(DIR);

	foreach my $infile ( @infiles )	
	{
		my $filepath = "$args->{inputdir}/$infile";
		next if not -f $filepath;
		next if not $infile =~ /\.gtf$/;

		$args->{inputfile} = $infile;
		($args->{feature}) = $infile;
		$args->{feature} =~ s/\.[^\.]{3,5}$//;
		$self->_gtfToGff($args);	
	}		
}

=head

	SUBROUTINE		_gtfToGff

	PURPOSE

		CONVERT A CHROMOSOME GTF FILE INTO GFF FORMAT

	NOTES

		TOP FOLDER OF inputdir PATH MUST BE THE NAME OF THE

		CHROMOSOME
=cut

method _gtfToGff ($args) {


	my $inputdir	=	$args->{inputdir};
	my $inputfile	=	$args->{inputfile};
	my $outputdir	=	$args->{outputdir};
	my $feature		=	$args->{feature};
	my $refseqfile	=	$args->{refseqfile};

	#### GET REFSEQ	RECORD FOR THIS CHROMOSOME
	my ($chromosome) = $inputdir =~ /([^\/]+)$/;
	my $refseq = $self->getRefseq($refseqfile, $chromosome);
	print "refseq:\n";
	print Dumper $refseq;

	return if not $refseq;

	my @args;
	push @args, "--inputfile $inputfile";
	push @args, "--inputdir $inputdir";
	push @args, "--outputdir $outputdir";
	push @args, "--feature $feature";
	push @args, "--refseqfile $refseqfile";

	print "Agua::JBrowse::_gtfToGff    Doing files in input directory: $inputdir\n";
	#### GET REFERENCE SEQUENCE INFO
	my $reference = $refseq->{name};
	my $start = $refseq->{start};
	my $end = $refseq->{end};

	#### CREATE OUTPUT SUBDIR IF NOT EXISTS
	my $output_subdir = "$outputdir/$reference";
	File::Path::mkpath($output_subdir) or die "Can't create output subdir: $output_subdir" if not -d $output_subdir;

	my $infile = "$inputdir/$inputfile";
	my $outfile = "$outputdir/$inputfile";
	$outfile =~ s/gtf$/gff/;

	#### OPEN OUTPUT FILE AND PRINT RUN COMMAND, DATE AND REFERENCE SEQUENCE LINE
	print "Agua::JBrowse::_gtfToGff    Can't find input file: $infile\n" and next if not -f $infile;
	open(OUTFILE, ">$outfile") or die "Can't open output file: $outfile\n";
	print OUTFILE "### $0 @args\n";
	print OUTFILE "####", Util::datetime(), "\n\n";
	print OUTFILE "$reference\trefseqfile\trefseqfile\t$start\t$end\t.\t.\t.\tName=$reference\n";

	#### OPEN INPUT FILE
	open(FILE, $infile) or die "Can't open input file: $infile\n";
	$/ = "\n";
	my $counter = 0;
	while ( <FILE> )
	{
		next if $_ =~ /^\s*$/;
		$counter++;
		if ( $counter % 10000 == 0 ) {	print "$counter\n";	}

		my ($start, $last) = $_ =~ /^(\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+)(.+)$/;
		print "Agua::JBrowse::_gtfToGff    last not defined in line: $_\n"
			and next if not defined $start;
		print "Agua::JBrowse::_gtfToGff    last not defined in line: $_\n"
			and next if not defined $last;


		#### PARSE END OF LINE INFO
		my @elements = split ";", $last;
		$last = '';
		foreach my $element ( @elements )
		{
			next if $element =~ /^\s*$/;
			$element =~ s/^\s+//;
			$element =~ s/\s+$//;
			my ($key, $value) = $element =~ /^(\S+)\s+(.+)$/;

			if ( $key eq "transcript_id" )
			{
				$value =~ s/"//g;
				$last = "Name=$value";
			}
		}
		$last =~ s/;$//;

		#### PARSE START OF LINE INFO
		my @fields = split " ", $start;

		#### SET FEATURE
		$fields[2] = $feature;

		#### CORRECT SCORE TO NO DECIMAL PLACES
		$fields[5] =~ s/\.[0]+//g;

		#### ADD LAST ENTRY TO FIELDS
		push(@fields, $last);
		my $line = join "\t", @fields;
		print OUTFILE "$line\n";
	}
	close(FILE);
	close(OUTFILE);
	print "$outfile\n";

}


method getRefseq ($refseqfile, $chromosome) {

	#### CREATE REFSEQ HASH
	my $skip_assign = "refSeqs =";
	my $refseqs = $self->readJson($refseqfile, $skip_assign);

	foreach my $entry ( @$refseqs )
	{
		return $entry if $entry->{name} eq $chromosome;
	}

	return;
}

#### GENERATE refSeqs.js
=head2

	SUBROUTINE		generateRefseq

	PURPOSE	

		CREATE A refSeqs.js JSON FILE CONTAINING ENTRIES FOR ALL CHROMOSOMES

		IN THE REFERENCE GENOME, E.G., IF NEED TO MODIFY SEQUENCE DIR OR

		CHUNK SIZE INDEPENDENTLY OF RUNNING JBROWSE'S prepare-refseqs.pl

    INPUT

			1. chromosome-sizes.txt FILE GENERATED BY chromosomeSizes.pl

			cat /nethome/bioinfo/data/sequence/chromosomes/rat/rn4/fasta/chromosome-sizes.txt

				chr1    0       273269105       273269105
				chr2    273269106       536640798       263371692
				chr3    536640799       711125402       174484603
				chr4    711125403       901993930       190868527
				chr5    901993931       1078552066      176558135

			2. CHUNK SIZE FOR GENERATING FEATURES

			3. OUTPUT DIRECTORY

    OUTPUT

        1. OUTPUT FILE refSeqs.js IN OUTPUT DIRECTORY

		refSeqs =
		[
		   {
			  "length" : 247249719,
			  "name" : "chr1",
			  "seqDir" : "data/seq/chr1",
			  "seqChunkSize" : 20000,
			  "end" : 247249719,
			  "start" : 0
		   }
		   ,
		   {
			  "length" : 242951149,
			  "name" : "chr2",
			  "seqDir" : "data/seq/chr2",
			  "seqChunkSize" : 20000,
			  "end" : 242951149,
			  "start" : 0
		   }
		   ,
		   ...
		]

	EXAMPLES

/nethome/bioinfo/apps/agua/0.5/bin/apps/jbrowseRefseq.pl \
--chromofile /nethome/bioinfo/data/sequence/chromosomes/rat/rn4/fasta/chromosome-sizes.txt \
--outputdir /nethome/syoung/base/pipeline/jbrowse/ucsc/0.5/rat/rn4 \
--chunk 20000



=cut

method generateRefseq ($outputdir, $chromofile, $chunksize) {

	#### GET VARIABLES IF NOT DEFINED
	$outputdir	=	$self->outputdir() if not defined $outputdir;
	$chromofile	=	$self->chromofile() if not defined $chromofile;
	$chunksize	=	$self->chunksize() if not defined $chunksize;

	#### CHECK INPUTS
	die "Agua::JBrowse::generateRefseq    outputdir not defined (Use --help for usage)\n" if not defined $outputdir;
	die "Agua::JBrowse::generateRefseq    chromofile not defined (Use --help for usage)\n" if not defined $chromofile;
	die "Agua::JBrowse::generateRefseq    chunksize not defined (Use --help for usage)\n" if not defined $chunksize;

	#### OPEN INPUT FILE
	open(FILE, $chromofile) or die "Can't open chromofile: $chromofile\n";
	my @lines = <FILE>;
	close(FILE);

	#### CREATE OUTPUT DIR IF NOT EXISTS	
	print "Agua::Cluster::generateRefseq    outputdir is a file: $outputdir\n" if -f $outputdir;
	File::Path::mkpath($outputdir) if not -d $outputdir;
	print "Can't create output dir: $outputdir\n" if not -d $outputdir;

	#### OPEN OUTPUT FILE
	my $outputfile = "$outputdir/refSeqs.js";
	open(OUTFILE, ">$outputfile") or die "Can't open outputfile: $outputfile\n";

	my $data;
	foreach my $line ( @lines )
	{
		next if $line =~ /^#/ or $line =~ /^\s*$/;

		#### FORMAT:
		#### chr1    0       273269105       273269105
		#### chr2    273269106       536640798       263371692
		my ($chromosome, $length) = $line =~ /^(\S+)\s+\S+\s+\S+\s+(\S+)$/;
		print "chromosome not defined in line: $line\n" and exit if not defined $chromosome;
		print "length not defined in line: $line\n" and exit if not defined $length;
		push @$data, {
			length 	=>	$length,
			name	=>	$chromosome,
			seqDir	=>	"data/seq/$chromosome",
			seqChunkSize	=>	$chunksize,
			end		=>	$length,
			start	=>	0
		};
	}

	#### WRITE TO JSON FILE WITH ASSIGNED NAME 'refSeqs'	
	my $callback = sub { return $data };
	$self->modifyJsonfile($outputfile, "refSeqs", $callback);
	print "Agua::JBrowse::generateRefseq    printed outputfile:\n";
	print "\n$outputfile\n\n";
}


#### GENERATE names.json
=head2

	SUBROUTINE		generateNames

	PURPOSE

		UPDATE OR CREATE NEW names.json FILE INCORPORATING

		ALL FEATURES IN THE trackInfo.js FILE

=cut

method generateNames () {
	my $inputdir	=	$self->inputdir();
	my $outputdir	=	$self->outputdir();
	my $refseqfile	=	$self->refseqfile();
	my $filetype	=	$self->filetype();
	my $jbrowse		=	$self->jbrowse();

#die "outputdir not defined (Use --help for usage)\n" if not defined $outputdir;
#die "inputdir not defined (Use --help for usage)\n" if not defined $inputdir;
#die "filetype not defined (Use --help for usage)\n" if not defined $filetype;
#die "filename not defined (Use --help for usage)\n" if not defined $filename;
#die "label not defined (Use --help for usage)\n" if not defined $label;
#die " not defined (Use --help for usage)\n" if not defined $key;
#die "refseqfile not defined (Use --help for usage)\n" if not defined $refseqfile;
#die "jbrowse not defined (Use --help for usage)\n" if not defined $jbrowse;

	print "Agua::JBrowse::generateNames    Agua::JBrowse::generateNames()\n";
	print "Agua::JBrowse::generateNames    inputdir: $inputdir\n";
	print "Agua::JBrowse::generateNames    outputdir: $outputdir\n";
	print "Agua::JBrowse::generateNames    refseqfile: $refseqfile\n";
	print "Agua::JBrowse::generateNames    jbrowse: $jbrowse\n";



	#### GET LIST OF REFERENCES
	my $references = $self->parseReferences();

	#### SET EXECUTABLE
	my $executable = "$jbrowse/generate-names.pl";

	#### DO ALL FILES SPECIFIED IN FEATURES HASH
	my $jobs = [];
	foreach my $reference ( @$references )
	{
		print "Agua::JBrowse::runGenerateNames    reference $reference\n";

		#### CONVERT TO GFF
		####/nethome/syoung/base/apps/aqwa/0.4/html/plugins/view/jbrowse/bin/generate-names.pl \
		####-v /nethome/syoung/base/pipeline/jbrowse/ucsc/reference/chr1/data/tracks/*/*/names.json

		my $command = "cd $outputdir/$reference; $executable -v $outputdir/$reference/data/tracks/*/*/names.json";
		print "command: $command\n";
		my $job = $self->setJob( [$command], "generateNames-$reference", $outputdir );
		push @$jobs, $job;
	}

	print "Agua::JBrowse::runGenerateNames    No. jobs: ", scalar(@$jobs), "\n";
	print "Agua::JBrowse::runGenerateNames    jobs[0]: \n";
	print Dumper $$jobs[0];
	exit;

	#### RUN COMMANDS IN PARALLEL
	$self->runJobs( $jobs, "generateNames" );	
}


method addTrackdata (Str $trackdatafile, Str $trackinfofile) {


	my $trackinfo = $self->getTrackinfo($trackinfofile);
	my $trackdata = $self->getTrackdata($trackdatafile);

	#### ADD TRACK 
	my $url = 'data/tracks/{refseq}/'
				. $trackdata->{key}
				. '/trackData.json';
	my $entry;
	$entry->{key} = $trackdata->{key};
	$entry->{label} = $trackdata->{label};
	$entry->{type} = $trackdata->{type};
	$entry->{url} = $url;

	#### QUIT IF FEATURE ALREADY EXISTS IN LIST
	print "Agua::JBrowse::addTrackdata    trackdata already found in trackinfo: \n"
	and print Dumper $entry and return 
	if $self->objectInArray($trackinfo, $trackdata, ["key", "label", "type", "url"]);

	#### OTHERWISE, ADD TO LIST
	push @$trackinfo, $entry;
	$self->writeJson($trackinfofile, $trackinfo, "trackInfo = ");	
}

method addTrackinfo (Str $sourceinfofile, Str $targetinfofile) {

	my $targetinfo = $self->getTrackinfo($targetinfofile);
	print "Agua::JBrowse::addTrackinfo    targetinfo not defined or empty for targetinfofile: $targetinfofile\n"
		and exit if not defined $targetinfo or not @$targetinfo;

	my $sourceinfo = $self->getTrackinfo($sourceinfofile);
	print "Agua::JBrowse::addTrackinfo    sourceinfo not defined or empty for sourceinfofile: $sourceinfofile\n"
		and exit if not defined $sourceinfo or not @$sourceinfo;

	#### QUIT IF FEATURE ALREADY EXISTS IN LIST
	if ( $self->objectInArray($targetinfo, $$sourceinfo[0], ["key", "label", "type", "url"]) )
	{
		return 0;
	}

	#### OTHERWISE, ADD TO LIST
	push @$targetinfo, $$sourceinfo[0];

	$self->writeJson($targetinfofile, $targetinfo, "trackInfo = ");	

	return 1;
}


method removeTrackinfo (HashRef $featureobject, Str $targetinfofile) {
	print "Agua::JBrowse::removeTrackinfo    Agua::JBrowse::removeTrackinfo(featureobject, targetinfofile)\n";	

	my $targetinfo = $self->getTrackinfo($targetinfofile);
	print "Agua::JBrowse::removeTrackinfo    targetinfo not defined or empty for targetinfofile: $targetinfofile\n"
		and exit if not defined $targetinfo or not @$targetinfo;

	#### QUIT IF FEATURE ALREADY EXISTS IN LIST
	my ($index) = $self->_indexInArray($targetinfo, $featureobject, ["key", "label"]);
	if ( not defined $index )
	{
		return ;
	}

	#### OTHERWISE, REMOVE FROM LIST
	splice(@$targetinfo, $index, 1);

	$self->writeJson($targetinfofile, $targetinfo, "trackInfo = ");
}


method removeTrackdata (Str $trackdatafile, Str $trackinfofile) {


	my $trackinfo = $self->getTrackinfo($trackinfofile);
	my $trackdata = $self->getTrackdata($trackdatafile);

	#### ADD TRACK 
	my $url = 'data/tracks/{refseq}/'
				. $trackdata->{key}
				. '/trackData.json';
	my $entry;
	$entry->{key} = $trackdata->{key};
	$entry->{label} = $trackdata->{label};
	$entry->{type} = $trackdata->{type};
	$entry->{url} = $url;
#exit;

	#### QUIT IF FEATURE ALREADY EXISTS IN LIST
	my $index = $self->objectInArray($trackinfo, $trackdata, ["key", "label", "type", "url"]);
	print "Agua::JBrowse::removeTrackdata    trackdata not found in trackinfo: \n"
	and print Dumper $entry and return 
	if not defined $index;
#exit;

	#### OTHERWISE, REMOVE FROM LIST
	splice(@$trackinfo, $index, 1);
	$self->writeJson($trackinfofile, $trackinfo, "trackInfo = ");	
}

method getTrackdata ($trackdatafile) {

	my $trackdata = $self->readJson($trackdatafile, 0);

	return $trackdata;
}

method getTrackinfo ($trackinfofile) {


	my $trackinfo = $self->readJson($trackinfofile, "trackInfo = ");

	return $trackinfo;
}

#### GENERATE FEATURES
=head2

	SUBROUTINE		generateFeatures

	PURPOSE

		DO ALL FEATURES OF ALL CHROMOSOMES IN PARALLEL:

			1. INPUT GFF FILES (E.G., DOWNLOADED FROM UCSC) IN

				EACH chr* DIRECTORY INSIDE THE referencedir

			2. USE flatfile-to-json.pl TO GENERATE FEATURES

	INPUTS

		1. DIRECTORY CONTAINING INPUT GFF OR BAM FILES

			OR

			A GFF OR BAM FILE

		2. refSeqs.json FILE CONTAINING 

		3. OUTPUT DIRECTORY TO PRINT data DIR

	NOTES

		1. CREATE DIRECTORIES FOR OUTPUT FILES

		2. COPY refSeqs.js TO EACH SUB DIRECTORY

		3. RUN flatfile-to-json.pl AS ARRAY JOB

=cut

method generateFeatures () {

	my $inputdir	=	$self->inputdir();
	my $featuresdir	=	$self->featuresdir();
	my $refseqfile	=	$self->refseqfile();
	my $configfile	=	$self->configfile();
	my $label		=	$self->label();
	my $key			=	$self->key();

	my $filetype	=	$self->filetype();
	my $jbrowse		=	$self->jbrowse();
	my $species		=	$self->species();
	my $build		=	$self->build();



	#### CHECK INPUTS
	die "inputdir not defined (Use --help for usage)\n" and exit if not defined $inputdir;
	die "featuresdir not defined (Use --help for usage)\n" and exit if not defined $featuresdir;
	die "refseqfile not defined (Use --help for usage)\n" and exit if not defined $refseqfile;
	die "filetype not defined (Use --help for usage)\n" and exit if not defined $filetype;
	die "jbrowse not defined (Use --help for usage)\n" and exit if not defined $jbrowse;
	die "species not defined (Use --help for usage)\n" and exit if not defined $species;
	die "build not defined (Use --help for usage)\n" and exit if not defined $build;

	#### PROCESS CONFIG FILE
	$self->processConfigfile();

	#### CHECK DIRS
	print "Can't find featuresdir: $featuresdir\n" and exit if not -d $featuresdir;

	#### 1. CREATE DIRECTORIES FOR OUTPUT FILES
	####  	AND COPY refSeqs.js TO EACH SUB DIRECTORY
	#$self->createFeaturedir();

	#### 2. RUN flatfile-to-json.pl IN PARALLEL IF MULTIPLE FILES PRESENT
	$self->runFlatfileToJson();

	print "Agua::JBrowse::generateFeatures    COMPLETED\n";
}


=head2

	SUBROUTINE	runFlatfileToJson

	PURPOSE

		1. CREATE DIRECTORIES FOR OUTPUT FILES

		2. COPY refSeqs.js TO EACH SUB DIRECTORY

=cut

method createFeaturedir (Str $feature) {

	#### GET OUTPUT DIR AND REFSEQ FILE
	my $featuresdir = $self->featuresdir();
	my $refseqfile = $self->refseqfile();

	#### CHECK INPUTS
	print "Agua::JBrowse::createFeaturedir    Can't find featuresdir: $featuresdir\n" and exit if not -d $featuresdir;
	print "Agua::JBrowse::createFeaturedir    Can't find refseqfile: $refseqfile\n" and exit if not -f $refseqfile;
	print "Agua::JBrowse::createFeaturedir    featuresdir already exists as a file: $featuresdir\n" and exit if -f $featuresdir;


	#### 1. CREATE OUTPUT DIRECTORY 
	my $featuredir = "$featuresdir/$feature";
	File::Path::mkpath($featuredir);
	print "Agua::JBrowse::createFeaturedir    Can't create featuredir: $featuredir\n" and exit if not -d $featuredir;

	#### 2. CREATE data DIRECTORY
	my $datadir = "$featuredir/data";
	print "Agua::JBrowse::createFeaturedir    datadir already exists as a file: $datadir\n" and exit if -f $datadir;
	File::Path::mkpath($datadir);
	print "Agua::JBrowse::createFeaturedir    Can't create datadir: $datadir\n" and exit if not -d $datadir;

	#### 3. COPY refSeqs.js FILE TO THE data DIRECTORY
	File::Copy::copy($refseqfile, "$datadir/refSeqs.js") or die "Can't copy refseqfile: $refseqfile\n";

	###### 4. CHANGE TO OUTPUT DIR
	#chdir($featuresdir) or die "Can't change to output dir: $featuresdir\n";

	return $featuredir;
}


=head2

	SUBROUTINE	runFlatfileToJson

	PURPOSE

		RUN flatfile-to-json FOR EVERY REFERENCE SEQUENCE

		AGAINST EVERY inputdir/reference/filename INPUT FILE

=cut

method runFlatfileToJson () {


	#### GET INPUT AND OUTPUT DIRS
	my $inputdir = $self->inputdir();
	my $refseqfile = $self->refseqfile();
	my $configfile = $self->configfile();

	#### GET FILETYPE, FILENAME, LABEL AND KEY
	my $filename = $self->filename();
	my $filetype = $self->filetype();
	my $label = $self->label();
	my $key = $self->key();

	#### SET EXECUTABLE
	my $jbrowse = $self->jbrowse();
	my $executable = "$jbrowse/flatfile-to-json.pl";

	#### GET LIST OF REFERENCES
	my $references = $self->parseReferences();
	print "Agua::JBrowse::runFlatfileToJson    references: @$references\n";

	#### DO ALL FILES SPECIFIED IN FEATURES HASH
	my $jobs = [];

	my $registered;
	if ( defined $filename )
	{
		my ($feature) = $filename =~ /^(.+)(\.[^\.]+)$/;
		print "Agua::JBrowse::runFlatfileToJson    filename: $filename\n";
		print "Agua::JBrowse::runFlatfileToJson    feature: $feature\n";
		$self->registerFeature($feature);

		foreach my $reference ( @$references )
		{
			next if $reference =~ /^\./;

			my $label = "generateFeatures-$reference-$feature";
			my $outdir = $self->createFeaturedir($feature);
			my $infile = "$inputdir/$reference/$filename";

			$jobs = $self->addJob($jobs, $label, $executable, $infile, $outdir);
		}
	}
	else
	{

		my $registered;
		foreach my $reference ( @$references )
		{
			next if $reference =~ /^\./;

			opendir(DIR, "$inputdir/$reference") or die "Can't open inputdir/reference: $inputdir/$reference\n";
			my @files = readdir(DIR);
			close(DIR);
			foreach my $file ( @files )
			{
				next if $file =~ /^\./;	
				print "Agua::JBrowse::runFlatfileToJson    file: $file\n";

				my $infile = "$inputdir/$reference/$file";
				print "Agua::JBrowse::runFlatfileToJson    infile: $infile\n";
				next if not -f $infile;

				#### REGISTER FEATURE
				my ($feature) = $file =~ /^(.+)(\.[^\.]+)$/;
				print "Agua::JBrowse::runFlatfileToJson    file: $file\n";
				print "Agua::JBrowse::runFlatfileToJson    feature: $feature\n";

				if ( not exists $registered->{$feature} )
				{
					$self->registerFeature($feature);
				}
				else
				{
					$registered->{$feature} = 1;
				}


				#### ADD JOB
				my $label = "generateFeatures-$reference-$feature";
				my $outdir = $self->createFeaturedir($feature);

				$jobs = $self->addJob($jobs, $label, $executable, $infile, $outdir);
			}
		}
	}
	print "Agua::JBrowse::runFlatfileToJson    No. jobs: ", scalar(@$jobs), "\n";

#exit;

	#### RUN COMMANDS IN PARALLEL
	$self->runJobs( $jobs, "generateFeatures" );	
}

method registerFeature($feature) {
	#### OVERRIDE THIS METHOD IN 	
}

=head2

	SUBROUTINE		addJob

	PURPOSE

		PUSH JBROWSE PROCESSING JOB ONTO THE jobs ARRAY

	NOTES

		flatfile-to-json.pl USAGE:

		USAGE: $0 [--gff <gff3 file> | --gff2 <gff2 file> | --bed <bed file>] [--out <output directory>] --tracklabel <track identifier> --key <human-readable track name> [--cssClass <CSS class for displaying features>] [--autocomplete none|label|alias|all] [--getType] [--getPhase] [--getSubs] [--getLabel] [--urltemplate "http://example.com/idlookup?id={id}"] [--extraData <attribute>] [--subfeatureClasses <JSON-syntax subfeature class map>] [--clientConfig <JSON-syntax extra configuration for FeatureTrack>]

			--out: defaults to "data"
			--cssClass: defaults to "feature"
			--autocomplete: make these features searchable by their "label", by their "alias"es, both ("all"), or "none" (default).
			--getType: include the type of the features in the json
			--getPhase: include the phase of the features in the json
			--getSubs:  include subfeatures in the json
			--getLabel: include a label for the features in the json
			--urltemplate: template for a URL that clicking on a feature will navigate to
			--arrowheadClass: CSS class for arrowheads
			--subfeatureClasses: CSS classes for each subfeature type, in JSON syntax
				e.g. '{"CDS": "transcript-CDS", "exon": "transcript-exon"}'
			--clientConfig: extra configuration for the client, in JSON syntax
				e.g. '{"featureCss": "background-color: #668; height: 8px;", "histScale": 2}'
			--type: only process features of the given type
			--nclChunk: NCList chunk size; if you get "json text or perl structure exceeds maximum nesting level" errors, try setting this lower (default: $nclChunk)
			--extraData: a map of feature attribute names to perl subs that extract information from the feature object
				e.g. '{"protein_id" : "sub {shift->attributes(\"protein_id\");} "}'
			--compress: compress the output (requires some web server configuration)
			--sortMem: the amount of memory in bytes to use for sorting

=cut

method addJob ($jobs, $label, $executable, $inputfile, $outputdir) {

	print "Agua::JBrowse::addJob    label: $label\n";
	print "Agua::JBrowse::addJob    inputfile: $inputfile\n";


	my $filetype = $self->filetype();	
	return if not -f $inputfile or $inputfile !~ /\.$filetype$/;

	my ($feature) = $inputfile =~ /([^\/]+)$/;
	$feature =~ s/\.[^\/]{3,5}$//;
	print "Agua::JBrowse::addJob    feature: $feature\n";
	my $config = $self->getFeatureConfig($feature);
#exit;

	my $command = "cd $outputdir; ";
	$command .= " $executable ";
	$command .= " --" . $self->filetype . " $inputfile ";
	$command .= " --tracklabel $feature ";
	$command .= " --key $feature ";
	$command .= " --autocomplete all ";
	$command .= " --getType ";
	$command .= " --nclChunk " 	. $self->chunksize() if $self->chunksize();
	$command .= " --compress " 	. $self->compress() if $self->compress();
	$command .= " --sortMem " 	. $self->sortmem() if $self->sortmem();

	if ( defined $config )
	{
		my $flags = ["--getPhase",  "--getSubs"];
		my $flag_slots = ["phase",  "subfeatures"];
		for ( my $i = 0; $i < @$flags; $i++ )
		{
			$command .= " $$flags[$i] "
				if defined $config->{$$flag_slots[$i]};
		}

		my $params = ["--cssClass", "--urltemplate", "--arrowheadClass"];
		my $param_slots  = [ "class", "urlTemplate", "arrowheadClass" ];
		for ( my $i = 0; $i < @$params; $i++ )
		{
			$command .= " $$params[$i] $config->{$$param_slots[$i]} "
				if defined $config->{$$param_slots[$i]};
		}

		my $objects = [ "subfeatureClasses",  "--clientConfig",  "--extraData"];
		my $object_slots = [ "subfeature_classes",  "clientConfig",  "extraData"];
		for ( my $i = 0; $i < @$objects; $i++ )
		{
			print "Agua::JBrowse::addJob    config->$$object_slots[$i]: \n" if defined $config->{$$object_slots[$i]};
			print Dumper $config->{$$object_slots[$i]} if defined $config->{$$object_slots[$i]};
			if ( defined $config->{$$object_slots[$i]} )
			{
				$command .= " $$objects[$i] \'";
				$command .= $self->parseJson(($config->{$$object_slots[$i]}));
				$command .= "\' ";
			}
		}	
	}

#exit;


	my $job = $self->setJob( [$command], $label, $outputdir );
	push @$jobs, $job;

	return $jobs;
}

#### JSON METHODS
method processConfigfile () {

	my $configfile	=	$self->configfile();
	return if not defined $configfile or not $configfile;

	print "Agua::JBrowse::processConfigfile    configfile: $configfile\n";
	my $config = $self->readJson($configfile, 0);
	print "Agua::JBrowse::processConfigfile    Can't process config from configfile: $configfile\n" and exit if not defined $config or not $config;


	$self->config($config);
}


method getFeatureConfig ($feature) {

	my $config	=	$self->config();


	return if not defined $config or not $config;

	foreach my $track ( @{$config->{tracks}} )
	{
		return $track if $track->{track} eq $feature;
	}

	return;
}

=head2

	SUBROUTINE	parseReferences

	PURPOSE

		PARSE REFERENCE SEQUENCE NAMES FROM refSeqs.js JSON FILE

=cut
method parseReferences () {


	my $refseqfile = $self->refseqfile();

	#### INITIALISE JSON PARSER
	my $jsonparser = JSON->new();

	#### GET .json FILE HASH
	open(ARGS, $refseqfile) or die "Can't open args file: $refseqfile\n";
	$/ = undef;
	my $contents = <ARGS>;
	close(ARGS);
	die "Json file is empty\n" if $contents =~ /^\s*$/;
	$contents =~ s/\n//g;
	$contents =~ s/^\s*refSeqs\s*=\s*//;

	my $refseq_hashes = $jsonparser->decode($contents);


	my $references = [];
	foreach my $refseq_hash ( @$refseq_hashes )
	{
		push @$references, $refseq_hash->{name} if defined $refseq_hash->{name};
	}

	return $references;	
}


method featurehash ($json, $inputfile) {

	my $featurehash;
	my @inputtypes = keys %{$json->{features}};
	@inputtypes = sort @inputtypes;
	foreach my $inputtype ( @inputtypes )
	{

		#use re 'eval';# EVALUATE $pattern AS REGULAR EXPRESSION
		if ( $inputfile =~ /$inputtype/ )
		{
			$featurehash = $json->{features}->{$inputtype};
			last;
		}
		#no re 'eval';# EVALUATE $pattern AS REGULAR EXPRESSION
	}

	return $featurehash;
}


method setJsonParser () {

	return $self->jsonparser() if $self->jsonparser();
	#my $jsonparser = JSON->new() ;
	my $jsonparser = JSON->new->allow_nonref;
	$self->jsonparser($jsonparser);

	return $self->jsonparser();
}


=head2

	SUBROUTINE		readJson

	PURPOSE

		PARSE DATA FROM JSON FILE

		(THIS SUBROUTINE AND FOLLOWING TWO ADAPTED FROM JBROWSE JsonGenerator.pm)

	INPUTS

		1. JSON FILE LOCATION

		2. DEFAULT CONTENT IF EMPTY

		3. SKIP ASSIGN IF NO VARIABLE NAME IN FIRST LINE

	OUTPUTS

		1. FILE CONTENTS MINUS VARIABLE NAME LINE IF skipAssign SPECIFIED

=cut

method readJson ($file, $skipAssign) {


	#### CHECK INPUT
	print "Agua::JBrowse::readJson    file not defined. Exiting\n" and exit if not defined $file;
    print "Agua::JBrowse::readJson    file is empty: $file\n" and return if not -s $file;

	my $jsonparser = $self->jsonparser();
	$self->jsonparser($jsonparser);

	my $OLDSEP = $/;
	my $fh = new IO::File $file, O_RDONLY or die "couldn't open $file: $!";
	flock $fh, LOCK_SH;
	undef $/;
	my $contents = <$fh>;
	$fh->close() or die "couldn't close $file: $!";
	$/ = $OLDSEP;

	#### SKIP TEXT IN skipAssign IF DEFINED
	$contents =~ s/$skipAssign// if $skipAssign;

	my $object = $jsonparser->jsonToObj($contents);

    return $jsonparser->jsonToObj($contents);
}

method parseJson ($object) {

	my $jsonparser = $self->jsonparser();

	return $jsonparser->objToJson($object);
}

=head2

	SUBROUTINE		writeJson

	PURPOSE

		1. PRINT DATA TO A JSON FILE

		2. PRINT IN 'PRETTY' FORMAT BY DEFAULT

=cut

method writeJson ($file, $toWrite, $assignation) {

	#### CHECK INPUT
	print "Agua::JBrowse::writeJson    file not defined. Exiting\n"
		and exit if not defined $file;

    #### CREATE JSON OBJECT
    my $jsonparser = $self->jsonparser();
	my $json = '';
	$json .= $assignation if defined $assignation;
	$json .= $jsonparser->objToJson($toWrite, {pretty => 1, indent => 2});

    #### WRITE JSON TO FILE
    my $fh = new IO::File $file, O_WRONLY | O_CREAT or die "couldn't open $file: $!";
    flock $fh, LOCK_EX;
    $fh->seek(0, SEEK_SET);
    $fh->truncate(0);
	$fh->print($json);
    $fh->close() or die "couldn't close $file: $!";

}

=head2

	SUBROUTINE		modifyJsonfile

	PURPOSE

		1. MODIFY THE DATA IN AN EXISTING JSON FILE OR WRITE A NEW FILE

		2. USE A CALLBACK SUBROUTINE WITH THE EXISTING DATA PROVIDED AS

			AN ARGUMENT

		3. ASSIGN VARIABLE NAME IN FIRST LINE OF FILE, E.G.: "refSeqs =\n"		

	INPUTS

		1. INPUT JSON FILE WITH VARIABLE NAME ASSIGNED IN FIRST LINE

		2. CALLBACK FUNCTION FOR MODIFYING EXISTING DATA, E.G., 'ADD

			AN ELEMENT TO A HASH ARRAY'

		3. VARIABLE NAME FOR THE DATA IN THE JSON FILE

	OUTPUTS

		1. JSON FILE WITH MODIFIED VARIABLE NAME AND DATA

=cut

method modifyJsonfile ($file, $varName, $callback) {
    my ($data, $assign);
    my $fh = new IO::File $file, O_RDWR | O_CREAT or die "couldn't open $file: $!";
    flock $fh, LOCK_EX;

    #### GET DATA AND REMOVE INITIAL VARIABLE LINE IF FILE IS NOT EMPTY
    if (($fh->stat())[7] > 0)
	{
        $assign = $fh->getline();
        my $jsonString = join("", $fh->getlines());
        $data = JSON::from_json($jsonString) if (length($jsonString) > 0);

        #### PREPARE TO WRITE OVER ANY EXISTING FILE DATA
        $fh->seek(0, SEEK_SET);
        $fh->truncate(0);
    }

    #### add assignment line
    $fh->print("$varName = ");

    #### modify data, write back
    $fh->print(JSON::to_json($callback->($data), {pretty => 1}));
    $fh->close() or die "couldn't close $file: $!";
}





}	#### END


