#!/usr/bin/perl -w
use strict;

$DEBUG = 1;

=head2

APPLICATION     archive

PURPOSE

    1. ARCHIVE AND PACKAGE A GIT REPOSITORY FOR DISTRIBUTION

INPUT

    1. BASE DIR OF git REPOSITORY

    2. OUTPUTDIR FOR PACKAGE


OUTPUT

    MYSQL DATABASE CONFIGURATION AND EDITED CONFIG FILE            

USAGE

sudo ./archive.pl \
    <--name String> <--basedir String> <--version String> \
    <--versionfile String> <--outputdir String> [--help]

--name      :   Name of application
--version   :	Version number (e.g., 0.6, 0.8)
--basedir   :	Name of basedir where the .git directory is located
--outputdir :	Create packages inside RELEASE dir in the outputdir
--help      :   Print help info

EXAMPLES

./archive.pl \
--name agua \
--basedir /home/syoung/0.6 \
--version 0.6 \
--versionfile bin/scripts/resources/VERSION \
--outputdir /home/syoung


=cut

#### FLUSH BUFFER
$| = 1;

#### USE LIB
use FindBin qw($Bin);

use lib "$Bin/../../lib";

#### EXTERNAL MODULES
use Getopt::Long;
use Data::Dumper;

#### INTERNAL MODULES
use Agua::Configure;
use Agua::DBaseFactory;
use Conf;

#### GET OPTIONS
my $name;
my $basedir;
my $version;
my $versionfile;
my $outputdir;
my $help;
GetOptions (
    'name=s'        => \$name,
    'basedir=s'     => \$basedir,
    'version=s'     => \$version,
    'versionfile=s' => \$versionfile,
    'outputdir=s'   => \$outputdir,
    'help'          => \$help
) or die "No options specified. Try '--help'\n";

print "archive.pl    name not defined\n" and exit if not defined $name;
print "archive.pl    basedir not defined\n" and exit if not defined $basedir;
print "archive.pl    version not defined\n" and exit if not defined $version;
print "archive.pl    versionfile not defined\n" and exit if not defined $versionfile;
print "archive.pl    outputdir not defined\n" and exit if not defined $outputdir;

#### 1. GET THE COMMIT COUNT
chdir($basedir) or die "Can't chdir to basedir: $basedir\n";
my $iteration = `git log --pretty=format:'' | wc -l`;
$iteration =~ s/\s+//g;
print "archive.pl    iteration: $iteration\n";

#### 2. GET THE SHORT SHA KEY AS THE BUILD ID
my $build = `git rev-parse --short HEAD`;
$build =~ s/\s+//g;
print "archive.pl    build: $build\n";

#### 3. CREATE THE RELEASE DIR AND VERSION SUBDIR
my $versiondir = "$outputdir/$version";
`mkdir -p $versiondir`;
print "archive.pl    Can't create versiondir: $versiondir\n"
    and exit if not -d $versiondir;

#### 4. UPDATE THE VERSION FILE
open(OUT, ">$versionfile") or die "Can't open versionfile: $versionfile\n";
print OUT "$version-$iteration-$build";
close(OUT);

#### 5. CREATE PACKAGE
my $archive = "git archive --format=tar --prefix=$version-$iteration-$build/ HEAD | gzip > $outputdir/$version/agua.$version-$iteration-$build.tar.gz";
print "$archive\n";
print `$archive`;

print "Completed $0\n";
