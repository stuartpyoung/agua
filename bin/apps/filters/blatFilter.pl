#!/usr/bin/perl

use strict;
use warnings;
#use Cwd;

##### ##### ##### ##### #####

use Getopt::Std;
use vars qw( $opt_a $opt_b $opt_c $opt_v $opt_w);

# Usage
my $usage = "

script - does something useful.

		      by
		Brian J. Knaus
		 October 2009

Copyright (c) 2009 Brian J. Knaus.
License is hereby granted for personal, academic, and non-profit use.
Commercial users should contact the author (http://brianknaus.com).

Great effort has been taken to make this software perform its said
task however, this software comes with ABSOLUTELY NO WARRANTY,
not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

Usage: perl script.pl options
 required:
  -a	fastq input file.
  -b	reference fasta file.
 optional:
  -c	inset in nucleotides, used to remove barcodes [default = 0].
  -v	verbose mode [T/F, default is F].
  -w	delete intermediate (blat) files [T/F, default is T]

";

# command line processing.
getopts('a:b:c:v:w:');
die $usage unless ($opt_a);
die $usage unless ($opt_b);

my ($inf, $ref, $inset, $verb, $clean);
$inf	= $opt_a if $opt_a;
$ref	= $opt_b if $opt_b;
$inset	= $opt_c ? $opt_c : 0;
$verb	= $opt_v ? $opt_v : "F";
$clean	= $opt_w ? $opt_w : "T";

##### ##### ##### ##### #####
# Subfunctions.

##### ##### ##### ##### #####
# Globals.

my ($temp, $in, $out, $outb, $outf, $outfb, $key, $value, $ns, $bs);
my @temp;
my %blats;

##### ##### ##### ##### #####
# Main.


# Make fasta file.

# Manage outfile name.
@temp = split(/\//, $inf);
$temp = $temp[$#temp];
@temp = split(/\./, $temp);
#$outf = join "", $temp[0], ".fa";
$outf = $temp[0];

$ns = 0; # 'N' read counter;

# Create a fasta file for input to blat.
open( $in,  "<",  $inf)  or die "Can't open $inf: $!";
open( $out,  ">",  $outf.".fa")  or die "Can't open $outf.fa: $!";

while (<$in>){
  chomp($temp[0] = $_);		# First line is an id.
  chomp($temp[1] = <$in>);	# Second line is a sequence.
  <$in>;			# Third line is an id.
  <$in>;			# Fourth line is quality.

  # Prune first char.
  $temp[0] = substr($temp[0], 1);

  # Substring to inset value.
  $temp[1] = substr($temp[1], $inset);

  if ($temp[1] =~ /[Nn\.]/){
    $ns++;
    next;
  }

  # Print to fasta file.
  print $out ">$temp[0]\n";
  print $out "$temp[1]\n";
}

close $in or die "$in: $!";
close $out or die "$out: $!";


# Call blat.

$temp = join("", "blat $ref $outf.fa ", $outf, "_blatout.psl -tileSize=6 -stepSize=1 -minScore=12 -noHead -fine -maxIntron=0");
system $temp;

$bs = 0; # Blat counter.

open( $in,  "<",  $outf.'_blatout.psl')  or die "Can't open ", $outf, "_blatout.psl: $!";

while (<$in>){
  @temp = split(/\t/);
  $blats{$temp[9]} = 0;
  $bs++;
}

close $in or die "$in: $!";

# Create sorted files.
open( $in, "<", $inf)  or die "Can't open $inf: $!";
open( $out, ">", $outf."_miss.fq")  or die "Can't open ", $outf, "_miss.fq: $!";
open( $outb, ">", $outf."_hit.fq")  or die "Can't open ", $outf, "_hit.fq: $!";

while (<$in>){
  chomp($temp[0] = $_);		# First line is an id.
  chomp($temp[1] = <$in>);	# Second line is a sequence.
  chomp($temp[2] = <$in>);	# Third line is an id.
  chomp($temp[3] = <$in>);	# Fourth line is quality.

  # Discard Ns.
  if ($temp[1] =~ /[Nn\.]/){
    next;
  }

  if (exists $blats{substr($temp[0], 1)}){
    # Print to fastq file.
    print $outb "$temp[0]\n";
    print $outb "$temp[1]\n";
    print $outb "$temp[2]\n";
    print $outb "$temp[3]\n";
  } else {
    # Print to fastq file.
    print $out "$temp[0]\n";
    print $out "$temp[1]\n";
    print $out "$temp[2]\n";
    print $out "$temp[3]\n";
  }
}

close $in or die "$in: $!";
close $out or die "$out: $!";
close $outb or die "$outb: $!";

if ($clean eq "T"){
  # Manage outfile name.
#  @temp = split(/\//, $inf);
#  $temp = $temp[$#temp];
#  @temp = split(/\./, $temp);
#  $outf = join "", $temp[0], ".fa";
#  unlink "blatout.psl", $outf;
  unlink $outf.".fa";
  unlink $outf."_blatout.psl";
}


##### ##### ##### ##### #####
# EOF.
