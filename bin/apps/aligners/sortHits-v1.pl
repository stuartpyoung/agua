#!/usr/bin/perl -w


#http://www.perlmonks.org/?node_id=53920

#constants
$record_size = 8;
$amp_index = 6;
$amp_format = 'S';
$time_length = 6;

$sorted = $big=$ARGV[0];
$sorted =~ s/(\.bin)$/_sorted$1/;

#load the data from $big
open DATA, $big;
binmode DATA;
my $amp_a = loadamp();
close DATA;

die "Could not read data from $big" unless @{$amp_a->[0]};

print "Ready with loading\n";

#sort on first item: amp
my @sortamp = sort{ $a->[0] <=> $b->[0] } @$amp_a; 

open OUT, ">$sorted";
binmode OUT;
for $amp( @sortamp ){
  print OUT $amp->[1].pack( $amp_format, $amp->[0]); 
}
close OUT;

# UNPACK THE DATA IN SORT RATHER THAN CARRYING AROUND UNPACKED DATA
#my @sortamp = sort{ unpack(..., substr(...$a)) <=> unpack(..., substr(...$b) } @$amp_a; 

# IF SORT KEY IS INTEGER AT END OF DATA
#my @sortamp = sort{ reverse($a) cmp reverse($b) } @$amp_a;


sub loadamp{
  my $str;
  my @amp_a;
  my $n = 1;
  do {
    $n = sysread( DATA, $str, $record_size);
    sysseek( DATA, $record_size, 1);
    (my $amp)=unpack( $amp_format, substr( $str, $amp_index, 2) );
    push(@amp_a, [$amp,substr( $str, 0, $time_length)]) if $n;
  } while ( $n);
  \@amp_a;
}
