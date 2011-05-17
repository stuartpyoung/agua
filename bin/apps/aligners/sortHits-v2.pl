
#http://www.perlmonks.org/?node_id=53920


my $maxrecs= 8*1024*1024; # Whatever you determine fits
my $recsize= 8;
my $sortsize= 2;
my $sortoff= 6;
# Here is the only memory hog:
my $sorton= " "x($maxrecs*$sortsize);
my $idx= 0;
# Note that I don't use sysread() here as I think the
# buffering offered by read() may improve speed:
while(  $idx < $maxrecs  &&  read(FILE,$rec,$recsize)  ) {
    substr( $sorton, $idx++*$sortsize, $sortsize )=
      substr( $rec, $sortoff, $sortsize );
}
my @idx= sort { substr($sorton,$a*$sortsize,$sortsize)
            cmp substr($sorton,$b*$sortsize,$sortsize)
         } 0..($idx-1);
for $idx (  @idx  ) {
    seek( FILE, $idx*$recsize, 0 );
    sysread( FILE, $rec, $recsize );
    print OUT, $rec; # or substr($rec,0,6)
}