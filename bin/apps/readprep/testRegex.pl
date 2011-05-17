#!/usr/bin/perl -w

my $input = "\@HWI-EAS185_6_20GVYAAXX_Jia_cDNA2_mtDNA_Total_Small_Medium_read2_JH:6:1:589:870/1";

if ( $input =~ /^(@[A-Z0-9\_\-\.]+):([0-9]+:[0-9]+:[0-9]+:[0-9]+)/i 
    or $input =~ /^(@[A-Z0-9\_\-\.]+):([0-9]+:[0-9]+:[0-9]+:[0-9]+#[0-9]+\/(1|2))\s*/i )
{
    print "matched\n";
}
else
{
    print "Not matched\n";
}
