#!/usr/bin/perl -w

# testHasTable.pl

use lib "../../../lib";
use HasTable;
use TestTrait;

my $hasTable = TestTrait->new();
$hasTable->foo();
print "Completed $0\n";