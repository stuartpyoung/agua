#!/usr/bin/perl -w

use FindBin qw($Bin);
use lib "$Bin/../../../lib";

use TestMooseDeclare;

my $checkingaccount = CheckingAccount->new();
print "balance: ";
print $checkingaccount->balance;
print "\n";
print "depositing...";
$checkingaccount->deposit(100);
print "\n";
print "balance: ";
print $checkingaccount->balance;
print "\n";
print "depositing...";
$checkingaccount->deposit(100);
print "\n";
print "balance: ";
print $checkingaccount->balance;
print "\n";