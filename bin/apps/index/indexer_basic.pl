#!/usr/local/bin/perl
require 5;
use DB_File;    # Access DB databases
use Fcntl;      # Needed for above...
use File::Find; # Directory searching
undef $/; # Don't obey line boundaries
$currentKey = 0;

#
# Single database version:
#
#  Stores file entries in index.db as <NULL><binary file number>
#  The leading NULL prevents any word entries from colliding.
#

############################################################################

# Delete old index.db and attach %indexdb to database
unlink("index.db");
tie(%indexdb,'DB_File',"index.db",
    O_RDWR | O_CREAT, 0644, $DB_File::DB_BTREE) ;

# Index all of the files
find(\&IndexFile,"/mnt/tmp/ddjcd6-0.4/articles/1998");

# Close databases
untie(%indexdb); # release database

###########################################################################

sub IndexFile {
    if(!-f) { return; }

    if(/\.html?$/) { # Handle HTML files
	print "$File::Find::name\n";
	open(HTML_FILE,$_);
	my($text) = <HTML_FILE>; # Read entire file
	$text =~ s/<[^>]*>//g; # Strip out all HTML tags
	# Index all the words under the current key
	my($wordsIndexed) = &IndexWords($text,$currentKey);
	# Map key to this filename
	$indexdb{pack"xn",$currentKey} = $File::Find::name;
	$currentKey++;
    }
}

###########################################################################

sub IndexWords {
    my($words, $fileKey) = @_;
    my(%worduniq); # for unique-ifying word list
    # Split text into Array of words
    my(@words) = split(/[^a-zA-Z0-9\xc0-\xff\+\/\_]+/, lc $words);
    @words = grep { $worduniq{$_}++ == 0 } # Remove duplicates
	     grep { s/^[^a-zA-Z0-9\xc0-\xff]+//; $_ } # Strip leading punct
             grep { length > 1 } # Must be longer than one character
             grep { /[a-zA-Z0-9\xc0-\xff]/ } # must have an alphanumeric
             @words;

    # For each word, add key to word database
    foreach (sort @words) {
	my($a) = $indexdb{$_};
	$a .= pack "n",$fileKey;
        $indexdb{$_} = $a;
    }

    # Return count of words indexed
    return scalar(@words);
}

###########################################################################
