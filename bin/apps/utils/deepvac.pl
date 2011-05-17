#!/usr/bin/perl -w

#### DEBUG

#### TIME
my $time = time();

=head2

    APPLICATION     deepvac

    PURPOSE

        PERFORM 'DEEP VACCUUM' DOWNLOADS OF FTP SITES

        BY DOWNLOADING BY wget ALL LINKED FILES:

            1. FROM HTML CONTENT /<A [^>]* href=(\S+)/

            2. FROM TEXT CONTENT

    INPUT

        1. REQUIRED: URL AND OUTPUT DIRECTORY (MUST ALREADY EXIST)

        TO DO:

        2. OPTIONAL: FILTER TERM AND REGEX TO FILTER FILES BEFORE DOWNLOADING

    OUTPUT

        1. FILES DOWNLOADED BY WGET TO OUTPUT DIRECTORYBLY 1

    USAGE

    ./deepvac.pl <--url String> <--outputdir String> [--filter String] [--regex String] [--help]

    --url               :   URL OF FTP SITE OR PAGE
    --outputdir       :   DIRECTORY TO DOWNLOAD FILES TO
    --filter            :   Download only files containing this text fragment
    --regex             :   Download only files matching this regular expression
    --help              :   PRINT HELP INFO

    EXAMPLES

./deepvac.pl --url ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/ --outputdir /home/syoung/base/pipeline/exome/ccds

    Run time: 00:00:37
    Completed ./deepvac.pl
    10:40PM, 6 January 2009
    ****************************************


/nethome/bioinfo/apps/agua/0.5/bin/apps/utils/deepvac.pl \
--url ftp://ftp.ncbi.nlm.nih.gov/sra/static/SRX016/SRX016231 \
--outputdir /scratch/syoung/base/pipeline/SRA/NA18507/SRP000239/SRX016231 \
--type text \
--regex "_(1|2)\.fastq"



	NOTES

		THIS EXPECTS FTP PAGES IN THE FORMAT:

			All the files and tables in this directory are freely usable for any purpose.
		</pre>
		<pre>
			<a href="?C=N;O=D">Name</a>                                   <a href="?C=M;O=A">Last modified</a>\
			<a href="?C=S;O=A">Size</a>  <a href="?C=D;O=A">Description</a><hr>      <a href="/goldenPath/hg19/">Pa\
			rent Directory</a>                                            -
			<a href="all_bacends.sql">all_bacends.sql</a>                        24-May-2009 11:40  2.1K
			<a href="all_bacends.txt.gz">all_bacends.txt.gz</a>                     24-May-2009 11:40   93M
			<a href="all_est.sql">all_est.sql</a>                            21-Feb-2010 14:49  2.3K
			<a href="all_est.txt.gz">all_est.txt.gz</a>                         21-Feb-2010 14:49  374M

=cut

use strict;

#### USE LIB
use FindBin qw($Bin);
use lib "$Bin/../../../lib";

#### INTERNAL MODULES
use Timer;
use Util;
use Conf::Agua;

#### EXTERNAL MODULES
use File::Path;
use LWP::Simple;
use Term::ANSIColor qw(:constants);
use Getopt::Long;

#### GET OPTIONS
my $url;
my $outputdir;
my $filter;
my $regex;
my $type;
my $help;
GetOptions (
	'url=s' => \$url,
	'outputdir=s' => \$outputdir,
	'filter=s' => \$filter,
	'regex=s' => \$regex,
	'type=s' => \$type,
	'help' => \$help
) or die "No options specified. Try '--help'\n";

#### PRINT HELP IF REQUESTED
if ( defined $help )	{	usage();	}

#### TYPE OF URL IS 'html' OR 'text'
print "Type not defined (option --url)\n" and exit if not defined $type;
print "Type must be 'html' or 'text'\n" and exit if $type !~ /^(html|text)$/i;

#### CHECK FOR REQUIRED URL
die "Url not defined (option --url)\n" if not defined $url;

#### CHECK FOR REQUIRED OUTPUT DIRECTORY
die "Output directory not defined (option --outputdir)\n" if not defined $outputdir;

#### CREATE DIRECTORY IF NOT EXISTS
File::Path::mkpath($outputdir) if not -d $outputdir;
die "Can't create output directory: $outputdir" if not -d $outputdir;

#### GET CONTENTS OF PAGE
$url =~ s/\/$//;
print "Doing get($url)\n";
my $contents = get($url);

#### PARSE OUT FILENAMES FROM PAGE
my $files;
if ( $type eq "text" )
{
    my @lines = split "\n", $contents;
    foreach my $line ( @lines )
    {
        if ( $line =~ /(\S+)\s*$/ )
        {
            push @$files, "$url/$1";
        }
    }
}
else
{
    my @lines = split "\n", $contents;
    foreach my $line ( @lines )
    {
		if ( $line =~ /<a\s*[^>]*\s*href\s*=\s*('|")(\S+)('|")/ )
		{

			push @$files, "$url/$2";
		}
    }
}

#### CHANGE TO OUTPUT DIRECTORY
chdir($outputdir) or die "Can't change to download directory: $outputdir\n";
foreach my $file ( @$files )
{
	next if defined $regex and not $file =~ /$regex/;

    my $wget = "wget $file";
    print "wget: $wget\n";

    `$wget`;
}

#### PRINT RUN TIME
my $runtime = Timer::runtime( $time, time() );
print "\nRun time: $runtime\n";
print "Completed $0\n";
print Util::datetime(), "\n";
print "****************************************\n\n\n";
exit;




###########################################################################################
####################                 S U B R O U T I N E S                ################# 
###########################################################################################


# PROCESS HTML FROM ONE PAGE
sub html2text
{
    my ($HTML_content, $symbol_pair, $ascii, $html);
    ($HTML_content) = @_;


    $HTML_content =~ s/&nbsp;/ /g;
    $HTML_content =~ s/\s\s*/ /g;
    $HTML_content =~ s/<p[^>]*>/\n\n/gi;   #<p>  -> \n\n
    $HTML_content =~ s/<br>|<\/*h[1-6][^>]*>|<li[^>]*>|<dt[^>]*>|<dd[^>]*>|<\/tr[^>]*>/\n/gi; 
    # <br> or <H*> or <li> or </tr> or <dt> or <dd> -> \n
    $HTML_content =~ s/(<[^>]*>)*//g;
    $HTML_content =~ s/\n\s*\n\s*/\n\n/g;
    $HTML_content =~ s/\n */\n/g;

    foreach $symbol_pair(&parseSymbols)
    {
        ($ascii, $html) = split(/\s\s*/,$symbol_pair);
        $HTML_content =~ s/$html/$ascii/g;
    }
    return $HTML_content;
}

# PARSE HTML SYMBOLS
sub parseSymbols
{
    return (
    "&	&amp;",
    "\"	&quot;",
    "<	&lt;",
    ">	&gt;",
    "©	&copy;",
    "®	&reg;",
    "Æ	&AElig;",
    "Á	&Aacute;",
    "Â	&Acirc;",
    "À	&Agrave;",
    "Å	&Aring;",
    "Ã	&Atilde;",
    "Ä	&Auml;",
    "Ç	&Ccedil;",
    "Ð	&ETH;",
    "É	&Eacute;",
    "Ê	&Ecirc;",
    "È	&Egrave;",
    "Ë	&Euml;",
    "Í	&Iacute;",
    "Î	&Icirc;",
    "Ì	&Igrave;",
    "Ï	&Iuml;",
    "Ñ	&Ntilde;",
    "Ó	&Oacute;",
    "Ô	&Ocirc;",
    "Ò	&Ograve;",
    "Ø	&Oslash;",
    "Õ	&Otilde;",
    "Ö	&Ouml;",
    "Þ	&THORN;",
    "Ú	&Uacute;",
    "Û	&Ucirc;",
    "Ù	&Ugrave;",
    "Ü	&Uuml;",
    "Ý	&Yacute;",
    "á	&aacute;",
    "â	&acirc;",
    "æ	&aelig;",
    "à	&agrave;",
    "å	&aring;",
    "ã	&atilde;",
    "ä	&auml;",
    "ç	&ccedil;",
    "é	&eacute;",
    "ê	&ecirc;",
    "è	&egrave;",
    "ð	&eth;",
    "ë	&euml;",
    "í	&iacute;",
    "î	&icirc;",
    "ì	&igrave;",
    "ï	&iuml;",
    "ñ	&ntilde;",
    "ó	&oacute;",
    "ô	&ocirc;",
    "ò	&ograve;",
    "ø	&oslash;",
    "õ	&otilde;",
    "ö	&ouml;",
    "ß	&szlig;",
    "þ	&thorn;",
    "ú	&uacute;",
    "û	&ucirc;",
    "ù	&ugrave;",
    "ü	&uuml;",
    "ý	&yacute;",
    "ÿ	&yuml;",
    " 	&#160;",
    "¡	&#161;",
    "¢	&#162;",
    "£	&#163;",
    "¥	&#165;",
    "¦	&#166;",
    "§	&#167;",
    "¨	&#168;",
    "©	&#169;",
    "ª	&#170;",
    "«	&#171;",
    "¬	&#172;",
    "­	&#173;",
    "®	&#174;",
    "¯	&#175;",
    "°	&#176;",
    "±	&#177;",
    "²	&#178;",
    "³	&#179;",
    "´	&#180;",
    "µ	&#181;",
    "¶	&#182;",
    "·	&#183;",
    "¸	&#184;",
    "¹	&#185;",
    "º	&#186;",
    "»	&#187;",
    "¼	&#188;",
    "½	&#189;",
    "¾	&#190;",
    "¿	&#191;",
    "×	&#215;",
    "Þ	&#222;",
    "÷	&#247;")
}


sub usage
{
	print GREEN <<"EOF";

    APPLICATION     deepvac

    PURPOSE

        PERFORM 'DEEP VACCUUM' DOWNLOADS OF FTP SITES

        BY DOWNLOADING BY wget ALL LINKED FILES:

            1. FROM HTML CONTENT /<A [^>]* href=(\\S+)/

            2. FROM TEXT CONTENT

    INPUT

        1. REQUIRED: URL AND OUTPUT DIRECTORY (MUST ALREADY EXIST)

        TO DO:

        2. OPTIONAL: FILTER TERM AND REGEX TO FILTER FILES BEFORE DOWNLOADING

    OUTPUT

        1. FILES DOWNLOADED BY WGET TO OUTPUT DIRECTORYBLY 1

    USAGE

    ./deepvac.pl <--url String> <--outputdir String> [--filter String] [--regex String] [--help]

    --url               :   URL OF FTP SITE OR PAGE
    --outputdir         :   DIRECTORY TO DOWNLOAD FILES TO
    --filter            :   DOWNLOAD ONLY FILENAMES CONTAINING THIS FILTER
    --regex             :   DOWNLOAD ONLY FILENAMES MATCHING THIS REGEX
    --help              :   PRINT HELP INFO

    EXAMPLES

 ./deepvac.pl --url ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/ --outputdir /home/syoung/base/pipeline/exome/ccds

    Run time: 00:00:37
    Completed ./deepvac.pl
    10:40PM, 6 January 2009
    ****************************************

=cut

EOF

	print RESET;

	exit;
}




