#!/usr/bin/perl

# snplookup.pl -- retrieves chromosome location information from Entrez

# Author: David Eccles (gringer), 2009 <programming@gringer.org>

# expects lines in the following form:
# rsXXXXXX <other data 1>
# rsXXXXXY <other data 2>
# rsXXXXXZ <other data 2>
# outputs lines in the following form form:
# rsXXXXXX <other data 1> 5:1204395
# rsXXXXXY <other data 2> 21:492940
# rsXXXXXZ <other data 2> 6:104030423
# [ie. chromosome:location]

use strict;
use warnings;
use LWP::UserAgent;
use HTTP::Request::Common;

sub usage {
  print("usage: ./snplookup.pl < <file name>\n");
  print("\n");
}

my $ua = LWP::UserAgent->
    new(env_proxy => 1, agent => "Mozilla/5.0 (compatible; Perl; ".
	"hacking\@gringer.dis.org.nz)", timeout => "10");

while(<>){
    my $line = $_;
    chomp($line);
    if(/^rs([0-9]+)/){
	my $address = "http://eutils.ncbi.nlm.nih.gov/".
	    "entrez/eutils/efetch.fcgi?db=snp&id=".
	    $1."&report=DocSet";
	my $res = ($ua->request(GET $address))->content;
	if($res =~ /CHROMOSOME BASE POSITION=(.*)$/m){
	    print $line." ".$1."\n";
	}
	else{
	    print $line." NOT FOUND -- rs".$1."\n";
	}
    }
}
