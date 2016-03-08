#!/usr/bin/perl

# calc_gtcounts.pl -- determines per-marker genotype counts for a
# simplegt formatted text file.

# Author: David Eccles (gringer), 2008 <programming@gringer.org>

use strict;
use warnings;

sub usage {
  print("usage: ./calc_gtcounts.pl <file name>\n");
  print("\n");
}

printf "%-12s %5s %5s %5s %5s %5s\n",
    "marker", "AA/TT", "AC/TG",
    "AG/TC", "CC/GG", "--/--";

while(<>){
    ($marker, @genotypes) = split(/\s/,$_);
    #printf "%12s", $marker;
    ## number of genotypes
    #print "(".(@genotypes)." genotypes):";
    %counts = ("AA/TT",0,"AC/TG",0,"AG/TC",0,"CC/GG",0,"--/--",0);
    foreach $gt (@genotypes){
	if ($gt =~ /(aa|tt)/i){
	    #print " aa";
	    $counts{"AA/TT"}++;
	}
	elsif ($gt =~ /(ac|ca|tg|gt)/i){
	    #print " ac";
	    $counts{"AC/TG"}++;
	}
	elsif ($gt =~ /(ag|ga|tc|ct)/i){
	    #print " ag";
	    $counts{"AG/TC"}++;
	}
	elsif ($gt =~ /(cc|gg)/i){
	    #print " cc";
	    $counts{"CC/GG"}++;
	}
	elsif ($gt =~ /../) {
	    #print " --";
	    $counts{"--/--"}++;
	}
	else{
	    print " -- ERROR -- ";
	}
    }

    printf "%-12s %-5s %-5s %-5s %-5s %-5s\n",
    $marker, $counts{"AA/TT"}, $counts{"AC/TG"},
    $counts{"AG/TC"}, $counts{"CC/GG"}, $counts{"--/--"};
}
