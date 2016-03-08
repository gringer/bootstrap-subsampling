#!/usr/bin/perl

# calc_gtfreq.pl -- determines per-marker genotype frequencies for a
# simplegt formatted text file.

# Author: David Eccles (gringer), 2008 <programming@gringer.org>

use strict;
use warnings;

sub usage {
  print("usage: ./calc_gtfreq.pl <file name>\n");
  print("\n");
}

printf "%-12s %5s %5s %5s %5s %5s %s\n",
    "marker", "AA/TT", "AC/TG",
    "AG/TC", "CC/GG", "--/--", "(n)";

while(<>){
    my ($marker, @genotypes) = split(/ /,$_);
    while($genotypes[0] !~ /[ACGT-][ACGT-]/){
        shift(@genotypes);
    }
    #printf "%12s", $marker;
    ## number of genotypes
    #print "(".(@genotypes)." genotypes):";
    my $numgenos = (@genotypes);
    my %counts = ("AA/TT",0,"AC/TG",0,"AG/TC",0,"CC/GG",0,"--/--",0);
    foreach my $gt (@genotypes){
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
    if($numgenos != $counts{"--/--"}){
        printf "%-12s %.3f %.3f %.3f %.3f %2.3f (%d)\n",
        $marker,
        ($counts{"AA/TT"}/($numgenos - $counts{"--/--"})),
        ($counts{"AC/TG"}/($numgenos - $counts{"--/--"})),
        ($counts{"AG/TC"}/($numgenos - $counts{"--/--"})),
        ($counts{"CC/GG"}/($numgenos - $counts{"--/--"})),
        ($counts{"--/--"}/$numgenos),
        ($numgenos - $counts{"--/--"});
    } else {
        printf "%-12s %.3f %.3f %.3f %.3f %2.3f (%d)\n",
        $marker,0,0,0,0,
        ($counts{"--/--"}/$numgenos),
        ($numgenos - $counts{"--/--"});
    }
}
