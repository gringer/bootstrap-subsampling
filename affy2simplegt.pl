#!/usr/bin/perl

# affy2simplegt.pl -- Converts data from the affymetrix chip format
# (marker, individual, genotype, QC) to simplegt format.

# Author: David Eccles (gringer), 2007 <programming@gringer.org>

sub usage {
  print("usage: ./affy2simplegt.pl <file name>\n");
  print("\n");
}

use strict;
use warnings;

my %data = ();

while(<>){
    my ($marker, $ind, $gt, $qc) = split("\\s");
    if(!$data{$marker}){
        $data{$marker} = {};
    }
    ${%data}{$marker}{$ind} = $gt;
#   print "$marker $ind $gt\n";
}

my $firstLine = 1; # true

foreach my $marker (sort keys %data){
    if($firstLine){
        print(join(" ",sort keys %{$data{$marker}})."\n");
        $firstLine = 0; # false
    }
    print $marker;
    foreach my $ind (sort keys %{$data{$marker}}){
        print " ".${%data}{$marker}{$ind};
    }
    print "\n";
}
