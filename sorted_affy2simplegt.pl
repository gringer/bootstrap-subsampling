#!/usr/bin/perl

# sorted_affy2simplegt.pl -- converts from affymetrix formatted file to
# simplegt formatted file.

# Author: David Eccles (gringer), 2007 <programming@gringer.org>

use strict;
use warnings;

sub usage {
  print("usage: ./sorted_affy2simplegt.pl <file name>\n");
  print("\n");
}

# assumes input is sorted by marker, then by individual (so once a
# different marker appears, the previous marker will never appear
# again in the file)

my %data = ();
my $oldmarker = 0; # false
my $firstLine = 1; # true
my ($marker, $ind, $gt, $qc);

while(<>){
    ($marker, $ind, $gt, $qc) = split("\\s");
    if($marker ne $oldmarker){
        unless (!$oldmarker){
            if($firstLine){
                print(join(" ",sort keys %data)."\n");
                $firstLine = 0; # false
            }
            print $oldmarker;
            foreach my $indiv (sort keys %data){
                print " ".$data{$indiv};
                $data{$indiv} = "--";
            }
            print "\n";
        }
        $oldmarker = $marker;
    }
    $data{$ind} = $gt;
}

print $oldmarker;
foreach my $indiv (sort keys %data){
    print " ".$data{$indiv};
}
print "\n";
