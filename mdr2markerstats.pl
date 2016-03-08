#!/usr/bin/perl

# mdr2markerstats.pl -- retrieves marker values from MDR formatted file
# (as produced by running MDR from the command line), and produces a
# list of SNPs with mean association values

# The code searches for the string "@Landscape" (this step can be
# diabled by the -nosearch option), and while there are no more lines
# beginning with "@", reads in marker names and values

# The design of the code is flexible enough that it can take as input
# a list of markers with a number of columns of numerical values, with
# the last value being the one associated with the marker. This allows
# the code to work with file formats other than the MDR format. The
# code will automatically detect the absence of an '@' on the first
# line (which should indicate the file is *not* MDR formatted), and
# adjust the landscapeReached parameter accordingly.

# Author: David Eccles (gringer), 2009 <programming@gringer.org>

use strict;
use warnings;

sub usage {
  print("usage: ./mdr2markerstats.pl <MDR output file>\n");
  print("\nOther Options:\n");
  print("-nosearch : disable searching for Landscape string\n");
  print("\n");
}

my $landscapeReached = 0; #false
my $firstLine = 1; # true

while ((@ARGV > 0) && !(-f $ARGV[0])){
    my $argval = shift(@ARGV);
    if($argval eq "-nosearch"){
        $landscapeReached = 1; #true
    } else {
        print("Warning: Unknown option '$argval'");
    }
}

my %markerValues = ();
my %markerCounts = ();

while(<>){
    if($firstLine){
        $firstLine = 0; # false
        if(/^[^\@]/ && !$landscapeReached){
            print STDERR "'\@' not found on first line... ".
                "assuming non-MDR file\n";
            $landscapeReached = 1; # true
        }
    }
    my @line = split(/\s/, $_);
    if(!$landscapeReached){
        $landscapeReached = (@line && ($line[0] =~ /^\@Landscape/));
    } else {
        if(((@line + 0) >= 2) && ($line[$#line] =~ /^[0-9\.]+$/)){
            for(@line){
                if(/^[^0-9]/){
                    $markerValues{$_} += $line[$#line];
                    if(!$markerCounts{$_}){
                        $markerCounts{$_} = 1;
                    } else {
                        $markerCounts{$_}++;
                    }
                }
            }
        }
        if(@line && ($line[0] =~ /^\@/)){
            last;
        }
    }
}

for(keys(%markerValues)){
    printf "%s %f\n", $_ , ($markerValues{$_} / $markerCounts{$_});
}
