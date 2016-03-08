#!/usr/bin/perl

# nwayfilter.pl -- Does an n-way comparison of SNPs, assuming a given
# 1-way complexity

# usage: ./snpchip2mdrfilter.pl [-c <complexity>] [-n <n-way>] <marker statistics file>

# Marker statistics file should be in the following format:
# <marker> ... ... ... <value>
# (i.e. some number of columns, first column is the marker name, last is the
# value associated with that marker)

# The output will be of the following form:
# (<marker>){n} <value>
# (i.e. n+1 columns, with n markers in the first n columns, followed
# by the value associated with that set of markers)

# Note: Due to the way that MDR currently calculates these statistics,
# n > ~4 is unlikely to return useful results.

# Author: David Eccles (gringer), 2009 <programming@gringer.org>

use strict;
use warnings;

sub usage {
  print("usage: ./nwayfilter.pl [options] <marker statistics file>\n");
  print("\nOther Options:\n");
  print("-c <integer> : equivalent complexity to this many one-way interactions\n");
  print("-n <integer> : number of interactions to consider\n");
  print("-v           : be verbose about what is being done\n");
  print("\n");
}

sub factorial {
	my $n = 1;
	$n *= $_ for 2..shift;
	return $n;
}

my $complexity = 500000;
my $nway = 2;
my %markerValues = ();
my $verbose = 0; # false

while ((@ARGV > 0) && !(-f $ARGV[0])){
    my $argval = shift(@ARGV);
    if($argval eq "-v"){
        $verbose = 1; #true
    } elsif($argval eq "-c"){
        if($ARGV[0]){
            $argval = shift(@ARGV);
        } else {
            $argval = 0; # false
        }
        if ($argval && ($argval =~ /^([0-9]+)$/)){
            $complexity = $1;
        } else {
            die("Error: Complexity setting used without valid value");
        }
    } elsif($argval eq "-n"){
        if($ARGV[0]){
            $argval = shift(@ARGV);
        } else {
            $argval = 0; # false
        }
        if ($argval && ($argval =~ /^([0-9]+)$/)){
            $nway = $1;
        } else {
            die("Error: N-way setting used without valid value");
        }
    } else {
        print("Warning: Unknown option '$argval'");
    }
}

if ($nway > 17){
    die("Error: Too many combinations to consider for factorial function" .
        "-- please choose smaller n");
}

# number of possible combinations = n!/(w!(n-w)!)
# n: number of markers
# w: nway

# from this, we can derive the geometric mean (gm) of the number of
# markers for a given complexity (c) as follows:
#    c = n! / (w!(n-w)!)
#      = ((n)*(n-1)*(n-2)*...*(n-w)) / (w!)
# => c*w! = gm((n),(n-1),(n-2),...,(n-w))^w
# => root(c*w!,w) = gm((n),(n-1),(n-2),...,(n-w))

# we take this mean as the upper value for the number of markers (as
# the given complexity is treated as a maximum value, don't want to
# exceed it)

my $numMarkers = int(($complexity * factorial($nway))**(1/$nway));

printf "A %d-way test doing around %d comparisons would require about %d " .
    "markers\n", $nway, $complexity, $numMarkers if $verbose;

$| = 1; # disable line buffering on standard out
print "Reading marker information file..." if $verbose;

while(<>){
    my @line = split(/\s/, $_);
    if(((@line + 0) >= 2) && ($line[$#line] =~ /^[0-9\.]+$/)){
        $markerValues{$line[0]} = $line[$#line];
    }
}

print " done. Read in information for ".(keys(%markerValues) + 0).
    " markers\n" if $verbose;

# sort, $b before $a to produce descending values
my @sortedMarkers = sort {$markerValues{$b} <=> $markerValues{$a}}
    keys %markerValues;

if((@sortedMarkers + 0) > $numMarkers){
    splice(@sortedMarkers, $numMarkers);
} else {
    print "Fewer than $numMarkers markers in the input list... " .
        "keeping all of them\n";
    $numMarkers = (@sortedMarkers + 0);
}

print "The following is a list of the top $numMarkers markers:\n" if
    $verbose;

for (@sortedMarkers) {
    printf "%s %f\n", $_, $markerValues{$_};
}


