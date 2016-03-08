#!/usr/bin/perl

# linecount.pl -- print out the line at which a marker has been seen a
# particular number of times. The marker is assumed to be in the first
# field of each line.

# usage: ./linecount.pl <sorted, ranked file> > linecount_output.txt

# Author: David Eccles (gringer), 2009 <programming@gringer.org>

# potential 'easy' improvements (as command-line options):
# * allow field separator to be changed
# * change number of times (for more bootstraps, or less strict
#   consistency)

sub usage {
    print("usage: ./linecount.pl <sorted, ranked file> ".
          "> linecount_output.txt\n");
}

use warnings;
use strict;

my %seen = ();

my $line = 0;

my $times = 100;

while(<>){
    $line++;
    my ($marker, @rest) = split(/,/);
    $seen{$marker}++;
    if($seen{$marker} == $times) {
        print($marker . "," . $line . "\n");
    }
}
