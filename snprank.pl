#!/usr/bin/perl

# snprank.pl -- ranks markers in a bootstrap-sorted file. Expected
# input is the output of bootstrap.r, sorted first by bootstrap, then
# by marker information statistic.

# usage: ./snprank.pl <bootstrap-sorted file> > output_ranked.txt

# Author: David Eccles (gringer), 2009 <programming@gringer.org>

# potential 'easy' improvements (as command-line options):
# * allow field separator to be changed
# * identify and ignore header lines

use strict;
use warnings;

sub usage {
    print("usage: ./snprank.pl <bootstrap-sorted file> > output_ranked.txt\n");
}

my $currentBootstrap = "";
my $rank = 0;

while(<>){
    chomp;
    my ($marker, $bootstrap, $value) = split(/,/);
    if($bootstrap ne $currentBootstrap){
        $currentBootstrap = $bootstrap;
        $rank = 1;
    } else {
        $rank++;
    }
    print ($marker.",".$rank."\n");
}
