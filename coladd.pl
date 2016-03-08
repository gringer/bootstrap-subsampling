#!/usr/bin/perl

# coladd.pl -- adds the columns of an input file
# (second column onwards, starting from the second line)

# Author: David Eccles (gringer), 2009 <programming@gringer.org>

use strict;
use warnings;

sub usage {
  print("usage: ./coladd.pl <file name>\n");
  print("\n");
}

#my $dummy = 0;
my $nw = '%-0.7f '; #number width number formatting
my $sw = '%-9s '; #number width string formatting

#my @colnames = ();

my ($dummy, @colnames) = split(/\s+/,<>);

my @colcounts = ();
for (my $i=0; $i < (@colnames); $i++){
    $colcounts[$i] = 0;
    printf $sw, $colnames[$i];
}

print "\n";

my $numrows = 0;
my @columns = ();

while(<>){
    ($dummy, @columns) = split(/\s+/,$_);
    for (my $i=0; $i < (@colnames); $i++){
        $colcounts[$i] += $columns[$i];
    }
    $numrows++;
}

for (my $i=0; $i < (@colnames); $i++){
    $colcounts[$i] = $colcounts[$i] / $numrows;
    printf $nw, $colcounts[$i];
}

print "\n";
