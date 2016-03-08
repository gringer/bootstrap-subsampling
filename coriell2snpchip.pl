#!/usr/bin/perl

# coriell2snpchip.pl -- read in (3+N)xM array (coriell .pre data file
# notation). Presumed file format is <study>[ -]<personID> 1 <GT1/1>
# <GT1/2> <GT2/1> <GT2/2> ...

# Author: David Eccles (gringer), 2009 <programming@gringer.org>

use strict;
use warnings;

sub usage {
  print("usage: ./coriell2snpchip.pl <.pre File> <.map File>\n");
  print("\n");
}


my @indarray = ();
my $numMarkers = 0;

if(scalar(@ARGV) < 2){
    print(STDERR "Error: Fewer than two input files specified\n");
    usage();
    exit(1);
}

open PREFILE, "< ".shift(@ARGV) or die("Error: could not open .pre file\n");
open MAPFILE, "< ".shift(@ARGV) or die("Error: could not open .map file\n");


while (<PREFILE>) {
    chomp;
    s/^[^\- ]+[\- ][^ ]+ [^ ]+ //; # data files seem to be inconstent...
    s/([^ ]+) ([^ ]+)/$1$2/g; # combine adjacent alleles into single genotype
    s/0/N/g; # replace null genotypes with 'N' (used in HapMap samples)
    my @tmp = split;
    # print join(", ", @tmp)."\n";
    if(!$numMarkers){
        $numMarkers = (@tmp + 0);
    }
    # FIXME: hrm... need to combine adjacent alleles
    die ("Error: .pre file has lines of differing lengths\n")
        if ($numMarkers != (@tmp + 0));
    push @indarray, [ @tmp ];
}

my $numInds = (@indarray + 0);


# read in 8xN array (coriell .map marker file notation)
# map file format is <CS> <location> <rsNum> <mut1> <mut2> \
# <freq1> <freq2> <?missingNo?>

my @markerNames = ();
my $markerCount = 0;
while (<MAPFILE>) {
    my @tmp = split;
    push @markerNames, $tmp[2];
    $markerCount++;
}

die ("Error: .pre file declares $numMarkers markers,\n".
     " while .map file declares $markerCount markers.\n".
     "(these numbers must be the same for the conversion to work)\n")
    if ($numMarkers != $markerCount);

# spit out (1+M)xN array

$markerCount = 0;

foreach my $marker (@markerNames) {
    printf "%-12s", $marker;
    for(my $i=0; $i<$numInds; $i++){
        print " $indarray[$i][$markerCount]";
    }
    $markerCount++;
    print "\n";
}
