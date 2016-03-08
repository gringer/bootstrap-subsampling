#!/usr/bin/perl

# 3_indiv_minimal.pl -- Make sure all individuals are different using
# a minimal set of information from the genotyped individuals.
# Practically, this means extracting additional lines from an input
# file until every individual appears different (or the total input
# file, whichever comes first)

# Author: David Eccles (gringer), 2007 <programming@gringer.org>

use strict;
use warnings;

sub usage {
  print("usage: ./3_individual.pl <individual columns> < <fileName>\n");
  print("\n");
}

# extract command line arguments
my @indivs = @ARGV;

# remove command line arguments (so while(<>) works)
foreach my $argnum (0 .. $#ARGV) {
    my $test=shift @ARGV;
}

my $nw = '%-0.5f '; #number width string formatting
my $sw = '%-7s '; #number width string formatting

my $nmarkers = 0;
my @diffcounts = 0;

for(my $c1=0;$c1<(@indivs-1);$c1++) {
    for(my $c2=$c1+1;$c2<(@indivs+0);$c2++) {
	$diffcounts[$c1*(@indivs+0)+$c2] = 0;
    }
}

while(<>){
    my ($marker, @genotypes) = split(/\s/,$_);
    $nmarkers++;

    # substitute complementary alleles to simplify comparisons later
    for(my $c1=0;$c1<(@indivs+0);$c1++) {
	$genotypes[$indivs[$c1]] =~ s/(a|t)/t/ig;
	$genotypes[$indivs[$c1]] =~ s/(g|c)/c/ig;
    }
    print $marker."\n";
    for(my $c1=0;$c1<(@indivs-1);$c1++) {
 	for(my $c2=$c1+1;$c2<(@indivs+0);$c2++) {
 	    my $gt1 = $genotypes[$indivs[$c1]];
 	    my $gt2 = $genotypes[$indivs[$c2]];
	    #same genotype (assuming a=t,c=g)
 	    if ($gt1 eq $gt2){
		$diffcounts[$c1*(@indivs+0)+$c2]++;
 	    }
	    #same genotype, but switched around
 	    elsif (($gt1 eq "ac" && $gt2 eq "ca") ||
		   ($gt1 eq "ca" && $gt2 eq "ac")){
		$diffcounts[$c1*(@indivs+0)+$c2]++;
 	    }
 	    # if genotypes are not similar, then don't consider similar
 	    else {
		# @diffcounts[$c1*(@indivs+0)+$c2] += 0;
 	    }
 	}
    }
    my $different = 1; # true
    for(my $c1=0;$c1<(@indivs-1);$c1++) {
 	for(my $c2=$c1+1;$c2<(@indivs+0);$c2++) {
	    if ($diffcounts[$c1*(@indivs+0)+$c2] == $nmarkers){
		$different = 0; # false
	    }
 	}
    }
    if ($different){
	last;
    }
}

# Header line 1, shows number of individuals
printf "Number of individuals: %d\n", (@indivs+0);

# Header line 2, shows number of markers
printf "Number of markers: %d\n", $nmarkers;

# Header line 3, shows first person comparison
printf "%-17s", "Pairwise: P1 vs.";

for(my $c1=0;$c1<(@indivs-1);$c1++) {
    printf $sw, $indivs[$c1];
    for(my $c2=$c1+1;$c2<(@indivs-1);$c2++) {
	if ($c2+1==$c1){
	    printf $sw, "*";
	}
	else{
	    printf $sw, "";
	}
    }
}

print "\n";

# Header line 4, shows second person comparison
printf "%-17s", "Pairwise: P2";

for(my $c1=0;$c1<(@indivs-1);$c1++) {
    for(my $c2=$c1+1;$c2<(@indivs+0);$c2++) {
	printf $sw, $indivs[$c2];
    }
}

print "\n";

printf "%-17s", "Pairwise: ";
for(my $c1=0;$c1<(@indivs-1);$c1++) {
    for(my $c2=$c1+1;$c2<(@indivs+0);$c2++) {
	printf $nw, ($diffcounts[$c1*(@indivs+0)+$c2] / $nmarkers);
    }
}
print "\n";
