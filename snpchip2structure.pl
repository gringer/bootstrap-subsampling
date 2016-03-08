#!/usr/bin/perl

# snpchip2structure.pl -- convert from simplegt formatted text file
# input file formatted for structure.

# Author: David Eccles (gringer), 2007 <programming@gringer.org>

# read in NxM array ('simple' HapMap notation)
# spit out (4+M)xN array (structure input file)
# 2008-05-27: added modification to allow labelling of individuals
# 2008-06-05: added modification to allow 4 different alleles
# 2008-08-22: fixed bugs identified by porting to R
#             - extract marker name *before* doing numerical conversions
#             - make sure failure on one allele = failure on both

use strict;
use warnings;

sub usage {
  print("usage: ./snpchip2structure.pl < <file name>\n");
  print("\nOther Options:\n");
  print("-nocombine : Don't combine complementary genotypes\n");
  print("\n");
}

my %genogroup = ();
my $numgts = 0;
my $combineComp = 1; # true

if((@ARGV + 0) > 0){
    $_ = shift @ARGV;
    if(/-nocombine/){
        $combineComp = 0; # false
    }
}

my @indlabels = ();
my @markerOrder = ();

while (<>){
    my $line = $_;
    if($line =~ "^#"){
        if($line =~ m%<Individual/Column IDs:\s*(.*?)\s*>%){
            @indlabels = split(/\s+/,$1);
        }
    }
    else{
        $line =~ s/\s+/ /g; #replace whitespace with spaces
        $line =~ s/^(.*?) (.*)$/$2/;
        my $marker = $1;
        push(@markerOrder, $marker);
        $line =~ s/a/1/ig;
        $line =~ s/c/2/ig;
        $line =~ s/g/3/ig;
        $line =~ s/t/4/ig;
        if($combineComp){
            $line =~ s/4/1/ig;
            $line =~ s/3/2/ig;
        }
        $line =~ s/\b21\b/12/ig;
        $line =~ s/[^1-4 ]/0/ig;
        $line =~ s/\b[^ ]?0[^ ]?\b/00/ig;
        my @genotypes = split(/ /,$line);
        @{$genogroup{$marker}} = @genotypes;
        if($numgts && ($numgts != (@genotypes))){
            die("Number of genotypes does not match");
        }
        else{
            $numgts = (@genotypes);
        }
    }
}

if(!@indlabels){
    @indlabels = (1 .. $numgts);
}
# print STDERR "(".join(",",@indlabels).")";

if((@indlabels + 0) != $numgts){
    die "Error: Number of labels (".(@indlabels + 0).") does not equal number of genotypes (".$numgts.")\n";
}

print join(" ",@markerOrder)."\n";

for(my $i=0; $i < $numgts; $i++){
    my @genotypes = ();
    foreach my $marker (@markerOrder){
        my $writeval = @{$genogroup{$marker}}[$i];
        $writeval =~ s/(.)(.)/$1 $2/;
        push(@genotypes, $writeval);
    }
    if($indlabels[$i] =~ /^[0-9]+$/){
        printf "%04d  ", $indlabels[$i];
    } else {
        printf "%4s  ", $indlabels[$i];
    }
    print join("  ",@genotypes)."\n";
}
