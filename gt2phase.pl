#!/usr/bin/perl

# gt2phase.pl -- convert from simplegt formatted text file input file
# formatted for PHASE or fastPHASE.

# (derived from snpchip2structure.pl)

# Author: David Eccles (gringer), 2008 <programming@gringer.org>

# read in NxM array ('simple' HapMap notation)
# spit out (4+M)xN array (PHASE input file)

use strict;
use warnings;

sub usage {
  print("usage: ./gt2phase.pl < <file name>\n");
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
        $line =~ s/\b[^ ]?0[^ ]?\b/\?\?/ig;
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

my @markers = ();
foreach my $marker (keys(%genogroup)){
    push(@markers, $marker);
}

# print join(" ",@markers)."\n"; # from snpchip2structure -- not needed

print $numgts."\n"; # no.individuals
print scalar(@markers)."\n"; # no.SNPsites

for(my $i=0; $i < $numgts; $i++){
    if($indlabels[$i] =~ /^[0-9]+$/){
        printf "# id %04d\n", $indlabels[$i];
    } else {
        printf "# id %4s\n", $indlabels[$i];
    }
    for(my $allele = 0; $allele < 2; $allele++){
        my @genotypes = ();
        foreach my $marker (keys(%genogroup)){
            my $writeval = substr(@{$genogroup{$marker}}[$i],$allele,1);
            push(@genotypes, $writeval);
        }
        print join("",@genotypes)."\n";
    }
}
