#!/usr/bin/perl

# snpchip2linkage.pl -- converts from a simplegt formatted file into a
# linkage formatted file.

use strict;
use warnings;

# Author: David Eccles (gringer), 2007 <programming@gringer.org>

# usage: ./snpchip2linkage.pl <marker location file> <genotype file>
# updated 2008 to accept individual label header

sub usage {
    print("usage: ./snpchip2linkage.pl <marker location file> ".
        "<genotype file>\n\n");
}

my %genogroup = ();

my $numgts = 0;
my %markers = ();
my @markerOrder = ();

if(!@ARGV){
    print STDERR "Error: no input files on command line\n";
    usage();
    exit(1);
}

if(-f $ARGV[0]){
    open MARKFILE, "< ".shift(@ARGV);
    while(<MARKFILE>){
        if(m/(rs[0-9]+)\s+([0-9]+)\b/){
            $markers{$1} = $2;
            push(@markerOrder,$1);
        }
    }
}

my @indLabels = ();
my @phenoVals = ();
my $phenoNum = 1;

while (<>){
    my $line = $_;
    if($line =~ /^##/){
        ## Determine individual labels
        ## This works even if more than one <ID> region is present in the
        ## header line, as might be the case in a 'join'ed file
        if($line =~ /IDs:\s+(.*?)\s*>/){
            @indLabels = ();
            @phenoVals = ();
        }
        while($line =~ /IDs:\s+(.*?)\s*>/){
            my @lineData = split(/\s+/, $1);
            push(@indLabels, @lineData);
            ## generate an array of (@lineData) copies of $phenoNum
            push(@phenoVals, (($phenoNum) x scalar(@lineData)) );
            $phenoNum++;
            $line =~ s/^.*?>//;
        }
    } else{
        $line =~ s/\s+/ /g; #replace whitespace with spaces
        $line =~ s/a/1/ig;
        $line =~ s/c/2/ig;
        $line =~ s/g/3/ig;
        $line =~ s/t/4/ig;
        $line =~ s/ ([34])([12])/ $2$1/ig; # order heterozygotes so
                                           # smallest number is first
        $line =~ s/ [^0-9][^0-9]/ 00/ig;
        my ($marker, @genotypes) = split(/ /,$line);
        if($markers{$marker}){
            @{$genogroup{$marker}} = @genotypes;
            if($numgts && ($numgts != (@genotypes))){
                die("Number of genotypes does not match");
            }
            else{
                $numgts = (@genotypes);
            }
        }
    }
}

if(@indLabels){
    if(scalar(@indLabels) != $numgts){
        die("Number of genotypes does not match labels in header line");
    }
} else {
    @indLabels = (1 .. $numgts);
    @phenoVals = ((0) x scalar(@indLabels));
}

for(my $i=0; $i < $numgts; $i++){
    my @genotypes = ();
    foreach my $marker (@markerOrder){
        my $writeval = @{$genogroup{$marker}}[$i];
        $writeval =~ s/(.)(.)/$1 $2/;
        push(@genotypes, $writeval);
    }
    printf("%s 1 0 0 0 %d  ", $indLabels[$i], $phenoVals[$i]);
    print join("  ",@genotypes)."\n";
}
