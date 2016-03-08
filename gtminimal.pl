#!/usr/bin/perl

# gtminimal.pl -- Determine minimal set sizes for an input data
# set. Data is processed in a linear fashion, taking successive lines
# and finding out how large an untrained selection of markers needs to
# be in order to distinguish populations

# Author: David Eccles (gringer), 2009 <programming@gringer.org>

use strict;
use warnings;

sub usage {
  print("usage: ./gtminimal.pl -(individual|population) [options] < <input file>\n");
  print("\nOther Options:\n");
  print("-help       : Only display this help message\n");
  print("-population : Calculate population-based sizes\n");
  print("-individual : Calculate individual-based sizes\n");
  print("-verbose    : Print out generated hash values during run\n");
  print("\n");
}

my $method = "none";
my $verbose = 0; # false
my @fileNames = ();

# extract command line arguments
while(@ARGV){
    my $argument = shift @ARGV;
    if(-f $argument){ # file existence check
        push(@fileNames, $argument);
    } else {
        if($argument eq "-help"){
            usage();
            exit(0);
        }
        elsif($argument eq "-population"){
            $method = "population";
        }
        elsif($argument eq "-individual"){
            $method = "individual";
        }
        elsif($argument eq "-verbose"){
            $verbose = 1; # true
        }
        else{
            printf(STDERR "Error: Unknown command-line option, '$argument'\n");
            usage();
            exit(4);
        }
    }
}

push(@ARGV, @fileNames);

if($method eq "none"){
    print(STDERR "Error: No method specified (population/individual)\n");
    usage();
    exit(1);
}

my $lineCount = 0;
my @indLabels = ();
my @phenoVals = ();
my $phenoNum = 1;
my $marker = "";
my @genotypes = ();

my @hashVals = ();
my $nextHash = 0;
my %hashGTMaps = ();

my %seenHashes = ();
my %phenoHashes = ();
my $curPheno = 0;
my $intersectionFound = 0; # false
my @markerSet = ();

while (<>){
    chomp;
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
        @hashVals = ((0) x scalar(@phenoVals));
    } else{
        $nextHash = 0;
        $intersectionFound = 0; # false
        ($marker, @genotypes) = split(/\s+/, $line);
        push(@markerSet, $marker);
        if(!@indLabels){
            @indLabels = (1 .. (@genotypes));
            @phenoVals = ((0) x scalar(@genotypes));
            @hashVals = ((0) x scalar(@phenoVals));
        }
        %hashGTMaps = ();
        %phenoHashes = ();
        %seenHashes = ();
        $curPheno = $phenoVals[0];
        for(my $i = 0; $i < scalar(@genotypes); $i++){
            my $nextCode = $hashVals[$i].$genotypes[$i];
            if($hashGTMaps{$nextCode}){
                $hashVals[$i] = $hashGTMaps{$nextCode};
            } else {
                $hashVals[$i] = $nextHash++;
                $hashGTMaps{$nextCode} = $hashVals[$i];
            }
            if($curPheno != $phenoVals[$i]){
                if($method eq "population"){
                    %seenHashes = (%seenHashes, %phenoHashes);
                    %phenoHashes = ();
                }
                $curPheno = $phenoVals[$i];
            }
            if((($method eq "population") && ($seenHashes{$hashVals[$i]})) ||
               (($method eq "individual") && ($phenoHashes{$hashVals[$i]}))){
                $intersectionFound = 1; # true
            } else {
                $phenoHashes{$hashVals[$i]} = 1; # true
            }
        }
        if($verbose){
            printf("%s %s\n", $marker, join(",", @hashVals));
        }
        if(!$intersectionFound){
            printf("%d %s\n", scalar(@markerSet), join(",", @markerSet));
            @markerSet = ();
            @hashVals = ((0) x scalar(@phenoVals));
        }
    }
}

if($intersectionFound){
    printf("%d END,%s\n", scalar(@markerSet), join(",", @markerSet));
}
