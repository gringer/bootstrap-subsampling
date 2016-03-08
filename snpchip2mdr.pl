#!/usr/bin/perl

# snpchip2mdr.pl -- convert from simplegt formatted text file to MDR
# formatted text file.

# Author: David Eccles (gringer), 2008 <programming@gringer.org>

use strict;
use warnings;

sub usage {
  print("usage: ./snpchip2mdr.pl [options] <file name>\n");
  print("\nOther Options:\n");
  print("-pMDR : format for pMDR instead of MDR\n");
  print("\n");
}

# read in Nx(M+1) array ('simple' HapMap notation)

my %genogroup = ();
my $startClassAt = 0;
my $pMDR = 0;
my $classAtStart = 0;

if($ARGV[0] && ($ARGV[0] eq "-pMDR")){
    shift(@ARGV);
    $pMDR = 1;
    $classAtStart = $pMDR;
}

if($ARGV[0]){
    $startClassAt = shift(@ARGV);
}

# spit out (M+1)xN array

my $numgts = 0;

while (<>){
    chomp;
    my $line = $_;
    $line =~ s/(a|t)/1/ig;
    $line =~ s/(g|c)/2/ig;
    $line =~ s/ 21/ 12/ig;
    $line =~ s/ [^0-9][^0-9]/ 00/ig;
    my ($marker, @genotypes) = split(/\s+/,$line);
    @{$genogroup{$marker}} = @genotypes;
    if($numgts && ($numgts != (@genotypes))){
        die("Number of genotypes does not match");
    }
    else{
        $numgts = (@genotypes);
    }
}

my @markers = ();
foreach my $marker (keys(%genogroup)){
    push(@markers, $marker);
}


print join("\t",@markers)."\tClass\n" unless $pMDR;

for(my $i=0; $i < $numgts; $i++){
    my @genotypes = ();
    foreach my $marker (keys(%genogroup)){
	my $writeval = @{$genogroup{$marker}}[$i];
	if ($writeval == "21"){
	    $writeval = "12"; # reorder heterozygotes
	}
	if ($writeval == "00"){
	    $writeval = "00"; # remove null values (not supported by MDR)
	}
	if($pMDR){# genotypes are not separated in MDR
	    #$writeval =~ s/(.)(.)/$1 $2/;
	    $writeval =~ s/11/0/; # convert Hom/Het/Hom to  0/1/2
	    $writeval =~ s/12/1/;
	    $writeval =~ s/22/2/;
	    $writeval =~ s/00/-1/g;
	}
	push(@genotypes, $writeval);
    }
#    printf "%04d\t", $i; #individual ID not used for MDR format
    my $classString = "0";
    if($i >= $startClassAt){
	$classString = "1";
    }
    if($classAtStart){
	print $classString."\t";
    }
    print join("\t",@genotypes);
    if(!$classAtStart){
	print "\t".$classString;
    }
    print "\n";
}
