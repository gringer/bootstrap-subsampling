#!/usr/bin/perl

# maf_stats.pl -- calculates allele frequencies, based on MAF of first
# population (derived from delta_stats.pl).

# Author: David Eccles (gringer), 2009 <programming@gringer.org>

use strict;
use warnings;

sub usage {
  print("usage: ./maf_stats.pl p<P1 name> <P1 cols> p<P2 name> <P2 cols> ...\n");
  print("\n");
}

sub max {
    my ($val1, $val2) = @_;
    if($val2 > $val1){
        return $val2;
    }
    else{
        return $val1;
    }
}

my %indivs = ();

# remove command line arguments (so while(<>) works)
my $popname = 0;
my $indnum = 0;
my $progname = $0;
my @pops = ();
my @files = ();

while ((@ARGV > 0) && !(-f $ARGV[0])){
    my $argval = shift(@ARGV);
    if ($argval =~ /^p(.*)$/i){
        $popname = $1;
        push(@pops,$popname);
				$indivs{$popname} = ();
    }
    elsif (($popname) && $argval =~ /^[0-9]+$/) {
				push (@{$indivs{$popname}}, $argval);
    }
}

my $nw = '%-0.3f    '; #number width float formatting
my $sw = '%-8s '; #number width string formatting
my $mw = '%-17s '; #marker width string formatting

my $ncounts = 0;

# print out Column Titles

printf $mw, "Marker";
printf $sw, "Allele";
foreach my $pop (@pops){
    printf $sw, $pop;
}
print "\n";

if ((@pops + 0)  < 1){
    print "\nERROR: At least one population must be defined\n";
    print "arguments: $progname p<P1 name> ".
	"<P1 cols> p<P2 name> <P2 cols> ... < <input>\n";
    exit(1);
}


#exit(0);

my $pA = 0;
my $pAC = 0;
my $pC = 0;
my $pT = 0;
my $pG = 0;
my @popgts = ();
my $gt;
my $delta;
my $nmarkers = 0;
my $badData;
my $MAFAllele = 0; # false

while(<>){
    if(/^#/){
        next;
    }
    $_ = lc; # lower case string
    my ($marker, @genotypes) = split(/\s/,$_);
    $badData = 0;
    $MAFAllele = 0; # false
    $nmarkers++;
    printf $mw, $marker;
    $delta = 0;
    foreach my $pop1 (@pops){
        $pA = 0; $pC = 0; $pG = 0; $pT = 0;
        $pAC = 0;
        my @popgts = @genotypes[@{$indivs{$pop1}}];
        $pAC = grep(/([at][gc]|[gc][at])/,@popgts); # heterozygote AC
        $pA = 2 * grep(/[at][at]/,@popgts) + $pAC; # homozygote AA
        $pC = 2 * grep(/[gc][gc]/,@popgts) + $pAC; # homozygote CC
        $pT = $pA + $pC;
        $pA = $pA / $pT; # frequencies based on observed counts, ignoring
        $pC = $pC / $pT; # null values
        if(($pA == 0) && ($pC == 0)){
            $badData = 1;
        }
        if(!$MAFAllele){
            if($pA <= $pC){
                $MAFAllele = "a";
            } else {
                $MAFAllele = "c";
            }
            printf $sw, $MAFAllele;
        }
        if($badData){
            printf $sw, "NA";
        } else {
            if($MAFAllele eq "a"){
                printf $nw, $pA;
            } else {
                printf $nw, $pC;
            }
        }
    }
    print "\n";
}
