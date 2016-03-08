#!/usr/bin/perl

# emp_fst.pl -- Calculates Fst values using pairwise comparisons
# calculations are carried out for between and within populations,
# with the resultant statistic being the following:
# Fst = (mean(between) - mean(mean(within1), mean(within2)) / mean(between))

# Author: David Eccles (gringer), 2008 <programming@gringer.org>

# Note: this is an incorrect definition of Fst. See snpchip2fst.r for
# a better definition.

use warnings;
use strict;

sub usage {
  print("usage: ./emp_fst.pl p<P1 name> <P1 cols> p<P2 name> <P2 cols> < <file>\n");
  print("\n");
}

sub gtCompare {
    my ($gt1, $gt2) = @_;
    $gt1 =~ s/(a|t)/a/ig;
    $gt1 =~ s/(c|g)/c/ig;
    $gt2 =~ s/(a|t)/a/ig;
    $gt2 =~ s/(c|g)/c/ig;
    my $compinc = 0;
    if (($gt1 eq $gt2) && ($gt1 =~ /[ac][ac]/)){ # must exclude nn/--/etc.
        $compinc = 2; #homozygote vs homozygote
    }
    elsif ((($gt1 eq "ac") || ($gt1 eq "ca")) &&
           (($gt2 eq "ac") || ($gt2 eq "ca"))){
        $compinc = 2; #heterozygote vs heterozygote
    }
    elsif((($gt1 eq "aa") || ($gt1 eq "cc")) &&
          (($gt2 eq "ac") || ($gt2 eq "ca"))){
        $compinc = 1; #homozygote vs heterozygote
    }
    elsif((($gt2 eq "aa") || ($gt2 eq "cc")) &&
          (($gt1 eq "ac") || ($gt1 eq "ca"))){
        $compinc = 1; #homozygote vs heterozygote
    }
#   print "checking $gt1 and $gt2: $compinc\n";
    return $compinc;
}

sub pwCompare {
    my ($list1Ref, $list2Ref, $gtRef) = @_;
    my @list1 = @$list1Ref;
    my @list2 = @$list2Ref;
#   print "Comparing ".join(",",@list1)." and ".join(",",@list2)."\n";
    my @gts = @$gtRef;
    my $pairwiseSum = 0;
    my $pairwiseCount = 0;
    foreach my $indiv1 (@list1){
        foreach my $indiv2 (@list2){
#            if($indiv1 != $indiv2){
                $pairwiseSum +=
                    gtCompare($gts[$indiv1],$gts[$indiv2]);
                $pairwiseCount++;
#            }
        }
    }
    return ($pairwiseSum / ($pairwiseCount * 2));
}

my %indivs = ();

# remove command line arguments (so while(<>) works)
my $pop1Name = 0; # false
my $pop2Name = 0; # false
my @pop1Indivs = ();
my @pop2Indivs = ();

while (@ARGV && !(-f $ARGV[0])) {
    my $argval = shift @ARGV;
    if ($argval =~ /^p(.*)$/i){
        if(!$pop1Name){
            $pop1Name = $1;
        }
        elsif(!$pop2Name){
            $pop2Name = $1;
        }
        else{
            die("Error: More than two populations specified");
        }
    }
    elsif ($argval =~ /^[0-9]+$/) {
        if($pop2Name){
            push(@pop2Indivs, $argval);
        }
        elsif($pop1Name){
            push(@pop1Indivs, $argval);
        }
        else{
            die("Error: No populations defined, but individual columns stated");
        }
    }
}

#print join(",",@pop1Indivs)."\n";
#print join(",",@pop2Indivs)."\n";

if(!$pop2Name){
    die("Error: Fewer than two populations specified");
}

my $nw = '%-0.5f '; #number width string formatting
my $sw = '%-7s '; #number width string formatting

while(<>){
    $_ = lc;
    my ($marker, @genotypes) = split(/\s+/);
    printf $sw, $marker;
    my $within1 = pwCompare(\@pop1Indivs,\@pop1Indivs, \@genotypes);
    my $within2 = pwCompare(\@pop2Indivs,\@pop2Indivs, \@genotypes);
    my $between = pwCompare(\@pop1Indivs,\@pop2Indivs, \@genotypes);
#    print "[".join(",",@genotypes[@pop1Indivs])."] [".
#           join(",",@genotypes[@pop2Indivs])."] : ";
    my $within = ($within1 + $within2) / 2;
    my $fst = 0;
    if($between != 1){
        $fst = ((1-$between) - (1-$within)) / (1-$between);
    }
    printf("w1: $nw, w2: $nw, wm: $nw, bt: $nw, fst: $nw\n",
           $within1, $within2, $within, $between, $fst);
}
