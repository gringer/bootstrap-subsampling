#!/usr/bin/perl

# gt2bayes.pl -- calculate bayesian probability for group assignment,
# assuming marker independence. Also calculates log(p(max)/p(min)) to
# determine reliability of group assignment

# This script currently assumes that the test group is the same as the
# training group

# Author: David Eccles (gringer), 2009 <programming@gringer.org>

use warnings;
use strict;

sub usage {
    print "usage: ./bayes.pl <simplegt file> p<P1 name> <P1 ids> p<P2 name> <P2 ids>\n";
}

sub sum {
  my ($sum, $elem);
  $sum = 0;
  foreach $elem (@_) {
    $sum += $elem;
  }
  return($sum);
}


sub countGTs {
    my @genotypes = @{\@_};
    my $aacount = scalar(grep(/[at][at]/i, @genotypes));
    my $account = scalar(grep(/([at][cg]|[cg][at])/i, @genotypes));
    my $cccount = scalar(grep(/[cg][cg]/i, @genotypes));
    $aacount = 0.5 unless ($aacount > 0);
    $account = 0.5 unless ($account > 0);
    $cccount = 0.5 unless ($cccount > 0);
    return ( $aacount, $account, $cccount );
}

my @filenames = ();
my %indivs = ();
my %allgenos = ();
my %popgtcounts = ();
my %gtcounts = ();
my %totalcounts = ();
my %poptotalcounts = ();
my @markers = ();
my $numindivs = 0;

my $debug = 0; # false

while(@ARGV){
    my $argval = shift(@ARGV);
    if(-e $argval){
    } else {
        if ($argval =~ /^p(.*)$/i){
            my $popname = $1;
            $indivs{$popname} = ();
            while((@ARGV) && ($ARGV[0] =~ /^[0-9]+$/)){
                push (@{$indivs{$popname}}, shift(@ARGV));
            }
        } elsif ($argval eq "-v"){
            $debug = 1; # true;
        } else {
            print(STDERR "Error: Unknown command line option %s\n", $_);
            usage();
            exit(1);
        }
    }
}

while(<>){
    chomp;
    if(!/^#/){
        my ($marker, @genotypes) = split(/[\s,]+/);
        $numindivs = scalar(@genotypes) unless ($numindivs > 0);
        push(@markers, $marker);
        if($debug) {
            printf("Marker = %s, ", $marker);
        }
        foreach my $pop (keys(%indivs)){
            my @countResult = countGTs(@genotypes[@{$indivs{$pop}}]);
            $popgtcounts{$pop}->{$marker} = [ @countResult ];
            $poptotalcounts{$pop}->{$marker} = sum(@countResult);
            if($debug){
                printf("Counts(%s) = %s, total = %f, ", $pop,
                       join(":",@countResult), sum(@countResult));
            }
        }
        $gtcounts{$marker} = [ countGTs(@genotypes) ];
        $totalcounts{$marker} = sum(@{$gtcounts{$marker}});
        if($debug){
            printf("TotalCounts = %s, total = %f\n",
                   join(":",@{$gtcounts{$marker}}), $totalcounts{$marker});
        }
        for(my $i = 0; $i < scalar(@genotypes); $i++){
            push(@{$allgenos{$i}}, $genotypes[$i]);
        }
    }
}

foreach my $ind (keys(%allgenos)){
    my @genotypes = @{$allgenos{$ind}};
    my $gtprob = 0;
    my %condprob = ();
    if($debug){
        printf("%4s (%s):\n", $ind, join(":", @genotypes));
    }
    foreach my $pop (keys(%indivs)){
        $condprob{$pop} = 0;
    }
    for(my $i = 0; $i < scalar(@genotypes); $i++){
        my $gtval = -1;
        if ($genotypes[$i] =~ /[at][at]/i){
            $gtval = 0;
        } elsif ($genotypes[$i] =~ /[cg][cg]/i){
            $gtval = 2;
        } elsif ($genotypes[$i] =~ /[atcg][atcg]/i){
            $gtval = 1;
        }
        unless($gtval < 0){
            $gtprob += log(@{$gtcounts{$markers[$i]}}[$gtval] / $totalcounts{$markers[$i]});
            if($debug){
                printf("GT(%s) = %f, ", $genotypes[$i], exp($gtprob));
            }
            foreach my $pop (keys(%indivs)){
                $condprob{$pop} += log(@{$popgtcounts{$pop}->{$markers[$i]}}[$gtval] /
                    $poptotalcounts{$pop}->{$markers[$i]});
                if($debug){
                    printf("Condprob(%s) = %f, ", $pop, exp($condprob{$pop}));
                }
            }
            if($debug){
                print("\n");
            }
        }
    }
    printf("%4s: ", $ind);
    my $maxprob = "";
    my $maxpop = "";
    my $minprob = "";
    my $minpop = "";
    foreach my $pop (sort(keys(%indivs))){
        my $prob = exp($condprob{$pop} - $gtprob) * (scalar(@{$indivs{$pop}}) / $numindivs);
        printf("p(%s, %d/%d) = %f, ", $pop, scalar(@{$indivs{$pop}}), $numindivs, $prob);
        if((!$maxprob) || ($prob > $maxprob)){
            $maxprob = $prob;
            $maxpop = $pop;
        }
        if((!$minprob) || ($prob < $minprob)){
            $minprob = $prob;
            $minpop = $pop;
        }
    }
    printf("%s (%f)\n",$maxpop, (log($maxprob) - log($minprob))/log(10));
}
