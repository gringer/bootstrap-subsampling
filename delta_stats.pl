#!/usr/bin/perl

# delta_stats.pl -- just produces delta for each marker

# Author: David Eccles (gringer), 2008 <programming@gringer.org>

# arguments:
# ./delta_stats.pl p<P1 name> <P1 cols> p<P2 name> <P2 cols>

# Note: this algorithm will only work on two populations

use strict;
use warnings;

sub usage {
    print("usage: /delta_stats.pl p<P1 name> <P1 cols> ".
          "p<P2 name> <P2 cols>\n");
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
my $showAlleles = 1; #true
my $showNull = 1; #true

while ((@ARGV > 0) && !(-f $ARGV[0])){
    my $argval = shift(@ARGV);
    if ($argval =~ /^p(.*)$/i){
        $popname = $1;
        push(@pops,$popname);
				$indivs{$popname} = ();
    }
    elsif ($argval =~ /^-noalleles/){
        $showAlleles = 0; #false
    }
    elsif ($argval =~ /^-nonull/){
        $showNull = 0; #false
    }
    elsif (($popname) && $argval =~ /^[0-9]+$/) {
				push (@{$indivs{$popname}}, $argval);
    }
}

my $nw = '%-0.6f '; #number width float formatting
my $sw = '%-8s '; #number width string formatting
my $mw = '%-17s '; #marker width string formatting

my $ncounts = 0;

# print out Column Titles

printf $mw, "Marker";
foreach my $pop1 (@pops){
    if($showAlleles){
        printf $sw, "p(A,$pop1)";
        printf $sw, "p(C,$pop1)";
        printf $sw, "MAF($pop1)";
    }
}
if($showNull){
    printf $sw, "missing";
}
printf $sw, "delta";
print "\n";

if ((@pops + 0)  != 2){
    print "\nERROR: Two populations have not been defined\n";
    print "arguments: $progname p<P1 name> ".
	"<P1 cols> p<P2 name> <P2 cols> < <input>\n";
    exit(1);
}


#exit(0);
# print "Total number of comparisons to do: $ncounts\n";

my $pA = 0;
my $pAC = 0;
my $pC = 0;
my $pNull = 0;
my %popScaleVal = ();
my @popgts = ();
my $gt;
my $delta;
my $nmarkers = 0;
my $badData;

foreach my $pop1 (@pops){
    $popScaleVal{$pop1} = 1 / ((@{$indivs{$pop1}} + 0)*2);
}

while(<>){
    $_ = lc; # lower case string
    my ($marker, @genotypes) = split(/\s/,$_);
    $badData = 0;
    $nmarkers++;
    printf $mw, $marker;
    $delta = 0;
    $pNull = 0;
    foreach my $pop1 (@pops){
        $pA = 0;
        $pAC = 0;
        $pC = 0;
        # retrieve genotype list for population individuals
        my @popgts = @genotypes[@{$indivs{$pop1}}];
        # grep returns the count of matching patterns
        $pAC = grep(/([at][gc]|[gc][at])/,@popgts); # heterozygote AC
        $pA = 2 * grep(/[at][at]/,@popgts) + $pAC; # homozygote AA
        $pC = 2 * grep(/[gc][gc]/,@popgts) + $pAC; # homozygote CC
        # determine frequencies using scale values
        $pA = $pA * $popScaleVal{$pop1};
        $pC = $pC * $popScaleVal{$pop1};
        $pNull += (1 - ($pA + $pC));
        if(($pA == 0) && ($pC == 0)){
            $badData = 1; # true
        }
        if($showAlleles){
            printf $nw, $pA;
            printf $nw, $pC;
        }
        if(!$badData){
            my $tA = $pA / ($pA + $pC);
            my $tC = $pC / ($pA + $pC);
            $pA = $tA;
            $pC = $tC;
        }
        if($showAlleles){
            if($pA < $pC){
                printf $nw, $pA;
            } else {
                printf $nw, $pC;
            }
        }
        if (!$delta){
            $delta = $pA;
        }
        else{
            $delta = abs($delta - $pA);
        }
    }
    if($badData){
        $delta = 0;
    }
    if($showNull){
        printf $nw, ($pNull / 2);
    }
    printf $nw, $delta;
    print "\n";
}
