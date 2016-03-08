#!/usr/bin/perl

# 9_predictive_power.pl -- Determine the predictive power of each
# marker, including power to distinguish between two populations
# outputs a list of markers, and probabilities associated with the
# distinguishing power of each marker

# Author: David Eccles (gringer), 2007 <programming@gringer.org>

sub usage {
    print("usage: /9_predictive_power.pl p<P1 name> <P1 cols> ".
          "p<P2 name> <P2 cols>\n");
    print("\nOther Options:\n");
    print("-summarise : Only create summary statistics (mean values)\n");
    print("\n");
}

# Note: this algorithm has not been tested on more than two
# populations. Final difference statistic will need to be changed.

use strict;
use warnings;

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
my $totargs = $#ARGV;
my $progname = $0;
my $summarise = 0; # false
foreach my $argnum (0 .. $totargs) {
    my $argval = shift @ARGV;
    if ($argval =~ /^p(.*)$/i){
				$popname = $1;
				$indivs{$popname} = ();
    }
    if ($argval =~ /^-summari[sz]e$/){
				$summarise = 1; #true
    }
    elsif (($popname) && $argval =~ /^[0-9]+$/) {
				push (@{$indivs{$popname}}, $argval);
    }
}

my $nw = '%-0.6f '; #number width float formatting
my $sw = '%-8s '; #number width string formatting
my $mw = '%-17s '; #marker width string formatting

my $ncounts = 0;
my @pops = ();

# print out Column Titles

printf $mw, "Marker";

foreach my $pop1 (keys(%indivs)){
    printf $sw, "p(A,$pop1)";
    printf $sw, "p(C,$pop1)";
    printf $sw, "hom($pop1)";
    push(@pops, $pop1);
}

printf $sw, "delta";
printf $sw, "pp(Sum)";
printf $sw, "pp(Mul)";
print "\n";

if ((@pops + 0)  == 0){
    print "\nERROR: No populations defined\n";
    print "arguments: $progname p<P1 name> ".
	"<P1 cols> p<P2 name> <P2 cols> < <input>\n";
    exit(1);
}

my @summaries;
my $spos;

foreach my $i (0 .. (@pops + 1)) { # + 2, -1 for indexing
    $summaries[$i] = 0;
}

#exit(0);
# print "Total number of comparisons to do: $ncounts\n";

my %pA = ();
my %pC = ();
my $gt;
my $pp;
my $delta;
my $endStat;
my $mulStat;
my $nmarkers = 0;
my $badData;

while(<>){
    $badData = 0;
    my ($marker, @genotypes) = split(/\s/,$_);
    $spos = 0;
    $nmarkers++;

    # substitute complementary alleles to simplify comparisons later
    for(my $c1=0;$c1<(@genotypes+0);$c1++) {
				$genotypes[$c1] =~ s/(a|t)/a/ig;
				$genotypes[$c1] =~ s/(c|g)/c/ig;
				$genotypes[$c1] =~ s/^.?[\-n].?$/NN/i; # any null value is considered null on both chromosomes
    }
    if(!$summarise){
        printf $mw, $marker;
    }
    $delta = 0;
    $endStat = 0;
    foreach my $pop1 (keys(%indivs)){
	$pA{$pop1} = 0;
	$pC{$pop1} = 0;
        for (@{$indivs{$pop1}}){
            $gt = $genotypes[$_];
            if ($gt eq 'aa'){
                $pA{$pop1} += 2;
            }
            elsif (($gt eq 'ac') || ($gt eq 'ca')){
                $pA{$pop1} += 1;
                $pC{$pop1} += 1;
            }
            elsif ($gt eq 'cc'){
                $pC{$pop1} += 2;
            }
	    elsif (($gt eq 'NN') || ($gt eq '--')){
		# null / not typed
                # $pA{$pop1} += 0;
                # $pC{$pop1} += 0;
	    }
            else{
                printf $sw, "ERR-'".$gt."'";
            }
	}
	$pA{$pop1} = $pA{$pop1} / ((@{$indivs{$pop1}} + 0)*2);
	$pC{$pop1} = $pC{$pop1} / ((@{$indivs{$pop1}} + 0)*2);
	if(($pA{$pop1} == 0) && ($pC{$pop1} == 0)){
	    $pp = 0;
	    $badData = 1;
	}
	else{
	    $pp =
		abs($pA{$pop1} - $pC{$pop1}) /
		max($pA{$pop1}, $pC{$pop1});
	}
	if($summarise){
	    $summaries[$spos++] += $pA{$pop1};
	    $summaries[$spos++] += $pC{$pop1};
	    $summaries[$spos++] += $pp;
	}
	else{
	    printf $nw, $pA{$pop1};
	    printf $nw, $pC{$pop1};
	    printf $nw, $pp;
	}
        $endStat += $pp;
        if (!$delta){
            $delta = $pA{$pop1};
        }
        else{
            $delta = abs($delta - $pA{$pop1});
        }
    }
    $endStat = $endStat / (@pops + 0);
    $mulStat = $endStat * $delta;
    $endStat = $endStat + $delta;
    if($badData){
	$delta = 0;
	$endStat = 0;
	$mulStat = 0;
    }
    if($summarise){
	    $summaries[$spos++] += $delta;
	    $summaries[$spos++] += $endStat;
    }
    else{
	printf $nw, $delta;
	printf $nw, $endStat;
	printf $nw."\n", $mulStat;
    }
}

if($summarise){
    printf $mw, "Summary:";
    foreach my $i (0 .. (@summaries - 1)) {
        printf $nw, ($summaries[$i]/ $nmarkers);
    }
    print "\n";
}
