#!/usr/bin/perl

# hardyweinberg.pl -- Determines genotype frequencies for each marker,
# as well as expected frequencies under HW equilibrium. The
# calculations are done on a per-population basis.

# Author: David Eccles (gringer), 2008 <programming@gringer.org>

# Note: This algorithm has not been tested on more than two
# populations. It might work...

use strict;
use warnings;

sub usage {
  print("usage: ./hardyweinberg.pl p<P1 name> <P1 cols> ".
        "p<P2 name> <P2 cols>\n");
  print("\nOther Options:\n");
  print("-summarise : Only generate mean summary statistics\n");
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
				$summarise = 1; # true
    }
    elsif (($popname) && $argval =~ /^[0-9]+$/) {
				push (@{$indivs{$popname}}, $argval);
    }
}

my $nw = '%-0.7f '; #number width float formatting
my $sw = '%-8s '; #number width string formatting
my $mw = '%-17s '; #marker width string formatting

my $ncounts = 0;
my @pops = ();

# print out Column Titles

printf $mw, "Marker";

foreach my $pop1 (keys(%indivs)){
    printf $sw, "o(AA|$pop1)";
    printf $sw, "o(AC|$pop1)";
    printf $sw, "o(CC|$pop1)";
    printf $sw, "e(AA|$pop1)";
    printf $sw, "e(AC|$pop1)";
    printf $sw, "e(CC|$pop1)";
    push(@pops, $pop1);
}

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

my %pAA = ();
my %pAC = ();
my %pCC = ();
my $fA = 0;
my $fC = 0;
my $gt;
my $pp;
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
				$genotypes[$c1] =~ s/^.?[\-n].?$/NN/i; # any null value is considered
    }                                          # null on both chromosomes
    if(!$summarise){
				printf $mw, $marker;
    }
    $endStat = 0;
    foreach my $pop1 (keys(%indivs)){
				$pAA{$pop1} = 0;
				$pAC{$pop1} = 0;
				$pCC{$pop1} = 0;
        for (@{$indivs{$pop1}}){
            $gt = $genotypes[$_];
            if ($gt eq 'aa'){
                $pAA{$pop1} += 1;
            }
            elsif (($gt eq 'ac') || ($gt eq 'ca')){
                $pAC{$pop1} += 1;
            }
            elsif ($gt eq 'cc'){
                $pCC{$pop1} += 1;
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
				if($summarise){
						$summaries[$spos++] += $pAA{$pop1};
						$summaries[$spos++] += $pAC{$pop1};
						$summaries[$spos++] += $pCC{$pop1};
						$pAA{$pop1} = $pAA{$pop1} / (@{$indivs{$pop1}} + 0);
						$pAC{$pop1} = $pAC{$pop1} / (@{$indivs{$pop1}} + 0);
						$pCC{$pop1} = $pCC{$pop1} / (@{$indivs{$pop1}} + 0);
						$fA = $pAA{$pop1} + (0.5 * $pAC{$pop1});
						$fC = $pCC{$pop1} + (0.5 * $pAC{$pop1});
						$summaries[$spos++] += ($fA * $fA);
						$summaries[$spos++] += (2*($fA * $fC));
						$summaries[$spos++] += ($fC * $fC);
				}
				else{
						printf $nw, $pAA{$pop1};
						printf $nw, $pAC{$pop1};
						printf $nw, $pCC{$pop1};
						$pAA{$pop1} = $pAA{$pop1} / (@{$indivs{$pop1}} + 0);
						$pAC{$pop1} = $pAC{$pop1} / (@{$indivs{$pop1}} + 0);
						$pCC{$pop1} = $pCC{$pop1} / (@{$indivs{$pop1}} + 0);
						$fA = $pAA{$pop1} + (0.5 * $pAC{$pop1});
						$fC = $pCC{$pop1} + (0.5 * $pAC{$pop1});
						printf $nw, ($fA * $fA * (@{$indivs{$pop1}} + 0));
						printf $nw, (2*($fA * $fC) * (@{$indivs{$pop1}} + 0));
						printf $nw, ($fC * $fC * (@{$indivs{$pop1}} + 0));
				}
    }
		if(!$summarise){
				print "\n";
		}
}

if($summarise){
    printf $mw, "Summary:";
    foreach my $i (0 .. (@summaries - 1)) {
				printf $nw, ($summaries[$i]/ $nmarkers);
    }
    print "\n";
}
