#!/usr/bin/perl

# 5_pop_minimal.pl -- Make sure all populations are different using
# a minimal set of information from the genotyped individuals
# practically, this means extracting additional lines from an input
# file until every individual within a population appears different to
# every individual from the other population(s) (or the total input
# file, whichever comes first)

# Author: David Eccles (gringer), 2007 <programming@gringer.org>

# arguments: ./5_pop_minimal.pl p<P1 name> <P1 ids> p<P2 name> <P2 ids> ...

use strict;
use warnings;

sub usage {
  print("usage: ./5_pop_minimal.pl p<P1 name> <P1 ids> ".
        "p<P2 name> <P2 ids> ...\n");
  print("\n");
}
my %indivs = ();

# remove command line arguments (so while(<>) works)
my $popname = 0;
my $indnum = 0;
my $totargs = $#ARGV;
foreach my $argnum (0 .. $totargs) {
    my $argval = shift @ARGV;
    if ($argval =~ /^p(.*)$/i){
	$popname = $1;
	$indivs{$popname} = ();
    }
    elsif (($popname) && $argval =~ /^[0-9]+$/) {
	push (@{$indivs{$popname}}, $argval);
    }
}

my $nw = '%-0.5f '; #number width string formatting
my $sw = '%-7s '; #number width string formatting

my $nmarkers = 0;
my $ncounts = 0;
my @diffcounts = ();
my @donepops = ();

# print out individual comparisons in each population, and set
# difference counts to zero


foreach my $pop1 (keys(%indivs)){
    # print $pop1.":\n"; # pop1 = first population to compare
    push (@donepops, $pop1);
    for (@{$indivs{$pop1}}){ # iterate through individuals in pop1
	# print $_." vs ( ";
	foreach my $pop2 (keys(%indivs)){  # pop2 = second population
	    # make sure not comparing to populations previously compared
	    if (!(grep { $_ eq $pop2 } @donepops)){
		# print $pop2.": ";
		for (@{$indivs{$pop2}}){
		    # print $_." ";
		    $diffcounts[$ncounts++] = 0;
		}
	    }
	}
	# print ")\n";
    }
    # print "\n";
}

# print "Total number of comparisons to do: $ncounts\n";

while(<>){
    my ($marker, @genotypes) = split(/\s/,$_);
    $nmarkers++;
    
    # substitute complementary alleles to simplify comparisons later
    for(my $c1=0;$c1<(@genotypes+0);$c1++) {
	$genotypes[$c1] =~ s/(a|t)/t/ig;
 	$genotypes[$c1] =~ s/(g|c)/c/ig;
    }
    print $marker."\n";
    $ncounts = 0;
    @donepops = ();
    foreach my $pop1 (keys(%indivs)){
	push (@donepops, $pop1);
	foreach my $c1 (@{$indivs{$pop1}}){ # c1 - column # for 1st comp
	    foreach my $pop2 (keys(%indivs)){
		if (!(grep { $_ eq $pop2 } @donepops)){
		    foreach my $c2 (@{$indivs{$pop2}}){
			# c2 - column # for 2nd comparison
			my $gt1 = $genotypes[$c1];
			my $gt2 = $genotypes[$c2];
			#same genotype (assuming a=t,c=g)
			if ($gt1 eq $gt2){
			    $diffcounts[$ncounts]++;
			}
			#same genotype, but switched around
			elsif (($gt1 eq "ac" && $gt2 eq "ca") ||
			       ($gt1 eq "ca" && $gt2 eq "ac")){
			    $diffcounts[$ncounts]++;
			}
			else {
			    # print $gt1." != ".$gt2;
			    # @diffcounts[$ncounts] += 0;
			}
			$ncounts++;
		    }
		}
	    }
	}
    }
    my $different = 1; # true
    for (@diffcounts) {
	if ($_ == $nmarkers){
	    $different = 0; # false
	}
    }
    if ($different){
	last;
    }
}

# Header line 1, shows number of markers
printf "Number of markers: %d\n", $nmarkers;

# Header line 2, shows second person comparison
printf "%-17s", "Pairwise: C1 vs.";

$ncounts = 0;
@donepops = ();

foreach my $pop1 (keys(%indivs)){
    push (@donepops, $pop1);
    foreach my $c1 (@{$indivs{$pop1}}){
	my $done = 0; #false
	foreach my $pop2 (keys(%indivs)){
	    if (!(grep { $_ eq $pop2 } @donepops)){
		for (@{$indivs{$pop2}}){
		    if (!$done){
			printf $sw, $c1;
			$done = 1; # true
		    }
		    else{
			printf $sw, "";
		    }
		}
	    }
	}
    }
}
print "\n";

# Header line 3, shows second person comparison
printf "%-17s", "Pairwise: C2";

$ncounts = 0;
@donepops = ();

foreach my $pop1 (keys(%indivs)){
    push (@donepops, $pop1);
    for (@{$indivs{$pop1}}){
	foreach my $pop2 (keys(%indivs)){
	    if (!(grep { $_ eq $pop2 } @donepops)){
		for (@{$indivs{$pop2}}){
		    printf $sw, $_;
		}
	    }
	}
    }
}
print "\n";

printf "%-17s", "Pairwise: ";

$ncounts = 0;
@donepops = ();

foreach my $pop1 (keys(%indivs)){
    push (@donepops, $pop1);
    for (@{$indivs{$pop1}}){
	foreach my $pop2 (keys(%indivs)){
	    if (!(grep { $_ eq $pop2 } @donepops)){
		for (@{$indivs{$pop2}}){
		    printf $nw, ($diffcounts[$ncounts++] / $nmarkers);
		}
	    }
	}
    }
}

print "\n";

# for(my $c1=0;$c1<(@indivs-1);$c1++) {
#     for(my $c2=$c1+1;$c2<(@indivs+0);$c2++) {
# 	printf $nw, ($diffcounts[$c1*(@indivs+0)+$c2] / $nmarkers);
#     }
# }
# print "\n";
