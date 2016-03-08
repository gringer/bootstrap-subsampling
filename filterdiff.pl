#!/usr/bin/perl

# filterdiff.pl -- Filters out similar genotype lines from a file.

# Author: David Eccles (gringer), 2009 <programming@gringer.org>


use warnings;
use strict;

sub usage {
  print("usage: ./filterdiff.pl [options]\n");
  print("\nFilters out similar genotype lines from a file\n");
  print("\nOther Options:\n");
  print("-help        : Only display this help message\n");
  print("-v <float>   : Threshold value for inclusion (Default: 0.5)\n");
  print("\n");
}

# Calculations for marker linkage (D):
# Consider two markers with alleles a/c for marker 1, and A/C for
# marker 2. It is assumed that a/A and c/C correspond to the major and
# minor alleles of the two mutations respectively. The correlation
# between these markers can be determined by weights in a 3x3 table:
#     AA  AC  CC
# aa |+1 |0  |-1
# ac |0  |0* |0
# cc |-1 |0  |+1

# These weights are summed up over all individuals, then divided by
# the number of [good] genotypes. This is a similar calculation to
# that of the kendall tau rank correlation coefficient. Markers that
# are correlated in all individuals (i.e. fully linked) will have a
# coefficient of 1. Unlinked markers have a coefficient of 0. A
# coefficient of less than 0 would indicate that major/minor allele
# allocation was incorrect.

# * The dual heterozygote poses a problem with unphased data, as two possible
# states exist, one where linkage is present, and one where it is absent:
#  a||c (1)     a||c (-1)
#  A||C         C||A
# [all other combinations of biallelic genotypes are unambiguous, even
# with unphased genotypes]
# In this case, both states are considered equally likely, so the
# dual-heterozygote is given a weighting of the mean of these values
# (0). Alternative "linked if possible" or "unlinked if possible"
# assumptions can be made (weighting of 1 or -1 respectively).

# Two modifications to D are suggested in the literature to remove the
# dependency on allele frequency :
# * D-prime = D/min(f(a)*f(C),f(c)*f(A))
# * r^2 = D^2 / f(a)*f(c)*f(A)*f(C)

# [see P.W. Hedrick and S. Kumar (2001). "Mutation and linkage disequilibrium in human mtDNA". Eur. J. Hum. Genet. 9: 969–972. doi:10.1038/sj.ejhg.5200735.]

# The code below uses the second of these two modifications

sub alleleLinkage {
    my %counts = ();
    # set up initial counts as <first SNP>,<second SNP>
    $counts{"aa,aa"} = 0;
    $counts{"ac,aa"} = 0;
    $counts{"cc,aa"} = 0;
    $counts{"aa,ac"} = 0;
    $counts{"ac,ac"} = 0;
    $counts{"cc,ac"} = 0;
    $counts{"aa,cc"} = 0;
    $counts{"ac,cc"} = 0;
    $counts{"cc,cc"} = 0;
    my ($lineData1, $lineData2) = @_;
    # determine minor/major alleles and adjust accordingly
    my $countA1 = ($lineData1 =~ tr/a//);
    my $countC1 = ($lineData1 =~ tr/c//);
    my $countA2 = ($lineData2 =~ tr/a//);
    my $countC2 = ($lineData2 =~ tr/c//);
    if((($countA1 > $countC1) && ($countA2 < $countC2)) ||
       (($countA1 < $countC1) && ($countA2 > $countC2))){
        $lineData2 =~ tr/ac/ca/;
    }
    # there's probably a quicker way to do this, basically lining up the
    # two strings and counting the differences, accounting for nulls
    my @genotypes1 = split(/\s+/, $lineData1);
    my @genotypes2 = split(/\s+/, $lineData2);
    my $numgts = scalar(@genotypes1);
    if($numgts == 0){
        die("not enough genotypes");
    }
    my $result = 0;
    my $goodCount = 0;
    for(my $i=0; $i < $numgts; $i++){
        my $gt1 = $genotypes1[$i];
        my $gt2 = $genotypes2[$i];
        if(($gt1 =~ /^[ac][ac]$/) &&
           ($gt2 =~ /^[ac][ac]$/)){
            $goodCount++;
            if($gt1 eq "aa"){
                $countA1 += 2;
            } elsif($gt1 eq "ac"){
                $countA1 += 1;
                $countC1 += 1;
            } else {
                $countC1 += 2;
            }
            if($gt2 eq "aa"){
                $countA2 += 2;
            } elsif($gt2 eq "ac"){
                $countA2 += 1;
                $countC2 += 1;
            } else {
                $countC2 += 2;
            }
        }
        $counts{$gt1 . "," . $gt2} ++;
    }
    my $sumCorrelation =
        $counts{"aa,aa"} - $counts{"aa,cc"} -
        $counts{"cc,aa"} + $counts{"cc,cc"};
    if($goodCount){
        $result = $sumCorrelation / $goodCount;
        # f(A1) = $countA1 / ($goodCount * 2); etc.
        $result = ($result * $result) /
            (($countA1   * $countC1   * $countA2   * $countC2)/
             ($goodCount * $goodCount * $goodCount * $goodCount * 16));
    }
    ## Older attempt -- 1 for linked, 0 for unlinked, 0.5 for heterozygotes
    # my $sumCorrelation =
    #     2 * $counts{"aa,aa"} + $counts{"aa,ac"} +
    #     $counts{"ac,aa"} + $counts{"ac,ac"} + $counts{"ac,cc"} +
    #     $counts{"cc,ac"} + 2 * $counts{"cc,cc"};
    # if($goodCount){
    #     $result = $differences / ($goodCount * 2);
    # }
}

my $threshold = 0.1;

while(@ARGV){
    if(-f $ARGV[0]){ # file existence check
        last;
    } else {
        my $argument = shift @ARGV;
        if($argument eq "-help"){
            usage();
            exit(0);
        }
        if($argument eq "-v"){
            $threshold = shift @ARGV;
        } else {
            print(STDERR "unknown argument: $argument\n");
            usage();
            exit(1);
        }
    }
}

print(STDERR "Rejecting markers with r^2 > $threshold for any marker in the set\n");

my %inputLines = ();
my %actualLines = ();
my @addedMarkers = ();

my $lineCount = 0;

print(STDERR "Reading in genotypes");

while(<>){
    my $actualLine = $_;
    tr/GgCTtAN/cccaaa-/; # translate to complementary base
    s/ca/ac/ig; # reorder heterozygotes
    my $marker1 = "";
    my $lineData = "";
    if(/^([^\s]*)\s+(.*)$/){
        $marker1 = $1;
        $lineData = $2;
    } else {
        last;
    }
    $lineCount++;
    if($lineCount % 100 == 0){
        print(STDERR ".");
    }
    my $rejected = 0; # false
    foreach(@addedMarkers){
        my $marker2 = $_;
        my $linkage = alleleLinkage($inputLines{$marker2}, $lineData);
        if($linkage > $threshold){
            $rejected = 1; # true
            last;
#            print(STDERR "rejecting $marker1: too similar to $marker2\n");
        } else {
#            print(STDERR "$marker1 vs $marker2 : $difference\n");
        }
    }
    if(!$rejected){
        push(@addedMarkers, $marker1);
        $inputLines{$marker1} = $lineData;
        $actualLines{$marker1} = $actualLine;
        print(STDERR "[including '$marker1']\n");
    }
}
print(STDERR "done!\n");

foreach(@addedMarkers){
    print($actualLines{$_});
#    print($_." ".$inputLines{$_}."\n");
}
