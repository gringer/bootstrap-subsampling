#!/usr/bin/perl

# snpimpute.pl -- infers unknown genotypes using a bayesian approach

# Author: David Eccles (gringer), 2009 <programming@gringer.org>

use warnings;
use strict;
my $progname = $0;

sub usage {
  print("usage: ./snpimpute.pl [options] < <file name>\n");
  print("\nOther Options:\n");
  print("-nohet    : don't do heterozygote calculations\n");
  print("-probs    : print prior/conditional probabilities\n");
  print("-inferall : infer all genotypes, rather than just unknown ones\n");
  print("\n");
}

# This ambitious program attempts to infer (impute) unknown genotypes
# based on the knowledge contained in other known genotypes. Using
# this program with a genome-wide set is probably not a good idea, as
# it will then use all SNPs to infer, which could take some time. The
# central idea of this process is Bayes' Theorem:
#
# p(B = b | A = a) = [p(A = a | B = b) * p(B = b)] / p(A = a)
#
# In this case, 'B = b' refers to the unknown genotype, 'A = a' refers
# to all the known SNP genotypes. If we assume that those other SNPs
# are independent, then the total conditional probability is the
# product sum of the conditional probabilities of each individual SNP
# genotype.
#
# The process can be carried out on both an allele-frequency level, as
# well as a genotype-frequency level. Note that both models assume
# bimorphic SNPs. With the allele-frequency model, a genotype is
# assumed to be the heterozygote when the predicted probabilities of
# both alleles are similar. With the genotype frequency model, the
# most likely genotype is chosen as the predicted genotype.
#
# This depends on the concept of a prior probability, p(B = b), which
# can be either population independent, or population dependent. In
# the population dependent case, the prior probability is generated
# based on known counts for members of a similar population (a
# potential extention is to contrast this with known counts for
# members of different population(s)). In the simpler population
# independent case, all individuals are used for this prior probability.
#
# NOTE: The current program determines prior probabilities based on
# the total dataset (population independent).

my $folddiff = 4; # how different A/C probabilities need to be to be
# considered not a heterozygote

# Read in command line arguments

# option defaults
my $doHet = 1; # true
my $printProbs = 0; # false
my $inferAllGTs = 0; # false

# read command line arguments until a file name is found
while($ARGV[0] && !(-e $ARGV[0])){
    $_ = shift(@ARGV);
    if($_ eq "-nohet"){ # don't do heterozygote calculations
        $doHet = 0;
    }
    if($_ eq "-probs"){ # print prior/conditional probabilities
        $printProbs = 1;
    }
    if($_ eq "-inferall"){ # infer all genotypes, rather than just unknown ones
        $inferAllGTs = 1;
    }
    else{
        print "ERROR: Unknown argument \"$_\" or file not found\n";
        usage();
        exit(1);
    }
}

# Read in lines of input file, storing genotypes for each individual

# the following is derived from snpchip2structure.pl
# it generates a list of genotypes, grouped by marker

my %genomarker = ();
my $numindivs = 0;

while (<>){
    chomp;
    my $line = $_;
    $line =~ s/(a|t)/1/ig;
    $line =~ s/(g|c)/2/ig;
    $line =~ s/ 21/ 12/ig;
    $line =~ s/ [^0-9][^0-9]/ 00/ig;
    my ($marker, @genotypes) = split(/\s/,$line);
    @{$genomarker{$marker}} = @genotypes;
    if($numindivs && ($numindivs != (@genotypes))){
        die("Number of genotypes does not match");
    }
    else{
        $numindivs = (@genotypes);
    }
}

# the following generates a list of genotypes, grouped by individual

my @markers = keys(%genomarker);
my @genoindivs = ();

for(my $i=0; $i < $numindivs; $i++){
    my %genotypes = ();
    foreach my $marker (@markers){
        ${$genoindivs[$i]}{$marker} = @{$genomarker{$marker}}[$i];
    }
}


# determine allele frequency counts for each SNP (this will be the
# prior probabilities)
my %afreq = ();
my %hfreq = (); # for heterozygotes
my %cfreq = ();
foreach my $marker (@markers){
    $afreq{$marker} = 0;
    $hfreq{$marker} = 0;
    $cfreq{$marker} = 0;
    for(my $i=0; $i < $numindivs; $i++){
        my $genotype = @{$genomarker{$marker}}[$i];
        if($genotype eq "11"){
            $afreq{$marker} += 2;
        }
        if($genotype eq "12"){
            if($doHet){
                $hfreq{$marker} += 2;
            }
            else{
                $afreq{$marker} += 1;
                $cfreq{$marker} += 1;
            }
        }
        if($genotype eq "22"){
            $cfreq{$marker} += 2;
        }
    }
    # assume not seen equates to a count of 0.5
    # (half minimum value)
    if($doHet){ # 1 when heterozygotes are considered
        $afreq{$marker} = 1 unless $afreq{$marker};
        $hfreq{$marker} = 1 unless $hfreq{$marker};
        $cfreq{$marker} = 1 unless $cfreq{$marker};
    }
    else{
        $afreq{$marker} = 0.5 unless $afreq{$marker};
        $cfreq{$marker} = 0.5 unless $cfreq{$marker};
        $hfreq{$marker} = 0;
    }
    # determine final frequencies
    my $totcount = $afreq{$marker} + $hfreq{$marker} + $cfreq{$marker};
    $afreq{$marker} = $afreq{$marker} / $totcount;
    $hfreq{$marker} = $hfreq{$marker} / $totcount;
    $cfreq{$marker} = $cfreq{$marker} / $totcount;
}

# for the bayesian stuff, need to determine counts, given the genotype
# at some location

# first refers to "removed" genotype, second refers to "current" genotype
my %aabayescount = ();
my %acbayescount = ();
my %cabayescount = ();
my %ccbayescount = ();

# the increased complication that heterozygote storage brings in...
my %ahbayescount = ();
my %chbayescount = ();
my %habayescount = ();
my %hhbayescount = ();
my %hcbayescount = ();

foreach my $marker1 (@markers){
    for(my $i=0; $i < $numindivs; $i++){
        my $gt1 = @{$genomarker{$marker1}}[$i]; # "removed" genotype
        if($gt1 eq "11"){
            foreach my $marker2 (@markers){          # "current" genotype
                my $gt2 = @{$genomarker{$marker2}}[$i];
                if($gt2 eq "11"){
                    ${$aabayescount{$marker1}}{$marker2} +=2;
                }
                if($gt2 eq "12"){
                    if($doHet){
                        ${$ahbayescount{$marker1}}{$marker2} +=2;
                    }
                    else{
                        ${$aabayescount{$marker1}}{$marker2} +=1;
                        ${$acbayescount{$marker1}}{$marker2} +=1;
                    }
                }
                if($gt2 eq "22"){
                    ${$acbayescount{$marker1}}{$marker2} +=2;
                }
            }
        }
        elsif($gt1 eq "12"){
            foreach my $marker2 (@markers){          # "current" genotype
                my $gt2 = @{$genomarker{$marker2}}[$i];
                if($gt2 eq "11"){
                    if($doHet){
                        ${$habayescount{$marker1}}{$marker2} +=2;
                    }
                    else{
                        ${$aabayescount{$marker1}}{$marker2} +=1;
                        ${$cabayescount{$marker1}}{$marker2} +=1;
                    }
                }
                if($gt2 eq "12"){
                    if($doHet){
                        ${$hhbayescount{$marker1}}{$marker2} +=2;
                    }
                    else{
                        ${$aabayescount{$marker1}}{$marker2} +=0.5;
                        ${$cabayescount{$marker1}}{$marker2} +=0.5;
                        ${$acbayescount{$marker1}}{$marker2} +=0.5;
                        ${$ccbayescount{$marker1}}{$marker2} +=0.5;
                    }
                }
                if($gt2 eq "22"){
                    if($doHet){
                        ${$hcbayescount{$marker1}}{$marker2} +=2;
                    }
                    else{
                        ${$acbayescount{$marker1}}{$marker2} +=1;
                        ${$ccbayescount{$marker1}}{$marker2} +=1;
                    }
                }
            }
        }
        elsif($gt1 eq "22"){
            foreach my $marker2 (@markers){          # "current" genotype
                my $gt2 = @{$genomarker{$marker2}}[$i];
                if($gt2 eq "11"){
                    ${$cabayescount{$marker1}}{$marker2} +=2;
                }
                if($gt2 eq "12"){
                    if($doHet){
                        ${$chbayescount{$marker1}}{$marker2} +=2;
                    }
                    else{
                        ${$cabayescount{$marker1}}{$marker2} +=1;
                        ${$ccbayescount{$marker1}}{$marker2} +=1;
                    }
                }
                if($gt2 eq "22"){
                    ${$ccbayescount{$marker1}}{$marker2} +=2;
                }
            }
        }
    }
}

if($printProbs){
    # display prior probabilities (rows should add up to 1)
    print "Prior probabilities:\n\n";
    if($doHet){
        printf("%-3s %-12s %7s %7s %7s\n", "No.", "Marker",
               "p(gt=A)", "p(gt=H)", "p(gt=C)");
    }
    else{
        printf("%-3s %-12s %7s %7s\n", "No.", "Marker",
               "p(gt=A)", "p(gt=C)");
    }
    my $num = 0;
    foreach my $marker1 (@markers){
        if($doHet){
            printf("%03d %-12s %0.5f %0.5f %0.5f\n", $num, $marker1,
                   $afreq{$marker1}, $hfreq{$marker1}, $cfreq{$marker1});
        }
        else{
            printf("%03d %-12s %0.5f %0.5f\n", $num, $marker1,
                   $afreq{$marker1}, $cfreq{$marker1});
        }
        $num++;
    }
}

# determine conditional probabilities, based on counts
# e.g. afreq = acount / (acount + ccount)

my %aabayesfreq = ();
my %acbayesfreq = ();
my %cabayesfreq = ();
my %ccbayesfreq = ();

# heterozygote stuff again...
my %ahbayesfreq = ();
my %chbayesfreq = ();
my %habayesfreq = ();
my %hhbayesfreq = ();
my %hcbayesfreq = ();

foreach my $marker1 (@markers){
    foreach my $marker2 (@markers){
        # half minimum value is 0.25, so choose that for unobserved instances
        if ($doHet){ # 1 when heterozygotes are considered
            ${$aabayescount{$marker1}}{$marker2} = 1 unless
                ${$aabayescount{$marker1}}{$marker2};
            ${$acbayescount{$marker1}}{$marker2} = 1 unless
                ${$acbayescount{$marker1}}{$marker2};
            ${$cabayescount{$marker1}}{$marker2} = 1 unless
                ${$cabayescount{$marker1}}{$marker2};
            ${$ccbayescount{$marker1}}{$marker2} = 1 unless
                ${$ccbayescount{$marker1}}{$marker2};
            ${$ahbayescount{$marker1}}{$marker2} = 1 unless
                ${$ahbayescount{$marker1}}{$marker2};
            ${$chbayescount{$marker1}}{$marker2} = 1 unless
                ${$chbayescount{$marker1}}{$marker2};
            ${$habayescount{$marker1}}{$marker2} = 1 unless
                ${$habayescount{$marker1}}{$marker2};
            ${$hhbayescount{$marker1}}{$marker2} = 1 unless
                ${$hhbayescount{$marker1}}{$marker2};
            ${$hcbayescount{$marker1}}{$marker2} = 1 unless
                ${$hcbayescount{$marker1}}{$marker2};
        }
        else{
            ${$aabayescount{$marker1}}{$marker2} = 0.25 unless
                ${$aabayescount{$marker1}}{$marker2};
            ${$acbayescount{$marker1}}{$marker2} = 0.25 unless
                ${$acbayescount{$marker1}}{$marker2};
            ${$cabayescount{$marker1}}{$marker2} = 0.25 unless
                ${$cabayescount{$marker1}}{$marker2};
            ${$ccbayescount{$marker1}}{$marker2} = 0.25 unless
                ${$ccbayescount{$marker1}}{$marker2};
            ${$ahbayescount{$marker1}}{$marker2} = 0.25;
            ${$chbayescount{$marker1}}{$marker2} = 0.25;
            ${$hhbayescount{$marker1}}{$marker2} = 0.25;
            ${$habayescount{$marker1}}{$marker2} = 0.25;
            ${$hcbayescount{$marker1}}{$marker2} = 0.25;
        }

        ${$aabayesfreq{$marker1}}{$marker2} =
            ${$aabayescount{$marker1}}{$marker2} /
            (${$aabayescount{$marker1}}{$marker2} +
             ${$ahbayescount{$marker1}}{$marker2} +
             ${$acbayescount{$marker1}}{$marker2});
        ${$ahbayesfreq{$marker1}}{$marker2} =
            ${$ahbayescount{$marker1}}{$marker2} /
            (${$aabayescount{$marker1}}{$marker2} +
             ${$ahbayescount{$marker1}}{$marker2} +
             ${$acbayescount{$marker1}}{$marker2});
        ${$acbayesfreq{$marker1}}{$marker2} =
            ${$acbayescount{$marker1}}{$marker2} /
            (${$aabayescount{$marker1}}{$marker2} +
             ${$ahbayescount{$marker1}}{$marker2} +
             ${$acbayescount{$marker1}}{$marker2});
        ${$habayesfreq{$marker1}}{$marker2} =
            ${$habayescount{$marker1}}{$marker2} /
            (${$habayescount{$marker1}}{$marker2} +
             ${$hhbayescount{$marker1}}{$marker2} +
             ${$hcbayescount{$marker1}}{$marker2});
        ${$hhbayesfreq{$marker1}}{$marker2} =
            ${$hhbayescount{$marker1}}{$marker2} /
            (${$habayescount{$marker1}}{$marker2} +
             ${$hhbayescount{$marker1}}{$marker2} +
             ${$hcbayescount{$marker1}}{$marker2});
        ${$hcbayesfreq{$marker1}}{$marker2} =
            ${$hcbayescount{$marker1}}{$marker2} /
            (${$habayescount{$marker1}}{$marker2} +
             ${$hhbayescount{$marker1}}{$marker2} +
             ${$hcbayescount{$marker1}}{$marker2});
        ${$cabayesfreq{$marker1}}{$marker2} =
            ${$cabayescount{$marker1}}{$marker2} /
            (${$cabayescount{$marker1}}{$marker2} +
             ${$chbayescount{$marker1}}{$marker2} +
             ${$ccbayescount{$marker1}}{$marker2});
        ${$chbayesfreq{$marker1}}{$marker2} =
            ${$chbayescount{$marker1}}{$marker2} /
            (${$cabayescount{$marker1}}{$marker2} +
             ${$chbayescount{$marker1}}{$marker2} +
             ${$ccbayescount{$marker1}}{$marker2});
        ${$ccbayesfreq{$marker1}}{$marker2} =
            ${$ccbayescount{$marker1}}{$marker2} /
            (${$cabayescount{$marker1}}{$marker2} +
             ${$chbayescount{$marker1}}{$marker2} +
             ${$ccbayescount{$marker1}}{$marker2});
    }
}

if($printProbs){
    # display condiditional probabilities (columns within each 2x2 box
    # should add up to 1)
    print "\nConditional probabilities, probability that m1 = X, given m2 = Y:\n\n";
    my $num1 = 0;
    my $num2 = 0;
    # header row
    print "  \\m2 ";
    foreach my $marker1 (@markers){
        if($doHet){
            printf("%03d  %03d  %03d   ", $num1, $num1, $num1);
        }
        else{
            printf("%03d  %03d   ", $num1, $num1);
        }
        $num1++;
    }
    print "\n";
    # print alleles/genotypes
    print " m1\\  ";
    foreach my $marker1 (@markers){
        if($doHet){
            printf("%-4s %-4s %-4s  ", "A","H","C");
        }
        else{
            printf("%-4s %-4s  ", "A", "C");
        }
        $num1++;
    }
    print "\n";
    # and finally... the probabilities
    $num1 = 0;
    foreach my $marker1 (@markers){
        printf("%03d A ", $num1);
        $num2 = 0;
        foreach my $marker2 (@markers){
            if($doHet){
                printf("%0.2f %0.2f %0.2f  ", ${$aabayesfreq{$marker1}}{$marker2},
                       ${$habayesfreq{$marker1}}{$marker2},
                       ${$cabayesfreq{$marker1}}{$marker2});
            }
            else{
                printf("%0.2f %0.2f  ", ${$aabayesfreq{$marker1}}{$marker2},
                       ${$cabayesfreq{$marker1}}{$marker2});
            }
            $num2++;
        }
        print "\n";
        if($doHet){
            printf("%03d H ", $num1);
            $num2 = 0;
            foreach my $marker2 (@markers){
                printf("%0.2f %0.2f %0.2f  ", ${$ahbayesfreq{$marker1}}{$marker2},
                       ${$hhbayesfreq{$marker1}}{$marker2},
                       ${$chbayesfreq{$marker1}}{$marker2});
                $num2++;
            }
            print "\n";
        }
        printf("%03d C ", $num1);
        $num2 = 0;
        foreach my $marker2 (@markers){
            if($doHet){
                printf("%0.2f %0.2f %0.2f  ", ${$acbayesfreq{$marker1}}{$marker2},
                       ${$hcbayesfreq{$marker1}}{$marker2},
                       ${$ccbayesfreq{$marker1}}{$marker2});
            }
            else{
                printf("%0.2f %0.2f  ", ${$acbayesfreq{$marker1}}{$marker2},
                       ${$ccbayesfreq{$marker1}}{$marker2});
            }
            $num2++;
        }
        print "\n";
        $num1++;
    }
}

# spit out most likely genotype probabilities for individuals
# this uses Bayes' Theorem
# p(B = b | A = a) = [p(A = a | B = b) * p(B = b)] / p(A = a)

foreach my $marker1 (@markers){
    printf("%s ",$marker1);
    for(my $i=0; $i < $numindivs; $i++){
        if(!$doHet){ # necessary because log(0) is undefined
            $hfreq{$marker1} = 1;
        }
        my $aprob = log($afreq{$marker1}); # the prior prob that gt1 is A
        my $hprob = log($hfreq{$marker1}); # the prior prob that gt1 is H
        my $cprob = log($cfreq{$marker1}); # the prior prob that gt1 is C
        my $gtprob = log(1); # the probability of the known genotypes
        my %gts = %{$genoindivs[$i]};
        foreach my $marker2 (@markers){
            if($marker2 ne $marker1){
                my $gt2 = $gts{$marker2};
                if($gt2 eq "11"){
                    $aprob = $aprob + log(${$aabayesfreq{$marker1}}{$marker2});
                    $hprob = $hprob + log(${$habayesfreq{$marker1}}{$marker2});
                    $cprob = $cprob + log(${$cabayesfreq{$marker1}}{$marker2});
                    $gtprob = $gtprob + log($afreq{$marker2});
                }
                elsif($gt2 eq "12"){
                    if($doHet){
                        $aprob = $aprob + log(${$ahbayesfreq{$marker1}}{$marker2});
                        $hprob = $hprob + log(${$hhbayesfreq{$marker1}}{$marker2});
                        $cprob = $cprob + log(${$chbayesfreq{$marker1}}{$marker2});
                        $gtprob = $gtprob + log($hfreq{$marker2});
                    }
                    else{ # a hack for heterozygotes... geometric mean / log average
                        $aprob = $aprob + (log(${$aabayesfreq{$marker1}}{$marker2})
                                           + log(${$acbayesfreq{$marker1}}{$marker2}))/2;
                        $cprob = $cprob + (log(${$cabayesfreq{$marker1}}{$marker2})
                                           + log(${$ccbayesfreq{$marker1}}{$marker2}))/2;
                        $gtprob = $gtprob +
                            (log($afreq{$marker2}) + log($cfreq{$marker2})) / 2;
                    }
                }
                elsif($gt2 eq "22"){
                    $aprob = $aprob + log(${$acbayesfreq{$marker1}}{$marker2});
                    $hprob = $hprob + log(${$hcbayesfreq{$marker1}}{$marker2});
                    $cprob = $cprob + log(${$ccbayesfreq{$marker1}}{$marker2});
                    $gtprob = $gtprob + log($cfreq{$marker2});
                }
                # unknown values don't get included
            }
        }
        my $gtInf = "NN";
        # divide by p(A = a)
        $aprob = $aprob - $gtprob;
        $hprob = $hprob - $gtprob;
        $cprob = $cprob - $gtprob;
        my $gtratio = exp($aprob - $cprob);
        if($doHet){
            if(($aprob > $cprob) && ($aprob > $hprob)){
                $gtratio = ($folddiff*2);
            }
            elsif(($cprob > $aprob) && ($cprob > $hprob)){
                $gtratio = 1/($folddiff*2);
            }
            else{ # all other cases [should] have heterozygote most likely
                $gtratio = 1;
            }
        }
        if($gtratio > $folddiff){ # by trial/error with IMS/HBM dataset
            $gtInf = "AA";
        }
        elsif($gtratio < (1/$folddiff)){
            $gtInf = "CC";
        }
        else{
            $gtInf = "AC";
        }
        my $gtOld = $gts{$marker1};
        $gtOld =~ s/0/N/g;
        $gtOld =~ s/1/A/g;
        $gtOld =~ s/2/C/g;
        # for debugging (prints out log(probability)
        # printf("(%0.1f/%0.1f/%0.1f)%2s/%-3s", $aprob,$hprob,$cprob,$gtInf,$gtOld);
        if($inferAllGTs){
            printf("%2s/%-3s",$gtInf,$gtOld);
        }
        else{
            if($gtOld =~ /N/){
                printf("%-3s",$gtInf);
            }
            else{
                printf("%-3s",$gtOld);
            }
        }
    }
    print("\n");
}
