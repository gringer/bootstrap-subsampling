#!/usr/bin/perl

# 8_modelfit.pl -- scans a marker pairwise comparison file (output
# like 2_maximal.pl) for the next marker that, in conjunction with
# previously discovered markers, best fits a reference model

# Author: David Eccles (gringer), 2007 <programming@gringer.org>

sub usage {
  print("usage: ./8_modelfit.pl <reference model> ".
        "<iteration set> <PW file>\n");
  print("\n");
}

use strict;
use warnings;

my $refmodel = shift @ARGV;
my $itset = shift @ARGV;


my $numcomparisons = 0;
my @rlabels = ();
my @rpwvalues = ();
my @seenmarkers = ();
open RFILE, "< $refmodel";
while(<RFILE>){
    if (!(/^Pairwise/ || /^Number/)){
        next;
    }
    if (/^Number.*: ([0-9]+)$/){
        $numcomparisons = $1;
    }
    elsif (/^Pairwise: P1 vs. (.*)$/){
        @rlabels = split(/\s+/,$1);
    }
    elsif (/^Pairwise: P2 +.* +([^ ]+) +$/){
        push (@rlabels,$1);
    }
    elsif (/^Pairwise:  +(.*)$/){
        @rpwvalues = split(/\s+/,$1);
    }
}
close (RFILE); 

my $nummarkers = 0;
my @ipwvalues = zeroes($numcomparisons); #create array, initialise at zero
open IFILE, "< $itset";
while(<IFILE>){
    if (!(/^rs/)){
        next;
    }
    my ($marker, @pwvalues) = split(/\s/,$_);
    push(@seenmarkers,$marker);
    matadd(\@ipwvalues,\@pwvalues); # array add
    $nummarkers++;
}
close (IFILE);

@rpwvalues *= ($nummarkers + 1); # multiply by scalar
@rpwvalues -= @ipwvalues;


my $bestmarker;
my @bestpwvalues;
my $besttest;

while (<>){
    if (!(/^rs/)){
        next;
    }
    my ($marker, @pwvalues) = split(/\s/,$_);
    if (!$bestmarker){
        $bestmarker = $marker;
        @bestpwvalues = @pwvalues;
        $besttest = sum(abs(@rpwvalues - @bestpwvalues));
        # array subtract, absolute values (arraywise), sum(of array)
    }
    else{
        my $newtest = sum(abs(@rpwvalues - @pwvalues));
        # array subtract, absolute values (arraywise), sum(of array)
        if($besttest == $newtest){
            print "$bestmarker and $marker give the same result";
        }
        elsif($newtest < $besttest){
            $bestmarker = $marker;
            @bestpwvalues = @pwvalues;
            $besttest = sum(abs(@rpwvalues - @bestpwvalues));
        }
    }
}

printf '%-17s', $bestmarker;
print join(" ",@bestpwvalues);
