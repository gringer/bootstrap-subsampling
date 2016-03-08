#!/usr/bin/perl

# addpops2structure.pl -- adds population identifier IDs to a
# structure file

# Author: David Eccles (gringer), 2008 <programming@gringer.org>

sub usage {
  print("usage: ./addpops2structure.pl <pop1Name><pop1Count>_".
        "<pop2Name><pop2Count>... < <file>\n");
  print("\n");
}


use strict;
use warnings;

my @poprefs = ();

# extract command line arguments
my $test = shift @ARGV;
$test =~ s/_/ /g;
$test =~ s/([^0-9])([0-9])/$1 $2/g;
@poprefs = split(/ /,$test);

my $curLine = 0;
my $curPop = 1;
my $changeAt = 1 + $curLine + @poprefs[($curPop-1)*2+1];

while (<>){
    if($curLine > 0){
        if($curLine >= $changeAt){
            $curPop++;
            $changeAt = $curLine + @poprefs[($curPop-1)*2+1];
        }
        s/^(....)/$1 $curPop /;
        print;
    }
    else{
        print;
    }
    $curLine++;
}
