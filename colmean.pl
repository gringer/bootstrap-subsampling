#!/usr/bin/perl

# colmean.pl -- averages the columns of an input file
# (second column onwards, starting from the second line).
# By default this groups by the first column in the file

# Author: David Eccles (gringer), 2009 <programming@gringer.org>

use strict;
use warnings;

sub usage {
  print("usage: ./colmean.pl <file name>\n");
  print("\n");
}

#my $dummy = 0;
my $nw = '%-0.6f '; #number width number formatting
my $sw = '%-17s '; #number width string formatting

my $firstLine = 1; # true

my $temprs = "";
my $headerRow = 0; # false
my @colcounts = ();
my @colnames = ();
my $numrows = 0;
my @columns = ();
my $rsnum = "";


while(<>){
    if($firstLine){
        ($rsnum, @colnames) = split(/\s+/);
        if ((@colnames > 1) && ($colnames[1] =~ /[a-zA-Z]/)){ # if there's a letter in the first column
            $headerRow = 1; # assume the input has a header row
        }
        for (my $i=0; $i < (@colnames); $i++){
            if($headerRow){
                $colcounts[$i] = 0;
                printf $sw, $colnames[$i];
            }
            else{
                $colcounts[$i] = $colnames[$i];
            }
        }
        if ($headerRow){
            print "\n";
        }
        if(!($headerRow)){
            $numrows = 1;
        }
        $firstLine = 0; # false
    }
    else{
        ($temprs, @columns) = split(/\s+/,$_);
        if (!($temprs eq $rsnum)){
            printf $sw, $rsnum;
            for (my $i=0; $i < (@colnames); $i++){
                $colcounts[$i] = $colcounts[$i] / $numrows;
                printf $nw, $colcounts[$i];
                $colcounts[$i] = 0;
            }
            print "\n";
            $numrows = 0;
            $rsnum = $temprs;
        }
        for (my $i=0; $i < (@colnames); $i++){
            $colcounts[$i] += $columns[$i];
        }
        $numrows++;
    }
}

printf $sw, $rsnum;
for (my $i=0; $i < (@colnames); $i++){
    $colcounts[$i] = $colcounts[$i] / $numrows;
    printf $nw, $colcounts[$i];
}
print "\n";
