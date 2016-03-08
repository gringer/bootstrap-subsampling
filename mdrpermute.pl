#!/usr/bin/perl

# mdrpermute.pl -- permutes the case (last) column of a MDR formatted
# text file

# Author: David Eccles (gringer), 2009 <programming@gringer.org>

use warnings;
use strict;
use List::Util 'shuffle';

sub usage {
  print("usage: ./mdrpermute.pl <input file> <number of permutations>\n");
  print("\n");
}

sub cumulSum{
    my @input = @_;
    my $total = 0;
    my @output = ();
    foreach(@input){
        $total += $_;
        push(@output, $total);
    }
    return(@output);
}

my $mdrFileName = "";
my $numPermutations = 1;

while(@ARGV){
    my $argument = shift @ARGV;
    if(-f $argument){ # file existence check
        if($mdrFileName){
            printf(STDERR "Error: More than one file specified\n");
            usage();
            exit(1);
        } else {
            $mdrFileName = $argument;
        }
    } else {
        if($argument eq "-help"){
            usage();
            exit(0);
        }
        elsif($argument =~ /^[0-9]+$/){
            $numPermutations = $argument;
        }
        else{
            printf(STDERR "Error: Unknown command-line option, '$argument'\n");
            usage();
            exit(2);
        }
    }
}

if(!$mdrFileName){
    print(STDERR "Error: No valid input file given\n");
    usage();
    exit(3);
}


my $lineNumber = 0;
my @indClasses = ();
my $currentClass = "";

# determine how many of each class there are. This is equivalent to
# the following piped shell command:

# $ tail -n +2 MDR-SampleData.txt | perl -pe 's/^.*\s([^\s])+/$1/' | \
#        sort | uniq -c

open(MDRFILE, "< $mdrFileName");

while(<MDRFILE>){
    my $line = $_;
    if($lineNumber == 0){
        # print("*".$line."*\n");
    } else {
        if($line =~ /\s([^\s]+)[\r\n]+$/){  # preserves DOS/UNIX format
            push(@indClasses, $1);
        }
    }
    $lineNumber++;
}

if(!@indClasses){
    die("Error: no classes defined, cannot continue");
}

# randomise classes and output to 'permute<num>_<inputFileName>'

for(my $permuteNum = 0; $permuteNum < $numPermutations; $permuteNum++){
    $lineNumber = 0;
    seek(MDRFILE,0,0);
    open(PERMUTEFILE, "> permute".sprintf("%05d",$permuteNum)."_$mdrFileName");
    # Permute the Class assignments
    my @tempClasses = shuffle(@indClasses);
#    print(join(":",@indClasses)."\n");
#    print(join(":",@tempClasses)."\n");
    while(<MDRFILE>){
        my $lineStart = $_;
        my $lineEnd = "";
        if($lineNumber != 0){
            if($lineStart =~ s/(\s)([^\s]+)([\r\n]+)$/$1/){ # preserves
                $lineEnd = $3;                              # DOS/UNIX format
                $lineStart .= pop(@tempClasses);
            }
        }
        print(PERMUTEFILE $lineStart.$lineEnd);
        $lineNumber++;
    }
    close(PERMUTEFILE);
    if(@tempClasses){
        warn("Warning: classes still left to be assigned...");
    }
}

close(MDRFILE);
