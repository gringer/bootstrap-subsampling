#!/usr/bin/perl

# firstn.pl -- extracts the first n lines with a given field repeated

# Author: David Eccles (gringer), 2009 <programming@gringer.org>

use warnings;
use strict;

sub usage {
  print("usage: ./firstn.pl [options] <filename>\n");
  print("\nOther Options:\n");
  print("-t <string>  : Define field separator\n");
  print("-f <integer> : Field to consider for repetitions\n");
  print("-n <integer> : Maximum count per identical field\n");
  print("\n");
}

my $field = 0;
my $countLimit = 1;
my $separator = "\\s+"; # defaults to whitespace

my @ARGFiles = ();
# extract command line arguments
while($_ = shift @ARGV){
    if(-f){
        push(@ARGFiles, $_);
    } else {
        if($_ eq "-t"){ # field separator
            $separator = shift @ARGV;
        }
        if($_ eq "-f"){ # field to consider
            $field = shift @ARGV;
        }
        if($_ eq "-n"){ # maximum count per identical field
            $countLimit = shift @ARGV;
        }
    }
}

@ARGV = @ARGFiles;

my $firstField = 1; # true
my $oldTest;
my $count;

while(<>){
    my @fields = split($separator);
    if(!$firstField && ($fields[$field] eq $oldTest)){
        if($count++ < $countLimit){
            print;
        }
    } else {
        $firstField = 0; # false
        $oldTest = $fields[$field];
        $count = 1;
        print;
    }
}
