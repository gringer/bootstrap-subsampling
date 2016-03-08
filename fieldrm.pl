#!/usr/bin/perl

# fieldrm.pl -- removes specified fields from a file

# Author: David Eccles (gringer), 2008 <programming@gringer.org>

use strict;
use warnings;

sub usage {
  print("usage: ./fieldrm.pl [f1 f2 f3... fn] [options] <file>\n");
  print("\nOther Options:\n");
  print("-startat [n] : Start removing fields from this field onwards\n");
  print("\n");
}

my @removecols = ();

# extract command line arguments
while(@ARGV && $ARGV[0] =~ /^[0-9]+$/){
    push(@removecols, shift @ARGV);
}

my $modifier = 0;
if(@ARGV && (($ARGV[0] eq "-startat") || ($ARGV[0] eq "-s"))){
    shift @ARGV;
    $modifier = shift(@ARGV);
}

my @line = ();

while (<>){
    @line = split(/\s+/);
    my $numcols = @line;
    for(@removecols){
        if(($_ + $modifier) >= $numcols){
            die(sprintf("Error: column %d exceeds number of columns (%d)",
                        ($_ + $modifier),$numcols));
        }
        $line[($_ + $modifier)] = "";
    }
    my $linetext = join(" ",@line); # removes extra blank element whitespace
    $linetext =~ s/^\s+//;
    @line = split(/\s+/, $linetext);
    print join(" ",@line)."\n";
}
