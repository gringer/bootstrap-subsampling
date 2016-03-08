#!/usr/bin/perl

# resort.pl -- sorts fields from a file, based on provided column
# numbers in arguments

# Warning: There is no error checking here, although column numbers
# outside the range of actual columns should result in errors by
# default

# Author: David Eccles (gringer), 2009 <programming@gringer.org>

use strict;
use warnings;

sub usage {
  print("usage: ./resort.pl <file name>\n");
  print("\n");
}

my @sortcols = ();

# extract command line arguments
while(@ARGV && $ARGV[0] =~ /^[0-9]+$/){
    push(@sortcols, shift @ARGV);
}

my $modifier = 0;
if(@ARGV && (($ARGV[0] eq "--startat") || ($ARGV[0] eq "-s"))){
    shift @ARGV;
    $modifier = shift(@ARGV);
}

my @line = ();

while (<>){
    @line = split(/\s+/); # separate line by whitespace
    my @newline = splice(@line,0,$modifier); # extract line up to start point
    foreach(@sortcols){
        push(@newline, $line[$_]); # add on sortcols; modifier not needed as
    }                              # start bit has already been spliced out
    print join(" ",@newline)."\n";
}
