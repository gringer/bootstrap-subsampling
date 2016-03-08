#!/usr/bin/perl

# flatten_tex.pl -- flattens a .tex file structure to remove include /
# input statements.

# Author: David Eccles (gringer), 2008 <programming@gringer.org>

use strict;
use warnings;

sub usage {
  print("usage: ./flatten_tex.pl <file name>\n");
  print("\n");
}

sub processFile {
    my ($inFile) = @_;
    if($inFile !~ /.tex$/){
        $inFile =~ s/$/.tex/;
    }
    local *TEXFILE; # allows the same FH to be used in a recursive method
    open(TEXFILE, "< $inFile") or die("File '$inFile' could not be opened");
    while(<TEXFILE>){
        if(m/\\input{(.*?)}/){
            my $nextFile = $1;
            processFile($nextFile);
            $_ = "";
        }
        print $_;
    }
    close(TEXFILE);
}

while(<>){
    if(/\\input{(.*?)}/){
        my $nextFile = $1;
        processFile($nextFile);
        $_ = "";
    }
    if(/\\include{(.*?)}/){
        my $nextFile = $1;
        processFile($nextFile);
        $_ = "";
    }
    print $_;
}
