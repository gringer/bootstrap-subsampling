#!/usr/bin/perl

# tfam_gender.pl -- modify gender data using sample data file

# Author: David Eccles (gringer), 2008 <programming@gringer.org>

# - modified Sep 2009 to use columns headings rather than assuming
#   [012] is gender. Also handles tab separated files (assuming
#   heading identifiers contain no whitespace).
# - Also can output phenotypes as pseudo-genotype data

use strict;
use warnings;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

sub usage {
  print("usage: ./tfam_gender.pl <sample file> [options] < <input file>\n");
  print("Note: sample file must have a header line including 'gender'");
  print("\nOther Options:\n");
  print("-help               : Only display this help message\n");
  print("-nowarn             : Don't warn about missing individuals in tfam file\n");
  print("-phenomarker <file> : Produce phenotypes as simplegt format\n");
  print("\n");
}

my %gender = ();
my %sampledata = ();
my %columndata = ();
my $sampleInFilename = "";
my $phenoFilename = "";
my @filenames = ();

my $genderWarn = 1; # true

# extract command line arguments
while(@ARGV){
    my $argument = shift @ARGV;
    if(-f $argument){ # file existence check
        if(!$sampleInFilename){
            print(STDERR "Assuming '$argument' contains gender information\n");
            $sampleInFilename = $argument;
        } else {
            push(@filenames, $argument);
        }
    } else {
        if($argument eq "-help"){
            usage();
            exit(0);
        }
        elsif($argument eq "-phenomarker"){
            $phenoFilename = shift @ARGV;
        }
        elsif($argument eq "-nowarn"){
            $genderWarn = 0; # false
        }
        else{
            printf(STDERR "Error: Unknown command-line option, '$argument'\n");
            usage();
            exit(4);
        }
    }
}

@ARGV = @filenames;

if(!$sampleInFilename){
    print(STDERR "Error: No valid sample data file given\n");
    usage();
    exit(1);
}

if($phenoFilename !~ /.tped$/){
    $phenoFilename =~ s/$/.tped/;
}

if(-f $phenoFilename){
    print(STDERR "\nError: Phenotype output file name ($phenoFilename) already exists.\n");
    print(STDERR "Please delete this file, or choose another file name.\n\n");
    usage();
    exit(2);
}


my $sampleInFile = new IO::Uncompress::Gunzip "$sampleInFilename" or
    die "Unable to open $sampleInFilename\n";

my $firstline = 1; # true
my @headings = ();

while(<$sampleInFile>){
    chomp;
    my $line = $_;
    my @data1 = split(/\s+/, $line);
    my @data2 = split(/\t/, $line);
    if($firstline){
        @headings = @data1;
        $firstline = 0; # false;
    } else {
        my @data = @data1;
        if((scalar(@data1) != scalar(@headings)) && (scalar(@data2) == scalar(@headings))){
            @data = @data2;
        }
        for(my $i = 1; $i < scalar(@headings); $i++ ){
            $sampledata{$data[0]}->{$headings[$i]} = $data[$i];
            $columndata{$headings[$i]}->{$data[0]} = $data[$i];
            if( $headings[$i] =~ /gender/i ){
                $gender{$data[0]} = $data[$i];
            }
#            printf(STDERR "%s = %s, ", $headings[$i], $data[$i]);
        }
#        print(STDERR "\n");
    }
}

my @indivOrder = ();

while(<>){
    if(/^([^\s]+)(.*\s)[012](\s[0-9]+)$/){
        push(@indivOrder, $1);
        if($gender{$1}){
            $_ = $1.$2.$gender{$1}.$3."\n";
        }
        elsif($genderWarn) {
            print(STDERR "Warning: No gender information for '$1'\n");
        }
    }
    if(!$phenoFilename){
        print $_;
    }
}

if($phenoFilename){
    open(PFILE, "> ".$phenoFilename) or
        die("Unable to create file $phenoFilename");
    print(STDERR "Now writing phenotypes to $phenoFilename...");
    foreach(keys(%columndata)){
        my $category = $_;
        my $line = "1 phd".uc($category)." 0 0 ";
        foreach (@indivOrder){
            my $indLabel = $_;
            if($gender{$indLabel}){ # work out if there is data for $indLabel
                my $data = ($columndata{$category}->{$indLabel});
                my $convertedData = $data;
                if($category =~ /^gender/i){
                    if    ($data == 1) { $convertedData = "AA"; }
                    elsif ($data == 2) { $convertedData = "CC"; }
                    else               { $convertedData = "NN"; }
                }
                if($category =~ /^age/i){
                    if    ($data eq ".") { $convertedData = "NN"; }
                    elsif ($data == 0)   { $convertedData = "NN"; }
                    elsif ($data < 5)    { $convertedData = "AA"; }
                    else                 { $convertedData = "CC"; }
                }
                $line = $line.$convertedData." ";
            } else {
                $line = $line."-- ";
            }
        }
        $line =~ s/\s+$//;
        print(PFILE $line."\n");
    }
    close(PFILE);
    print(STDERR "Done!\n");
}
