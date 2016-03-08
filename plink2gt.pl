#!/usr/bin/perl

# plink2gt.pl -- Convert from plink rotated input files to simplegt-formatted
# files

# tped files are fairly close to the simplegt format. This script
# strips the genetic position and chromosome number, then removes
# spaces between the alleles for each person.

# Author: David Eccles (gringer), 2009 <programming@gringer.org>

use strict;
use warnings;

# use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

sub usage {
  print("usage: ./plink2gt.pl <.tped file name>\n");
  print("\n");
}

my $currentFile = "";
my @genotypes = ();
my $marker = "";
my $line = "";
my $lineCount = 0;
my $outFileName = "out_simplegt.txt";

while(<>){
    $line = $_;
    if($ARGV ne $currentFile){
        close(OUTFILE);
        $currentFile = $ARGV;
        printf(STDERR "Currently processing file: %s\n",$currentFile);
        if($currentFile =~ /^(.*)\.tped$/){
            my $famFile = $1.".tfam";
            if(!(-e $famFile)){
                $famFile = $1.".fam";
            }
            $outFileName = $1."_simplegt.txt";
            if(-e $outFileName){
                die("Error: $outFileName already exists, ".
                    "refusing to overwrite.");
            }
            printf(STDERR "Writing result to: %s\n",$outFileName);
            open(OUTFILE, ">", $outFileName);
            if(-e $famFile){
                printf(STDERR "Getting column IDs from: %s\n",$famFile);
                print(OUTFILE "## <Individual/Column IDs:");
                open(FAMFILE, "<", $famFile);
                while(<FAMFILE>){
                    my @famData = split(/\s+/);
                    my $famName = shift(@famData);
                    my $indName = shift(@famData);
                    print(OUTFILE " ".$famName."_".$indName." ");
                }
                close(FAMFILE);
                print(OUTFILE " > ##\n");
            }
        }
        if($currentFile =~ /^(.*)\.ped$/){
            my $mapFile = $1.".map";
            $outFileName = $1."_simplegt.txt";
            die("Error: cannot currently handle non-transposed ".
                "linkage files");
        }
    }
    @genotypes = split(/\s+/, $line);
    $lineCount++;
    if(($lineCount % 1000) == 0){
        print(STDERR ".");
    }
    splice(@genotypes,2,2); # remove genetic, physical position
    splice(@genotypes,0,1); # remove chromosome number
    $marker = shift(@genotypes);
    # put spaces between every second genotype and print result
    my $a = 0;
    print OUTFILE $marker." ".
        join("",map({if(($a++ + 1) % 2){$_}else{$_." "}} @genotypes))."\n";
 }
close(OUTFILE);
print(STDERR "done!\n");
