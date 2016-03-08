#!/usr/bin/perl

# ged2linkage.pl -- creates a linkage formatted pedigree file based on
# family definitions in a GEDCOM file.

# Author: David Eccles (gringer), 2008 <programming@gringer.org>

use warnings;
use strict;

sub usage {
  print("usage: ./ged2linkage.pl <file name>\n");
  print("\n");
}

my $inFamily = 0; # false
my $wife = "";
my $husband = "";
my $child = "";
my $indiv = "";
my %asChild = (); # stores if a child has been seen
my %idNums = (); # stores numbers associated with each individual
my %active = (); # stores if an individual has been seen
my %names = ();
my %sex = ();
my %trios = ();
my $indCounter = 1;
my $unknownCounter = 1;
my $doHeader = 0; # generate header
my $asLinkage = 0; # use more strict linkage format
my $intCoded = 1; # should the program re-encode individuals as integers
my $unknownSex = 0;
if($asLinkage){
    $unknownSex = 0;
}
if($doHeader){
    print "FamilyID\nIndividualID\nFather\nMother\nGender\nAffected\n\n";
}
my $t = [ 1,2,3 ];
while(<>){
    chomp;
    if(/^0 /){
        $inFamily = 0; # false
        $wife = "";
        $husband = "";
        $indiv = "";
    }
    if(/0 (.*) INDI\s*$/){
        $indiv = $1; # true
    }
    if($indiv && (/^1 NAME (.*?)\s*$/)){
        $names{$indiv} = $1;
    }
    if(/^1 SEX (.)/){
        $sex{$indiv} = $1;
    }
    if(/FAM\s*$/){
        $inFamily = 1; # true
    } elsif($inFamily){
        if(/1 HUSB (.*)$/){
            $husband = $1;
            $husband =~ s/\s+$//g;
        }
        if(/1 WIFE (.*)$/){
            $wife = $1;
            $wife =~ s/\s+$//g;
        }
        if(/1 CHIL (.*)$/){
            $child = $1;
            $child =~ s/\s+$//g;
            if(!$idNums{$child}){
                $idNums{$child} = $indCounter++;
            }
            if($wife eq ""){
                $wife = "UNK".sprintf("%04d",$unknownCounter++);
                $sex{$wife} = "F";
            }
            if($husband eq ""){
                $husband = "UNK".sprintf("%04d",$unknownCounter++);
                $sex{$husband} = "M";
            }
            if(!$sex{$child}){
                $sex{$child} = "U";
            }
            if(!$idNums{$husband}){
                $idNums{$husband} = $indCounter++; # true
            }
            if(!$idNums{$wife}){
                $idNums{$wife} = $indCounter++; # true
            }
            if(!$sex{$wife} || ($sex{$wife} eq "U")){
                warn("warning: sex of #".$idNums{$wife}." (parent of $child) unspecified, setting to female");
                $sex{$wife} = "F";
            }
            if(!$sex{$husband} || ($sex{$husband} eq "U")){
                warn("warning: sex of #".$idNums{$husband}." (parent of $child) unspecified, setting to male");
                $sex{$husband} = "M";
            }
                if($sex{$wife} eq "M"){
                warn("warning: sex of $wife is male, recorded as wife");
            }
                if($sex{$husband} eq "F"){
                warn("warning: sex of $husband is female, recorded as husband");
            }
            $asChild{$child} = 1; # true
            # only add husband/wife once a child is added
#            $active{$child} = 1; # true;
            $active{$husband} = 1; # true;
            $active{$wife} = 1; # true;
            if($intCoded){
                $trios{$child} = [$idNums{$child}, $idNums{$husband},
                                  $idNums{$wife}];
            } else {
                $trios{$child} = [$child, $husband, $wife];
            }
        }
    }
}

while(my ($indiv, $indid) =  each(%idNums)){
    if(!$asChild{$indiv}){
        if(!$sex{$indiv}){
            $sex{$indiv} = "U";
        }
        if($intCoded){
            $trios{$indiv} = [$idNums{$indiv},0,0];
        } else {
            $trios{$indiv} = [$indiv,0,0];
        }
    }
}

open(NAMEFILE, ">indiv_names.txt");
printf NAMEFILE ("%5s %7s %s\n", "IDNUM", "GEDID","Name");


foreach my $indiv (sort(keys(%idNums))){
    my @tmpTrio = @{ $trios{$indiv} };
    my $tmpSex = $sex{$indiv};
    $tmpSex =~ s/M/1/;
    $tmpSex =~ s/F/2/;
    $tmpSex =~ s/U/$unknownSex/;
    if(($tmpTrio[1]) || ($tmpTrio[2]) || $active{$indiv}){
        printf("1 %8s %8s %8s  %2d 1\n",$tmpTrio[0],
               $tmpTrio[1],$tmpTrio[2], $tmpSex);
        if($names{$indiv}){
            printf NAMEFILE ("%5d %7s %s\n", $idNums{$indiv},
                   $indiv, $names{$indiv});
        }
    }
}

close(NAMEFILE);
