#!/usr/bin/perl

# gt2plink.pl -- Convert from simplegt-formatted file to plink rotated
# input files

# Author: David Eccles (gringer), 2008 <programming@gringer.org>

# derived from rsfilter.pl / gt2plink.r

use strict;
use warnings;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

sub usage {
  print("usage: ./gt2plink.pl <input file> <map file> [options]\n");
  print("\nOther Options:\n");
  print("-help               : Only display this help message\n");
  print("-output             : output base file name\n");
  print("-t <character>      : map file separator character\n");
  print("\n");
}

my %markers = ();

my $genotypesInFilename = 0; # false
my $mapInFilename = 0; # false
my $mapSep = "\\s+"; # whitespace separated

my $outputBaseName = 0; # false

# extract command line arguments
while(@ARGV){
    my $argument = shift @ARGV;
    if(-f $argument){ # file existence check
        if($genotypesInFilename){
            $mapInFilename = $argument;
        } else {
            $genotypesInFilename = $argument;
        }
    } else {
        if($argument eq "-help"){
            usage();
            exit(0);
        }
        elsif($argument eq "-output"){
            $outputBaseName = shift @ARGV;
        }
        elsif($argument eq "-t"){
            $mapSep = shift @ARGV;
        }
        else{
            printf(STDERR "Error: Unknown command-line option, '$argument'\n");
            usage();
            exit(4);
        }
    }
}

if(!$genotypesInFilename){
    print(STDERR "Error: No valid genotype input file given\n");
    usage();
    exit(1);
}

if(!$mapInFilename){
    print(STDERR "Error: No valid map file given\n");
    usage();
    exit(2);
}

if(!$outputBaseName){
    if($genotypesInFilename =~ m/^(.*?)\./){
        # truncate at location of first '.'
        $outputBaseName = $1;
    } else {
        ## just in case the genotype input file doesn't contain '.'
        $outputBaseName = $genotypesInFilename;
    }
}

if((-f $outputBaseName.".tped") || (-f $outputBaseName.".tfam")){
    print(STDERR "Error: output tPED and/or output tFAM file already exists\n".
          "Please delete these before running this conversion program:\n".
          $outputBaseName.".tped\n".
          $outputBaseName.".tfam\n");
    usage();
    exit(3);
}


my $genotypesInFile = 0;
my $mapInFile = 0;

print(STDERR "Now reading entire map file ($mapInFilename)...");
# Note: this works with both gzipped files and plain text files files
$mapInFile = new IO::Uncompress::Gunzip "$mapInFilename" or
    die "Unable to open $mapInFilename\n";
my $hasHeader = 0; # false
# default assumption: plink .map file
my %colNums = (
    "Chromosome" => 0,
    "Marker", => 1,
    "Map" => 2,
    "Location" => 3
    );
my $foundMarkerCol = 0; # false

# arrays to store temporary "marker" sets until the marker column is known
my %markerLines = ();
my %csLines = ();
my %mapLines = ();
my %locLines = ();

my $linePos = 0;
my @lineData = ();
my @mapData = (); # array containing all data from the map file
                  # (excluding header)
while(<$mapInFile>){
#    print $_;
    if(!$linePos){ # if this is the first line
        $linePos++;
        push(@mapData, $_);
        @lineData = split(/$mapSep/, $_);
        if(!grep(/\d+/,@lineData)){
            # first line contains no numbers, assume it's a header
            pop(@mapData); # header, so remove from data array
            if(grep(/(Chromosome|Marker|ID|Location|Position)/,@lineData) < 3){
                print(STDERR "Error: File header does not contain ".
                      "necessary columns\n");
                print(STDERR "Please make sure it has at least ".
                      "Chromosome, Marker/ID, and Location/Position\n");
                usage();
                exit(4);
            }
            if(grep(/ID/,@lineData) > 1){
                print(STDERR "Error: File header has more than one ".
                      "ID column\n");
                print(STDERR "Please modify the file so that only one ".
                      "ID column is present\n");
                usage();
                exit(5);
            }
            foreach(@lineData){
            # any 'ID' columns will be renamed 'Marker'
                $_ =~ s/^.*ID.*$/Marker/;
            # any columns containing 'Map' will be renamed to just 'Map'
                $_ =~ s/^.*Map.*$/Map/;
            # any 'Position' columns will be renamed 'Location'
                $_ =~ s/^.*Position.*$/Location/;
            }
            %colNums = ();
            for(my $colNum = 0; $colNum < (@lineData); $colNum++){
                $colNums{$lineData[$colNum]} = $colNum;
            }
            $hasHeader = 1; # true
            $foundMarkerCol = 1; # true
        } else { # digits in the first line, assume it's data
            # make best guess at what columns are
            if((@lineData) == 3){
                # only 3 columns, so remove "Map"
                %colNums = (
                    "Chromosome" => 0,
                    "Marker" => 1,
                    "Location" => 2
                    );
            }
            # Note: if there are more than 4 columns, the subsequent
            # ones will not be considered as Marker / Chromosome /
            # Location candidates
            $markerLines{$lineData[$colNums{"Marker"}]} = $linePos - 1;
            $csLines{$lineData[$colNums{"Chromosome"}]} = $linePos - 1;
            $locLines{$lineData[$colNums{"Location"}]} = $linePos - 1;
            if($colNums{"Map"}){
                $mapLines{$lineData[$colNums{"Map"}]} = $linePos - 1;
            }
        }
    } else { # not the first line
        $linePos++;
        push(@mapData, $_);
        @lineData = split(/$mapSep/, $_);
        if($foundMarkerCol){
            if($markerLines{$lineData[$colNums{"Marker"}]}){
                printf(STDERR "Error: duplicate data (%s) found in ".
                       "marker column (%d). Unable to continue.\n",
                       $lineData[$colNums{"Marker"}], $colNums{"Marker"});
                usage();
                exit(6);
            } else{
                $markerLines{$lineData[$colNums{"Marker"}]} = $linePos - 1;
            }
        } else {
            # not sure which column is the marker, so carry out a few checks
            @lineData = split(/$mapSep/, $_);
            my @markerCandidates = ();
            # check for unique values in candidate columns. If a
            # column has a value that has been seen before, it can't
            # be the marker column. Also, marker columns need some
            # text in them, so if a pure number is found, it can't be the
            # marker column
            if(%markerLines && (!$markerLines{$lineData[$colNums{"Marker"}]})
               && ($lineData[$colNums{"Marker"}] !~ /^(\d+\.?\d*|\.\d+)$/)){
                $markerLines{$lineData[$colNums{"Marker"}]} = $linePos - 1;
                push (@markerCandidates, "Marker");
            }
            if(%csLines && (!$csLines{$lineData[$colNums{"Chromosome"}]})
               && ($lineData[$colNums{"Chromosome"}] !~ /^(\d+\.?\d*|\.\d+)$/)){
                $csLines{$lineData[$colNums{"Chromosome"}]} = $linePos - 1;
                push (@markerCandidates, "Chromosome");
            }
            if(%locLines && (!$locLines{$lineData[$colNums{"Location"}]})
               && ($lineData[$colNums{"Location"}] !~ /^(\d+\.?\d*|\.\d+)$/)){
                $locLines{$lineData[$colNums{"Location"}]} = $linePos - 1;
                push (@markerCandidates, "Location");
            }
            if(%mapLines && (!$mapLines{$lineData[$colNums{"Map"}]})
               && ($lineData[$colNums{"Map"}] !~ /^(\d+\.?\d*|\.\d+)$/)){
                $mapLines{$lineData[$colNums{"Map"}]} = $linePos - 1;
                push (@markerCandidates, "Map");
            }
            # delete hash arrays for columns that aren't candidates
            %markerLines = () if(!grep(/Marker/,@markerCandidates));
            %csLines = () if(!grep(/Chromosome/,@markerCandidates));
            %locLines = () if(!grep(/Loc/,@markerCandidates));
            %mapLines = () if(!grep(/Map/,@markerCandidates));
            if((@markerCandidates) == 0){
                print(STDERR "Error: duplicate data found in all ".
                      "potential marker columns. Unable to continue.\n");
                usage();
                exit(6);
            }
            if((@markerCandidates) == 1){
                # only one candidate, so call that the marker column,
                # swapping the column number with the other one
                $foundMarkerCol = 1; # true
                if($markerCandidates[0] eq "Chromosome"){
                    my $tmpVal = $colNums{"Marker"};
                    $colNums{"Marker"} = $colNums{"Chromosome"};
                    $colNums{"Chromosome"} = $tmpVal;
                    %markerLines = %csLines;
                    %csLines = ();
                }
                if($markerCandidates[0] eq "Location"){
                    my $tmpVal = $colNums{"Marker"};
                    $colNums{"Marker"} = $colNums{"Location"};
                    $colNums{"Location"} = $tmpVal;
                    %markerLines = %locLines;
                    %locLines = ();
                }
                if($markerCandidates[0] eq "Map"){
                    my $tmpVal = $colNums{"Marker"};
                    $colNums{"Marker"} = $colNums{"Map"};
                    $colNums{"Map"} = $tmpVal;
                    %markerLines = %mapLines;
                    %mapLines = ();
                }
                # if the initial assumption of the Marker column is
                # correct, don't need to do anything
            }
        }
    }
}
close($mapInFile);
printf (STDERR "done (%d lines read)!\n", scalar(@mapData));

$genotypesInFile = new IO::Uncompress::Gunzip "$genotypesInFilename" or
    die "Unable to open $genotypesInFilename\n";

my $gtDataLength = 0;

printf (STDERR "Writing output .tped file (%s.tped)...", $outputBaseName);
open(TPEDFILE, "> ".$outputBaseName.".tped");
my $lineCount = 0;
my $writtenLines = 0;
my $marker = "";
my @indLabels = ();
my @phenoVals = ();
my $phenoNum = 1;
my %seenMarker = ();
while (<$genotypesInFile>){
    $lineCount++;
    if($lineCount % 1000 == 0){
        print (STDERR ".");
    }
    my $line = $_;
    if($line =~ /^##/){
        ## Determine individual labels
        ## This works even if more than one <ID> region is present in the
        ## header line, as might be the case in a 'join'ed file
        if($line =~ /IDs:\s+(.*?)\s*>/){
            @indLabels = ();
            @phenoVals = ();
        }
        while($line =~ /IDs:\s+(.*?)\s*>/){
            @lineData = split(/\s+/, $1);
            push(@indLabels, @lineData);
            ## generate an array of (@lineData) copies of $phenoNum
            push(@phenoVals, (($phenoNum) x scalar(@lineData)) );
            $phenoNum++;
            $line =~ s/^.*?>//;
        }
    } else{
        if($line =~ /^(.*?)[\s,]/){
            $marker = $1;
            if($seenMarker{$marker}){
                printf(STDERR "\nWarning: Already seen data for marker %s, found again on line %d. ".
                       "This repeated line will be ignored\n", $marker, $lineCount);
            } else {
                $line =~ s/^.*?[\s,]+//;
                chomp($line);
                if(!$gtDataLength){
                    $gtDataLength = length($line);
                }
                # checks for data corruption
                if(length($line) != $gtDataLength){
                    # carry out a basic check for input data oddness
                    printf(STDERR "\nWarning: Length of line %d of input (%d) doesn't match the expected length (%d). ".
                           "This line will be ignored\n", $lineCount, length($line), $gtDataLength);
                } else {
                    if(!@indLabels){
                        @lineData = split(/\s+/, $line);
                        @indLabels = (1 .. (@lineData));
                        @phenoVals = ((0) x scalar(@indLabels));
                    }
                    $seenMarker{$marker} = 1; # true
                    if (defined $markerLines{$marker}){
                        $writtenLines++;
                        @lineData = split(/$mapSep/, $mapData[$markerLines{$marker}]);
                        print(TPEDFILE $lineData[$colNums{"Chromosome"}]." ");
                        print(TPEDFILE $lineData[$colNums{"Marker"}]." ");
                        if($colNums{"Map"}){
                            print(TPEDFILE $lineData[$colNums{"Map"}]." ");
                        } else {
                            print(TPEDFILE "0 ");
                        }
                        print(TPEDFILE $lineData[$colNums{"Location"}]." ");
                        # add spaces between alleles in genotypes
                        $line =~ s/([^\s,])[\s,]/ $1 /g;
                        # add spaces between alleles in last genotype of line
                        $line =~ s/([^\s,])([^\s,])$/$1 $2/;
                        # anything that isn't AGCT (or 1234) gets recoded as N
                        $line =~ s/[^AGCT1234 ]/N/g;
                        print(TPEDFILE $line."\n");
                    } else {
                        print(STDERR "\nWarning: mutation $marker not found in mapfile, ".
                              "will not be included in output file\n");
                    }
                }
            }
        }
    }
}
close(TFAMFILE);
printf (STDERR "done (%d lines written)!\n", $writtenLines);

printf (STDERR "Writing output .tfam file (%s.tfam)...", $outputBaseName);
open(TFAMFILE, "> ".$outputBaseName.".tfam");
for(my $i = 0; $i < (@indLabels); $i++){
    printf(TFAMFILE "%s 1 0 0 0 %d\n", $indLabels[$i], $phenoVals[$i]);
}
close(TFAMFILE);
printf (STDERR "done (%d lines written)!\n", scalar(@indLabels));
