#!/usr/bin/Rscript

## gt2plink.r -- Convert from simplegt-formatted file to plink rotated
## input files

## Author: David Eccles (gringer), 2008 <programming@gringer.org>

## simplegt-formatted input file
genotypes.inFile = FALSE;
map.inFile = FALSE;
output.baseName = FALSE;
map.sep = "";

usage <- function(){
  cat("usage: ./gt2plink.r",
      "<input file> <map file> [options]\n");
  cat("\nOther Options:\n");
  cat("-help               : Only display this help message\n");
  cat("-output             : output base file name\n");
  cat("-t <character>      : map file separator character\n");
  cat("\n");
}

argLoc <- grep("--args",commandArgs()) + 1; # hack to get around R v2.4
                                            # issue stopping
                                            # commandArgs(TRUE) from
                                            # working
while(!is.na(commandArgs()[argLoc])){
  if(file.exists(commandArgs()[argLoc])){ # file existence check
    if(genotypes.inFile != FALSE){
      map.inFile <- commandArgs()[argLoc];
    } else {
      genotypes.inFile <- commandArgs()[argLoc];
    }
  } else {
    if(commandArgs()[argLoc] == "-help"){
      usage();
      quit(save = "no", status=0);
    }
    if(commandArgs()[argLoc] == "-output"){
      output.baseName <- commandArgs()[argLoc+1];
      argLoc <- argLoc + 1;
    }
    if(commandArgs()[argLoc] == "-t"){
      map.sep <- commandArgs()[argLoc+1];
      argLoc <- argLoc + 1;
    }
  }
  argLoc <- argLoc + 1;
}

if((genotypes.inFile == FALSE) || (!file.exists(genotypes.inFile))){
  cat("Error: No valid genotype input file given\n\n");
  usage();
  quit(save = "no", status = 1);
}

if((map.inFile == FALSE) || (!file.exists(map.inFile))){
  cat("Error: No valid map file given\n\n");
  usage();
  quit(save = "no", status = 2);
}

if(output.baseName == FALSE){
  if(length(grep("\\.",genotypes.inFile)) > 0){
    ## truncate at location of first '.'
    output.baseName <- substring(genotypes.inFile, 1,
                                 min(regexpr("\\.",genotypes.inFile))-1);
  } else {
    ## just in case the genotype input file doesn't contain '.'
    output.baseName <- genotypes.inFile;
  }
}

output.tped <- paste(output.baseName,"tped",sep=".");
output.tfam <- paste(output.baseName,"tfam",sep=".");

if(file.exists(output.tped) || file.exists(output.tfam)){
  cat("Error: output tPED and/or output tFAM file already exists\n",
      "Please delete these before running this conversion program:\n",
      output.tped,"\n",
      output.tfam,"\n",
      "\n");
  usage();
  quit(save = "no", status = 3);
}

has.mutation <- FALSE;
map.skip = 0;
## read 2 lines to check for header, formatting, etc.
cat("Reading first 2 lines of map file...", file = stderr());
map.df <- read.table(gzfile(map.inFile), stringsAsFactors = FALSE,
                     sep = map.sep, nrows = 2);
if(length(grep("/",map.df)) > 0){
  has.mutation <- TRUE;
}
if(length(grep("[^0-9]",map.df$Location[1])) > 0){
  ## file probably has a header. Ignore it
  cat("Deciding to skip first line due to header\n", file = stderr());
  map.skip <- 1;
}
cat("done!\n", file = stderr());
## the real read procedure
cat("Now reading entire map file...", file = stderr());
map.df <- read.table(gzfile(map.inFile), sep = map.sep, skip = map.skip);
if(dim(map.df)[2] == 3) {
  colnames(map.df) <- c("Chromosome","Marker","Location");
  map.df$MapDist <- "0";
} else {
  if(has.mutation){
    colnames(map.df) <- c("Chromosome","Mutation","Marker","Location");
    map.df$MapDist <- "0";
  } else {
    ## assume 4 columns, third is the map distance
    colnames(map.df) <- c("Chromosome","Marker","MapDist","Location");
  }
}
map.df$Marker <- factor(map.df$Marker); # so levels(map.df$Marker) will work
## Check to make sure the correct decision was made about
## Chromosome/Marker ordering in map file. Incorrect ordering is
## suspected if there are 26 or fewer unique values in the "mutation"
## column *and* none of those values have more than 2 characters
## [the second condition protects against small map files]
if((length(levels(map.df$Marker)) <= 26) &&
   (max(nchar(levels(map.df$Marker))) < 6)){ ## allow for chrXX notation
  ## Marker/Chromosome labels are incorrect, so swap them
  colnames(map.df)[c(which(colnames(map.df) == "Chromosome"),
                     which(colnames(map.df) == "Marker"))] <-
                       c("Marker","Chromosome");
}
levels(map.df$Chromosome) <- sub("^chr","",levels(map.df$Chromosome));
if(length(levels(map.df$Marker)) != dim(map.df)[1]){
  cat("Error: mutation names are not unique, cannot continue\n");
  usage();
  quit(save = "no", status = 4);
}
## reference rows by mutation -- allows quick lookup later on
rownames(map.df) <- map.df$Marker;
cat("done!\n", file = stderr());

indSeq <- FALSE;
indLabels <- FALSE;
numIndivs <- FALSE;
phenoVals <- numeric(0);
phenoNum <- 1;
linesDone <- 0;

cat("Writing output .tped file...", file = stderr());
output.tped.file <- file(output.tped);
open(output.tped.file, open = "w+");
genotypes.gzFile <- gzfile(genotypes.inFile);
open(genotypes.gzFile);
z <- "dummy value";
while(length(z) > 0){
  z <- scan(genotypes.gzFile, nlines = 1, what = character(), quiet = TRUE);
  linesDone <- linesDone + 1;
  ## every 1000 lines, output another '.'
  if((linesDone %% 1000) == 0){
    cat(".", file = stderr());
  }
  if((length(z) > 0) && (z[1] == "##")){
    indPos <- grep("IDs:",z);
    endPos <- grep(">",z);
    ## Determine individual labels
    ## This works even if more than one <ID> region is present in the
    ## header line, as might be the case in a 'join'ed file
    if(length(indPos) > 0){
      indSeq <- integer(0);
      indLabels <- character(0);
      for(x in 1:length(indPos)){
        indSeq <- (indPos[x] + 1):(endPos[x] - 1);
        indLabels <- c(indLabels,z[indSeq]);
        phenoVals <- c(phenoVals,rep(phenoNum,length(indSeq)));
        phenoNum <- phenoNum + 1;
      }
      numIndivs <- length(indLabels);
    }
  } else {
    if(length(z) > 0){
      if(numIndivs == FALSE){
        numIndivs <- length(z) - 1;
      }
      mapLine <-
        as.matrix(map.df[z[1],c("Chromosome","Marker","MapDist","Location")]);
      if(is.na(mapLine[1])){
        cat("Warning: mutation '",z[1],
            "' not found in mapfile, will not be included in output file\n",
            sep="", file = stderr());
      } else {
        z <- sub("^(.)","\\1 ",z[-1]);
        cat(mapLine, "",file = output.tped.file); ## output marker details
        cat(z, file = output.tped.file); ## output genotypes
        cat("\n", file = output.tped.file);
      }
    }
  }
}
close(genotypes.gzFile);
close(output.tped.file);
cat("done!\n", file = stderr());
## In the absence of individual labels, label as 1..numIndivs and
## treat all phenotypes as unknown
if(indLabels[1] == FALSE){
  indLabels <- 1:numIndivs;
  phenoVals <- rep(0,numIndivs);
}

cat("Writing output .tfam file...", file = stderr());
output.tfam.file <- file(output.tfam);
open(output.tfam.file, open = "w+");
out.data <- apply(cbind(indLabels,1,0,0,0,phenoVals),
                  1, paste, collapse = " ");
cat(out.data, sep = "\n", file = output.tfam.file);
close(output.tfam.file);
cat("done!\n", file = stderr());
