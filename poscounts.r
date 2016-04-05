#!/usr/bin/Rscript

source("/itsshared/phd/common.r");

## Author: David Eccles (gringer), 2009 <programming@gringer.org>

## poscounts.r -- determine clinical parameters for structure file outputs
## TP: True positive (count of "positive" results that are clinically positive)
## FN: False negative (count of "negative" results that are clinically positive)
## TN: True negative (count of "negative" results that are clinically negative)
## FP: False positive (count of "positive" results that are clinically negative)

usage <- function(){
  cat("usage: ./poscounts.r <structure \"q\" file> -s <split point> [options]\n");
  cat("\nOther Options:\n");
  cat("  -r <value>       : cutoff resolution (default 0.1)\n");
  cat("  -min <value>     : minimum cutoff value (default 0)\n");
  cat("  -max <value>     : maximum cutoff value (default 1)\n");
  cat("  -flip            : invert positive/negative labels\n");
 ## Should probably have 'flop' to invert case/control labels as well
  cat("  -label <lA> <lO> : case/control labels\n");
  cat("\n");
  }

infile.name <- FALSE;
splitPoint <- FALSE; # location at which positives change to negatives
minCO <- 0; # minimum cutoff value
maxCO <- 1; # maximum cutoff value
COResolution <- 0.1; # cutoff resolution
flip <- FALSE; # whether to flip positive/negative labels
ccLabels <- NULL; # case/control labels

argLoc <- 1;
while(!is.na(commandArgs(TRUE)[argLoc])){
  if(file.exists(commandArgs(TRUE)[argLoc])){ # file existence check
    if(infile.name == FALSE){
      infile.name <- commandArgs(TRUE)[argLoc];
    } else{
      cat("Error: More than one input file specified\n");
      usage();
      quit(save = "no", status=1);
    }
  }
  if(commandArgs(TRUE)[argLoc] == "-flip"){
      flip <- TRUE;
  }
  if(commandArgs(TRUE)[argLoc] == "-label"){
      ccLabels <- c(commandArgs(TRUE)[argLoc+1],commandArgs(TRUE)[argLoc+2]);
      argLoc <- argLoc+2;
  }
  if(commandArgs(TRUE)[argLoc] == "-s"){
      splitPoint <- as.numeric(commandArgs(TRUE)[argLoc+1]);
      argLoc <- argLoc+1;
  }
  if(commandArgs(TRUE)[argLoc] == "-r"){
      COResolution <- as.numeric(commandArgs(TRUE)[argLoc+1]);
      argLoc <- argLoc+1;
  }
  if(commandArgs(TRUE)[argLoc] == "-min"){
      minCO <- as.numeric(commandArgs(TRUE)[argLoc+1]);
      argLoc <- argLoc+1;
  }
  if(commandArgs(TRUE)[argLoc] == "-max"){
      maxCO <- as.numeric(commandArgs(TRUE)[argLoc+1]);
      argLoc <- argLoc+1;
  }
  if(length(grep("-h",commandArgs(TRUE)[argLoc])) > 0){
    usage();
    quit(save = "no", status = 0);
  }
  argLoc <- argLoc+1;
}

if(splitPoint == FALSE){
  cat("Error: No split point defined (use -s)\n");
  usage();
  quit(save = "no", status=1);
}

a <- read.table(infile.name, row.names = 1);

getValues <- function(cutoff){
  if(length(cutoff)>1){
    return(data.frame(t(sapply(unique(cutoff),getValues))));
  }
  TP <- length(which(a$V2[1:splitPoint] >= cutoff));
  FN <- length(which(a$V2[1:splitPoint] < cutoff));
  FP <- length(which(a$V2[-(1:splitPoint)] >= cutoff));
  TN <- length(which(a$V2[-(1:splitPoint)] < cutoff));
  if(!flip){
    return (data.frame(cutoff = cutoff, TP = TP, TN = TN, FP = FP, FN = FN));
  } else {
    return (data.frame(cutoff = cutoff, TP = FN, TN = FP, FP = TN, FN = TP));
  }
}

cat("\n");

cat("TP: Cases >= cutoff value\n");
cat("FP: Controls >= cutoff value\n");
cat("TN: Cases < cutoff value\n");
cat("FN: Controls < cutoff value\n");

cat("\n");

getValues(seq(minCO,maxCO,COResolution));

pdf("output_QHist.pdf", paper = "a4r", width = 11, height = 8);

if(length(ccLabels) == 0){
  ccLabels <- c("CASE","CONTROL");
}

hist(a$V2[1:splitPoint], breaks = 20, col = "red", xlab = "Q value",
     main = paste("Histogram of Q values for ",ccLabels[1],"",sep=""));
hist(a$V2[-(1:splitPoint)], breaks = 20, col = "red", xlab = "Q value",
     main = paste("Histogram of Q values for ",ccLabels[2],"",sep=""));

dummy <- dev.off();
