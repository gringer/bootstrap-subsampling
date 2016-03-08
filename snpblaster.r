#!/usr/bin/Rscript

## Author: David Eccles (gringer), 2007 <programming@gringer.org>

usage <- function(){
  cat("Snpblaster -- Calculates location differences to determine\n",
      "              which markers can be removed without much loss.\n",
      "usage: ./snpblaster.r <file> [-size <windowSize>]\n",
      "Expects a CSV file with headings: [Marker],Delta,Mutation,Chromosome,Location\n",
      sep="");
}

infile = FALSE;

windowSize = 10^6;
argLoc <- 1;

while(!is.na(commandArgs(TRUE)[argLoc])){
  if(file.exists(commandArgs(TRUE)[argLoc])){ # file existence check
    if(infile == FALSE){
      infile <- commandArgs(TRUE)[argLoc];
    } else{
      cat("Error: More than one input file specified\n");
      usage();
      quit(save = "no", status=1);
    }
  } else {
    if(commandArgs(TRUE)[argLoc] == "-window"){
      windowSize <- eval(parse(text = commandArgs(TRUE)[argLoc+1]));
      cat("Setting window size to ",windowSize,"\n",sep="");
      argLoc <- argLoc + 1;
    }
  }
  argLoc <- argLoc + 1;
}

if(infile == FALSE){
  cat("Error: No input file specified\n");
  usage();
  quit(save = "no", status=1);
}

meanstats.location <- read.csv(infile, row.names = 1);

num.rowcol <- dim(meanstats.location)[1];
snpblast.loc <- matrix(meanstats.location$Location, num.rowcol, num.rowcol);
snpblast.cs <- matrix(meanstats.location$Chromosome,num.rowcol, num.rowcol);
snpblast.delta <- matrix(meanstats.location$Delta,num.rowcol, num.rowcol);
snpblast.names <- matrix(rownames(meanstats.location),num.rowcol,
                        num.rowcol);
snpblast.mat <- abs(snpblast.loc-t(snpblast.loc));
snpblast.mat[snpblast.cs != t(snpblast.cs)] <- NA;
snpblast.mat[lower.tri(snpblast.mat, diag = TRUE)] <- NA;

near.markers <- which(snpblast.mat < windowSize);
snpblast.df <- data.frame(Marker1 = snpblast.names[near.markers], Marker2 = t(snpblast.names)[near.markers],
                         Chromosome = as.numeric(snpblast.cs[near.markers]),
                         Location1 = snpblast.loc[near.markers], Location2 = t(snpblast.loc)[near.markers],
                         Delta1 = snpblast.delta[near.markers], Delta2 = t(snpblast.delta)[near.markers],
                         Distance = snpblast.mat[near.markers]);
snpblast.df <- snpblast.df[order(as.numeric(snpblast.df$Chromosome),snpblast.df$Distance),];
snpblast.df$Remove <- as.character(snpblast.df$Marker2);
snpblast.df$Remove[snpblast.df$Delta1 < snpblast.df$Delta2] <-
  as.character(snpblast.df$Marker1[snpblast.df$Delta1 < snpblast.df$Delta2]);

write.csv(snpblast.df, "", row.names = FALSE);
