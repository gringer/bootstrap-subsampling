#!/usr/bin/Rscript

# data2plot.r -- takes as input a space-separated text file, converts
# it into a plot

## Author: David Eccles (gringer), 2009 <programming@gringer.org>

popNames <- NULL;
popLimits <- NULL;

argLoc <- 1;
while(!is.na(commandArgs(TRUE)[argLoc])){
  if(substr(commandArgs(TRUE)[argLoc],1,1) == "p"){
    inName <- substring(commandArgs(TRUE)[argLoc],2);
    popNames <- append(popNames, sub("_"," ",inName));
    popLimits <- cbind(popLimits,c(as.numeric(commandArgs(TRUE)[argLoc+1]),
                                   as.numeric(commandArgs(TRUE)[argLoc+2])));
    argLoc <- argLoc + 2;
  }
  argLoc <- argLoc + 1;
}


data.df <- read.table(file("stdin"));
pdf("output.pdf", paper = "a4r", height = 8, width = 11);
plot(data.df$V2);
dummy <- dev.off();

dm <- abs(mean(data.df$V2[popLimits[1,1]:popLimits[2,1]]) -
          mean(data.df$V2[popLimits[1,2]:popLimits[2,2]]));

if(dim(popLimits)[2] == 2){
  cat(sprintf("Difference of means: %f (adjusted = %f)\n",
              dm, dm / diff(range(data.df$V2))));
}
