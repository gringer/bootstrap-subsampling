#!/usr/bin/Rscript

## bs2meancalc.r -- calculates mean / SD from bootstrap summary results

## Author: David Eccles (gringer), 2009 <programming@gringer.org>


usage <- function(){
  cat("usage: ./bs2meanvar.r <file> [options]\n");
  cat("\nOther Options:\n");
  cat("-help            : Only display this help message\n");
## ignore string... not yet implemented
##  cat("-ignore <string> : consider <string> to indicate missing values");
  cat("\n");
}

means.df <- NULL;
infile.name <- FALSE;
ignore.string <- FALSE;

argLoc <- 1;

# print(commandArgs());

argLoc <- grep("--args",commandArgs()) + 1; # hack to get around R v2.4
                                            # issue stopping
                                            # commandArgs(TRUE) from
                                            # working

while(!is.na(commandArgs()[argLoc])){
  if(file.exists(commandArgs()[argLoc])){ # file existence check
    if(infile.name == FALSE){
      infile.name <- commandArgs()[argLoc];
    } else{
      cat("Error: More than one input file specified\n");
      usage();
      quit(save = "no", status=1);
    }
  } else {
    if(commandArgs()[argLoc] == "-help"){
      usage();
      quit(save = "no", status=0);
    }
    if(commandArgs()[argLoc] == "-ignore"){
      ignore.string = commandArgs()[argLoc+1];
      argLoc <- argLoc + 1;
    }
  }
  argLoc <- argLoc + 1;
}

if(infile.name == FALSE){
  cat("Error: no input file specified\n");
  usage();
  quit(save = "no", status=2);
}

con <- gzfile(infile.name);
open(con);

## determine number of bootstraps
cat("Checking first marker to determine number of bootstraps...",file = stderr());
numBootstraps <- FALSE;
marker.df <- NULL;
line.df <- read.csv(con, header = TRUE, nrows = 1, colClasses =
                    c("character", "integer", "numeric"), col.names
                    = c("marker","bs","value"));
first.name <- line.df$marker[1];
marker.name <- line.df$marker[1];
while(marker.name == first.name){
  marker.df <- rbind(marker.df, line.df);
  numBootstraps <- numBootstraps + 1;
  line.df <- read.csv(con, header = FALSE, nrows = 1, colClasses =
                      c("character", "integer", "numeric"), col.names
                      = c("marker","bs","value"));
  marker.name <- line.df$marker[1];
}

## analyse first marker
marker.mean <- mean(marker.df$value);
marker.sd <- sd(marker.df$value);
write.table(data.frame(name = marker.name, mean = marker.mean, sd = marker.sd), row.names = FALSE, col.names = TRUE, sep = ",");

cat("Done, found", numBootstraps, "bootstraps\n",file = stderr());
firstLine <- TRUE;
lineCount <- 1;

cat("Calculating mean/SD (one '.' per 1000 markers)",file = stderr());
while(dim(marker.df)[1] > 0){
  if(lineCount %% 1000 == 0){
    cat(".",file = stderr());
  }
  if(firstLine){ # picks up remainder from bootstrap number check
    marker.df <-
      read.csv(con, header = FALSE, nrows = numBootstraps - 1,
               colClasses = c("character", "integer", "numeric"),
               col.names = c("marker","bs","value"));
    marker.df <- rbind(line.df, marker.df);
    firstLine <- FALSE;
  } else {
    marker.df <-
      read.csv(con, header = FALSE, nrows = numBootstraps,
               colClasses = c("character", "integer", "numeric"),
               col.names = c("marker","bs","value"));
  }
#  cat("Read",dim(marker.df)[1],"lines\n");
  if(dim(marker.df)[1] == numBootstraps){
    marker.name <- marker.df$marker[1];
    marker.mean <- mean(marker.df$value);
    marker.sd <- sd(marker.df$value);
    if(marker.df[numBootstraps,"marker"] != marker.df[1,"marker"]){
      stop("Error: first marker is not the same as the last marker");
    }
    write.table(data.frame(name = marker.name, mean = marker.mean, sd = marker.sd), row.names = FALSE, col.names = FALSE, sep = ",");
  }
  lineCount <- lineCount + 1;
}
close(con);

cat("done!\n",file = stderr());
