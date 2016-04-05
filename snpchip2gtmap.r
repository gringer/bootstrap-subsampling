#!/usr/bin/Rscript

## snpchip2gtmap.r -- Creates a chromosome location diagram for a sequence of SNPs, with
## a genotype summary below

## Author: David Eccles (gringer), 2009 <programming@gringer.org>

usage <- function(){
  cat("usage: ./snpchip2gtmap.r <genotype file> <SNP location file>\n");
}

argLoc <- 1;
gtFile <- FALSE;
locFile <- FALSE;
useLocation <- TRUE;
featureNames <- NULL;
featureLocs <- NULL;
popLimits <- NULL;
popNames <- NULL;

while(!is.na(commandArgs(TRUE)[argLoc])){
  if(file.exists(commandArgs(TRUE)[argLoc])){ # file existence check
    if(gtFile == FALSE){
      gtFile <- commandArgs(TRUE)[argLoc];
    } else {
      if(locFile == FALSE){
        locFile <- commandArgs(TRUE)[argLoc];
      } else{
        cat("Error: More than two input files specified\n");
        cat("usage: ./snpchip2gtmap.r <genotype file> <SNP location file>\n");
        quit(save="no", status=2);
      }
    }
  }
  else{ # other command line arguments
    if(commandArgs(TRUE)[argLoc] == "-feature"){
      featureNames <- append(featureNames,
                             commandArgs(TRUE)[argLoc+1]);
      featureLocs <- append(featureLocs,
                            as.numeric(commandArgs(TRUE)[argLoc+2]));
      argLoc <- argLoc + 2;
    }
    if(commandArgs(TRUE)[argLoc] == "-nolocation"){
      useLocation <- FALSE;
    }
    if(substr(commandArgs(TRUE)[argLoc],1,1) == "p"){
      popNames <- append(popNames, substring(commandArgs(TRUE)[argLoc],2));
      popLimits <- cbind(popLimits,c(as.numeric(commandArgs(TRUE)[argLoc+1]),
                                     as.numeric(commandArgs(TRUE)[argLoc+2])));
      argLoc <- argLoc + 2;
    }
  }
  argLoc <- argLoc + 1;
}

if(locFile == FALSE){
  cat("Error: Fewer than two [valid] input files specified\n");
  cat("usage: ./snpchip2gtmap.r <genotype file> <SNP location file>\n");
  quit(save="no", status=2);
}

gt.df <-
  read.table(gtFile, colClasses = "character");
numIndivs <- dim(gt.df)[2] - 1;
colnames(gt.df) <- c("Marker",1:numIndivs);
rownames(gt.df) <- gt.df[,"Marker"];

loc.df <- read.table(locFile, nrows = 1);
if(length(grep("^rs",loc.df[1,1])) == 0){
  loc.df <- read.table(locFile, header = TRUE);
} else {
  loc.df <- read.table(locFile);
}
colnames(loc.df) <- c("Marker","Mutation","Chromosome","Location");

gt.df <- merge(gt.df,loc.df);
rownames(gt.df) <- gt.df[,"Marker"];
gt.df <- gt.df[,-1]; # remove marker, now that it's added to rownames

for(x in 1:numIndivs){
  ## substitute complementary bases
  gt.df[,x] <- gsub("T","A",gt.df[,x]);
  gt.df[,x] <- gsub("G","C",gt.df[,x]);
  if(length(which(gt.df[,x] == "CA"))>0){
    ## not including this check complains gives the following error:
    ## In max(i) : no non-missing arguments to max; returning -Inf
    gt.df[which(gt.df[,x] == "CA"),x] <- "AC";
  }
  if(length(grep("[^AC]",gt.df[,x]))>0){
    ## not including this check complains gives the following error:
    ## In max(i) : no non-missing arguments to max; returning -Inf
    gt.df[grep("[^AC]",gt.df[,x]),x] <- NA;
  }
}
## replace text with "C-purity" value
gt.df[,1:numIndivs] <- (function(x){
  ifelse(is.na(x),NA,ifelse(x=="AA",0,ifelse(x=="CC",1,0.5)));
})(gt.df[,1:numIndivs])

##print(gt.df[1:5,c(1:5,(dim(gt.df)[2]-5):dim(gt.df)[2])]);

if(useLocation){
  gt.df[,"Location"] <- gt.df[,"Location"] / 1000000;
  gt.df <- gt.df[order(gt.df[,"Location"]),];
} else {
  gt.df[,"Location"] <- 1:length(gt.df[,"Location"]);
}

gt.numeric <- gt.df[,1:numIndivs];
gt.numeric[is.na(gt.numeric)] <- -2;
gt.df[,1:numIndivs] <- gt.numeric;

ref.means <- colMeans(gt.df[,1:numIndivs], na.rm = TRUE);
hcc <- hclust(dist(t(gt.df[,1:numIndivs])));
ddc <- reorder(as.dendrogram(hcc), ref.means);
hcc <- order.dendrogram(ddc);

gt.numeric <- gt.df[,hcc];
gt.numeric[gt.numeric == -2] <- NA;
colnames(gt.numeric) <- hcc;
gt.numeric[,numIndivs+1] <- NA;
colnames(gt.numeric)[numIndivs+1] <- "";
if(!is.null(popNames)){
  for(x in 1:length(popNames)){
    gt.numeric[,numIndivs+x+1] <-
      rowMeans(gt.numeric[,popLimits[1,x]:popLimits[2,x]]);
    colnames(gt.numeric)[numIndivs+x+1] <- popNames[x];
  }
}

pdf("output_gtmap.pdf",width=11,height=8,paper="a4r");
library(gtools);
library(gdata);
library(gplots, warn.conflicts = FALSE);
if(useLocation){
  xLabel <- "Location (Mbp)";
} else {
  xLabel <- "SNP Number";
}
image(gt.df[,"Location"], 1:(numIndivs + length(popNames)+1),
      as.matrix(gt.numeric),
      col = rgb(
        c(seq(1,1,length.out=5),seq(1,0,length.out=5)),
        c(seq(0,1,length.out=5),seq(1,1,length.out=5)),0),
       axes=FALSE, xlab = xLabel, ylab = "",cex.lab = 2);

axis(1, las = 2, line = 0.5, cex.axis = 1);

axis(4, at = 1:(numIndivs + length(popNames)+1), labels=colnames(gt.numeric),
     las = 2, line = -0.5, tick =
     0, cex.axis =  min(0.5, 30 / (numIndivs + length(popNames)+1)) );

if(!is.null(featureNames)){
  for(x in 1:length(featureNames)){
    segments(x0=featureLocs[x]/1000000, y0=0,
             x1=featureLocs[x]/1000000, y1=numIndivs + length(popNames)+2);
    rect(xleft=featureLocs[x]/1000000 - 0.125,
         ybottom=(numIndivs + length(popNames)+2)/2 - 10,
         xright=featureLocs[x]/1000000 + 0.125,
         ytop=(numIndivs + length(popNames)+2)/2 + 10,
         col="white", border = "black");
    text(srt=90, x=featureLocs[x]/1000000,
         y=(numIndivs + length(popNames)+2)/2, featureNames[x]);
  }
}
dummy <- dev.off();
