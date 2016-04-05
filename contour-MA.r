#!/usr/bin/Rscript

## contour-MA.r -- converts a delta summary file into a Minor
## Allele Frequency contour graph.

## Author: David Eccles (gringer), 2014 <bioinformatics@gringene.org>

usage <- function(){
  cat("usage: ./contour_MA.r <file>\n");
}

fileName <- "data/stats_HBM30_CEU90_lowMAFfiltered.txt.gz";

if(length(commandArgs(TRUE)) > 0){
    fileName <- commandArgs(TRUE)[1];
}
inFile <- gzfile(fileName);
open(inFile);
cat("Reading in file...", file = stderr());
table.delta <- read.table(inFile, header = TRUE, row.names = 1);
close(inFile);
cat("done!\n", file = stderr());

popNames <- gsub(".","",sub("^p.AA?.","",colnames(table.delta)[grep("^p.AA?",colnames(table.delta))]), fixed = TRUE);
popNames <- unique(popNames);
cat("Populations:",popNames,"\n");
cat("Column names:",colnames(table.delta),"\n");


cat("Calculating prior statistics...", file = stderr());
cat(" mean allele frequencies...", file = stderr());
if(length(grep("^p.AC",colnames(table.delta))) > 0){
  # contains genotype rather than allele counts
  table.delta$meanA <- rowMeans(table.delta[,paste("p.AA.",popNames,".",sep="")] + table.delta[,paste("p.AC.",popNames,".",sep="")] / 2);
  table.delta$meanC <- rowMeans(table.delta[,paste("p.CC.",popNames,".",sep="")] + table.delta[,paste("p.AC.",popNames,".",sep="")] / 2);
  table.delta$meanBoth <- table.delta$meanA + table.delta$meanC;
} else {
  table.delta$meanA <- rowMeans(table.delta[,paste("p.A.",popNames,".",sep="")]);
  table.delta$meanC <- rowMeans(table.delta[,paste("p.C.",popNames,".",sep="")]);
  table.delta$meanBoth <- table.delta$meanA + table.delta$meanC;
}

cat(" minor allele...", file = stderr());
table.delta$mAllele <- c("A","C")[2 - (table.delta$meanA < table.delta$meanC)];

cat(" minor allele frequencies...", file = stderr());
table.delta[,paste("MAF.",popNames,".",sep="")] <- NA;
if(length(grep("^p.AC",colnames(table.delta))) > 0){
  table.delta[table.delta$mAllele == "A",paste("MAF.",popNames,".",sep="")] <-
      (table.delta[table.delta$mAllele == "A",
                   paste("p.AA.",popNames,".",sep="")] +
     table.delta[table.delta$mAllele == "A",
                 paste("p.AC.",popNames,".",sep="")]/2);
  table.delta[table.delta$mAllele == "C",
              paste("MAF.",popNames,".",sep="")] <-
    (table.delta[table.delta$mAllele == "C",
                 paste("p.CC.",popNames,".",sep="")] +
     table.delta[table.delta$mAllele == "C",
                 paste("p.AC.",popNames,".",sep="")]/2);
} else {
    table.delta[,paste("MAF.",popNames,".",sep="")] <- NA;
    table.delta[table.delta$mAllele == "A",
                paste("MAF.",popNames,".",sep="")] <-
                    table.delta[table.delta$mAllele == "A",
                                paste("p.A.",popNames,".",sep="")];
    table.delta[table.delta$mAllele == "C",
                paste("MAF.",popNames,".",sep="")] <-
                    table.delta[table.delta$mAllele == "C",
                                paste("p.C.",popNames,".",sep="")];
}
flipRows <- table.delta[,paste("MAF.",popNames[1],".",sep="")] > 0.5;
table.delta[flipRows,paste("MAF.",popNames,".",sep="")] <-
    1 - table.delta[flipRows,paste("MAF.",popNames,".",sep="")];
cat("done!\n", file = stderr());

head(table.delta);

## table.delta <- table.delta[,c(paste("MAF.",popNames,".",sep=""),"missing","mAllele")];

cat("Producing graph...", file = stderr());

pdf(paste("contour_MAF_",paste(popNames, collapse = "_"),".pdf",sep=""));
makeContour <- function(){
    par(mar = c(5,5,1,1));
    smoothScatter(table.delta[,paste("MAF.",popNames[1],".",sep="")],
                  table.delta[,paste("MAF.",popNames[2],".",sep="")],
                  colramp = colorRampPalette(c("red","orange","yellow",
                      "green","cyan","blue","violet")),
                  transformation = function(x){log(pmax(1,x))/log(10)},
                  nrpoints = 0, xlab = paste("MAF(",popNames[1],")",sep=""),
                  ylab = paste("AF(",popNames[2],")",sep=""),
                  xlim = range(table.delta[,paste("MAF.",popNames[1],".",sep="")]),
                  ylim = range(table.delta[,paste("MAF.",popNames[2],".",sep="")]),
                  xaxs = "i", yaxs = "i",
                  las = 1, cex.lab = 1.5);
    abline(0,1,lty = "dotted", lwd = 3);
}
makeContour();
dummy <- dev.off();
cat("done!\n", file = stderr());


#pdf(paste("contour_MAF_",paste(popNames, collapse = "_"),".pdf",sep=""))
#for(x in 1:(length(popNames) - 1)){
#  for(y in (x+1):length(popNames)){
#    pop1Str <- paste("MAF.",popNames[x],".",sep="");
#    pop2Str <- paste("MAF.",popNames[y],".",sep="");
#    hb <- hexbin(table.delta[,pop1Str],table.delta[,pop2Str], xbins = 60);
#    hb@xlab <- popNames[x];
#    hb@ylab <- popNames[y];
#    plot(hb, trans = function(x){log(pmax(1,x))/log(10)}, inv = function(x){10^x},
#         colramp = colorRampPalette(c("blue","green","yellow","orange","red")),
#         mincnt = 0);
#  }
#}
#dummy <- dev.off();
#cat("done!\n", file = stderr());
