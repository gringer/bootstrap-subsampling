#!/usr/bin/Rscript

## Author: David Eccles (gringer), 2008 <programming@gringer.org>

## stats2contour.r -- converts a delta summary file into a Minor Allele Frequency contour graph

source('/itsshared/phd/common.r')
library(cairoDevice)

usage <- function(){
  cat("usage: ./stats2contour.r <file>\n");
}

table.delta <- read.table(commandArgs(TRUE)[1],
                          header = TRUE, row.names = 1)

A.names <- names(table.delta)[grep("p.A",names(table.delta))];

## determine minor allele frequencies, and comparative SNP in other population
table.delta$MAF.p1 <- table.delta[A.names[1]];
names(table.delta$MAF.p1) <- "MAF.p1";
table.delta$comp.p2 <- table.delta[A.names[2]];
names(table.delta$comp.p2) <- "comp.p2";
##make sure this is minor allele frequency, compare with same allele
##in second population
major.locs <- which(table.delta$MAF.p1 > 0.5);
table.delta$MAF.p1[major.locs,] =
  1 - (table.delta$MAF.p1[major.locs,])
table.delta$comp.p2[major.locs,] =
  1 - (table.delta$comp.p2[major.locs,])
# population names
pn <- sub("\\.","",sub("p.A.","",A.names));

##pdf(paste("histdelta_",pn[1],"_",pn[2],".pdf", paper = "a4");
##hist(table.delta$delta, main = paste(
##     "Histogram of Delta values for each SNP\n(",
##                                  pn[1]," vs ",pn[2],")",sep=""), col =
##     "red", xlab = "Delta", xlim = c(0,1));

## dev.off()


for(dx in seq(0,1,by=0.1)){
  cat("Number of SNPs with delta >= ",dx,": ",length(which(abs(table.delta$MAF.p1 - table.delta$comp.p2) >= dx)),"\n", sep="");
}

library(MASS); # for kde2d
dd <- kde2d(unlist(table.delta$MAF.p1), unlist(table.delta$comp.p2),
            lims = c(0,0.5,0,1));

pdf(paste("contour_MAF_",pn[1],"_",pn[2],".pdf",sep=""))
#Cairo_svg(paste("contour_MAF_",pn[1],"_",pn[2],".svg",sep=""))
par(cex.axis = 1.5, cex.lab = 1.5, mar = c(4.5,4.7,2,1));
filled.contour(dd, nlevels = 40,
               xlab = paste("MAF(",pn[1],")",sep=""),
               ylab = paste("f(",pn[2],")",sep=""),
##               main = paste("Contour graph of ",pn[1]," Minor Allele ",
##               "Frequencies", "\n(",pn[2]," vs ",pn[1],")",sep=""),
               levels = c(exp(pretty(log(c(1,max(dd$z)+1)),39))-1),
               col = rainbow(40), cex = 1.5, cex.main = 1.5, cex.axis=1.5,
               ylim = c(0,1), xlim = c(0,0.5))
dummy <- dev.off();
