#!/usr/bin/Rscript

## pwsummary2pdf.r -- Converts a pairwise summary (such as that generated from
## pwsummary.txt in 2_maximal.pl) into a heatmap PDF file

## Author: David Eccles (gringer), 2008 <programming@gringer.org>

usage <- function(){
  cat("usage: ./pwsummary2pdf.r <top/right file> [bottom/left file]\n");
  cat("[creates output_pwheatmap.(pdf|svg)]\n");
}

argLoc <- 1;
topFile <- FALSE;
bottomFile <- FALSE;
doSVG <- FALSE;
asImage <- TRUE;
scaleBottom <- TRUE; # should the bottom SNPs be scaled?
preserveGroups <- FALSE; # preserve known groupings when clustering
out.filename <- "output_pwheatmap.pdf";
popLimits <- NULL;
popNames <- NULL;
ignoreIndivs <- NULL;
indLabelSize <- NULL;
showLabels <- TRUE;
startIgnore <- FALSE;
makeDendrogram <- FALSE;

while(!is.na(commandArgs(TRUE)[argLoc])){
  if(file.exists(commandArgs(TRUE)[argLoc])){ # file existence check
    print(paste("found file:",commandArgs(TRUE)[argLoc]));
    if(topFile == FALSE){
      topFile <- commandArgs(TRUE)[argLoc];
    } else {
      if(bottomFile == FALSE){
        bottomFile <- commandArgs(TRUE)[argLoc];
      } else{
        cat("Error: More than two input files specified\n");
        cat("usage: ./pwsummary2pdf.r <top/right file> [bottom/left file]\n");
        quit(save="no", status=2);
      }
    }
  }
  else{ # other command line arguments
    if(commandArgs(TRUE)[argLoc] == "-ignore"){
      startIgnore <- TRUE;
    } else {
      if(startIgnore && length(grep("^[0-9]",commandArgs(TRUE)[argLoc])) != 0){
        ignoreIndivs <-
          append(ignoreIndivs,as.numeric(commandArgs(TRUE)[argLoc]));
        cat("Ignoring individual ",
            as.numeric(commandArgs(TRUE)[argLoc]),"\n",sep="");
      }
      else{
        startIgnore <- FALSE;
      }
    }
    if(commandArgs(TRUE)[argLoc] == "-svg"){
      library(cairoDevice);
      doSVG <- TRUE;
    }
    if(commandArgs(TRUE)[argLoc] == "-image"){
      asImage <- TRUE;
    }
    if(commandArgs(TRUE)[argLoc] == "-nolabels"){
      showLabels <- FALSE;
    }
    if(commandArgs(TRUE)[argLoc] == "-dendrogram"){
      makeDendrogram <- TRUE;
    }
    if(commandArgs(TRUE)[argLoc] == "-heatmap"){
      asImage <- FALSE;
    }
    if(commandArgs(TRUE)[argLoc] == "-noscale"){
      scaleBottom <- FALSE;
    }
    if(commandArgs(TRUE)[argLoc] == "-size"){
      indLabelSize <- as.numeric(commandArgs(TRUE)[argLoc+1]);
      cat("Changing label size to ", as.numeric(commandArgs(TRUE)[argLoc+1]),"\n",sep="");
      argLoc <- argLoc + 1;
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
if(doSVG){
  out.filename <- sub(".pdf", ".svg", out.filename);
}

if(topFile == FALSE){
  cat("Error: No [valid] input files specified\n");
  cat("usage: ./pwsummary2pdf.r <top/right file> [bottom/left file]\n");
  quit(save="no", status=2);
}

if(bottomFile == FALSE){
  bottomFile <- topFile;
}

#if(makeDendrogram && (!is.null(popNames))){
#  cat("Error: Dendrogram is incompatible with defined populations\n");
#  cat("usage: ./pwsummary2pdf.r <top/right file> [bottom/left file]\n");
#  quit(save="no", status=2);
#}

##load library for heatmap clusters and sorts
library(gtools, warn.conflicts = FALSE);
library(gdata, warn.conflicts = FALSE);
library(gplots, warn.conflicts = FALSE);
##load top/right comparisons
line = "dummy value"
pw.file = gzfile(topFile);
open(pw.file);
while(length(line)>0){
  line = unlist(readLines(pw.file, n=1));
  if(length(line)>0){
    line = unlist(line);
    if(length(grep("^Pairwise: +[0-9]", line))>0){
      top.pwvalues = as.numeric(unlist(strsplit(line," +"))[-1])
    }
  }
}
close(pw.file);
## determine number of rows
## derived from triangle(x) = (x + x^2) / 2
top.pw.rows <- (sqrt(1+8*length(top.pwvalues)) - 1) / 2 + 1;
## load bottom/left comparisons
line = "dummy value"
pw.file = file(bottomFile, "r", blocking = FALSE)
while(length(line)>0){
  line = unlist(readLines(pw.file, n=1));
  if(length(line)>0){
    line = unlist(line);
    if(length(grep("^Pairwise: +[0-9]", line))>0){
      bottom.pwvalues = as.numeric(unlist(strsplit(line," +"))[-1])
    }
  }
}
close(pw.file)
## determine number of rows
bottom.pw.rows <- (sqrt(1+8*length(bottom.pwvalues)) - 1) / 2 + 1;
if(top.pw.rows != bottom.pw.rows){
  cat("Error: Number of rows in input files differ\n");
  cat("Please choose data sets that have the same cardinality.\n");
  cat("usage: ./pwsummary2pdf.r <top/right file> [bottom/left file]\n");
  quit(save="no", status=2);
}
## scale lower pairwise values to make top colours look nicer
## (makes pw values fit into same range as ref values)
if(scaleBottom){
  top.min = min(top.pwvalues, na.rm = TRUE);
  top.max = max(top.pwvalues, na.rm = TRUE);
  bottom.min = min(bottom.pwvalues, na.rm = TRUE);
  bottom.max = max(bottom.pwvalues, na.rm = TRUE);
  bottom.scale = (bottom.max - bottom.min);
  top.scale = (top.max - top.min);
  bottom.pwvalues = ((bottom.pwvalues - bottom.min)/bottom.scale) *
    top.scale + top.min;
}
## convert values to matrix
top.matrix <- matrix(NA,top.pw.rows,top.pw.rows);
top.matrix[which(lower.tri(top.matrix))] <- top.pwvalues;
top.matrix <- t(top.matrix);
top.matrix[which(lower.tri(top.matrix))] <- top.pwvalues;
bottom.matrix <- matrix(NA,top.pw.rows,top.pw.rows);
bottom.matrix[which(lower.tri(bottom.matrix))] <- bottom.pwvalues;
bottom.matrix <- t(bottom.matrix);
bottom.matrix[which(lower.tri(bottom.matrix))] <- bottom.pwvalues;
## set up reference matrix for dendrogram
ref.matrix <- top.matrix;
if(!is.null(popLimits)){
  cat("Mean, SD of similarity for named population groups:\n");
  for(x in 1:(dim(popLimits)[2])){
    for(y in x:(dim(popLimits)[2])){
      cat(popNames[x]," vs ",popNames[y],": mean = ",
          mean(ref.matrix[popLimits[1,x]:popLimits[2,x],
                         popLimits[1,y]:popLimits[2,y]], na.rm = TRUE),
          ", sd = ",
          sd(as.vector(ref.matrix[popLimits[1,x]:popLimits[2,x],
                                 popLimits[1,y]:popLimits[2,y]]), na.rm = TRUE),
          "\n",sep="");
    }
  }
}
## remove data for intermediate values for dendrogram creation
cluster.matrix <- ref.matrix;
if(!is.null(popLimits)){
  for(x in 1:(dim(popLimits)[2])){
    for(y in 1:(dim(popLimits)[2])){
      if(x != y){
##        cluster.matrix[popLimits[1,x]:popLimits[2,x],
##                       popLimits[1,y]:popLimits[2,y]] <- 0;
      }
    }
  }
}
ref.matrix <- cluster.matrix;
## reorder individuals based on hierarchical ranking
## if groupings should be preserved, then do things on a per-group basis
## (and to keep things simple, only re-order individuals in specifed groups)
if(!is.null(popLimits) && preserveGroups){
  popOrder <- NULL;
  for(x in 1:(dim(popLimits)[2])){
    ## copy subset of table for this population
    popref.matrix <-
      ref.matrix[popLimits[1,x]:popLimits[2,x],
                popLimits[1,x]:popLimits[2,x]];
    ref.means <- colMeans(popref.matrix, na.rm = TRUE);
    hcc <- hclust(dist(popref.matrix));
    ddc <- reorder(as.dendrogram(hcc), ref.means);
    hcc <- rev(order.dendrogram(ddc)) + popLimits[1,x] - 1;
    popOrder <- append(popOrder,hcc);
  }
  if(length(popOrder) == top.pw.rows){
    hcc <- rev(popOrder);
    if(!is.null(ignoreIndivs)){
      ref.matrix[ignoreIndivs,] <- NA;
      ref.matrix[,ignoreIndivs] <- NA;
    }
    dend.matrix <- ref.matrix[hcc,hcc];
  } else{
    top.pw.rows <- length(popOrder);
    bottom.pw.rows <- length(popOrder);
    cat("Warning: Populations defined, but do not cover all individuals.\n");
    cat("Warning: [The resultant table will be smaller,",
        "and may not display correctly]\n");
    hcc <- rev(popOrder);
    if(!is.null(ignoreIndivs)){
      ref.matrix[ignoreIndivs,] <- NA;
      ref.matrix[,ignoreIndivs] <- NA;
    }
    dend.matrix <- ref.matrix[hcc,hcc];
  }
} else {
  ref.means <- colMeans(ref.matrix, na.rm = TRUE);
  hcc <- hclust(dist(ref.matrix));
  ddc <- reorder(as.dendrogram(hcc), ref.means);
  hcc <- rev(order.dendrogram(ddc));
  if(!is.null(ignoreIndivs)){
    ref.matrix[ignoreIndivs,] <- NA;
    ref.matrix[,ignoreIndivs] <- NA;
  }
  dend.matrix <- ref.matrix[hcc,hcc];
}
## put reordered bottom values into matrix
if(!is.null(ignoreIndivs)){
  bottom.matrix[ignoreIndivs,] <- NA;
  bottom.matrix[,ignoreIndivs] <- NA;
}
pw.matrix <- bottom.matrix[hcc,hcc];
rownames(pw.matrix) <- rev(hcc);
colnames(pw.matrix) <- rev(hcc);
## put reordered top values into matrix
pw.matrix[which(upper.tri(pw.matrix))] =
  dend.matrix[which(upper.tri(dend.matrix))];
##calculate range for data for heatmap z limits
zRange <- c(min(pw.matrix,na.rm=TRUE),max(pw.matrix,na.rm=TRUE));
## make key
if(doSVG){
  Cairo_svg(paste("key_",out.filename,sep=""), width = 4.8, height = 2.7);
} else {
  pdf(paste("key_",out.filename,sep=""), width = 4.8, height = 2.7);
}
oldPar <- par(no.readonly = TRUE);
par(mar = c(2,1,2.5,1), mgp = c(2,0.75,0), bg = "white");
image(x=seq(zRange[1],zRange[2],length.out=100),y=0.5,
      z = cbind(seq(zRange[1],zRange[2],length.out=100)),
      col = topo.colors(100), main = "Mean IBS Value",
      yaxt = "n", ylab = "", xlab = "", cex.axis = 1.5, cex.main = 2);
box("figure",lwd = 2);
dummy <- dev.off();
par(oldPar);
##initiate PDF generation
fileWidth <- 8;
if(makeDendrogram && asImage){
  fileWidth <- 11;
}
fileHeight <- 8;
if(doSVG){
    Cairo_svg(out.filename, width = fileWidth, height = fileHeight);
} else {
  pdf(out.filename, width = fileWidth, height = fileHeight);
}
##make heatmap
if(is.null(indLabelSize)){
  ## a rough guess on the best size
  ## based on 0.4 being good for 120 individuals]
  indLabelSize <- ( 50 / top.pw.rows );
}
if(asImage){
  if(makeDendrogram){
    layout(rbind(c(1,2)), widths = c(3/11,8/11));
    par(mar=c(2,0,0.5,0)); # bottom, left, top, right
    plot(rev(ddc), horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none");
    if(!is.null(popLimits)){
      for(x in 1:(dim(popLimits)[2])){
        popNames[x] <- gsub("\\\\n","\n",popNames[x]);
          text(0.45,top.pw.rows - (popLimits[1,x] + popLimits[2,x]) / 2 + 1,
               popNames[x], cex = 2);
      }
    }
  }
  par(mar=c(2,0.5,0.5,2)); # bottom, left, top, right
  image(1:top.pw.rows, 1:top.pw.rows, pw.matrix[ncol(pw.matrix):1,], col =
        topo.colors(100), axes=FALSE, xlab = "", ylab = "", zlim = zRange);
  axis(1, at = top.pw.rows:1, labels=hcc, las = 2, line = -0.5, tick =
       0, cex.axis = indLabelSize);
  axis(4, at = 1:top.pw.rows, labels=hcc, las = 2, line = -0.5, tick =
       0, cex.axis = indLabelSize);
  if(!makeDendrogram){
    if(!is.null(popLimits)){
      for(x in 1:(dim(popLimits)[2])){
        for(y in x:(dim(popLimits)[2])){
          rect(popLimits[1,y] - 0.5, top.pw.rows - popLimits[2,x] + 0.5,
               popLimits[2,y] + 0.5, top.pw.rows - popLimits[1,x] + 1.5, lwd = 2);
          if(showLabels){
            if(x == y){
              text(popLimits[1,y] + (popLimits[2,y] - popLimits[1,y]) * 0.75,
                   top.pw.rows -
                   (popLimits[1,x] + (popLimits[2,x] - popLimits[1,x]) * 0.25),
                   popNames[x], cex = 2, adj = c(0.5,0.5),
                   srt = -45, offset = 0);
            } else {
              text(popLimits[1,y] + (popLimits[2,y] - popLimits[1,y]) * 0.5,
                   top.pw.rows -
                   (popLimits[1,x] + (popLimits[2,x] - popLimits[1,x]) * 0.5),
                   paste(popNames[x],"vs",popNames[y]), cex = 2,
                   adj = c(0.5,0.5), offset = 0);
            }
          }
        }
      }
    }
  }
} else {
  if(makeDendrogram){
    heatmap.2(as.matrix(pw.matrix), symm = TRUE, trace = 'none', Rowv =
              TRUE, key = TRUE, dendrogram = 'row', density.info =
              'none', col = topo.colors(100), main =
              paste("Pairwise similarity heatmap"), xlab = "Row ID",
              zlim = zRange,
              keysize = 1.2, cexRow = indLabelSize, cexCol = indLabelSize);
  } else {
    heatmap.2(as.matrix(pw.matrix), symm = TRUE, trace = 'none', Rowv =
              FALSE, key = TRUE, dendrogram = 'none', density.info =
              'none', col = topo.colors(100), main =
              paste("Pairwise similarity heatmap"), xlab = "Row ID",
              zlim = zRange,
              keysize = 1.2, cexRow = indLabelSize, cexCol = indLabelSize);
  }
}
dummy <- dev.off();
