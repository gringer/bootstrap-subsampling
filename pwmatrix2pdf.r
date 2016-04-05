#!/usr/bin/Rscript

## pwmatrix2pdf.r -- Converts a pairwise matrix (such as
## that generated from gt2pw.r) into a heatmap PDF file.
##
## usage: ./pwmatrix2pdf.r <top/right file> [bottom/left file]
## [creates output_pwheatmap.(pdf|svg)]

## Author: David Eccles (gringer), 2008 <programming@gringer.org>

usage <- function(){
  cat("usage: ./pwmatrix2pdf.r <top/right file> [bottom/left file]\n");
  cat("\nOther Options:\n");
  cat("-help               : Only display this help message\n");
  cat("-ignore (value)*    : Ignore particular individuals\n");
  cat("-svg                : Output to an SVG file\n");
  cat("-bitmap             : Output to a bitmap (XPM) file\n");
  cat("-image              : Create a more customisable image\n");
  cat("-invert <value>     : Create a similarity matrix (value is no. of alleles)\n");
  cat("-nolabels           : Don't output population labels\n");
  cat("-norects            : Don't draw squares around populations\n");
  cat("-dendrogram         : Sort, and draw a dendrogram / tree\n");
  cat("-heatmap            : Use the built-in heatmap function\n");
  cat("-noscale            : Don't scale bottom half to range of top half\n");
  cat("-outliers <value>   : Ignore values outside specified probability\n");
  cat("-size <value>       : Change label size\n");
  cat("p<label> <s> <e>    : Population group, starting at s, ending at e\n");
  cat("\n");
}
argLoc <- 1;
topFile <- FALSE;
bottomFile <- FALSE;
doSVG <- FALSE;
asImage <- TRUE;
asBitmap <- FALSE;
scaleBottom <- TRUE; # should the bottom SNPs be scaled?
preserveGroups <- FALSE; # preserve known groupings when clustering
out.filename <- "output_pwheatmap.pdf";
out.bmname <- FALSE;
popLimits <- NULL;
popNames <- NULL;
ignoreIndivs <- NULL;
indLabelSize <- NULL;
textSize <- NULL;
invertValues <- FALSE;
showLabels <- TRUE;
outlierProb <- FALSE;
drawRects <- TRUE;
startIgnore <- FALSE;
makeDendrogram <- FALSE;

maxColour <- NULL; # set to red if individuals are identical

while(!is.na(commandArgs(TRUE)[argLoc])){
  if(file.exists(commandArgs(TRUE)[argLoc])){ # file existence check
    cat(sprintf("found file: %s\n",commandArgs(TRUE)[argLoc]), file = stderr());
    if(topFile == FALSE){
      topFile <- commandArgs(TRUE)[argLoc];
    } else {
      if(bottomFile == FALSE){
        bottomFile <- commandArgs(TRUE)[argLoc];
      } else{
        cat("Error: More than two input files specified\n", file = stderr());
        usage();
        quit(save="no", status = 1);
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
            as.numeric(commandArgs(TRUE)[argLoc]),"\n",sep="", file = stderr());
      }
      else{
        startIgnore <- FALSE;
      }
    }
    if(commandArgs(TRUE)[argLoc] == "-help"){
      usage();
      quit(save="no", status = 0);
    }
    if(commandArgs(TRUE)[argLoc] == "-svg"){
      doSVG <- TRUE;
    }
    if(commandArgs(TRUE)[argLoc] == "-image"){
      asImage <- TRUE;
    }
    if(commandArgs(TRUE)[argLoc] == "-bitmap"){
      asBitmap <- TRUE;
    }
    if(commandArgs(TRUE)[argLoc] == "-nolabels"){
      showLabels <- FALSE;
    }
    if(commandArgs(TRUE)[argLoc] == "-norects"){
      drawRects <- FALSE;
    }
    if(commandArgs(TRUE)[argLoc] == "-invert"){
      invertValues <- as.numeric(commandArgs(TRUE)[argLoc+1]);
      cat("Inverting using max value of ",
          as.numeric(commandArgs(TRUE)[argLoc+1]),
          " alleles\n", sep="", file = stderr());
      argLoc <- argLoc + 1;
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
    if(commandArgs(TRUE)[argLoc] == "-outliers"){
      outlierProb <- as.numeric(commandArgs(TRUE)[argLoc+1]);
      argLoc <- argLoc + 1;
    }
    if(commandArgs(TRUE)[argLoc] == "-size"){
      textSize <- as.numeric(commandArgs(TRUE)[argLoc+1]);
      cat("Changing text size to ", as.numeric(commandArgs(TRUE)[argLoc+1]),
          "\n", sep="", file = stderr());
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

if(asBitmap){
  out.bmname <- sub("\\.(pdf|svg)",".xpm",out.filename);
}


if(topFile == FALSE){
  cat("Error: No [valid] input files specified\n", file = stderr());
  usage();
  quit(save="no", status = 2);
}

##load top/right comparisons
cat("Reading in file(s)...", file = stderr());
top.matrix <- as.matrix(read.csv(gzfile(topFile), row.names = 1));
if(bottomFile != FALSE){
  ##load bottom/left comparisons
  bottom.matrix <- as.matrix(read.csv(gzfile(bottomFile), row.names = 1));
} else {
  bottom.matrix <- top.matrix;
}
cat("done!\n", file = stderr());

## remove individuals to ignore
if(!is.null(ignoreIndivs)){
  top.matrix[ignoreIndivs,] <- NA;
  top.matrix[,ignoreIndivs] <- NA;
  bottom.matrix[ignoreIndivs,] <- NA;
  bottom.matrix[,ignoreIndivs] <- NA;
}

top.pw.rows <- dim(top.matrix)[1];
bottom.pw.rows <- dim(bottom.matrix)[1];
if(top.pw.rows != bottom.pw.rows){
  cat("Error: Number of rows in input files differ\n");
  cat("Please choose data sets that have the same cardinality.\n");
  cat("usage: ./pwsummary2pdf.r <top/right file> [bottom/left file]\n");
  quit(save="no", status = 3);
}

# If there are no defined populations, treat them all as a population
if(is.null(popLimits)){
      popNames <- "All";
      popLimits <- cbind(NULL,c(1,top.pw.rows));
}

## determine clusters, assuming top matrix is the most informative one
popClusters <- list();
popOrderList <- list();
curPos <- 1;
for(pop in 1:(dim(popLimits)[2])){
  popOrderList <- append(popOrderList, list(popLimits[1,pop]:popLimits[2,pop]));
  popLimits[,pop] <- popLimits[,pop] - min(popLimits[,pop]) + curPos;
  curPos <- curPos + abs(diff(popLimits[,pop])) + 1;
}
popOrder <- unlist(popOrderList);
top.matrix <- top.matrix[popOrder,popOrder];
bottom.matrix <- bottom.matrix[popOrder,popOrder];
top.pw.rows <- dim(top.matrix)[1];
popSizes <- top.pw.rows;
popOrderLocs <- popOrder;
if(makeDendrogram){
  cat("Sorting...", file = stderr());
  popOrder <- NULL;
  popSizes <- NULL;
  for(pop in 1:(dim(popLimits)[2])){
    tmpClust <- hclust(dist(top.matrix[popLimits[1,pop]:popLimits[2,pop],
                                       popLimits[1,pop]:popLimits[2,pop]]));
    popClusters <- append(popClusters, list(tmpClust));
    popSizes <- c(popSizes, length(tmpClust$order));
    popOrder <- c(popOrder, tmpClust$order + min(popLimits[,pop]) - 1);
    cat(".", file = stderr());
  }
  top.matrix <- top.matrix[popOrder,popOrder];
  bottom.matrix <- bottom.matrix[popOrder,popOrder];
  cat("done!\n", file = stderr());
}
if(length(popOrder) != top.pw.rows){
  cat("Warning: Population sizes do not match original pairwise matrix\n",
      file = stderr());
}

diag(top.matrix) <- NA;
diag(bottom.matrix) <- NA;

if(min(top.matrix, na.rm = TRUE) == 0){
  maxColour <- "red";
}

## scale lower pairwise values to make top colours look nicer
## (makes pw values fit into same range as ref values)
if(scaleBottom){
  top.min <- min(top.matrix, na.rm = TRUE);
  top.max <- max(top.matrix, na.rm = TRUE);
  top.mean <- mean(top.matrix, na.rm = TRUE);
  top.sd <- sd(as.vector(top.matrix), na.rm = TRUE);
  top.minOutlier <- qnorm(outlierProb / 2, mean = top.mean, sd = top.sd);
  top.maxOutlier <- qnorm(1 - (outlierProb / 2), mean = top.mean, sd = top.sd);
  if(outlierProb != FALSE){
    cat("Removing outliers...", file = stderr());
    cat(length(which(top.matrix < top.minOutlier))," too low...",sep="");
    top.matrix[which(top.matrix < top.minOutlier)] <- top.minOutlier;
    cat(length(which(top.matrix > top.maxOutlier))," too high...",sep="");
    top.matrix[which(top.matrix > top.maxOutlier)] <- top.maxOutlier;
    if(bottomFile == FALSE){
      bottom.matrix <- top.matrix;
    }
    top.min = min(top.matrix, na.rm = TRUE);
    top.max = max(top.matrix, na.rm = TRUE);
    top.mean = mean(top.matrix, na.rm = TRUE);
    top.sd = sd(top.matrix, na.rm = TRUE);
    cat("done\n", file = stderr());
  }
  if(bottomFile != FALSE){
    bottom.min = min(bottom.matrix, na.rm = TRUE);
    bottom.max = max(bottom.matrix, na.rm = TRUE);
    bottom.mean = mean(bottom.matrix, na.rm = TRUE);
    bottom.sd = sd(as.vector(bottom.matrix), na.rm = TRUE);
    bottom.scale = (bottom.max - bottom.min);
    top.scale = (top.max - top.min);
    ### straight transformation
    ## bottom.matrix = ((bottom.matrix - bottom.min)/bottom.scale) *
    ##  top.scale + top.min;
    ### normal transformation
    bottom.matrix = (bottom.matrix - bottom.mean) *
      top.sd /bottom.sd + top.mean;
  }
}

## put bottom matrix into its place
pw.matrix <- top.matrix;
pw.matrix[lower.tri(pw.matrix)] <- bottom.matrix[lower.tri(bottom.matrix)];
## blank out diagonals
diag(pw.matrix) <- NA;
if(invertValues != FALSE){
  ## invert values if necessary / requested
  pw.matrix = (invertValues - pw.matrix) / invertValues;
}

##calculate range for data for heatmap z limits
zRange <- c(min(pw.matrix, na.rm = TRUE),max(pw.matrix,na.rm = TRUE));


## make key
if(doSVG){
  ## load graphics library
  library(cairoDevice);
  Cairo_svg(file = paste("key_", out.filename,sep = ""),
            width = 4.8, height = 2.7);
} else {
  pdf(file = paste("key_", out.filename,sep = ""),
      width = 4.8, height = 2.7);
}
oldPar <- par(no.readonly = TRUE);
par(mar = c(2,1,2.5,1), mgp = c(2,0.75,0), bg = "white");
image(x=seq(zRange[1],zRange[2],length.out=100),y=0.5,
      z = cbind(seq(zRange[1],zRange[2],length.out=100)),
      col = c(topo.colors(100),maxColour), main = "Mean IBS Value",
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
  indLabelSize <- ( 80 / top.pw.rows );
}
if(is.null(textSize)){
  ## a rough guess on the best size
  ## based on 0.4 being good for 120 individuals]
  textSize <- 2;
}
if(asBitmap){
  ## create an XPM bitmap image
  ## The image is broken into 91 colours (more would require two
  ## characters per colour), converted to a character array, then
  ## output to the file
  ## Note: this is the raw matrix data, no dendrograms are attached to this
  par(mar=c(0,0,0,0)); # bottom, left, top, right
  cat("Determining ranges / counts...", file = stderr());
  imageVals <- cut(pw.matrix,91, dig.lab = 6); # 91 "usable" ascii characters
  cutCounts <- table(imageVals);
  cutSum <- cumsum(cutCounts);
  cat("done!\n", file = stderr());
  levels(imageVals) <- intToUtf8(as.integer((35:126)[-(92-34)]),
                                 multiple = TRUE); # exclude '"'
  image.mat <- as.character(imageVals);
  image.mat[is.na(image.mat)] <- " "; # space for NA values
  dim(image.mat) <- c(top.pw.rows,top.pw.rows);
  outFile <- file(out.bmname, open = "w");
  cat("/* XPM */\n", file = outFile);
  cat("static char * ",sub("\\.","_",out.bmname),"[] = {\n", sep = "",
      file = outFile);
  cat("\"",paste(top.pw.rows, top.pw.rows,92,1, collapse = " "),"\",\n",
      sep = "", file = outFile);
  if(invertValues != FALSE){
    cat("/* Values indicate similarity in range 0..1*/", file = outFile);
    cat(sprintf("/* using max number of shared alleles = %d */\n",
                invertValues), file = outFile);
  }
  cat(sprintf("\"%-6s %s %7s\", /* NA: %d */\n"," ","c","#FFFFFF",
              length(which(is.na(pw.matrix)))),
      file = outFile);
  cat(sprintf("\"%-6s %s %7s\", /* %s: %d */\n",levels(imageVals),"c",
              sub("FF$","",topo.colors(91)), names(cutCounts), cutCounts),
      sep = "", file = outFile);
  cat("\"", file = outFile);
  image.text <- apply(image.mat, 1, paste, collapse = "");
  cat(paste(image.text, collapse = "\",\n\""), sep = "", file = outFile);
  cat("\"};\n", file = outFile);
  close(outFile);
}
if(asImage || asBitmap){
  if(makeDendrogram){
    layout(cbind(1,2:(dim(popLimits)[2] + 1),dim(popLimits)[2] + 2),
           widths = c(0.5/11,1.5/11,8/11), heights = popSizes);
    library(stats);
    par(mar=c(0,0,0,0)); # bottom, left, top, right
    plot.new();
    for(x in 1:(dim(popLimits)[2])){
      popNames[x] <- gsub("\\\\n","\n",popNames[x]);
      text(labels = popNames[x], x = 0.5,
           y = 1 - mean(popLimits[,x]) / top.pw.rows,
           adj = c(0.5,0.5), srt = 90, cex = textSize);
    }
    for(x in 1:(dim(popLimits)[2])){
      par(mar=c(
            if(x == dim(popLimits)[2]) 2 else 0, # bottom
            0,                                   # left
            if(x == 1) 0.5 else 0,               # top
            0));                                 # right
      plot(as.dendrogram(popClusters[[x]]), horiz = TRUE, axes = FALSE,
           yaxs = "i", leaflab = "none");
    }
  }
  par(mar=c(2,0.5,0.5,2)); # bottom, left, top, right
  if(!asBitmap){
    image(1:top.pw.rows, 1:top.pw.rows, t(pw.matrix[ncol(pw.matrix):1,]), col =
          c(topo.colors(100),maxColour), axes=FALSE, xlab = "", ylab = "", zlim = zRange);
  } else {
    image(top.pw.rows/2, top.pw.rows/2, matrix(1,1,1), col =
          "red", axes=FALSE, xlab = "", ylab = "",
          xlim = c(0.5,top.pw.rows+0.5),
          ylim = c(0.5,top.pw.rows+0.5), zlim = zRange);
  }
  axis(1, at = 1:top.pw.rows, labels = popOrder, las = 2, line = -0.5,
       tick = 0, cex.axis = indLabelSize);
  axis(4, at = top.pw.rows:1, labels = popOrder, las = 2, line = -0.5,
       tick = 0, cex.axis = indLabelSize);
  if(drawRects){
    for(x in 1:(dim(popLimits)[2])){
      for(y in x:(dim(popLimits)[2])){
        rect(popLimits[1,y] - 0.5, top.pw.rows - popLimits[2,x] + 0.5,
             popLimits[2,y] + 0.5, top.pw.rows - popLimits[1,x] + 1.5, lwd = 2);
        if(showLabels){
          if(x == y){
            text(popLimits[1,y] + (popLimits[2,y] - popLimits[1,y]) * 0.75,
                 top.pw.rows -
                 (popLimits[1,x] + (popLimits[2,x] - popLimits[1,x]) * 0.25),
                 popNames[x], cex = textSize, adj = c(0.5,0.5),
                 srt = -45, offset = 0);
            cat("mean for ", popNames[x], ": ",
                mean(pw.matrix[popLimits[1,x]:popLimits[2,x],
                               popLimits[1,y]:popLimits[2,y]],
                     na.rm = TRUE), "\n", sep="", file=stderr());
            cat("SD for ", popNames[x], ": ",
                sd(as.vector(pw.matrix[popLimits[1,x]:popLimits[2,x],
                                       popLimits[1,y]:popLimits[2,y]]),
                   na.rm = TRUE), "\n", sep="", file=stderr());
          } else {
            text(popLimits[1,y] + (popLimits[2,y] - popLimits[1,y]) * 0.5,
                 top.pw.rows -
                 (popLimits[1,x] + (popLimits[2,x] - popLimits[1,x]) * 0.5),
                 paste(popNames[x],"vs",popNames[y]), cex = textSize,
                 adj = c(0.5,0.5), offset = 0);
            cat("mean for ", popNames[x], " vs ", popNames[y], ": ",
                mean(pw.matrix[popLimits[1,x]:popLimits[2,x],
                               popLimits[1,y]:popLimits[2,y]],
                     na.rm = TRUE), "\n", sep="", file=stderr());
            cat("SD for ", popNames[x], " vs ", popNames[y], ": ",
                sd(as.vector(pw.matrix[popLimits[1,x]:popLimits[2,x],
                                       popLimits[1,y]:popLimits[2,y]]),
                   na.rm = TRUE), "\n", sep="", file=stderr());
          }
        }
      }
    }
  }
} else {
  ##load library for heatmap clusters and sorts
  library(gtools, warn.conflicts = FALSE);
  library(gdata, warn.conflicts = FALSE);
  library(gplots, warn.conflicts = FALSE);
  if(makeDendrogram){
    heatmap.2(as.matrix(pw.matrix), symm = TRUE, trace = 'none', Rowv =
              TRUE, key = TRUE, dendrogram = 'row', density.info =
              'none', col = c(topo.colors(100), maxColour), main =
              paste("Pairwise similarity heatmap"), xlab = "Row ID",
              zlim = zRange, labCol = popOrder, labRow = popOrder,
              keysize = 1.2, cexRow = indLabelSize, cexCol = indLabelSize);
  } else {
    heatmap.2(as.matrix(pw.matrix), symm = TRUE, trace = 'none', Rowv =
              FALSE, key = TRUE, dendrogram = 'none', density.info =
              'none', col = c(topo.colors(100), maxColour), main =
              paste("Pairwise similarity heatmap"), xlab = "Row ID",
              zlim = zRange, labCol = popOrder, labRow = popOrder,
              keysize = 1.2, cexRow = indLabelSize, cexCol = indLabelSize);
  }
}
dummy <- dev.off();
