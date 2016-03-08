#!/usr/bin/Rscript

## Author: David Eccles (gringer), 2007 <programming@gringer.org>

source('/itsshared/phd/common.r');

#usage ./structure2pdf.r <file> <name1> <break1> <name2> <break2> ... <nameN>

usage <- function(){
  cat("usage: ./structure2pdf.r",
      "<file>( p<name> <range from> <range to>)* [options]\n");
  cat("\nOther Options:\n");
  cat("-help               : Only display this help message\n");
  cat("-gdionly   [K=2]    : Only calculate Genome Diagnostic Index\n");
  cat("-line <value>       : Draw a horizontal line at <value>\n");
  cat("-sort               : Sort individuals by Q values\n");
  cat("-basicsort          : Sort individuals strictly by Q values\n");
  cat("-barplot   [K=2]    : Always do a barplot (rather than scatterplot)\n");
  cat("-error     [k=2]    : Show error bars from \"_f\" file\n");
  cat("-mean      [K=2]    : Draw mean and SE lines for each popualtion\n");
  cat("-noshade   [K=2]    : Don't draw background coloured stripes\n");
  cat("-nostack   [K>2]    : Don't put highest Q on the bottom\n");
  cat("-stacksort [K>2]    : Sort each individual's bar vertically\n");
  cat("-halfheight         : Output a PDF file half the usual height\n");
  cat("-svg                : Output to an SVG file (instead of PDF)\n");
  cat("-pointsize <value>  : Scaling factor for point size\n");
  cat("-flip      [K=2]    : swap Q values for first and second clusters\n");
  cat("-rotatelabels       : Make population labels display vertically\n");
  cat("-labelaxis          : Place population labels on the axis\n");
  cat("c<name>             : Set the colour for the next population to <name>\n");
  cat("\n");
}

invertOrder <- function(x){
  order(x, decreasing = TRUE)[order(x, decreasing = TRUE)];
}

orderHighest <- function(x){
  highOrder <- order(x, decreasing = TRUE);
  highestNum <- highOrder[1];
  swapOrder <- 1:length(x);
  swapOrder[highestNum] <- 1;
  swapOrder[1] <- highestNum;
  return(swapOrder);
}

sortHighest <- function(x){
  return(x[orderHighest(x)]);
}


popNames <- character(0);
popLimits <- NULL;
breakstyle <- integer(0);
infile.name <- FALSE;
prob.df <- FALSE;
flip <- FALSE;
doSVG <- FALSE;
sortByQ <- FALSE;
basicSort <- FALSE;
labelAxis <- FALSE;
rotateLabels <- FALSE;
halfHeight <- FALSE;
pointSize <- 1;
popColours <- FALSE;
shadedBG <- TRUE;
onlyBarPlot <- FALSE;
showError <- FALSE;
drawMeans <- FALSE;
sortStacks <- TRUE;
stackFlip <- TRUE;
gdiOnly <- FALSE;
horizLines <- NULL;
qlow.df <- NULL;
qhigh.df <- NULL;

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
  } else {
    if(commandArgs(TRUE)[argLoc] == "-help"){
      usage();
      quit(save = "no", status=0);
    }
    if(substr(commandArgs(TRUE)[argLoc],1,1) == "p"){
      inName <- substring(commandArgs(TRUE)[argLoc],2);
      if(substr(inName,1,1) == '_'){
        breakstyle <- c(breakstyle,1);
        inName <- substr(inName, 2, nchar(inName));
      }
      else{
        breakstyle <- c(breakstyle,2);
      }
      popNames <- append(popNames, sub("_"," ",inName));
      popLimits <- cbind(popLimits,c(as.numeric(commandArgs(TRUE)[argLoc+1]),
                                     as.numeric(commandArgs(TRUE)[argLoc+2])));
      argLoc <- argLoc + 2;
    }
    if(commandArgs(TRUE)[argLoc] == "-gdionly"){
      gdiOnly <- TRUE;
    }
    if(commandArgs(TRUE)[argLoc] == "-line"){
      horizLines = c(horizLines, as.numeric(commandArgs(TRUE)[argLoc+1]));
      argLoc <- argLoc + 1;
    }
    if(commandArgs(TRUE)[argLoc] == "-barplot"){
      onlyBarPlot <- TRUE;
    }
    if(commandArgs(TRUE)[argLoc] == "-error"){
      showError <- TRUE;
    }
    if(commandArgs(TRUE)[argLoc] == "-sort"){
      sortByQ <- TRUE;
    }
    if(commandArgs(TRUE)[argLoc] == "-basicsort"){
      sortByQ <- TRUE;
      basicSort <- TRUE;
    }
    if(commandArgs(TRUE)[argLoc] == "-mean"){
      drawMeans <- TRUE;
    }
    if(commandArgs(TRUE)[argLoc] == "-noshade"){
      shadedBG <- FALSE;
    }
    if(commandArgs(TRUE)[argLoc] == "-nostack"){
      sortStacks <- FALSE;
    }
    if(commandArgs(TRUE)[argLoc] == "-stacksort"){
      sortStacks <- TRUE;
      stackFlip <- FALSE;
    }
    if(commandArgs(TRUE)[argLoc] == "-halfheight"){
      halfHeight <- TRUE;
    }
    if(commandArgs(TRUE)[argLoc] == "-svg"){
      library(cairoDevice);
      doSVG <- TRUE;
    }
    if(commandArgs(TRUE)[argLoc] == "-pointsize"){
      pointSize <- as.numeric(commandArgs(TRUE)[argLoc+1]);
      cat("Setting point size to ",pointSize,"\n",sep="");
      argLoc <- argLoc + 1;
    }
    if(commandArgs(TRUE)[argLoc] == "-flip"){
      flip <- TRUE;
    }
    if(commandArgs(TRUE)[argLoc] == "-rotatelabels"){
      labelAxis <- TRUE;
      rotateLabels <- TRUE;
    }
    if(commandArgs(TRUE)[argLoc] == "-labelaxis"){
      labelAxis <- TRUE;
    }
    if(substr(commandArgs(TRUE)[argLoc],1,1) == "c"){
      inName <- substring(commandArgs(TRUE)[argLoc],2);
      if(popColours[1] == FALSE){
        popColours <- inName;
      }
      else{
        popColours <- c(popColours,inName);
      }
    }
  }
  argLoc <- argLoc + 1;
}

if(infile.name == FALSE){
  cat("Error: No valid file given\n\n");
  usage();
  quit(save = "no", status=1);
} else {
  tmpWarn <- getOption("warn");
  infile.con <- file(infile.name);
  open(infile.con);
  options(warn = -1);
  tmpLines <- readLines(infile.con);
  options(warn = tmpWarn); rm(tmpWarn);
  idx <- grep("([0-9] individuals|Inferred clusters)",tmpLines);
  if(length(idx) == 2){
    ## Assume "_f"-style file [standard Structure output]
    numIndivs <- as.numeric(gsub("[^0-9]","",tmpLines[idx[1]]));
    seek(infile.con, 0);
    tmp.df <- read.table(infile.con, nrows = numIndivs, skip = idx[2],
                         row.names = 2, stringsAsFactors = FALSE);
    colnames(tmp.df)[c(1,2,3)] <- c("Line","Missing%","Separator");
    tmp.df$"Missing%" = as.numeric(gsub("[^0-9]","",tmp.df$"Missing%"))
    tmp.df$"Missing%" = NULL; ## simplifies conversion later
    tmp.df$Line = NULL
    tmp.df$Separator = NULL;
    idx <- grep("\\(",tmp.df[1,]);
    if(length(idx) > 0){
      ## probability intervals exist, so we put them in prob.df
      infile.df <- tmp.df[1:(idx[1]-1)]; # copy over data table
      tmp.df <- tmp.df[-(1:(idx[1]-1))]; # remove Q values from temp table
      tmp.names <- rownames(tmp.df);
      tmp.df <- data.frame(lapply(tmp.df,gsub,pattern="[^0-9,\\.]",
                                  replacement=""), stringsAsFactors = FALSE);
      # extract low/high limits, assumes #.###,#.### (i.e. 3 d.p)
      qlow.df <- data.frame(lapply(tmp.df,substring,1,5),
                            stringsAsFactors = FALSE);
      # convert character to numeric value
      qlow.df <- data.frame(lapply(qlow.df,as.numeric));
      rownames(qlow.df) <- tmp.names;
      colnames(qlow.df) <- paste("Q",1:(dim(qlow.df)[2]), sep = "")
      qhigh.df <- data.frame(lapply(tmp.df,substring,7,11),
                             stringsAsFactors = FALSE);
      # convert character to numeric value
      qhigh.df <- data.frame(lapply(qhigh.df,as.numeric));
      rownames(qhigh.df) <- tmp.names;
      colnames(qhigh.df) <- paste("Q",1:(dim(qhigh.df)[2]), sep = "")
    } else {
      infile.df <- tmp.df; # no probability intervals, so we're done
    }
    rm(tmp.df); # clean up
    rm(numIndivs);
  } else{
    ## Assume "_q"-style file
    ## NOTE: "_q" files seem to have 1 d.p. more precision in results
    seek(infile.con, 0);
    infile.df <- read.table(infile.con, row.names = 1);
  }
  rm(tmpLines);
  close(infile.con);
}

cat("Input file has data for ",dim(infile.df)[1]," individuals\n",sep="");

if(length(popNames) == 0){
  cat("Warning: no populations defined.\n",
      "         All individuals will be treated as the same population\n"
      , sep="");
  popNames <- "All";
  popLimits <- rbind(1,dim(infile.df)[1]);
  breakstyle <- c(breakstyle,1);
} else {
  totLength <- 0;
  for(x in 1:dim(popLimits)[2]){
    totLength <- totLength + (popLimits[2,x] - popLimits[1,x] + 1);
    if((popLimits[2,x] > dim(infile.df)[1]) ||
        (popLimits[1,x] < 1)){
      cat("Error: range of population", popNames[x] ,
          " outside range of individuals",
          "in the file\n");
      usage();
      quit(save = "no", status = 2);
    }
  }
  if(totLength < dim(infile.df)[1]){
    cat("Warning: fewer individuals specified than present in the file\n");
  }
  if(totLength > dim(infile.df)[1]){
    cat("Warning: more individuals specified than present in the file\n");
  }
}

nextPos <- 1;
tmpfile.df <- NULL;
tmpLimits <- NULL;
tmplow.df <- NULL;
tmphigh.df <- NULL;
for(x in 1:dim(popLimits)[2]){
  tmpLimits <- cbind(tmpLimits,c(nextPos,
                                 nextPos +
                                 (popLimits[2,x] - popLimits[1,x])));
  cat("extracting population '",popNames[x],
      "' (",popLimits[2,x] - popLimits[1,x] + 1," individuals)\n",sep="");
  tmpfile.df <- rbind(tmpfile.df,
                      infile.df[popLimits[1,x]:popLimits[2,x],]);
  tmplow.df <- rbind(tmplow.df,
                      qlow.df[popLimits[1,x]:popLimits[2,x],]);
  tmphigh.df <- rbind(tmphigh.df,
                      qhigh.df[popLimits[1,x]:popLimits[2,x],]);
  nextPos <- tmpLimits[2,x] + 1;
}

infile.df <- data.frame(tmpfile.df);
qlow.df <- data.frame(tmplow.df);
qhigh.df <- data.frame(tmphigh.df);
popLimits <- tmpLimits;

colnames(infile.df) <- paste("Q",1:(dim(infile.df)[2]), sep = "")
if(flip){ # flip order of the first and second population
  tmp <- infile.df$Q1;
  infile.df$Q1 <- infile.df$Q2;
  infile.df$Q2 <- tmp;
  qlow.df <- 1 - qlow.df;
  qhigh.df <- 1 - qhigh.df;
}

if(sortByQ){ # sort based on max Q within each population
  cat("sorting...");
  maxPopFormat <- paste("%0",nchar(dim(infile.df)[2]),"d",sep="");
  for(x in 1:dim(popLimits)[2]){
    tmp.df <- infile.df[popLimits[1,x]:popLimits[2,x],];
    # work out position of maximum Q per individual
    ## Sort 0 -- sort by dominant Q for the population
    popMaxQ <- unlist(tmp.df[which(colSums(tmp.df)
                                   == max(colSums(tmp.df)))]);
    ## Sort 1 -- Q > 0.5
    indLargeQpos <- apply((tmp.df > 0.5) *
                          rep(seq(1,dim(tmp.df)[2]),
                              each = dim(tmp.df[1])),1,max);
    ## Sort 2 -- greater than mean Q
    indMeanQpos <- apply((tmp.df > apply(tmp.df,1,mean))*
                        rep(seq(1,dim(tmp.df)[2]),
                            each = dim(tmp.df[1])),1,max);
    ## Sort 3 -- max Q
    indMaxQpos <- apply((tmp.df == apply(tmp.df,1,max))*
                        rep(seq(1,dim(tmp.df)[2]),
                            each = dim(tmp.df[1])),1,max);
    maxCounts <- # how many times a population is max for each individual
      colSums((tmp.df) == apply(tmp.df,1,max));
    ## Make sure populations with the same counts don't end up alternating
    maxCounts[duplicated(maxCounts)] <-
      seq(0.1,0.4, length.out = length(which(duplicated(maxCounts))));
    ## Ordering based on Q ranking of populations for each individual (not currently used)
    indQOrder <- order(apply(apply(apply(tmp.df,1,order),1,sprintf,fmt=maxPopFormat),1,paste,collapse=""));
    # work out value of maximum Q per individual
    indMaxQ <- apply(tmp.df * (tmp.df == apply(tmp.df,1,max)),1,max);
    indMeanQ <- apply(tmp.df * (tmp.df == apply(tmp.df,1,mean)),1,max);
    indMinQ <- apply(tmp.df * (tmp.df == apply(tmp.df,1,min)),1,max);
    ## older sort -- by maximum Q for all individuals
    #sumRank <- apply(infile.df[popLimits[1,x]:popLimits[2,x],],2,sum);
    #sortRank <- infile.df[popLimits[1,x]:popLimits[2,x],
    #                      sumRank == max(sumRank)];
    #infile.df[popLimits[1,x]:popLimits[2,x],] <-
    #  infile.df[order(sortRank)+popLimits[1,x]-1,];
    # Sort by max Q position, then by max Q
    if(basicSort){
      infile.df[popLimits[1,x]:popLimits[2,x],] <-
        infile.df[order(popMaxQ,decreasing = TRUE)+popLimits[1,x]-1,];
    } else {
      infile.df[popLimits[1,x]:popLimits[2,x],] <-
        infile.df[order(maxCounts[indMaxQpos],indMaxQ,decreasing = TRUE)+popLimits[1,x]-1,];
    }
    if(!is.null(qlow.df)){
      if(basicSort){
        qlow.df[popLimits[1,x]:popLimits[2,x],] <-
          qlow.df[order(popMaxQ,decreasing = TRUE)+popLimits[1,x]-1,];
        qhigh.df[popLimits[1,x]:popLimits[2,x],] <-
          qhigh.df[order(popMaxQ,decreasing = TRUE)+popLimits[1,x]-1,];
      } else {
        qlow.df[popLimits[1,x]:popLimits[2,x],] <-
          qlow.df[order(maxCounts[indMaxQpos],indMaxQ,decreasing = TRUE)+popLimits[1,x]-1,];
        qhigh.df[popLimits[1,x]:popLimits[2,x],] <-
          qhigh.df[order(maxCounts[indMaxQpos],indMaxQ,decreasing = TRUE)+popLimits[1,x]-1,];
      }
    }
  }
  cat("done\n");
}

if(doSVG){
  if(halfHeight){
    Cairo_svg("output_Q.svg", width = 11, height = 4);
  } else {
    Cairo_svg("output_Q.svg", width = 11, height = 8);
  }
} else {
 if(halfHeight){
   pdf("output_Q.pdf", width = 11, height = 4);
 } else {
   pdf("output_Q.pdf", width = 11, height = 8, paper = 'a4r');
 }
}
par(mar=c(5,6,1,1));
if((onlyBarPlot) || (length(colnames(infile.df)) > 2)){
  if(popColours[1] == FALSE){
    popColours <- rainbow(length(colnames(infile.df)));
  }
  if(sortStacks){
    if(stackFlip){
      popColours <- popColours[apply(infile.df,1,orderHighest)];
      barMat <- apply(infile.df,1,sortHighest);
    }
    else{
      popColours <- popColours[apply(infile.df,1,order, decreasing = TRUE)];
      barMat <- apply(infile.df,1,sort, decreasing = TRUE);
    }
  }
  else{
    barMat <- t(as.matrix(infile.df[,1:length(colnames(infile.df))]));
  }
  barMat <- (barMat*1.001) / colSums(barMat);
    ## set up plot, but don't draw bars/blocks
    barCentres <-
      barplot(barMat, col = NA, border
              = NA, space = 0, ylab = "Q", xlab =
              "Individual", cex.lab = 2, cex.axis = 1.5, las = 2, xaxt = "n",
              ylim = c(0,1), mgp = c(4,1,0));
    rect(xleft = rep(barCentres - 0.5, each = dim(barMat)[1]),
         xright = rep(barCentres + 0.5, each = dim(barMat)[1]),
         ybottom = rbind(0,apply(t(barMat),1,cumsum)[1:(dim(barMat)[1]-1),]),
         ytop =    rbind(  apply(t(barMat),1,cumsum)[2:(dim(barMat)[1]),],1),
         col = popColours, border = popColours, lwd = 0.5);
  if(length(horizLines)>0){
    abline(h = horizLines, lwd = 3, col = "grey");
  }
  for(x in 1:dim(popLimits)[2]){
    if(x > 1){
      lines(x=rep(barCentres[popLimits[1,x]]-0.5,2),y=c(0,1),
            lty = breakstyle[x], lwd = 1);
    }
    if(!labelAxis){
      srtLabel = 0;
      if(rotateLabels){
        srtLabel = 90;
      }
      text((barCentres[popLimits[1,x]] +
            barCentres[popLimits[2,x]])/2, 0.5, popNames[x],
           cex = 2, srt = srtLabel);
    }
  }
  if(labelAxis){
    axLas <- 0;
    if(rotateLabels){
      axLas <- 2;
    }
    axis(side = 1,labels = popNames, at = (barCentres[popLimits[1,]] +
                    barCentres[popLimits[2,]])/2, col = NA, cex.axis =
         1.5, las = axLas);
  }
} else {
  if(!gdiOnly){
    if(halfHeight){
      plot(infile.df$Q1,
           ylab = "Q",
           pch = NA, xlab = "Individual", main = "", las = 2,
           cex.lab = 2, cex.axis = 1.5, mgp=c(4,1,0), ylim = c(0,1));
    } else {
      plot(infile.df$Q1,
           ylab = "Estimated Ancestral Fraction (Q)",
           pch = NA, xlab = "Individual", main = "", las = 2,
           cex.lab = 2, cex.axis = 1.5, mgp=c(4,1,0), ylim = c(0,1));
    }
    if(length(horizLines)>0){
      abline(h = horizLines, lwd = 3, col = "grey");
    }
    for(x in 1:dim(popLimits)[2]){
      if(shadedBG){
        bgShade(popLimits[1,x]-0.5,popLimits[2,x]+0.5,50, density = 25,
                cl = c(rgb(0.5+seq1(25)*0.5,0.5+seq1(25)*0.5,1),
                  rgb(1-seq1(25)*0.5,1,1-seq1(25)*0.5)));
        ## redraw cutoff lines over shading
        if(length(horizLines)>0){
          abline(h = horizLines, lwd = 3, col = "grey");
        }
      }
      if(x > 1){
        lines(x=rep(popLimits[1,x]-0.5,2),y=c(0,1),
              lty = breakstyle[x], lwd = 1);
      }
      text((popLimits[1,x] + popLimits[2,x])/2,0.5,popNames[x],cex = 2);
    }
    par(new=TRUE)
    if(showError && (!is.null(qlow.df))){
      arrows(x0 = 1:(dim(infile.df)[1]), x1 = 1:(dim(infile.df)[1]),
             y0 = qlow.df$Q1, y1 = qhigh.df$Q1, angle = 90, code = 3,
             length = 0.03, lwd = 2 * pointSize, col = "#EF6464");
    }
    points(x=1:dim(infile.df)[1],y=infile.df$Q1, col="#000000",
           pch = 16, cex = 1 * pointSize);
  }
  for(x in 1:dim(popLimits)[2]){
    if(drawMeans){
      qdata <- infile.df[popLimits[1,x]:popLimits[2,x],1];
      sDev <- sd(qdata);
      cat(sprintf("Mean Q for population %d: %f (SD = %f)\n",x, mean(qdata), sDev));
      if(!gdiOnly){
        lines(x=c(popLimits[1,x]-0.5,popLimits[2,x]+0.5),y=rep(mean(qdata) + sDev,2), lwd = 1, col = "red");
        lines(x=c(popLimits[1,x]-0.5,popLimits[2,x]+0.5),y=rep(mean(qdata),2), lwd = 2, col = "grey");
        lines(x=c(popLimits[1,x]-0.5,popLimits[2,x]+0.5),y=rep(mean(qdata) - sDev,2), lwd = 1, col = "red");
      }
    }
  }
  if(dim(popLimits)[2] == 2){
    cat(sprintf("Genome Diagnostic Index (GDI): %f\n",
                abs(mean(infile.df[popLimits[1,1]:popLimits[2,1],1]) -
                    mean(infile.df[popLimits[1,2]:popLimits[2,2],1]))));
    cat(sprintf("Mean (SD) for populations: %f (%f), %f (%f)\n",
        mean(infile.df[popLimits[1,1]:popLimits[2,1],1]),
        sd(infile.df[popLimits[1,1]:popLimits[2,1],1]),
        mean(infile.df[popLimits[1,2]:popLimits[2,2],1]),
        sd(infile.df[popLimits[1,2]:popLimits[2,2],1])));
    if(flip){
      print(t.test(1 - infile.df[popLimits[1,1]:popLimits[2,1],1],
                   1 - infile.df[popLimits[1,2]:popLimits[2,2],1]));
    } else {
      print(t.test(infile.df[popLimits[1,1]:popLimits[2,1],1],
                   infile.df[popLimits[1,2]:popLimits[2,2],1]));
    }
  }
}
dummy <- dev.off();
