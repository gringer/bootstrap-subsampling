#!/usr/bin/Rscript

## Author: David Eccles (gringer), 2007 <programming@gringer.org>

## infocred.r -- carries out a difference-of-means test for difference between two populations

source("/itsshared/phd/common.r");

usage <- function(){
  cat("usage: ./infocred.r",
      "<file> <split point> [options]\n");
  cat("\nOther Options:\n");
  cat("-help               : Only display this help message\n");
  cat("notitle, hist, credible, showgreen, yrange, pointsize,svg, showdata, qrange,stdev");
  cat("\n");
}

if(commandArgs(TRUE)[1] == "-help"){
  usage();
  quit(save = "no", status=0);
}

cred.df <- read.table(commandArgs(TRUE)[1], header = TRUE)
splitPoint <- as.numeric(commandArgs(TRUE)[2]);

argLoc <- 3;
doGreen <- FALSE;
minYRange <- 0;
maxYRange <- 1;
pointSize <- 1;
doHist <- FALSE;
doSVG <- FALSE;
doTitle <- FALSE;
showData <- FALSE;
doSD <- TRUE;
doCredible <- FALSE;

while(!is.na(commandArgs(TRUE)[argLoc])){
  if(commandArgs(TRUE)[argLoc] == "-help"){
    usage();
    quit(save = "no", status=0);
  }
  if(commandArgs(TRUE)[argLoc] == "-notitle"){
    doTitle <- FALSE;
  }
  if(commandArgs(TRUE)[argLoc] == "-hist"){
    doHist <- TRUE;
  }
  if(commandArgs(TRUE)[argLoc] == "-credible"){
    doCredible <- TRUE;
  }
  if(commandArgs(TRUE)[argLoc] == "-showgreen"){
    doGreen <- TRUE;
  }
  if(commandArgs(TRUE)[argLoc] == "-yrange"){
    minYRange <- as.numeric(commandArgs(TRUE)[argLoc+1]);
    cat("Setting minYRange to",minYRange,"\n");
    maxYRange <- as.numeric(commandArgs(TRUE)[argLoc+2]);
    cat("Setting maxYRange to",maxYRange,"\n");
    argLoc <- argLoc + 2;
  }
  if(commandArgs(TRUE)[argLoc] == "-pointsize"){
    pointSize <- as.numeric(commandArgs(TRUE)[argLoc+1]);
    argLoc <- argLoc + 1;
  }
  if(commandArgs(TRUE)[argLoc] == "-svg"){
    library(cairoDevice);
    doSVG <- TRUE;
  }
  if(commandArgs(TRUE)[argLoc] == "-showdata"){
    showData <- TRUE;
  }
  if(commandArgs(TRUE)[argLoc] == "-qrange"){
    doSD <- FALSE;
  }
  if(commandArgs(TRUE)[argLoc] == "-stdev"){
    doSD <- TRUE;
  }
  argLoc <- argLoc + 1;
}

cred.df$Qrange <- abs(cred.df$Qlow - cred.df$Qhigh);
Q.matrix <- xtabs(Q ~ SNPs + ID,data = cred.df);
dataSize <- dim(Q.matrix)[2];
numSNPs <- dim(Q.matrix)[1];

p1Cluster <- (apply(Q.matrix[,1:splitPoint],1,mean) > 0.5) * 1;
p2Cluster <- (apply(Q.matrix[,(splitPoint+1):dataSize],1,mean) > 0.5) * 1;
Q.matrix.pops <-
t(rbind(apply(matrix(p1Cluster,1,length(p1Cluster)),2,rep, splitPoint),
        apply(matrix(p2Cluster,1,length(p2Cluster)),2,rep, dataSize - splitPoint)));
Q.matrix.diff <- abs(Q.matrix - Q.matrix.pops);

gs.diff <- abs(apply(Q.matrix.diff,1,function(x){x - Q.matrix.diff[numSNPs,]}));


Qrange.matrix <- xtabs(Qrange ~ SNPs + ID,data = cred.df);

Q.means <- rbind(apply(Q.matrix[,1:splitPoint],1,mean),
                 apply(Q.matrix[,(splitPoint+1):dataSize],1,mean));
Q.sd <- rbind(apply(Q.matrix[,1:splitPoint],1,sd),
              apply(Q.matrix[,(splitPoint+1):dataSize],1,sd));
Q.n <- rbind(apply(Q.matrix[,1:splitPoint],1,length),
             apply(Q.matrix[,(splitPoint+1):dataSize],1,length));
Qdiff.means <- abs(Q.means[1,] - Q.means[2,]);
Qdiff.sd <- sqrt(((Q.sd[1,]^2)/Q.n[1,]) + ((Q.sd[2,]^2)/Q.n[2]));

Qrange.means <- apply(Qrange.matrix,1,mean)

if(showData){
  cat("Difference of means:\n");
  print(Qdiff.means);
  print(Qdiff.means[Qdiff.means == max(Qdiff.means)]);
  cat("SD of difference:\n");
  print(Qdiff.sd);
  print(Qdiff.sd[Qdiff.sd == min(Qdiff.sd)]);
  print(round((1 - Qdiff.sd[2:length(Qdiff.sd)] /
         Qdiff.sd[1:(length(Qdiff.sd)-1)])*100),4);
  print(round(100 * Qdiff.sd[length(Qdiff.sd)] / Qdiff.sd),4);
  cat("Mean credible range:\n");
  print(Qrange.means);
  print(Qrange.means[Qrange.means == min(Qrange.means)]);
  print(round((1 - Qrange.means[2:length(Qrange.means)] / Qrange.means[1:(length(Qrange.means)-1)])*100),4);
  print(round(100 * Qrange.means[length(Qrange.means)] / Qrange.means),4);
  cat("mean, SD for",length(Qdiff.means),"SNPs:\n");
  print(Q.means[,length(Qdiff.means)]);
  print(Q.sd[,length(Qdiff.means)]);
}

if(doHist){
  if(doSVG){
    Cairo_svg("output_infohist.svg", width = 11, height = 8);
  } else {
    pdf("output_infohist.pdf", width = 11, height = 8, paper = "a4r");
  }
  par(mar=c(5,5,3,1));
  histcurve(Qdiff.means,
            main = paste("Histogram of Group Membership Difference Across",
              length(Q.means),"Samplings"), col = "red");
  histcurve(Qdiff.sd,
            main = paste("Histogram of Standard Deviation of Difference Across",
              length(Qdiff.sd),"Samplings"), col = "red");
  histcurve(Qrange.means,
            main = paste("Histogram of Mean Credible Region Size Across",
              length(Qrange.means),"Samplings"), col = "red");
  dummy <- dev.off();
}
if(doCredible){
  if(doSVG){
    Cairo_svg("output_credible.svg", width = 11, height = 8);
  } else {
    pdf("output_credible.pdf", width = 11, height = 8, paper = "a4r");
  }
  if(!doTitle){
    par(mar=c(5,5,1,1));
  } else {
    par(mar=c(5,5,3,1));
  }
  tmpcred <- subset(cred.df, No <= splitPoint);
  credmean.pop1 <- 1 - tapply(tmpcred$Qrange, tmpcred$SNPs, mean);
  tmpcred <- subset(cred.df, No > splitPoint);
  credmean.pop2 <- 1 - tapply(tmpcred$Qrange, tmpcred$SNPs, mean);
  tmpcred <- cred.df;
  credmean.all <- 1 - tapply(tmpcred$Qrange, tmpcred$SNPs, mean);
  plot(x = as.numeric(names(credmean.all)), y = credmean.all, ylim =
       c(minYRange,maxYRange), xlab = "Number of SNPs", main =
       "Q Value Accuracy (by credible range) for top n SNPs",
       ylab = "Mean Q Value Accuracy (%)", pch = NA, cex =
       0.5 * pointSize, col = "blue", cex.lab = 2, cex.axis = 1.5,
       frame.plot = FALSE);
  points(x = as.numeric(names(credmean.pop1)), y = credmean.pop1, pch =
         16, cex = 0.5 * pointSize, col = "green");
  points(x = as.numeric(names(credmean.pop2)), y = credmean.pop2, pch =
         16, cex = 0.5 * pointSize, col = "blue");
  points(x = as.numeric(names(credmean.all)), y = credmean.all, pch =
         16, cex = 0.5 * pointSize, col = "black");
  credmean.pop1 <- abs(1 - apply(Q.matrix.diff[,1:splitPoint], 1, mean));
  credmean.pop2 <- abs(1 - apply(Q.matrix.diff[,(splitPoint+1):dataSize], 1, mean));
  credmean.all <- 1 - apply(Q.matrix.diff, 1, mean);
  plot(x = as.numeric(names(credmean.all)), y = credmean.all, ylim =
       c(minYRange,maxYRange), xlab = "Number of SNPs", main =
       "Q Value Accuracy (by difference from 0/1) for top n SNPs",
       ylab = "Mean Q Value Accuracy (%)", pch = NA, cex =
       0.5 * pointSize, col = "blue", cex.lab = 2, cex.axis = 1.5,
       frame.plot = FALSE);
  points(x = as.numeric(names(credmean.pop1)), y = credmean.pop1, pch =
         16, cex = 0.5 * pointSize, col = "green");
  points(x = as.numeric(names(credmean.pop2)), y = credmean.pop2, pch =
         16, cex = 0.5 * pointSize, col = "blue");
  points(x = as.numeric(names(credmean.all)), y = credmean.all, pch =
         16, cex = 0.5 * pointSize, col = "black");
  credmean.all <- 1 - apply(Q.matrix.diff, 1, mean);
  credsd.all <- tapply(cred.df$Qrange, cred.df$SNPs, sd);
  credse.all <- sqrt(credsd.all / dataSize);
  plot(x = as.numeric(names(credmean.all)), y = credmean.all, ylim =
       c(minYRange,maxYRange), xlab = "Number of SNPs", main =
       "Q Value Accuracy (0/1 difference, with error as SE for credible range) for top n SNPs",
       ylab = "Mean Q Value Accuracy (%)", pch = NA, cex =
       0.5 * pointSize, col = "blue", cex.lab = 2, cex.axis = 1.5,
       frame.plot = FALSE);
  arrows(x0=as.numeric(names(credmean.all)),
         x1=as.numeric(names(credmean.all)),
         y0=apply(cbind(0,credmean.all - (credse.all / 2)),1,max),
         y1=apply(cbind(1,credmean.all + (credse.all / 2)),1,min),
         angle = 90, code = 3, length = 0.03*sqrt(pointSize), lwd = 1
         * pointSize);
  points(x = as.numeric(names(credmean.all)), y = credmean.all, pch =
         16, cex = 0.5 * pointSize, col = "black");
  mainLabel = paste("Q Value Accuracy (difference from ",numSNPs,"SNP result) for top n SNPs in pop1");
  if(!doTitle){mainLabel = ""};
  boxplot((1 - gs.diff[1:splitPoint,])*100, col = "green", main = mainLabel,
          xlab = "Number of SNPs", ylab = "Q Value Accuracy (%)", cex.lab = 2, cex.axis = 1.5);
  boxplot(1 - gs.diff[(splitPoint+1):dataSize,], col = "blue", main =
          paste("Q Value Accuracy (difference from ",numSNPs,"SNP result) for top n SNPs in pop2"),
          xlab = "Number of SNPs", ylab = "Q Value Accuracy (%)", cex.lab = 2, cex.axis = 1.5);
  boxplot(1 - gs.diff, col = "white", main =
          paste("Q Value Accuracy (difference from ",numSNPs,"SNP result) for top n SNPs"),
          xlab = "Number of SNPs", ylab = "Q Value Accuracy (%)", cex.lab = 2, cex.axis = 1.5);
  dummy <- dev.off();
} else {
  if(doSVG){
    Cairo_svg("output_diffpower.svg", width = 11, height = 8);
  } else {
    pdf("output_diffpower.pdf", width = 11, height = 8, paper = "a4r");
  }
  if(!doTitle){
    par(mar=c(5,5,1,1));
  } else {
    par(mar=c(5,5,3,1));
  }
  if(doSD){
    mainLabel <- paste("Group membership Difference for top n SNPs",
                       "(showing standard deviation)",
                       sep="\n");
  } else {
    mainLabel <- paste("Group membership Difference for top n SNPs",
                       "(with mean credible range across all individuals)",
                       sep="\n");
  }
  if(!doTitle){
    mainLabel <- "";
  }
  plot(x = as.numeric(names(Qdiff.means)), y = Qdiff.means, ylim =
       c(minYRange,maxYRange), xlab = "Number of SNPs", main =
       mainLabel, ylab = "Difference of mean Q", pch = NA, cex =
       0.5 * pointSize, col = "blue", cex.lab = 2, cex.axis = 1.5,
       frame.plot = FALSE);
  if(doSD){
    arrows(x0=as.numeric(names(Qdiff.means)),
           x1=as.numeric(names(Qdiff.means)),
           y0=apply(cbind(0,Qdiff.means - (Qdiff.sd / 2)),1,max),
           y1=apply(cbind(1,Qdiff.means + (Qdiff.sd / 2)),1,min),
           angle = 90, code = 3, length = 0.03*sqrt(pointSize), lwd = 1
           * pointSize);
  } else {
    arrows(x0=as.numeric(names(Qdiff.means)),
           x1=as.numeric(names(Qdiff.means)),
           y0=apply(cbind(0,Qdiff.means - (Qrange.means / 2)),1,max),
           y1=apply(cbind(1,Qdiff.means + (Qrange.means / 2)),1,min),
           angle = 90, code = 3, length = 0.03*sqrt(pointSize), lwd = 1
           * pointSize);
  }
  points(x = as.numeric(names(Qdiff.means)), y = Qdiff.means, pch =
       16, cex = 0.5 * pointSize, col = "blue");
  if(doGreen){
    if(doSD){
      points(x = as.numeric(names(Qrange.means)),
             y = Qrange.means, col = "green", pch = 25);
    } else {
      points(x = as.numeric(names(Qrange.means)),
             y = Qrange.means, col = "green", pch = 25);
    }
  }
  dummy <- dev.off();
}
