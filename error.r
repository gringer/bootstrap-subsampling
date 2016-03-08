#!/usr/bin/Rscript

## Author: David Eccles (gringer), 2009 <programming@gringer.org>

usage <- function(){
  cat("usage: ./error.r",
      "<file> <Case/Control Split location> [-flip]]\n");
  cat("Other options:\n");
  cat("-help          : show only this screen\n");
  cat("-flip          : invert positive/negative classification\n");
  cat("-vline <float> : place a vertical line at <location>\n");
  cat("-hline <float> : place a horizontal line at <location>\n");
  cat("-batch         : batch mode output\n");
  cat("-header        : display batch mode header\n");
  cat("-roc           : ROC analysis, with AUC calculation\n");
}

infile.df <- FALSE;
splitLoc <- FALSE;
flip <- FALSE;
doROC <- FALSE;
batchMode <- FALSE;
batchHeader <- FALSE;

vlineLocs <- NULL;
hlineLocs <- NULL;

argLoc <- 1;
while(!is.na(commandArgs(TRUE)[argLoc])){
  if(file.exists(commandArgs(TRUE)[argLoc])){ # file existence check
    if(infile.df == FALSE){
      infile.df <- read.table(commandArgs(TRUE)[argLoc], row.names = 1);
    } else{
      cat("Error: More than one input file specified\n");
      usage();
      quit(save = "no", status=1);
    }
  } else {
    if(length(grep("^[0-9]+$",commandArgs(TRUE)[argLoc])) > 0){
      splitLoc <- as.numeric(commandArgs(TRUE)[argLoc]);
    }
    if(commandArgs(TRUE)[argLoc] == "-flip"){
      flip <- TRUE;
    }
    if(commandArgs(TRUE)[argLoc] == "-roc"){
      doROC <- TRUE;
    }
    if(commandArgs(TRUE)[argLoc] == "-batch"){
      batchMode <- TRUE;
    }
    if(commandArgs(TRUE)[argLoc] == "-header"){
      batchHeader <- TRUE;
      batchMode <- TRUE;
    }
    if(commandArgs(TRUE)[argLoc] == "-vline"){
      vlineLocs <- c(vlineLocs, as.numeric(commandArgs(TRUE)[argLoc+1]));
      argLoc <- argLoc + 1;
    }
    if(commandArgs(TRUE)[argLoc] == "-hline"){
      hlineLocs <- c(hlineLocs, as.numeric(commandArgs(TRUE)[argLoc+1]));
      argLoc <- argLoc + 1;
    }
    if(commandArgs(TRUE)[argLoc] == "-help"){
      usage();
      quit(save = "no", status=0);
    }
  }
  argLoc <- argLoc + 1;
}

if((infile.df == FALSE) || (splitLoc == FALSE)){
  usage();
  quit(save = "no", status=2);
}

if(flip){
  infile.df$V2 <- 1 - infile.df$V2;
}

error.cutoff <- seq(0,1,length.out = 101);
if(doROC){
  error.cutoff <- c(0,sort(unique(infile.df[,"V2"])),1);
}
error.res <- data.frame(Cutoff = error.cutoff,
                        TP = tapply(error.cutoff, 1:length(error.cutoff),
                          function(x){length(which(infile.df[(1:splitLoc),"V2"] >= x))}) /
                        length(infile.df[(1:splitLoc),"V2"]),
                        FP = tapply(error.cutoff, 1:length(error.cutoff),
                          function(x){length(which(infile.df[-(1:splitLoc),"V2"] >= x))}) /
                        length(infile.df[-(1:splitLoc),"V2"]),
                        TN = tapply(error.cutoff, 1:length(error.cutoff),
                          function(x){length(which(infile.df[-(1:splitLoc),"V2"] < x))}) /
                        length(infile.df[-(1:splitLoc),"V2"]),
                        FN = tapply(error.cutoff, 1:length(error.cutoff),
                          function(x){length(which(infile.df[(1:splitLoc),"V2"] < x))}) /
                        length(infile.df[(1:splitLoc),"V2"]));

#error.res$PPV = error.res$TP / (error.res$TP + error.res$FP);
#error.res$NPV = error.res$TN / (error.res$TN + error.res$FN);

pdf("output_error.pdf",width = 8, height = 8);
if(doROC){
  plot(0.5,0.5,xlim = c(0,1), ylim = c(0,1), type = "n",
       xlab = "False Positive Rate (1 - specificity)",
       ylab = "True Positive Rate (sensitivity)");
  tpr <- error.res$TP / (error.res$TP+error.res$FN); # sensitivity
  fpr <- error.res$FP / (error.res$FP+error.res$TN); # 1 - specificity
  plotAUC <- -1 * sum((tpr[-1] - (diff(tpr)/2)) * diff(fpr));
  if(plotAUC < 0.5){
    tpr = 1 - tpr;
    fpr = 1 - fpr;
    plotAUC = 1 - plotAUC;
  }
  cornerDistance <- sqrt((1-tpr)^2 + (0-fpr)^2);
  closestPoint <- which(cornerDistance == min(cornerDistance))[1];
  points(fpr, tpr, type = "l", col = "black");
  points(fpr[closestPoint], tpr[closestPoint],
         type = "p", cex = 2, col = "red");
  ## points(fpr[round(length(fpr)/4)], tpr[round(length(tpr)/4)],
  ##        type = "p", cex = 2, col = "red");
  ## points(fpr[round(1 + length(fpr)/2)], tpr[round(1 + length(tpr)/2)],
  ##        type = "p", cex = 2, col = "red");
  ## points(fpr[round(length(fpr)*3/4)], tpr[round(length(tpr)*3/4)],
  ##        type = "p", cex = 2, col = "red");
  abline(h = hlineLocs, v = vlineLocs);
  ## AUC, assumes curve is broken into trapezoids, with tpr[1] being highest point on curve
  ## text(fpr[round(1 + length(fpr)/2)], tpr[round(1 + length(tpr)/2)],
  ##      paste("Median value = ",error.cutoff[round(1 + length(fpr)/2)]));
  text(fpr[closestPoint], tpr[closestPoint], round(error.cutoff[closestPoint],3), adj = c(-0.2,2));
  ## text(fpr[round(length(fpr)/4)], tpr[round(length(tpr)/4)], round(error.cutoff[round(length(tpr)/4)],3), adj = c(-0.2,2));
  ## text(fpr[round(1+length(fpr)/2)], tpr[round(1+length(tpr)/2)], round(error.cutoff[round(1+length(tpr)/2)],3), adj = c(-0.2,2));
  ## text(fpr[round(length(fpr)*3/4)], tpr[round(length(tpr)*3/4)], round(error.cutoff[round(length(tpr)*3/4)],3), adj = c(-0.2,2));
  text(0.75, 0.25, paste("AUC =",round(plotAUC,4)));
  lines(c(0,1), c(0,1), lty = 2);
  if(!batchMode){
    cat("AUC =",round(plotAUC,4),"; ");
  }
  if(length(vlineLocs) > 0){
    cutoffLoc <- which(abs(fpr - vlineLocs[1]) == min(abs(fpr - vlineLocs[1])))[1];
    cat("Value closest to FPR of ", vlineLocs[1], " = ", error.cutoff[cutoffLoc],
        " (TPR = ",tpr[cutoffLoc],", FPR = ",
        fpr[cutoffLoc],") [counts TP=",error.res$TP[cutoffLoc], ", FP=",error.res$FP[cutoffLoc],
        ", TN=",error.res$TN[cutoffLoc], ", FN=",error.res$FN[cutoffLoc], ", RR = ",error.res$TP[cutoffLoc]/error.res$FP[cutoffLoc],"]\n",sep="");
  }
  if (length(hlineLocs) > 0){
    cutoffLoc <- which(abs(tpr - hlineLocs[1]) == min(abs(tpr - hlineLocs[1])))[1];
    cat("Value closest to TPR of ", hlineLocs[1], " = ", error.cutoff[cutoffLoc],
        " (TPR = ",tpr[cutoffLoc],", FPR = ",
        fpr[cutoffLoc],") [counts TP=",error.res$TP[cutoffLoc], ", FP=",error.res$FP[cutoffLoc],
        ", TN=",error.res$TN[cutoffLoc], ", FN=",error.res$FN[cutoffLoc], ", RR = ",error.res$TP[cutoffLoc]/error.res$FP[cutoffLoc],"]\n",sep="");
  }
  if(batchMode){
    if(batchHeader){
      cat("AUC,Cutoff,TPR,FPR\n");
    }
    cat(round(plotAUC,4),",",error.cutoff[closestPoint],",",tpr[closestPoint],",",fpr[closestPoint],"\n",sep="");
  } else {
    cat("Closest cutoff to TPR 1, FPR 0 = ",error.cutoff[closestPoint],
        " (TPR = ",tpr[closestPoint],", FPR = ",
        fpr[closestPoint],")\n",sep="");
  }
} else {
  plot(0.5,0.5,xlim = c(0,1), ylim = c(0,1), type = "n", xlab = "Cutoff value", ylab = "Rate (as proportion of population size)");
  points(error.res$Cutoff,error.res$TP, type = "l", col = "green");
  points(error.res$Cutoff,error.res$FP, type = "l", col = "red");
  points(error.res$Cutoff,error.res$TN, type = "l", col = "blue");
  points(error.res$Cutoff,error.res$FN, type = "l", col = "grey");
  legend("top", legend = c("TP","FP","TN","FN"), fill = c("green","grey","blue","red"), inset = 0.01, bg = "white");
}
dummy <- dev.off();
