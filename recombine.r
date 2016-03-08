#!/usr/bin/Rscript

#usage ./recombine.r

## Author: David Eccles (gringer), 2008 <programming@gringer.org>

usage <- function(){
  cat("usage: ./recombine.r\n");
  cat("\nOther Options:\n");
  cat("-help               : Only display this help message\n");
  cat("\n");
}

argLoc <- 1;
while(!is.na(commandArgs(TRUE)[argLoc])){
  if(commandArgs(TRUE)[argLoc] == "-help"){
    usage();
    quit(save = "no", status=0);
  }
  argLoc <- argLoc + 1;
}

makeUniformChromosome <- function(colour = "NA"){
  if(colour == "NA"){
    colour = sample(colors(),1);
  }
  return(data.frame(Block = colour, Points = 0));
}

drawChromosome <- function(ChrX, xPos = 0.5, yPos = 0, scale = 1,
                           yScale = 0.5, new = TRUE){
  xLeft = xPos;
  yBottom = yPos;
  xRight = xPos + 0.05 * scale;
  if(new){
    plot.new();
  }
  for(i in 1:(dim(ChrX)[1])){
    if(i == dim(ChrX)[1]){
      rect(xLeft,yBottom+ChrX[i,2]*scale*yScale,xRight,
           yBottom+1*scale*yScale,
           col=as.character(ChrX[i,1]));
    } else {
      rect(xLeft,yBottom+ChrX[i,2]*scale*yScale,xRight,
           yBottom+ChrX[i+1,2]*scale*yScale,
           col=as.character(ChrX[i,1]));
    }
  }
}

recombine <- function(ChrX, ChrY, minPoints = 3, maxPoints = 4){
  result <- NULL;
  recombPoints <- sort(runif(sample(minPoints:maxPoints,1)));
  currentChr <- sample(c(0,1),1);
  startPos <- 0;
  newPoints <- NULL;
  for(endPos in c(recombPoints,1)){
    if(currentChr == 0){
      stopX <- c(ChrX$Points[-1],1);
      # list of recombination points on ChrX in the correct interval
      curPoints = which((ChrX$Points<=endPos) &
        (stopX >= startPos));
      newPoints <- ChrX[curPoints,];
    } else {
      stopY <- c(ChrY$Points[-1],1);
      # list of recombination points on ChrY in the correct interval
      curPoints = which((ChrY$Points<=endPos) &
        (stopY >= startPos));
      newPoints <- ChrY[curPoints,];
    }
    # set start recombination point to be where the recombination took place
    newPoints$Points[1] <- startPos;
    result <- rbind(result,newPoints);
    startPos <- endPos;
    currentChr <- (currentChr + 1) %% 2;
  }
  return(result);
}

# merges regions of the same colour
chrClean <- function(ChrX){
  result <- ChrX[1,];
  oldColour <- ChrX$Block[1];
  for(pos in 2:(dim(ChrX)[1])){
    if(oldColour != ChrX$Block[pos]){
##      cat(oldColour, "does not equal", ChrX$Block[pos],"\n");
      result <- rbind(result,ChrX[pos,]);
      oldColour <- ChrX$Block[pos];
    }
  }
  return(result);
}

# combines two chromosomes, so that regions of target colour common to
# both input chromosomes are in the result chromosome
chrUnion <- function(colourName, ChrX, ChrY, ChrOther = NULL, ...){
  if(is.null(ChrOther)){
    result <- rbind(ChrX, ChrY);
    result <- result[order(result$Points),];
    for(pos in 1:(dim(result)[1])){
      blockX <- ChrX$Block[max(which(ChrX$Points <= result$Points[pos]))];
      blockY <- ChrY$Block[max(which(ChrY$Points <= result$Points[pos]))];
      if(blockX == colourName){
        result$Block[pos] <- blockX;
      }
      if(blockY == colourName){
        result$Block[pos] <- blockY;
      }
    }
    result <- chrClean(result);
    return(result);
  } else {
    return(chrUnion(colourName, chrUnion(colourName, ChrX, ChrY), ChrOther, ...));
  }
}

## recombination descendancy, haplotype dilution
oc <- makeUniformChromosome("yellow");
g1 <- list(makeUniformChromosome("red"));
g2 <- list(
           chrClean(recombine(g1[[1]],oc)),
           chrClean(recombine(g1[[1]],oc))
           );
g2u<- chrUnion("red",g2[[1]],g2[[2]]);
g3 <- list(
           chrClean(recombine(g2[[1]],oc)),
           chrClean(recombine(g2[[1]],oc)),
           chrClean(recombine(g2[[2]],oc)),
           chrClean(recombine(g2[[2]],oc))
           );
g3u<- chrUnion("red",g3[[1]],g3[[2]],g3[[3]],g3[[4]]);
g4 <- list(
           chrClean(recombine(g3[[1]],oc)),
           chrClean(recombine(g3[[1]],oc)),
           chrClean(recombine(g3[[2]],oc)),
           chrClean(recombine(g3[[2]],oc)),
           chrClean(recombine(g3[[3]],oc)),
           chrClean(recombine(g3[[3]],oc)),
           chrClean(recombine(g3[[4]],oc)),
           chrClean(recombine(g3[[4]],oc))
           );
g4u<- chrUnion("red",g4[[1]],g4[[2]],g4[[3]],g4[[4]],
                g4[[5]],g4[[6]],g4[[7]],g4[[8]]);
g5 <- list(
           chrClean(recombine(g4[[1]],oc)),
           chrClean(recombine(g4[[1]],oc)),
           chrClean(recombine(g4[[2]],oc)),
           chrClean(recombine(g4[[2]],oc)),
           chrClean(recombine(g4[[3]],oc)),
           chrClean(recombine(g4[[3]],oc)),
           chrClean(recombine(g4[[4]],oc)),
           chrClean(recombine(g4[[4]],oc)),
           chrClean(recombine(g4[[5]],oc)),
           chrClean(recombine(g4[[5]],oc)),
           chrClean(recombine(g4[[6]],oc)),
           chrClean(recombine(g4[[6]],oc)),
           chrClean(recombine(g4[[7]],oc)),
           chrClean(recombine(g4[[7]],oc)),
           chrClean(recombine(g4[[8]],oc)),
           chrClean(recombine(g4[[8]],oc))
           );
g5u <- chrUnion("red",g5[[1]],g5[[2]],g5[[3]],g5[[4]],
                g5[[5]],g5[[6]],g5[[7]],g5[[8]],
                g5[[9]],g5[[10]],g5[[11]],g5[[12]],
                g5[[13]],g5[[14]],g5[[15]],g5[[16]]);

par(mar=c(0,0,0,0));
plot.new();
plot.window(xlim=c(0,1),ylim=c(0.2,1));
for(i in 1:8){
  drawChromosome(g5[[i*2-1]],xPos = (i-1)*0.125, yPos = 0.8,
                 scale = 0.4, new = FALSE);
  drawChromosome(g5[[i*2]],xPos = (i-1)*0.125+0.05, yPos = 0.8,
                 scale = 0.4, new = FALSE);
}
drawChromosome(g5u,xPos = (7)*0.125+0.125, yPos = 0.8, scale = 0.4, new = FALSE);
for(i in 1:4){
  drawChromosome(g4[[i*2-1]],xPos = (i-1)*0.25+0.0625, yPos = 0.58,
                 scale = 0.4, new = FALSE);
  drawChromosome(g4[[i*2]],xPos = (i-1)*0.25+0.0625+0.05, yPos = 0.58,
                 scale = 0.4, new = FALSE);
}
drawChromosome(g4u,xPos = (7)*0.125+0.095, yPos = 0.58, scale = 0.4, new = FALSE);
for(i in 1:2){
  drawChromosome(g3[[i*2-1]],xPos = (i-1)*0.5+0.1875, yPos = 0.4,
                 scale = 0.4, new = FALSE);
  drawChromosome(g3[[i*2]],xPos = (i-1)*0.5+0.1875+0.05, yPos = 0.4,
                 scale = 0.4, new = FALSE);
}
drawChromosome(g3u,xPos = (7)*0.125+0.065, yPos = 0.4, scale = 0.4, new = FALSE);
drawChromosome(g2[[1]],xPos = 0.4375, yPos = 0.22,
               scale = 0.4, new = FALSE);
drawChromosome(g2[[2]],xPos = 0.4375+0.05, yPos = 0.22,
               scale = 0.4, new = FALSE);
drawChromosome(g2u,xPos = (7)*0.125+0.035, yPos = 0.22, scale = 0.4, new = FALSE);
drawChromosome(g1[[1]],xPos = 0.0875, yPos = 0.22, scale = 0.4, new = FALSE);
arrows(x1=seq(0.035,0.91,length.out=8),
       x0=seq(0.035,0.91,length.out=8) + rep(c(0.015,-0.015),4),
       y1=0.79, y0=0.68, length = 0.1);
arrows(x1=seq(0.0975,0.8475,length.out=4),
       x0=seq(0.0975,0.8475,length.out=4) + rep(c(0.08,-0.08),2),
       y1=0.57, y0=0.5, length = 0.1);
arrows(x1=seq(0.2225,0.7225,length.out=2),
       x0=seq(0.2225,0.7225,length.out=2) + rep(c(0.2,-0.2),2),
       y1=0.39, y0=0.32, length = 0.1);
arrows(x1 = 0.42, x0 = 0.1275, y1 = 0.30, y0 = 0.30, length = 0.1);
arrows(x1 = 0.975, x0 = 0.875, y1 = 1.01, y0 = 0.21, length = 0);

## recombination ancestry, unbalanced haplotypes
g1 <- list(
           makeUniformChromosome("black"),
           makeUniformChromosome("red"),
           makeUniformChromosome("magenta"),
           makeUniformChromosome("blue"),
           makeUniformChromosome("cyan"),
           makeUniformChromosome("green"),
           makeUniformChromosome("yellow"),
           makeUniformChromosome("white"));
g2 <- list(
           recombine(g1[[1]],g1[[1]]),
           recombine(g1[[2]],g1[[2]]),
           recombine(g1[[3]],g1[[3]]),
           recombine(g1[[4]],g1[[4]]),
           recombine(g1[[5]],g1[[5]]),
           recombine(g1[[6]],g1[[6]]),
           recombine(g1[[7]],g1[[7]]),
           recombine(g1[[8]],g1[[8]]));
g3 <- list(
           recombine(g2[[1]],g2[[2]]),
           recombine(g2[[3]],g2[[4]]),
           recombine(g2[[5]],g2[[6]]),
           recombine(g2[[7]],g2[[8]]));
g4 <- list(
           recombine(g3[[1]],g3[[2]]),
           recombine(g3[[3]],g3[[4]]));
g5 <- list(
           recombine(g4[[1]],g4[[2]]));

pdf("recombination_ancestry.pdf",width=8,height=8);
par(mar=c(0,0,0,0));
plot.new();
plot.window(xlim=c(0,1),ylim=c(0.2,1));
for(i in 1:8){
  drawChromosome(g1[[i]],xPos = (i-1)*0.125, yPos = 0.8,
                 scale = 0.4, new = FALSE);
  drawChromosome(g1[[i]],xPos = (i-1)*0.125+0.05, yPos = 0.8,
                 scale = 0.4, new = FALSE);
}
for(i in 1:4){
  drawChromosome(g2[[i*2-1]],xPos = (i-1)*0.25+0.0625, yPos = 0.58,
                 scale = 0.4, new = FALSE);
  drawChromosome(g2[[i*2]],xPos = (i-1)*0.25+0.0625+0.05, yPos = 0.58,
                 scale = 0.4, new = FALSE);
}
for(i in 1:2){
  drawChromosome(g3[[i*2-1]],xPos = (i-1)*0.5+0.1875, yPos = 0.4,
                 scale = 0.4, new = FALSE);
  drawChromosome(g3[[i*2]],xPos = (i-1)*0.5+0.1875+0.05, yPos = 0.4,
                 scale = 0.4, new = FALSE);
}
drawChromosome(g4[[1]],xPos = 0.4375, yPos = 0.22,
               scale = 0.4, new = FALSE);
drawChromosome(g4[[2]],xPos = 0.4375+0.05, yPos = 0.22,
               scale = 0.4, new = FALSE);
drawChromosome(g5[[1]],xPos = 0.8375, yPos = 0.22, scale = 0.4, new = FALSE);
arrows(x0=seq(0.035,0.91,length.out=8),
       x1=seq(0.035,0.91,length.out=8) + rep(c(0.015,-0.015),4),
       y0=0.79, y1=0.68, length = 0.1);
arrows(x0=seq(0.0975,0.8475,length.out=4),
       x1=seq(0.0975,0.8475,length.out=4) + rep(c(0.08,-0.08),2),
       y0=0.57, y1=0.5, length = 0.1);
arrows(x0=seq(0.2225,0.7225,length.out=2),
       x1=seq(0.2225,0.7225,length.out=2) + rep(c(0.2,-0.2),2),
       y0=0.39, y1=0.32, length = 0.1);
arrows(x0 = 0.52, x1 = 0.82, y0 = 0.30, y1 = 0.30, length = 0.1);
dummy <- dev.off();
