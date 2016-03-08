#!/usr/bin/Rscript

## icdiff.r -- Calculate information content difference for
## fastphase-formatted genotype file

## Author: David Eccles (gringer), 2008 <programming@gringer.org>

usage <- function(){
  cat("usage: ./icdiff.r","<input file> <output file>\n");
  cat("\nOther Options:\n");
  cat("-help               : Only display this help message\n");
  cat("-hapmap3            : Hapmap phase 3 formatted Files\n");
  cat("\n");
}

infile.name <- FALSE;
outfile.name <- FALSE;
hapmap3 <- FALSE;

argLoc <- 1;
while(!is.na(commandArgs(TRUE)[argLoc])){
  if(file.exists(commandArgs(TRUE)[argLoc])){ # file existence check
    if(infile.name == FALSE){
      infile.name <- commandArgs(TRUE)[argLoc];
    } else {
      cat("Error: More than one existing [input] file specified\n");
      usage();
      quit(save = "no", status=1);
    }
  } else {
    if(commandArgs(TRUE)[argLoc] == "-help"){
      usage();
      quit(save = "no", status=0);
    } else if(commandArgs(TRUE)[argLoc] == "-hapmap3"){
      hapmap3 <- TRUE;
    } else {
      outfile.name <- commandArgs(TRUE)[argLoc];
    }
  }
  argLoc <- argLoc + 1;
}

if((infile.name == FALSE) || (outfile.name == FALSE)){
  cat("Error: No valid input/output file name given\n\n");
  usage();
  quit(save = "no", status=1);
}

## max(table(apply(head(gt.df,33), 2, paste, collapse = ""))) == 1;

## # show heatmap of first 33 genotypes (at start of chromosome 21)
## heatmap(t(head(gt.df,33)), scale = "none", Colv = NA);

## dist.ind <- as.matrix(dist(t(gt.df), method = "manhattan"))
## diag(dist.ind) <- NA;

## ## because of recombination, bits of one parent will be scattered
## ## through each chromosome, so the differences will even out over
## ## large stretches of the chromosome, so differences between each
## ## homologous chromosome and chromosomes from other people will
## ## balance out.
## heatmap(as.matrix(dist.ind), col = topo.colors(100), symm = TRUE)


## information content of some split with counts A1,A2,An:
## let P(An) = An / Sum(A1..An)
## Ic = Sum(-P(An) * log(An,2))
## [Russell, pages 659-660]

## What I want is the difference between the maximum split and the
## actual split. Sum that, and it should give an indication of the
## prevalence of particular haplotypes over others (i.e. extended
## haplotype homozygosity).

## Need to make sure it continues for at least log(n) iterations,
## because until then, we don't know if a haplotype is more prevalent
## than expected in the random case.

## this spreads out in both directions, averaging one choosing to go
## left first and one choosing to go right first

icDiffSum <- function(genotype.matrix, location, constrain = TRUE){
  if(is.null(dim(genotype.matrix))){
    dim(genotype.matrix) <- c(length(genotype.matrix),1);
  }
  ic.target <- log(dim(genotype.matrix)[1],2);
  num.haplotypes <- dim(genotype.matrix)[1];
  ic.diffs.left <- NULL;
  lengths.left <- NULL;
  ic.diffs.right <- NULL;
  lengths.right <- NULL;
  for(leftBias in c(TRUE, FALSE)){
    goLeft <- leftBias;
    diffpos <- 1;
    ic.value <- 0;
    loc.left <- location;
    loc.right <- location;
    while(ic.value < ic.target){
      if(constrain){
        ## constrain information content to maximum possible with
        ## given number of haplotypes
        ic.max <- min(((loc.right - loc.left) + 1), ic.target);
      } else {
        ic.max <- (loc.right - loc.left) + 1;
      }
      ## warning: null genotypes appear as genotype "NA"
      if(loc.right != loc.left){
        prob.hap <- table(apply(genotype.matrix[,loc.left:loc.right],
                                1, paste, collapse = "")) / num.haplotypes;
      } else {
        prob.hap <- table(genotype.matrix[,loc.left]) / num.haplotypes;
      }
      ic.value <- sum(-prob.hap * log(prob.hap,2));
      ##print(paste(leftBias, " -> ", ic.value));
      if(leftBias == TRUE){
        ic.diffs.left[diffpos] <- (ic.max - ic.value) / ic.max;
        lengths.left[diffpos] <- (loc.right-loc.left) + 1;
      } else {
        ic.diffs.right[diffpos] <- (ic.max - ic.value) / ic.max;
        lengths.right[diffpos] <- (loc.right-loc.left) + 1;
      }
      diffpos <- diffpos + 1;
      ## Extra logic to make sure it doesn't duplicate results by
      ## hanging in the same place for two iterations. These two
      ## 'goLeft' cases could be combined, but it makes the code less
      ## clear.
      if(goLeft){
        goLeft <- !goLeft; # i.e. goLeft <- FALSE
        if(loc.left > 1){
          loc.left <- max(1,loc.left - 1);
        } else {
          loc.right <- min(dim(genotype.matrix)[2],loc.right + 1);
        }
      } else {
        goLeft <- !goLeft; # i.e. goLeft <- TRUE
        if(loc.right < dim(genotype.matrix)[2]){
          loc.right <- min(dim(genotype.matrix)[2],loc.right + 1);
        } else {
          loc.left <- max(1,loc.left - 1);
        }
      }
      if((loc.left == 1) && (loc.right == dim(genotype.matrix)[2])){
        ## if the entire chromosome is covered, do no more
        ic.value <- ic.target;
      }
    }
  }
  ## print(rbind(round(lengths.left,5),round(ic.diffs.left,5)));
  ## print(rbind(round(lengths.right,5),round(ic.diffs.right,5)));
  return((sum(ic.diffs.left) + sum(ic.diffs.right))/2);
}

oldicDiffSum <- function(genotype.matrix){
  ic.diffs <- NULL;
  ic.diff <- 1;
  ic.value <- 0;
  x <- 1;
  if(is.null(dim(genotype.matrix))){
    dim(genotype.matrix) <- c(1,length(genotype.matrix));
  }
  ic.target <- log(dim(genotype.matrix)[2],2);
  while((ic.value < ic.target) && (x <= dim(genotype.matrix)[1])){
    ic.max <- min(x,ic.target);
    counts.hap <- table(apply(head(genotype.matrix,x), 2, paste, collapse = ""));
    prob.hap <- counts.hap / sum(counts.hap);
    ic.value <- sum(-prob.hap * log(prob.hap,2));
    ic.diff <- ((ic.max - ic.value) / ic.max);
    ic.diffs <- c(ic.diffs,ic.diff);
    x <- x+1;
  }
  return(sum(ic.diffs));
}

gtval <- function(x){
  if((is.matrix(x)) || (is.data.frame(x))){
    apply(x,2,function(data){
      as.numeric(sub("A|T",1,sub("G|C",2,data)));
    });
  } else if(is.vector(x)) {
    as.numeric(sub("A|T",1,sub("G|C",2,x)));
  } else {
    if((x = "A") || (x = "T")){
      1;
    } else if((x = "C") || (x = "G")){
      2;
    } else {
      0;
    }
  }
}

## ## The following test suggests that this algorithm is tractable and
## ## shouldn't need further optimisations:

## ## 1) It is linear in complexity for number of markers [O(n)]

## system.time(for(x in 1:10){oldicDiffSum(gt.df[-(1:x),])})
## ##   user  system elapsed
## ##  2.120   0.028   2.267
## system.time(for(x in 1:100){oldicDiffSum(gt.df[-(1:x),])})
## ##   user  system elapsed
## ## 17.393   0.240  17.688

## ## 2) Estimated runtime for chromosome 21 is less than half an
## ##    hour (actually, 16 minutes, 5435 markers)

## ## > (5435/100 * 17.393) / 60
## ## [1] 15.75516

## ## [Estimated runtime for chromosome 1 is about 1.5 hours, 23062
## ## markers]


cat("Retrieving data for chromosome ...");
inFile <- gzfile(infile.name);
open(inFile);
lineRead <- "";
if(!hapmap3){
  while(lineRead != "BEGINGENOTYPES"){ # read until genotype start
    lineRead <- paste(scan(inFile, nlines = 1,
                           what = character(0),
                           quiet = TRUE)[1:2], collapse = "");
  }
  lineRead <- "";
  person.id <- "";
  gt.df <- NULL;
  while(lineRead != "END GENOTYPES"){ # read until genotype end
    lineRead <- readLines(inFile, n = 1, warn = FALSE);
    if(substr(lineRead,3,4) == "id"){
      person.id <- substr(lineRead,6,nchar(lineRead));
      hap.1 <- scan(inFile, nlines = 1, quiet = TRUE);
      hap.2 <- scan(inFile, nlines = 1, quiet = TRUE);
      new.df <- cbind(hap.1, hap.2);
      colnames(new.df) <- paste(person.id,c(1,2), sep = ".");
      gt.df <- cbind(gt.df,new.df);
    }
  }
  colnames(gt.df) <- sub("0+([0-9][0-9]\\.)","\\1",
                         colnames(gt.df), perl=TRUE);
  gt.df <- t(gt.df);
} else { # if hapmap3-formatted input (but no header)
  gt.df <- read.table(inFile, row.names = 1, header = FALSE);
  gt.df <- gt.df[,-1];
  marker.names <- rownames(gt.df);
  gt.df <- data.frame(t(gt.df));
}
close(inFile);
cat("done\n");
cat("Calculating distance matrix to eliminate identical chromosomes...");
if(hapmap3){
  dist.ind <- as.matrix(dist(gtval(gt.df), method = "manhattan"));
} else {
  dist.ind <- as.matrix(dist(gt.df), method = "manhattan");
}
remove.columns <- NULL;
for(haplotype in 1:(dim(gt.df)[1])){
  ## if fewer than 20 differences, treat as the same chromosome
  ## this allows for some genotyping error
  if(length(which(which(dist.ind[haplotype,] < 20) > haplotype) > 0)){
    remove.columns <- c(remove.columns, (which(dist.ind[haplotype,] < 20))
                        [which(dist.ind[haplotype,] < 20) > haplotype]);
  }
}
if(length(remove.columns) > 0){
  gt.df <- gt.df[,-remove.columns];
}
cat("done (",length(unique(remove.columns))," columns removed)\n",sep="");
outFile <- file(outfile.name);
open(outFile, open = "w");
for(marker in 1:(dim(gt.df)[2])){
  cat(colnames(gt.df)[marker]," ",icDiffSum(gt.df, marker),
      "\n", file = outFile, sep = "");
}
close(outFile);


