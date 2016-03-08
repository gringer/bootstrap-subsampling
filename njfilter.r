#!/usr/bin/Rscript

## njfilter -- iteratively adds SNPs that don't have a genotype
## profile similar to SNPs already added to a set of informative SNPs

## The most informative SNP is always selected, and each additional
## considered SNP will have the highest delta. So, hopefully,
## maximally informative SNPs that are not similar in profile to
## already selected SNPs will be selected. Cutoff values are defined
## as a probability based on a normal distribution of pairwise
## similarity values (a good first-guess at this value would be 0.5).
## Lower values correspond to a set with many SNPs (lower cutoff...
## the probability of SNPs being considered close is low), while
## higher values will result in a set with fewer SNPs.

## Author: David Eccles (gringer), 2009 <programming@gringer.org>

usage <- function(){
  cat("usage: ./njfilter.r <Genotype file> -s <split point> [<SNP ranking file>] [options]\n");
  cat("\nOther Options:\n");
  cat("  -r   probability resolution (default 0.01)\n");
  cat("  -min minimum probability (default 0)\n");
  cat("  -max maximum probability (default 1)\n");
  }

infile.name <- FALSE;
SNPrankfile.name <- FALSE;
splitPoint <- FALSE;
probResolution <- 0.01;
minProb <- 0;
maxProb <- 1;

argLoc <- 1;
while(!is.na(commandArgs(TRUE)[argLoc])){
  if(file.exists(commandArgs(TRUE)[argLoc])){ # file existence check
    if(infile.name == FALSE){
      infile.name <- commandArgs(TRUE)[argLoc];
    } else{
      if(SNPrankfile.name == FALSE){
        SNPrankfile.name <- commandArgs(TRUE)[argLoc];
      }
      else{
        cat("Error: More than two input files specified\n");
        usage();
        quit(save = "no", status=1);
      }
    }
  }
  if(commandArgs(TRUE)[argLoc] == "-s"){
      splitPoint <- as.numeric(commandArgs(TRUE)[argLoc+1]);
      argLoc <- argLoc+1;
  }
  if(commandArgs(TRUE)[argLoc] == "-r"){
      probResolution <- as.numeric(commandArgs(TRUE)[argLoc+1]);
      argLoc <- argLoc+1;
  }
  if(commandArgs(TRUE)[argLoc] == "-min"){
      minProb <- as.numeric(commandArgs(TRUE)[argLoc+1]);
      argLoc <- argLoc+1;
  }
  if(commandArgs(TRUE)[argLoc] == "-max"){
      maxProb <- as.numeric(commandArgs(TRUE)[argLoc+1]);
      argLoc <- argLoc+1;
  }
  if(length(grep("-h",commandArgs(TRUE)[argLoc])) > 0){
    usage();
    quit(save = "no", status = 0);
  }
  argLoc <- argLoc+1;
}

if(infile.name == FALSE){
  cat("Error: No input file specified\n");
  usage();
  quit(save = "no", status=1);
}

if(splitPoint == FALSE){
  cat("Error: No split point (-s <value>) specified\n");
  usage();
  quit(save = "no", status=2);
}

gt2Num <- function(x){
  if(is.data.frame(x)){
    tmpRownames <- rownames(x);
    out.df <- data.frame(lapply(x,gt2Num));
    rownames(out.df) <- tmpRownames;
    return(out.df);
  } else if(length(x)>1){
    return(mapply(gt2Num,x));
  }
  x <- toupper(x);
  if((x == "AT") | (x == "TA")){
#    cat(file="STDERR","Warning: complementary mutation found");
    return(NA);
  }
  if((x == "CG") | (x == "GC")){
#    cat(file="STDERR","Warning: complementary mutation found");
    return(NA);
  }
  if((x == "AA") | (x == "TT")){
    return(0);
  } else if ((x == "AC") | (x == "CA") |
             (x == "AG") | (x == "GA") |
             (x == "TC") | (x == "CT") |
             (x == "TG") | (x == "GT")){
    return(1);
  } else if ((x == "CC") | (x == "GG")){
    return(2);
  } else {
    return(NA);
  }
}

## Determines pairwise differences between each SNP (or each
## individual if the matrix is rotated...). The value returned in each
## cell is the number of allelic differences across all individuals.
## For a simple implementation of this (which occurs when
## checkComplement = FALSE), there's an assumption that the A/T allele
## of one SNP is linked to the A/T allele of the other SNP, which is
## not always correct. To highlight this, consider a case where one
## SNP is homozygous AA for all individuals, and another SNP is
## homozygous CC for all individuals. The simple implementation will
## consider these two SNPs to have maximal difference (i.e. a value of
## 200 for 100 individuals), where they should really be treated as
## having no difference. A better comparison (which occurs when
## checkComplement = TRUE) involves choosing the minimum value of two
## comparisons, one assuming A/T = A/T, and the other assuming A/T =
## C/G. This minimum check is probably a bit slower, so the option to
## turn it off has been provided. Also, if you're comparing
## indivduals, you'd *want* the alleles to be consistent, because in
## that case the compared SNPs are the same.

pairwiseDiff <- function(in.mat, compareRows, na.rm = TRUE,
                         checkComplement = TRUE){
  matSize <- length(compareRows);
  if(checkComplement){
    inv.mat <- 2 - in.mat;
    res.mat <- matrix(NA,matSize,matSize*2);
    dim(res.mat) <- c(matSize,matSize,2);
  } else {
    res.mat <- matrix(NA,matSize,matSize);
    dim(res.mat) <- c(matSize,matSize,1);
  }
  ## For loops are *slow* in R... it would be nice to be able to
  ## vectorise this completely.
  for(x in 1:matSize){
    if(x %% 10 == 0){
      cat(".");
    }
    res.mat[x,,1] <- colSums(abs((in.mat[compareRows[x],] %*%
                                  t(rep(1,matSize))) -
                                 t(in.mat[compareRows,])),
                             na.rm = TRUE);
    if(checkComplement){
      res.mat[x,,2] <- colSums(abs((in.mat[compareRows[x],] %*%
                                    t(rep(1,matSize))) -
                                   t(inv.mat[compareRows,])),
                               na.rm = TRUE);
    }
  }
  if(checkComplement){
    res.mat <- apply(res.mat,c(1,2),min); # calculate minimum difference
  } else {
    res.mat <- res.mat[,,1];
  }
  rownames(res.mat) <- compareRows;
  colnames(res.mat) <- compareRows;
  return(res.mat);
}

cat("Reading in file...");
con <- file(infile.name);
open(con);
gt.df <- read.table(con, colClasses = "character", row.names = 1);
close(con);
cat("done!\n");

## takes about 50s on assimilis
cat("Converting genotype data into matrix...");
gt.mat <- as.matrix(gt2Num(gt.df));
cat("done!\n");

gt.deltaStats <-
  data.frame(
             T1D = rowSums(gt.mat[,1:splitPoint], na.rm = TRUE) /
             (splitPoint * 2),
             NBS = rowSums(gt.mat[,-(1:splitPoint)], na.rm = TRUE) /
             ((dim(gt.mat)[2] - splitPoint) * 2),
             T1DAdj = rowSums(gt.mat[,1:splitPoint], na.rm = TRUE) /
             (apply(!is.na(gt.mat[,1:splitPoint]),1,
                   function(x){length(which(x))}) * 2),
             NBSAdj = rowSums(gt.mat[,-(1:splitPoint)], na.rm = TRUE) /
             (apply(!is.na(gt.mat[,-(1:splitPoint)]),1,
                   function(x){length(which(x))}) * 2));
gt.deltaStats$Delta <- abs(gt.deltaStats$T1DAdj - gt.deltaStats$NBSAdj);

if(!(SNPrankfile.name == FALSE)){
  cat("Reading in SNP ranking file...");
  gt.deltaStats <- read.table(SNPrankfile.name, row.names = 1);
  if(is.null(colnames(gt.deltaStats))){
    colnames(gt.deltaStats) <- 1:(dim(gt.deltaStats)[2]);
  }
  colnames(gt.deltaStats)[dim(gt.deltaStats)[2]] <- "Delta";
  cat("done!\n");
}


## takes about 2min on assimilis
cat("Calculating pairwise SNP vs SNP matrix...");
pairwise.mat <- pairwiseDiff(gt.mat,rownames(gt.mat));
cat("done!\n");

## njFilterSNP -- removes SNPs that are close neighbours to SNPs
## already in the selection set. The most informative SNP is always
## selected, and each additional considered SNP will have the highest
## delta. So, hopefully, maximally informative SNPs that are not
## similar in profile to already selected SNPs will be selected. The
## cutoff here refers to a probability based on a normal distribution
## of pairwise similarity values (a good first-guess at this value
## would be 0.5). Lower values include more SNPs (lower cutoff), while
## higher values include fewer SNPs.
njFilterSNP <- function(cutoff){
  cutoff <- qnorm(cutoff,
                  mean = mean(pairwise.mat),
                  sd = sd(as.vector(pairwise.mat)));
  deltaSortedSNPs <- rownames(gt.deltaStats[rev(order(gt.deltaStats$Delta)),
                                            c("Delta","Delta")]);
  limit <- length(deltaSortedSNPs);
  SNPset <- deltaSortedSNPs[1];
  if(limit>1){
    for(x in 2:limit){
      testSNP <- deltaSortedSNPs[x];
      if(min(pairwise.mat[testSNP,SNPset]) > cutoff){
        SNPset <- c(SNPset,testSNP);
      }
    }
  }
  return(SNPset);
}

setSize <- NULL;
for(x in seq(minProb,maxProb,by=probResolution)){
  cat("Number of SNPs at p>",x,": ",sep="");
  cat(length(njFilterSNP(x)));
  cat("\n",sep="");
  cat("SNPs at p>",x,": ",sep="");
  cat(njFilterSNP(x));
  cat("\n",sep="");
}
