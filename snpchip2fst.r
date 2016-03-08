#!/usr/bin/Rscript

#usage ./snpchip2fst.r <file> (p<name> <range from> <range to>)+ [options]

## Author: David Eccles (gringer), 2008 <programming@gringer.org>

usage <- function(){
  cat("usage: ./snpchip2fst.r",
      "<file> (p<name> <range from> <range to>)+ [options]\n");
}

argLoc <- 1;
popNames <- character(0);
popLimits <- NULL;
infile <- FALSE;
infile.df <- FALSE;
flip <- FALSE;
doSVG <- FALSE;
sortByQ <- FALSE;
labelAxis <- FALSE;
rotateLabels <- FALSE;
halfHeight <- FALSE;
tailOnly <- TRUE;
pointSize <- 1;
numIndivs <- 0;
linesIn <- 0;

while(!is.na(commandArgs(TRUE)[argLoc])){
  if(file.exists(commandArgs(TRUE)[argLoc])){ # file existence check
    if(infile == FALSE){
      numIndivs <- length(scan(commandArgs(TRUE)[argLoc], nlines = 1, what = character(0), quiet = TRUE)) - 1;
      infile <- file(commandArgs(TRUE)[argLoc]);
    } else{
      cat("Error: More than one input file specified\n");
      usage();
      quit(save = "no", status=1);
    }
  } else {
    if(substr(commandArgs(TRUE)[argLoc],1,1) == "p"){
      inName <- substring(commandArgs(TRUE)[argLoc],2);
      popNames <- append(popNames, inName);
      popLimits <- cbind(popLimits,c(as.numeric(commandArgs(TRUE)[argLoc+1]),
                                     as.numeric(commandArgs(TRUE)[argLoc+2])));
      cat("Adding population: ",inName,"(",
          commandArgs(TRUE)[argLoc+1],"-",
          commandArgs(TRUE)[argLoc+2],", ",
          as.numeric(commandArgs(TRUE)[argLoc+2]) -
          as.numeric(commandArgs(TRUE)[argLoc+1])+1," individuals)\n",sep="");
      argLoc <- argLoc + 2;
    }
    if(commandArgs(TRUE)[argLoc] == "-lines"){
      linesIn <- as.numeric(commandArgs(TRUE)[argLoc+1]);
      cat("Setting line read limit to ",linesIn,"\n",sep="");
      argLoc <- argLoc + 1;
    }
    if(commandArgs(TRUE)[argLoc] == "-alldata"){
      tailOnly <- FALSE;
      cat("Displaying all data\n");
      argLoc <- argLoc + 1;
    }
  }
  argLoc <- argLoc + 1;
}

cat("Input file has data for ", numIndivs, " individuals\n", sep="");

if(numIndivs == 0){
  numIndivs <- sum((popLimits[2,] - popLimits[1,]) + 1);
}

if(length(popNames) == 0){
  cat("Warning: no populations defined.\n",
      "         All individuals will be treated as the same population\n"
      , sep="");
  popNames <- "All";
  popLimits <- rbind(1,numIndivs);
} else {
  totLength <- 0;
  for(x in 1:dim(popLimits)[2]){
    totLength <- totLength + (popLimits[2,x] - popLimits[1,x] + 1);
    if((popLimits[2,x] > numIndivs) ||
        (popLimits[1,x] < 1)){
      cat("Error: range of population", popNames[x] ,
          " outside range of individuals",
          "in the file\n");
      usage();
      quit(save = "no", status = 2);
    }
  }
  if(totLength < numIndivs){
    cat("Warning: fewer individuals specified than present in the file\n");
  }
  if(totLength > numIndivs){
    cat("Warning: more individuals specified than present in the file\n");
  }
}
if(!infile){
  infile = "";
}

cat("reading file...");
a <- proc.time();
infile.list <- scan(infile, what = as.list(c("marker",character(numIndivs))), nlines = linesIn, quiet = TRUE);
cat("done\n");
##print(proc.time() - a);
infile.rownames <- infile.list[[1]];
cat("substituting...");
a <- proc.time();
infile.df <-
  data.frame(sapply(infile.list[-1],
         function(x){
           x <- gsub("(\\w)","\\L\\1", x, perl = TRUE); # convert to lower case
           x <- gsub("t","a", x, fixed = TRUE);
           x <- gsub("g","c", x, fixed = TRUE);
           x <- gsub("ca","ac", x, fixed = TRUE);
         }));
rownames(infile.df) <- infile.rownames;
cat("done\n");
##print(proc.time() - a);
numMarkers <- length(infile.rownames);

## Begin Fst algorithm, taken from Weir & Cockerham (1984)

fst.n <- (popLimits[2,] - popLimits[1,]) + 1; # population sizes

pop.hets <- matrix(NA,nrow = numMarkers,ncol = dim(popLimits)[2]);
pop.afreqs <- matrix(NA,nrow = numMarkers,ncol = dim(popLimits)[2]);
rownames(pop.hets) <- infile.rownames;
rownames(pop.afreqs) <- infile.rownames;
cat("calculating allele statistics...");
for(x in 1:dim(popLimits)[2]){
  cat(x,"/het...",sep="");
  pop.hets[,x] <- (rowSums(1 * (infile.df[,popLimits[1,x]:popLimits[2,x]] == "ac"))
        / (fst.n[x]));
  cat(x,"/frq...",sep="");
  pop.afreqs[,x] <- (rowSums(2 * (infile.df[,popLimits[1,x]:popLimits[2,x]] == "aa") +
                             1 * (infile.df[,popLimits[1,x]:popLimits[2,x]] == "ac"))
                     / (fst.n[x] * 2));
}
cat("done\n");
##print(proc.time() - a);

fst.r <- dim(popLimits)[2]; # number of populations
fst.nhat <- mean(fst.n); # average sample size
fst.nc <- ((fst.r * fst.nhat) - sum((fst.n^2) / (fst.r*fst.nhat))) /
          (fst.r - 1) # squared coefficient of variation of sample sizes
fst.phat <- colSums((fst.n * t(pop.afreqs)) /
                    (fst.r * fst.nhat)); # average sample frequency of allele A
fst.s2 <- colSums((fst.n*t(pop.afreqs - fst.phat)^2) /
                  ((fst.r-1)*fst.nhat)); # sample variance of allele A
                                         # frequencies over populations
fst.hhat <- colSums((fst.n * t(pop.hets)) /
             (fst.r * fst.nhat)); # average heterozygote frequency for allele A

#print(fst.n);
#print(fst.r);
#print(head(pop.hets));
#print(head(pop.afreqs));
#print(head(fst.phat));
#print(head(fst.s2));
#print(head(fst.hhat));


fst.a = (fst.nhat/fst.nc) * (fst.s2 - 1/(fst.nhat-1) * (fst.phat*(1-fst.phat)
            - ((fst.r-1)/fst.r)*fst.s2 - 1/4*fst.hhat));
fst.b <- fst.nhat/(fst.nhat-1) *
  (fst.phat*(1-fst.phat) - (fst.r-1)/fst.r*fst.s2
   - (2*fst.nhat-1)/(4*fst.nhat) * fst.hhat);
fst.c <- 1/2 * fst.hhat;

fst.fst <- fst.a / (fst.a + fst.b + fst.c);
pop.delta <- abs(pop.afreqs[,1] - pop.afreqs[,2]);

fst.mean <- data.frame(Marker = "mean.fst", fst = mean(fst.fst, na.rm = TRUE), delta = mean(pop.delta), meanhet = mean(fst.hhat));
rownames(fst.mean) <- "mean.fst";
fst.sd <- data.frame(Marker = "sd.fst", fst = sd(fst.fst, na.rm = TRUE), delta = sd(pop.delta), meanhet = sd(fst.hhat));
rownames(fst.sd) <- "sd.fst";
fst.df <- rbind(data.frame(Marker = infile.rownames, fst = fst.fst, delta = pop.delta, meanhet = fst.hhat), fst.mean, fst.sd);

library(gdata);

#fst.df[,2] <- round(fst.df[,2],3);
#fst.df[,3] <- round(fst.df[,3],3);
#fst.df[,4] <- round(fst.df[,4],3);

##write.fwf(fst.df, quote = FALSE, rownames = TRUE);
if(tailOnly){
  write.fwf(tail(fst.df), quote = FALSE);
  print(cor.test(fst.fst,pop.delta));
} else {
  write.fwf(fst.df, quote = FALSE);
}
