#!/usr/bin/Rscript

## bootstrap.r -- runs a bootstrap process across a genotype file

## Author: David Eccles (gringer), 2009 <programming@gringer.org>

## simplegt-formatted input file
genotypes.inFile = "";
## columns of genotype file to consider cases/controls
## (first value is ignored, so can be used as a group identifier)
casecontrolColumns.inFile = "cases_controls.txt";

## number of subsamples to use for bootstrapping
replicates.proportion = 0.50;
replicates.cases = NULL; # if NULL, will be replicates.proportion or 5 less, whichever is smaller
replicates.controls = replicates.cases;
rep.case.outFile = "";
rep.control.outFile = "";
createReplicates = TRUE; # if false, use existing replicate data
## number of bootstraps to carry out
bootstrap.count = 100;

## file to place bootstrapped results into
bootstraps.outFile = "bootstrap_results.csv";

## should the output values be sorted within each marker?
sortValues <- FALSE;

## controls are on the first line of this file (normally expect cases first)
controlsFirst <- FALSE;

## method to use for calculating results
## current options are Adelta, GTdelta, Chisqmax
processMethod <- "Adelta";

## consider complementary alleles to be the same (allows different
## genotyping methods for the same mutation to be combined)
combineComplementary <- TRUE;

## keep zero counts separate when calculating chi^2 (produces invalid
## results when zero counts are observed for any genotype)
strictChi <- FALSE;

usage <- function(){
  cat("usage: ./bootstrap.r",
      "<case/control column file> [options]\n");
  cat("\nOther Options:\n");
  cat("-help               : Only display this help message\n");
  cat("-input              : File containing genotype data\n");
  cat("-controlfile        : File containing column data for cases/controls\n");
  cat("-repfiles <f1> <f2> : Files containing replicate columns (for repeat experiments)\n");
  cat("-controlsfirst      : case/control file has controls as first line\n");
  cat("-count              : Number of bootstrapts to carry out\n");
  cat("-casereps           : case subpopulation size for boostraps (overrides proportion)\n");
  cat("-controlreps        : control subpopulation size for boostraps (overrides proportion)\n");
  cat("-proportion         : proportion of individuals for boostraps (currently ", replicates.proportion,")\n", sep="");
  cat("-sort               : sort bootstrap results by value\n");
  cat("-strictGT           : Keep complementary alleles separate (don't combine)\n");
  cat("-strictChi          : Respect zero counts in chi^2 table, creates null results\n");
  cat("-output             : output file for results\n");
  cat("-method             : method to use for calculating results\n");
  cat("                      (Adelta, GTdelta, gChisq, Chisqmax, ChisqmaxAll, ShowValues)\n");
  cat("\n");
}

argLoc <- 1;

argLoc <- grep("--args",commandArgs()) + 1; # hack to get around R v2.4
                                            # issue stopping
                                            # commandArgs(TRUE) from
                                            # working
while(!is.na(commandArgs()[argLoc])){
  if(file.exists(commandArgs()[argLoc])){ # file existence check
    casecontrolColumns.inFile <- commandArgs()[argLoc];
  } else {
    if(commandArgs()[argLoc] == "-help"){
      usage();
      quit(save = "no", status=0);
    }
    else if(commandArgs()[argLoc] == "-input"){
      genotypes.inFile <- commandArgs()[argLoc+1];
      argLoc <- argLoc + 1;
    }
    else if(commandArgs()[argLoc] == "-controlfile"){
      casecontrolColumns.inFile <- commandArgs()[argLoc+1];
      argLoc <- argLoc + 1;
    }
    else if(commandArgs()[argLoc] == "-repfiles"){
      rep.case.outFile <- commandArgs()[argLoc+1];
      rep.control.outFile  <- commandArgs()[argLoc+2];
      createReplicates <- FALSE;
      argLoc <- argLoc + 2;
    }
    else if(commandArgs()[argLoc] == "-controlsfirst"){
      controlsFirst <- TRUE;
    }
    else if(commandArgs()[argLoc] == "-sort"){
      sortValues <- TRUE;
    }
    else if(commandArgs()[argLoc] == "-strictGT"){
      combineComplementary <- FALSE;
    }
    else if(commandArgs()[argLoc] == "-strictChi"){
      strictChi <- TRUE;
    }
    else if(commandArgs()[argLoc] == "-proportion"){
      replicates.proportion <- as.numeric(commandArgs()[argLoc+1]);
      cat(file=stderr(), "setting proportion to", replicates.proportion,
          "(",replicates.proportion*100,"% )\n");
      argLoc <- argLoc + 1;
    }
    else if(commandArgs()[argLoc] == "-casereps"){
      replicates.cases <- as.numeric(commandArgs()[argLoc+1]);
      argLoc <- argLoc + 1;
    }
    else if(commandArgs()[argLoc] == "-count"){
      bootstrap.count <- as.numeric(commandArgs()[argLoc+1]);
      argLoc <- argLoc + 1;
    }
    else if(commandArgs()[argLoc] == "-controlreps"){
      replicates.controls <- as.numeric(commandArgs()[argLoc+1]);
      argLoc <- argLoc + 1;
    }
    else if(commandArgs()[argLoc] == "-output"){
      bootstraps.outFile <- commandArgs()[argLoc+1];
      argLoc <- argLoc + 1;
    }
    else if(commandArgs()[argLoc] == "-method"){
      processMethod <- commandArgs()[argLoc+1];
      argLoc <- argLoc + 1;
    }
    else {
      cat("Error: Argument '",commandArgs()[argLoc],
          "' is not understood by this program\n\n", sep="");
      usage();
      quit(save = "no", status=0);
    }
  }
  argLoc <- argLoc + 1;
}

if(!file.exists(casecontrolColumns.inFile)){
  cat("Error: No valid case/control column file given\n\n");
  usage();
  quit(save = "no", status=1);
}

if((bootstraps.outFile != "") && (file.exists(bootstraps.outFile))){
  cat("Error: Output file (", bootstraps.outFile,
      ") exists, please delete it before running this program\n\n",sep="");
  usage();
  quit(save = "no", status=1);
}

if(!createReplicates && (!file.exists(rep.case.outFile) || !file.exists(rep.control.outFile))){
  cat("Error: Replicate files do not exist\n\n");
  usage();
  quit(save = "no", status=1);
}

## carries out a chisquare test of a vector, assuming entries 1..n/2
## are observed, n/2+1..n are expected counts.
vector.chisq <- function(in.vector){
  if(!is.vector(in.vector)){
    stop("Can only be called on a vector");
  }
  if(length(in.vector) %% 2 == 1){
    stop("observed, expected must be same length");
  }
  obs <- in.vector[seq(1,length(in.vector)/2)];
  exp <- in.vector[seq(length(in.vector)/2+1,length(in.vector))];
  ## last section taken (and simplified) from chisq.test code
  ## rescale expected values to same totals as observed values
  E = sum(obs) * exp / sum(exp);
  chisq.values <- (obs-E)^2/E;
  ## !is.nan removes NaN values from result, so cells with zero
  ## !expected counts are ignored
  if(!strictChi){
    return(sum(chisq.values[!is.nan(chisq.values)]));
  }  else {
    return(sum((obs-E)^2/E));
  }
  ## don't warn about counts < 5 -- it will flood the output with
  ## warnings for markers in a whole genome
}

GTCalc <- function(in.genotypes, columns.pop1, columns.pop2, method = "Adelta"){
  num.reps <- dim(columns.pop1)[2];
  ## make everything upper case (simplifies search expressions)
  in.genotypes <- toupper(in.genotypes);
  if(combineComplementary){
    ## substitute complementary alleles
    in.genotypes <- gsub("T","A",in.genotypes);
    in.genotypes <- gsub("G","C",in.genotypes);
    in.genotypes <- gsub("CA","AC",in.genotypes);
  }
  ## calculate major/minor alleles
  gt.counts <- table(in.genotypes);
  if(sum(gt.counts) == 0){ # no individuals
    return(rep(NA, num.reps));
  }
  allele.names <- unique(sub("(A|C|G|T)","",names(gt.counts)));
  homozygous.counts <- sort(gt.counts[paste(allele.names, allele.names, sep="")]);
  if(sum(homozygous.counts) == 0){ # no homozygous individuals
    return(rep(NA, num.reps));
  }
  minor.allele <- substr(names(homozygous.counts[1]),1,1);
  major.allele <- substr(names(homozygous.counts[length(homozygous.counts)]),1,1);
  if(major.allele == minor.allele){ # no minor allele (from homozygous counts)
    het.potentials <- names(gt.counts[grep(major.allele,names(gt.counts))]);
    if(length(het.potentials) == 1){ # no heterozygotes found, so no minor allele
      minor.allele <- "";
#      return(rep(NA, num.reps));
    } else {
      het.name <- het.potentials[het.potentials !=
                                 paste(major.allele, major.allele, sep="")];
      minor.allele <- sub(major.allele, "", het.name);
    }
  }
  ## recode tables as Major/minor (as plink says it *should* be doing)
  ## by substituting alleles for M/m
  ## Note: everything that isn't M/m is considered an invalid
  ## genotype, so if mutation is trimorphic or tetramorphic, then only
  ## the most frequent homozygotes and least frequent homozygotes will
  ## be counted
  if(minor.allele != ""){
    in.genotypes <- gsub(major.allele,"M",gsub(minor.allele,"m",in.genotypes));
  } else {
    in.genotypes <- gsub(major.allele,"M",in.genotypes);
  }
  in.genotypes <- sub("Mm","mM",in.genotypes); # make heterozygotes consistently mM
  in.genotypes <- factor(in.genotypes, levels = c("MM","mM","mm", NA));
  ## generate table based on genotype frequencies
  table1 <- apply(columns.pop1, 2, function(x){table(in.genotypes[x], exclude = NULL)})
  table2 <- apply(columns.pop2, 2, function(x){table(in.genotypes[x], exclude = NULL)})
  ## replace <NA> row name with "XX"
  rownames(table1)[is.na(rownames(table1))] <- "XX";
  rownames(table2)[is.na(rownames(table2))] <- "XX";
  ## return NA if there are no [non-null] genotypes for either population
  if((sum(table1[c("MM","mM","mm"),]) == 0) ||
     (sum(table2[c("MM","mM","mm"),]) == 0)){
    return(rep(NA, num.reps));
  }
  ## at this point, table1/2 will look like the following:
  ##    [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
  ## MM   12   11   10   10   11   12   12   12   11    11
  ## mM   41   42   44   43   42   43   42   43   44    43
  ## mm   32   32   31   32   32   30   31   30   30    31
  ## XX    0    0    0    0    0    0    0    0    0     0
  ## table1: cases, table2: controls
  ## Rows are genotypes, columns are bootstrap number. MM,mM,mm are
  ## the raw genotype counts, XX is a null/bad genotype.
  ## Next step, calculate appropriate statistic
  if(method == "GTdelta") { # genotype freq based delta
    ## remove null values, scale for number of individuals
    table1a <- table1[c("MM","mM","mm"),]/
      colSums(as.matrix(table1[c("MM","mM","mm"),]));
    table2a <- table2[c("MM","mM","mm"),]/
      colSums(as.matrix(table2[c("MM","mM","mm"),]));
    retVal <-
      colSums(abs(table1a - table2a))/2;
  } else if(method == "Adelta") { # allele freq based delta
    ## convert to allele frequencies, scale for number of [non-null]
    ## individuals
    table1a <- (table1["MM",]+table1["mM",]/2) /
      colSums(as.matrix(table1[c("MM","mM","mm"),]));
    table2a <- (table2["MM",]+table2["mM",]/2) /
      colSums(as.matrix(table2[c("MM","mM","mm"),]));
    ## note: if only two alleles, then comparing one allele is the
    ## same as comparing both and dividing by two
    retVal <-
      abs(table1a - table2a);
  } else if((method == "Chisqmax")|| # Max chi square
            (method == "ChisqmaxAll")||
            (method == "gChisq")||
            (method == "ShowValues")) {

    ## calculations for observed, expected derived from plink code
    obs.1 <- colSums(as.matrix(table1[c("mm","mM","MM"),]));
    obs.2 <- colSums(as.matrix(table2[c("mm","mM","MM"),]));
    obs.all <- obs.1+obs.2;
    obs.mm <- (table1["mm",] + table2["mm",]);
    obs.mM <- (table1["mM",] + table2["mM",]);
    obs.MM <- (table1["MM",] + table2["MM",]);
    obs.m <- 2*obs.mm + obs.mM;
    obs.M <- 2*obs.MM + obs.mM;
    exp.1.m <- (obs.1 * obs.m) / obs.all;
    exp.1.M <- (obs.1 * obs.M) / obs.all;
    exp.2.m <- (obs.2 * obs.m) / obs.all;
    exp.2.M <- (obs.2 * obs.M) / obs.all;
    exp.1.mm <- (obs.1 * obs.mm) / obs.all;
    exp.1.mM <- (obs.1 * obs.mM) / obs.all;
    exp.1.MM <- (obs.1 * obs.MM) / obs.all;
    exp.2.mm <- (obs.2 * obs.mm) / obs.all;
    exp.2.mM <- (obs.2 * obs.mM) / obs.all;
    exp.2.MM <- (obs.2 * obs.MM) / obs.all;

    ## calculate genotype statistics for each homozygote vs others

    ## Dominant -DOM.CHISQ (testing major allele genotypes vs other)
    ## i.e. minor allele is dominant
    dom.test <- rbind(table1["mm",]+table1["mM",],
                      table1["MM",],
                      table2["mm",]+table2["mM",],
                      table2["MM",],
                      exp.1.mm+exp.1.mM,
                      exp.1.MM,
                      exp.2.mm+exp.2.mM,
                      exp.2.MM);
                                        #    cat("dom.chisq =",mean(apply(dom.test,2,vector.chisq)));

    ## Recessive -REC.CHISQ (testing minor allele genotypes vs other)
    ## i.e. minor allele is recessive
    rec.test <- rbind(table1["mm",],
                      table1["MM",]+table1["mM",],
                      table2["mm",],
                      table2["MM",]+table2["mM",],
                      exp.1.mm,
                      exp.1.MM+exp.1.mM,
                      exp.2.mm,
                      exp.2.MM+exp.2.mM);
                                        #    cat(", rec.chisq =",mean(apply(rec.test,2,vector.chisq)));

    ## CA trend test -TREND.CHISQ
    ## Note: this assumes allele 1: minor, allele 2: Major
    CA <- ((obs.2/obs.all * table1["mM",]) - (obs.1/obs.all * table2["mM",])) +
      2 * ((obs.2/obs.all * table1["MM",]) - (obs.1/obs.all * table2["MM",]));
    var.CA <- obs.1 * obs.2 *
      (( obs.all * ( obs.mM + 4 * obs.MM )
        - ( obs.mM + 2 * obs.MM )^2 )
       / (obs.all^3));
    CA.chisq <- (CA^2) / var.CA;
                                        #    cat(", CA.chisq =",mean(CA.chisq));
    if ((method == "ChisqmaxAll") || (method == "ShowValues") || (method == "gChisq")) {
      geno.test <- rbind(table1["mm",],
                         table1["mM",],
                         table1["MM",],
                         table2["mm",],
                         table2["mM",],
                         table2["MM",],
                         exp.1.mm,
                         exp.1.mM,
                         exp.1.MM,
                         exp.2.mm,
                         exp.2.mM,
                         exp.2.MM);
                                        #      cat(", geno.chisq =",mean(apply(geno.test,2,vector.chisq)));

      ## Multiplicative / Allelic (actually testing allele quantities of m vs M)
      mult.test <- rbind(2*table1["mm",]+table1["mM",],
                         2*table1["MM",]+table1["mM",],
                         2*table2["mm",]+table2["mM",],
                         2*table2["MM",]+table2["mM",],
                         exp.1.m,
                         exp.1.M,
                         exp.2.m,
                         exp.2.M);
                                        #      cat(", mult.chisq =",mean(apply(mult.test,2,vector.chisq)),"\n");
    }
    ## -1 needed in apply to avoid warnings
    if(method == "Chisqmax"){
      retVal <-
        apply(rbind(
                    CA.chisq,
                    apply(dom.test,2,vector.chisq),
                    apply(rec.test,2,vector.chisq),
                    -Inf),2,max, na.rm = TRUE);
    } else if (method == "gChisq") {
      retVal <-
        apply(geno.test,2,vector.chisq);
    } else if (method == "ChisqmaxAll") {
      retVal <-
        apply(rbind(
                    CA.chisq,
                    apply(geno.test,2,vector.chisq),
                    apply(mult.test,2,vector.chisq),
                    apply(dom.test,2,vector.chisq),
                    apply(rec.test,2,vector.chisq),
                    -Inf),2,max, na.rm = TRUE);
    } else if (method == "ShowValues") {
      retVal <-
        apply(rbind(
                    CA.chisq,
                    apply(geno.test,2,vector.chisq),
                    apply(mult.test,2,vector.chisq),
                    apply(dom.test,2,vector.chisq),
                    apply(rec.test,2,vector.chisq)),2,paste, collapse = ",");
    }
  } else {
    retVal <- NULL;
  }
  ## retVal should be a vector of values, with length equal to the
  ## number of bootstraps
  return(retVal);
}

## matmap -- maps a vector onto a matrix of indexes to the vector
matmap <- function(vector.in, matrix.indices){
  if(max(matrix.indices) > length(vector.in)){
    cat("Error: maximum index greater than vector length\n", file = stderr());
  }
  res <- vector.in[matrix.indices];
  if(is.null(dim(matrix.indices))){
    dim(res) <- c(length(matrix.indices),1);
  } else {
    dim(res) <- dim(matrix.indices);
  }
  return(res);
}


## carries out allele frequency based delta on a single line of simplegt input
processLine <-
  function(marker.name, genotypes, popSamples1, popSamples2, output.file = "",
           append = TRUE, GTmethod = "Adelta", sortbyValue = FALSE){
#  cat(dim(popSamples1), "\n", file=stderr());
  if(length(marker.name) > 1){
    stop("Cannot process more than one line at a time");
  }
  if((length(marker.name) == 0) || (marker.name == "")){
    stop("Empty input: input is not valid");
  }
  bs.results <-
    GTCalc(genotypes,popSamples1,popSamples2,GTmethod);
  ## output 1 line per result
  if(sortValues){
    bs.order <- rev(order(bs.results));
  } else {
    bs.order <- 1:length(bs.results);
  }
  if((GTmethod == "ShowValues") && (append == FALSE)){
    cat("marker,bs.run,CA.chisq,geno.test,mult.test,dom.test,rec.test\n",
        file = output.file, append = FALSE);
    write.table(data.frame(marker = marker.name, bs.run = bs.order,
                           bs.value = bs.results[bs.order]),
                output.file, sep = ",", quote = FALSE, append = TRUE,
                col.names = FALSE, row.names = FALSE);
  } else {
    write.table(data.frame(marker = marker.name, bs.run = bs.order,
                           bs.value = bs.results[bs.order]),
                output.file, sep = ",", quote = FALSE, append = append,
                col.names = !append, row.names = FALSE);
  }
}

casecontrolColumns.con <- file(casecontrolColumns.inFile);
open(casecontrolColumns.con);

if(controlsFirst){
  controls.columns <-
    as.numeric((strsplit(readLines(casecontrolColumns.con, n = 1)," ")[[1]])[-1]);
  cases.columns <-
    as.numeric((strsplit(readLines(casecontrolColumns.con, n = 1)," ")[[1]])[-1]);
} else {
  cases.columns <-
    as.numeric((strsplit(readLines(casecontrolColumns.con, n = 1)," ")[[1]])[-1]);
  controls.columns <-
    as.numeric((strsplit(readLines(casecontrolColumns.con, n = 1)," ")[[1]])[-1]);
}
close(casecontrolColumns.con);

if(is.null(replicates.cases)){
  replicates.cases <- min(length(cases.columns) - 5,
                          trunc(length(cases.columns)*replicates.proportion));
}

if(is.null(replicates.controls)){
  replicates.controls <- min(length(controls.columns) - 5,
                          trunc(length(controls.columns)*replicates.proportion));
}

## generate bootstrap replicates first. This improves speed and
## maintains a consistent population subsample across multiple SNPs.
## note: this is sampling with replacement
## file to place subsample columns into
cases.samples <- NULL;
controls.samples <- NULL;
if(createReplicates){
  rep.case.outFile = paste("caseReplicates",bootstraps.outFile,sep="_");
  rep.control.outFile = paste("controlReplicates",bootstraps.outFile,sep="_");
  cases.samples <- replicate(bootstrap.count,
                             sample(cases.columns,replicates.cases));
  controls.samples <- replicate(bootstrap.count,
                                sample(controls.columns,replicates.controls));
  ## write out replicates to a file, in case validation is needed
  write.table(t(cases.samples),rep.case.outFile, quote = FALSE,
              col.names = FALSE);
  write.table(t(controls.samples),rep.control.outFile, quote = FALSE,
              col.names = FALSE);
} else {
  cases.samples <- t(read.table(rep.case.outFile, row.names = 1));
  controls.samples <- t(read.table(rep.control.outFile, row.names = 1));
  dcc <- c(dim(cases.samples), dim(controls.samples));
  cat(sprintf("Retrieved data for %d bootstraps (cases %d x %d, controls %d x %d)\n",
             dcc[2], dcc[1], dcc[2], dcc[3], dcc[4]), file = stderr());
}

if(genotypes.inFile != ""){
  genotypes.con <- gzfile(genotypes.inFile);
  open(genotypes.con)
} else {
  genotypes.con <- file("stdin");
  open(genotypes.con)
}

input.line <- scan(genotypes.con, what = character(), nlines = 1, quiet = TRUE);
while((length(input.line)>0) && (substr(input.line[1],1,1) == "#")){
  input.line <- scan(genotypes.con, what = character(), nlines = 1, quiet = TRUE);
}

num.indivs <- length(input.line[-1]);
if(num.indivs != (length(cases.columns) + length(controls.columns))){
  cat(sprintf("Warning: number of individuals detected on first line (%d) is not the same as number of cases (%d) + number of controls (%d)\n", num.indivs, length(cases.columns), length(controls.columns)), file = stderr());
}

## process line *without* append (creates header)
processLine(input.line[1], input.line[-1], cases.samples, controls.samples,
            append = FALSE, output.file = bootstraps.outFile,
            GTmethod = processMethod, sortbyValue = sortValues);
input.line <- scan(genotypes.con, what = character(), nlines = 1, quiet = TRUE);
while((length(input.line)>0) && (substr(input.line[1],1,1) == "#")){
  input.line <- scan(genotypes.con, what = character(), nlines = 1, quiet = TRUE);
}
line.number <- 1;
while(length(input.line) > 0){
  ## process line *with* append (no header)
  processLine(input.line[1], input.line[-1], cases.samples, controls.samples,
              append = TRUE, output.file = bootstraps.outFile,
              GTmethod = processMethod, sortbyValue = sortValues);
  line.number <- line.number+1;
  if(line.number %% 1000 == 0){
    cat(".", file = stderr());
  }
  input.line <- scan(genotypes.con, what = character(), nlines = 1, quiet = TRUE);
  while((length(input.line)>0) && (substr(input.line[1],1,1) == "#")){
    input.line <- scan(genotypes.con, what = character(), nlines = 1, quiet = TRUE);
  }
}
if(genotypes.inFile != ""){
  close(genotypes.con);
}

if(line.number > 1){
  cat("(",line.number," lines processed)\n", file = stderr(), sep = "");
} else {
  cat("(",line.number," line processed)\n", file = stderr(), sep = "");
}
