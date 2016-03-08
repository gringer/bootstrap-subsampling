#!/usr/bin/Rscript

## agrf2simplegt.r -- Convert a three-column Excel file from AGRF into
## a simplegenotype file (similar to those used in the HapMap samples)

## Author: David Eccles (gringer), 2007 <programming@gringer.org>

usage <- function(){
  cat("usage: ./agrf2simplegt.r <input file>\n");
  cat("\n");
}

library(gdata)

if(is.na(commandArgs(TRUE)[1])){
  cat("ERROR: no input file given\n\n");
  usage();
  quit(save = "no", status = 1);
} else {
  xl.filename = commandArgs(TRUE)[1];
}

# Note: the read.xls function prints the file name to STDOUT
xldata.df <- read.xls(xl.filename, sheet = 2);


# Convert lower to upper, Call -> Sample_ID
# [change in AGRF result format 2008-MAY-05]
colnames(xldata.df) <- sub("\\.","_",toupper(colnames(xldata.df)));
colnames(xldata.df) <- sub("^CALL$","GENOTYPE_ID",colnames(xldata.df));

xldata.df$SAMPLE_ID <- factor(xldata.df$SAMPLE_ID);

## temp variable to improve readability
gid <- xldata.df$GENOTYPE_ID;
## replace "Fail" with "NN", like HapMap
levels(gid)[which(levels(gid) == "Fail")] <- "NN";
## reorder alleles
levels(gid)[which(levels(gid) == "CA")] <- "AC";
levels(gid)[which(levels(gid) == "GA")] <- "AG";
levels(gid)[which(levels(gid) == "GC")] <- "CG";
levels(gid)[which(levels(gid) == "TC")] <- "CT";
levels(gid)[which(levels(gid) == "TG")] <- "GT";
levels(gid)[which(levels(gid) == "TA")] <- "AT";
## convert haploid to diploid
levels(gid)[which(levels(gid) == "A")] <- "AA";
levels(gid)[which(levels(gid) == "C")] <- "CC";
levels(gid)[which(levels(gid) == "G")] <- "GG";
levels(gid)[which(levels(gid) == "T")] <- "TT";
## put back into the dataframe
xldata.df$GENOTYPE_ID <- gid;

xldata.wide <- reshape(xldata.df, idvar="ASSAY_ID",
                     timevar="SAMPLE_ID", direction="wide");
colnames(xldata.wide) <- sub("GENOTYPE_ID.","",colnames(xldata.wide));
## rename rows to markers
rownames(xldata.wide) <- xldata.wide$ASSAY_ID;
## remove (now) redundant marker column
xldata.wide <- xldata.wide[,-1];

## output individual IDs
cat("## <Individual/Column IDs: ",colnames(xldata.wide)," > ##\n")

## output result
if(is.null(
           apply(rbind(sprintf("%-12s",rownames(xldata.wide)),
                       apply(xldata.wide,
                             1,paste)),
                 2,cat, "\n"))){
  cat(""); ## stops 'NULL' from being output on last line
}
