#!/usr/bin/Rscript

## Author: David Eccles (gringer), 2009 <programming@gringer.org>

if(length(commandArgs(TRUE)) < 1){
  cat("Error: No file specified on the command line");
  quit(save="no",status = 1);
}
a <- as.matrix(read.table(commandArgs(TRUE)[1], stringsAsFactors = FALSE,
                          row.names = 1));

a <- gsub(a,pattern = "G", replacement = "C", ignore.case = TRUE);
a <- gsub(a,pattern = "T", replacement = "A", ignore.case = TRUE);
a <- gsub(a,pattern = "CA", replacement = "AC", ignore.case = TRUE);
a <- gsub(a,pattern = "AA", replacement = "0", ignore.case = TRUE);
a <- gsub(a,pattern = "AC", replacement = "1", ignore.case = TRUE);
a <- gsub(a,pattern = "CC", replacement = "2", ignore.case = TRUE);
a <- gsub(a,pattern = "..", replacement = NA);
rnames <- rownames(a);
a <- matrix(as.numeric(a),dim(a)[1],dim(a)[2]);
rownames(a) <- rnames;

z <- dist(a);
write.table(sort(colMeans(as.matrix(dist(a)), na.rm = TRUE),
                 decreasing = TRUE), quote = FALSE, col.names = FALSE);
