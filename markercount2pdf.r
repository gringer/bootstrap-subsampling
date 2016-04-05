#!/usr/bin/Rscript

## Author: David Eccles (gringer), 2008 <programming@gringer.org>

## markercount2pdf.r -- generates plot of marker counts versus number of bootstraps

usage <- function(){
  cat("usage: ./markercount2pdf.r <file>\n");
}


marker.counts <- read.table(commandArgs(TRUE)[1], row.names = 2);

counts.df <- as.data.frame(table(marker.counts));
counts.df$marker.counts <-
  as.numeric(levels(counts.df$marker.counts))[counts.df$marker.counts];
counts.df$marker.counts <-
  counts.df$marker.counts / max(counts.df$marker.counts);

par(mar=c(5,4.5,0.1,0.1));
pdf("output_marker_counts.pdf", paper = "a4r", width = 11, height = 8);
plot(counts.df$marker.counts, counts.df$Freq, cex.axis = 1.5, cex.lab = 1.5,
     xlab = "Proportion of bootstraps", ylab = "Number of Markers");
dummy <- dev.off();
pdf("output_marker_counts_log.pdf", paper = "a4r", width = 11, height = 8);
plot(counts.df$marker.counts, log(counts.df$Freq)/log(10),
     cex.axis = 1.5, cex.lab = 1.5,
     xlab = "Proportion of bootstraps", ylab = "log10(Number of Markers)");
dummy <- dev.off();
