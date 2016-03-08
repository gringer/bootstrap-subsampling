gdi.data <- read.csv("/media/Linux_backup/snpchip/structure/results/wtccc/GDI_DISC_RAND001_100T1D_100NBS_level3_all.csv");
pdf("gdi_test_RAND001_100T1D_100NBS_level3_all.pdf",paper="a4r",
    width=11,height=8);
plot(gdi.data$SNPs,gdi.data$GDI, ylab = "Genome Diagnostic Index",
     xlab = "Number of SNPs");
plot(gdi.data$SNPs,gdi.data$GDI - c(0,gdi.data$GDI[-(dim(gdi.data)[1])]),
     ylab = "GDI difference", xlab = "Number of SNPs",
     col = sapply(gdi.data$GDI - c(0,gdi.data$GDI[-(dim(gdi.data)[1])]),
       function(x){
       if(x<0)
         return ("red");
       if(x==0)
         return("yellow");
       if(x>0)
         return("green");
     }));
abline(h=0);
dummy <- dev.off();
