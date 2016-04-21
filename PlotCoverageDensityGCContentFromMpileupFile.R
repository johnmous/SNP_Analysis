### GC content in relation to Depth 

fileA="/zfs/datastore0/group_root/MAD-RBAB/03_RBAB-internal/MAD1000-P012-variant_calling/MAD1000-P012-E002_2014_Filtering_imousta1/Scratch/compare_mapping/preparePresentationExperiments/YFPauregWT01REfCov.mpileup"
dataA <- read.delim(fileA, header=F)
depthA <- dataA$V2
GCContentA <- vector()
for (i in seq(0, 796, by=4)) {
  bases <- dataA[ depthA>i & depthA<i+4 , 1 ]
  totalBases <- length(bases)
  table <- table(bases)
  GCsum <- table[names(table)=="G"] + table[names(table)=="C"]
  GCContentA <- c(GCContentA, GCsum/totalBases)
}

fileB="/zfs/datastore0/group_root/MAD-RBAB/03_RBAB-internal/MAD1000-P012-variant_calling/MAD1000-P012-E002_2014_Filtering_imousta1/Scratch/compare_mapping/preparePresentationExperiments/NdEcoliBC01RefCov.mpileup"
dataB <- read.delim(fileB, header=F)
depthB <- dataB$V2
GCContentB <- vector()
for (i in seq(0, 396, by=4)) {
  bases <- dataB[ depthB>i & depthB<i+4 , 1 ]
  totalBases <- length(bases)
  table <- table(bases)
  GCsum <- table[names(table)=="G"] + table[names(table)=="C"]
  GCContentB <- c(GCContentB, GCsum/totalBases)
}

par(mar = c(5, 4, 4, 4) + 0.3) 
plot(density(depthA [depthA<600]), col= "cyan3", xlim = c(0, 500), ylim = c(0, 0.015),  xlab = "Coverage", ylab="Fraction of Genome", 
     main="Fraction of Genome/GC Content vs Coverage", lwd=3)
par(new=T)
plot(density(depthB), col= "red", xlim = c(0, 500), ylim = c(0, 0.015), xlab = "", ylab="", main="", lwd = 4)
par(new=T)
plot(seq(0, 596, by=4), GCContentA[1:150], axes=F ,  bty = "n", xlab = "", ylab = "", 
     col = "cyan3", xlim = c(0,500), ylim= c(0,1), pch=20)
par(new=T)
plot(seq(0, 396, by=4), GCContentB, axes=F ,  bty = "n", xlab = "", ylab = "", 
     col = "red",xlim = c(0,500), ylim= c(0,1), pch=20)

axis(side=4)
mtext("G+C content", side=4, line=3)
grid(lwd=4)
abline(h=0.5, lty= 2)
legend(x="topright", legend = c("EColi", "PAureginosa"), col = c("red", "cyan3" ), pch =19, cex=1.2)


