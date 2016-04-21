# title: Graphs of collapsed fastq counts
# Author: Ioannis Moustakas, i.moustakas@uva.nl

numberOfReads <- vector()
label <- vector()
countsOfReadsInGroup <- list()
cumsumCountsOfReadsInGroup <- list()

for (i in 1:4){
path="/zfs/datastore0/group_root/MAD-RBAB/03_RBAB-internal/MAD1000-P003-ZF_egg_RNA_inventory/MAD1000-P003-E005_2014_sRNA_gel_isolation_mlocati1/Scratch/ioannis/"
data <- read.delim(paste0(path, "collapsedCount00", i ,".txt"), header=F)
vector <- data[[1]]

numberOfReads <- c(numberOfReads, sum(vector))
countsAboveFive <- vector[vector>5]
countsOfReadsInGroup[[i]] <- countsAboveFive
cumsumCountsAboveFive <- cumsum(countsAboveFive)
cumsumCountsOfReadsInGroup[[i]] <-  cumsumCountsAboveFive
label <- c(label, paste0("BC", i))
}
colors <- c("red", "blue", "yellow", "green")

png(paste0(path, "plots.png"),width = 800, height = 800)
plot(1:length(countsOfReadsInGroup[[1]]), countsOfReadsInGroup[[1]], log="y", 
     ylab = "Number of Reads in Group (log Scale)", xlab = "Number of Group", col=colors[1], cex=.5)
for (i in 2:4) {
  points(1:length(countsOfReadsInGroup[[i]]), countsOfReadsInGroup[[i]], col=colors[i], cex=.5)
}
legend("bottomright", legend= label, fill=colors)
grid(lwd = 2, lty = 1)
dev.off()

png(paste0(path, "plots1.png"),width = 800, height = 800)
plot(1:length(cumsumCountsOfReadsInGroup[[1]]), cumsumCountsOfReadsInGroup[[1]], 
     ylab = "Cumulative Number of Reads in Groups", xlab = "Number of Group", col=colors[1], cex=.5, ylim=c(1, max(numberOfReads)))   
for (i in 2:4) {
  points(1:length(cumsumCountsOfReadsInGroup[[i]]), cumsumCountsOfReadsInGroup[[i]], col=colors[i], cex=.5)
}
legend("bottomright", legend= label, fill=colors)
axis(4, at = numberOfReads, label, lwd.ticks =4)
grid(lwd = 2, lty = 1)
dev.off()

png(paste0(path, "plots1.png"), width = 800, height = 800)
plot(1:length(cumsumCountsAboveFive), cumsumCountsAboveFive, ylim=c(1, sum),
     ylab = "Cumulative Number of Reads in Groups", xlab = "Number of Group")#, log="y", ylim=c(0, sum))

dev.off()

png(paste0(path, "test.png"), width = 800, height = 800)
plot(cumsumCountsOfReadsInGroup[[1]], countsOfReadsInGroup[[1]], xlim=c(1, sum), log="y", col=colors[1], cex=.5,
     xlab = "Cumulative Number of Reads in Groups", ylab="Number of Reads in Group (log Scale)")
for (i in 2:4) {
  points(cumsumCountsOfReadsInGroup[[i]], countsOfReadsInGroup[[i]], col=colors[i], cex=.5)
}
axis(3, at = numberOfReads, label, lwd.ticks =4)
legend("bottomright", legend= label, fill=colors)
grid(lwd = 2, lty = 1)
dev.off()

