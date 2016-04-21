# Author: Ioannis Moustakas, i.moustakas@uva.nl 
# Title: Coverage analysis of BAM files

library(GenomicFeatures)
library(GenomicRanges)
library(Biostrings)
library(ShortRead)
library(gridExtra)
library(S4Vectors)
library(seqinr)

path="/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1101-Hamoen/MAD1101-P001-Leendert/MAD1101-P001-E008_2015_gDNA_Bacillus_svleeuw1/Results/alignment/coverageDistribution/"

# get the names of all BAM files in the directory
listOfFiles = grep(".bam$", dir(path), value=TRUE)
# An emtpy data.frame to store statistics on genome Coverage. 
# Each row of this data frame is a different experiment
statistics <- data.frame(
                  a = NULL,
                  b = NULL,
                  c = NULL,
                  d = NULL,
                  e = NULL,
                  g = NULL, 
                  h = NULL,
                  i = NULL
)

coverageVsCumulativeLengthList <- list()
coverageVsLengthList <- list()
step = 5
i = 1
for (fileName in listOfFiles){
  #window.size <- 1001 # configure the window size
  x <- readGAlignments(paste(path, fileName, sep=""))
  cvg <- coverage(x)
  
  coverage <- vector()
  lengths <- vector()
  for (chrom in cvg@listData) {
    coverage<- c(coverage, chrom@values)
    lengths <- c(lengths, chrom@lengths)
  }
  genomeLength <- sum(lengths)
  
  uniqueDepths <- sort(unique(coverage))
  depthSubsetIndexes <- seq(1, length(uniqueDepths), by = step+1)
  
  sumOfLengths <- vector()
  meanDepth <- vector()
  for ( depthIndex in depthSubsetIndexes ) {
              depth <- uniqueDepths[ depthIndex : (depthIndex + step) ]
              # print(depth)
              indexes <- which(coverage%in%depth)
              sumOfLengths <- c(sumOfLengths, sum(lengths[indexes], na.rm = T))
              meanDepth <- c(meanDepth, mean(depth, na.rm = T))
            }
  
  genomeLength <- sum(sumOfLengths)
  coverageVsCumulativeLength <- data.frame( depth=meanDepth, PercentageOfGenome= cumsum(sumOfLengths)/genomeLength)
  coverageVsCumulativeLengthList[[i]] <- coverageVsCumulativeLength
  # plot(coverageVsCumulativeLength, xlim=c(0, 1000))
  coverageVsLength <- data.frame( depth=meanDepth, PercentageOfGenome=sumOfLengths/genomeLength)
  coverageVsLengthList[[i]] <- coverageVsLength
  # plot(coverageVsLength, xlim=c(0, 1000))
  statistics <- rbind(statistics, data.frame(
                  a = coverageVsCumulativeLength[ coverageVsCumulativeLength$depth>=20, ][1,]$PercentageOfGenome,
                  b = coverageVsCumulativeLength[ coverageVsCumulativeLength$depth>=50, ][1,]$PercentageOfGenome,
                  c = coverageVsCumulativeLength[ coverageVsCumulativeLength$depth>=100, ][1,]$PercentageOfGenome,
                  d = coverageVsCumulativeLength[ coverageVsCumulativeLength$depth>=200, ][1,]$PercentageOfGenome,
                  e = coverageVsCumulativeLength[ coverageVsCumulativeLength$depth>=350, ][1,]$PercentageOfGenome,
                  h = coverageVsCumulativeLength[ coverageVsCumulativeLength$depth>=500, ][1,]$PercentageOfGenome,
                  i = coverageVsCumulativeLength[ coverageVsCumulativeLength$depth>=1000, ][1,]$PercentageOfGenome,
                  row.names= fileName
                  ))
i = i+1
}

names(statistics) <- c("Below 20", "Below 50", "Below 100", "Below 200",
                      "Below 350", "Below 500", "Below 1000")
statistics <- round(statistics, digits=3)

# print the report in a PDF
colors= rainbow(length(coverageVsLengthList))
pdf(paste(path, "CovearageAnalysis.pdf", sep=""))

plot(coverageVsLengthList[[1]], xlim=c(0, 1000), ylim=c(0, 0.08), col=colors[1], cex=0.7, 
     ylab="Fraction of Genome", xlab="Coverage")
title(main = "Distribution of Coverage Across Genome")
for (i in 2:length(coverageVsLengthList)){
  points(coverageVsLengthList[[i]], col=colors[i], cex=0.7)
}
legend(x="topright", legend = listOfFiles, col = colors[1:i], pch =19, cex=0.7)

plot(coverageVsCumulativeLengthList[[1]], xlim=c(0, 1000), ylim=c(0,1), col=colors[1], cex=0.7, 
     ylab="Cumulative Fraction of Genome", xlab="Coverage")
title(main = "Cumulative Distribution of Coverage Across Genome")

for (i in 2:length(coverageVsCumulativeLengthList)){
  points(coverageVsCumulativeLengthList[[i]], col= colors[i], cex=0.7)
}
legend(x="bottomright", legend = listOfFiles, col = colors[1:i], pch=19, cex = 0.7)

grid.newpage()
grid.table(statistics, gp=gpar(fontsize=8))
title(main = "Fraction of Genome With Coverage:")

dev.off()



######## Reference Analysis ########

window.size <- 1001 # configure the window size
x <- readGAlignments(paste(path, "IrisGCrichSorted_reads_001.bam", sep=""))
cvg <- coverage(x)

fastaFile<- read.fasta(paste(path, "ps_LCT_PA102.fasta", sep=""), seqtype= "DNA")
sequence <- unlist(fastaFile)
windowSize = 50
GCContent <- vector()
breakPoints <- seq(windowSize, length(sequence), by = windowSize)
for (coord in breakPoints) {
  GCContent <- c(GCContent, GC(sequence[(coord-windowSize):coord]))                    
  }

length(GCContent)
length(breakPoints)
plot(breakPoints, GCContent, xlab="Genome" )

# compute gc content in a sliding window (as a fraction, if you want % : *100):
gc <- rowSums(letterFrequencyInSlidingView(fasta, 20, c("G","C")))/20
length(gc)
gc[(length(gc)+1):length(cvg$'gi|556503834|ref|NC_000913.3|' )]<-gc[length(gc)]


# the coverage vector can be shorter than the sequence. 
# The length is equal to the largest end position of an alignment in the data:
pad.right <- function (rle, value=0, len) {
  if (length(rle) > len) stop("length mismatch: len must be >= length(rle)!")
  if (len > length(rle)) 
    c(rle, Rle(value, len-length(rle))) 
  else rle
}
cvg <- pad.right(cvg$'gi|556503834|ref|NC_000913.3|', len=length(gc))


# we have to process the coverage vector to have the same window size:
cvgw <- runmean(x=cvg, k=window.size, endrule = c("constant"))

#?runmean
# endrule is important to get a vector of same size, see ?runmean



pdf(paste(path, "All reads GC content histogram.pdf", sep=""))
hist(gc, breaks = 100)
mean(gc)
dev.off()


