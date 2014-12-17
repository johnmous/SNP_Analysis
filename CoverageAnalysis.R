# Author: Ioannis Moustakas, i.moustakas@uva.nl 
# Title: Coverage analysis of BAM files

library(GenomicFeatures)
library(GenomicRanges)
library(Biostrings)
library(ShortRead)
library(gridExtra)


path="/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1117-Ter_Kuile/MAD1117-P001-gDNA_sequencing/MAD1117-P001-E001_2014_32x_gDNA_Yanfang_svleeuw1/Results/BAM_files/"

# get the names of all BAM files in the directory
listOfFiles = grep(".bam$", dir(path), value=TRUE)
statistics <- data.frame(
                  a = NULL,
                  b = NULL,
                  c = NULL,
                  d = NULL,
                  e = NULL,
                  g = NULL
)

coverageVsCumulativeLengthList <- list()
coverageVsLengthList <- list()
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
  depthSubsetIndexes <- seq(1, length(uniqueDepths), by = 11)
  
  sumOfLengths <- vector()
  meanDepth <- vector()
  for ( depthIndex in depthSubsetIndexes ) {
              depth <- uniqueDepths[ depthIndex : (depthIndex + 10) ]
              # print(depth)
              indexes <- which(coverage%in%depth)
              sumOfLengths <- c(sumOfLengths, sum(lengths[indexes], na.rm = T))
              meanDepth <- c(meanDepth, mean(depth, na.rm = T))
            }
  
  genomeLength <- sum(sumOfLengths)
  coverageVsCumulativeLength <- data.frame( depth=meanDepth, PercentageOfGenome=1-(cumsum(sumOfLengths)/genomeLength))
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
                  e = coverageVsCumulativeLength[ coverageVsCumulativeLength$depth>=500, ][1,]$PercentageOfGenome,
                  g = coverageVsCumulativeLength[ coverageVsCumulativeLength$depth>=1000, ][1,]$PercentageOfGenome,
                  row.names= fileName
                  ))
i = i+1
}

names(statistics) <- c("Over 20", "Over 50", "Over 100", 
                       "Over 200", "Over 500", "Over 1000")
statistics <- round(statistics, digits=3)

pdf(paste(path, "data_output.pdf", sep=""))
plot(coverageVsCumulativeLengthList[[1]], xlim=c(0, 1000))
for (i in 2:length(coverageVsCumulativeLengthList)){
  points(coverageVsCumulativeLengthList[[i]], col=14*i)
}
legend(500, 1, legend = listOfFiles, col = c(1,22,44))
plot(coverageVsLengthList[[1]], xlim=c(0, 1000))
for (i in 2:length(coverageVsLengthList)){
  points(coverageVsLengthList[[i]], col=14*i)
}
grid.newpage()
grid.table(statistics, gp=gpar(fontsize=10))
dev.off()


fasta2<-readDNAStringSet(paste(path, "ps_LCT_PA102.noCommas.fasta", sep=""), format="fasta"
                         ,nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
fasta<-DNAString(toString(fasta2))


# compute gc content in a sliding window (as a fraction, if you want % : *100):
gc <- rowSums(letterFrequencyInSlidingView(fasta, 20, c("G","C")))/20
length(gc)
gc[(length(gc)+1):length(cvg$`gi|32141095|ref|NC_003888.3|`)]<-gc[length(gc)]


# the coverage vector can be shorter than the sequence. 
# The length is equal to the largest end position of an alignment in the data:
pad.right <- function (rle, value=0, len) {
  if (length(rle) > len) stop("length mismatch: len must be >= length(rle)!")
  if (len > length(rle)) 
    c(rle, Rle(value, len-length(rle))) 
  else rle
}
cvg <- pad.right(cvg$`gi|32141095|ref|NC_003888.3|`, len=length(gc))


# we have to process the coverage vector to have the same window size:
cvgw <- runmean(x=cvg, k=window.size, endrule = c("constant"))

#?runmean
# endrule is important to get a vector of same size, see ?runmean



pdf("All reads GC content histogram.pdf")
hist(gc)
dev.off()


