# Title: compare the sets of SNPs output resulting from different TVC settings
# Author: Ioannis Moustakas, i.moustakas@uva.nl

# a function to produce few extra columns with variant statistics with TVC input
extraStatsTVC <- function(tableOfSNPs) {
  statistics <- strsplit(as.character(tableOfSNPs[["INFO"]]), ";")
  allelicFreq <- sapply(statistics, function(row) { strsplit(row[1], "=")[[1]][2]})
  readFreq <- sapply(statistics, function(row) { strsplit(row[3], "=")[[1]][2]})
  tableOfSNPs$SNPFreq <- allelicFreq
  tableOfSNPs$Depth <- readFreq
  return(tableOfSNPs)
}

path = "/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1117-Ter_Kuile/MAD1117-P001-gDNA_sequencing/MAD1117-P001-E001_2014_32x_gDNA_Yanfang_svleeuw1/Scratch/reRun/VCFCompiled/"
listOfFiles <- grep("vcf$", dir(path), value=TRUE)

# 1. Read all variant files
FetchFilesInTable <- function(listOfFiles, path){
    
  # Make a data.frame variable to store all samples
  tableOfSNPs <- data.frame()
  
  # read in and process all samples (vcf files) 
  for (i in 1:length(listOfFiles)) {
    #catch the error in the case the vcf file contains no lines
    possibleError <- tryCatch(read.table(paste(path, listOfFiles[i], sep = ""), sep="\t"),
                              error=function(e) e) 
    
    if (inherits(possibleError, "error")){
      warning(paste("This file is emtpy: ", listOfFiles[i]), immediate.= TRUE)
    }
    
    # if no error message was caught, continue
    if (!inherits(possibleError, "error")){
      # read file and store
      Sample <- read.table(paste(path, listOfFiles[i], sep = ""), sep="\t", head=T, comment.char = "#")
      
      # change header
      colnames(Sample)<-c("CHROM","POS","TYPE","REF","ALT","QUAL","FILTER","INFO","FORMAT","2-1")
      
      # assign filename as ID to samples
      Sample$ID <- paste0("V",i)
      
      # add aditonal variable to filter out tri-allelic variants
      Sample$POS2<-paste(Sample$POS,Sample$REF,Sample$ALT,sep="_")
      
      # add variant type (SNP/DEL/INS)
      Sample$TYPE<-with(Sample,ifelse(Sample$REF!="A" & Sample$REF!="T" & Sample$REF!="G"& Sample$REF!="C",c("DEL"),
                                      ifelse(Sample$ALT=="A"|Sample$ALT=="T"|Sample$ALT=="C"|Sample$ALT=="G",c("SNP"),c("INS"))))
      
      # put the sample in the data.frame 
      tableOfSNPs <- rbind(tableOfSNPs, Sample, deparse.level = 1)
    }
  }
  return(tableOfSNPs)
}


tableOfSNPs <- FetchFilesInTable(listOfFiles, path)
tableOfNoCalls <- FetchFilesInTable(listOfFiles, path)
nrow(tableOfSNPs)

# compile a vector with a unique identifier for each SNP
uniqueSNPID <- paste(tableOfSNPs[,1], tableOfSNPs[,2], tableOfSNPs[,3], sep="_")
uniqueVariantsInTable <- unique(uniqueSNPID)
# number of unique variants
length(uniqueVariantsInTable)
# number of variants per sample
uniqueSampleIDs <- unique(tableOfSNPs$ID)
sampleSizePerID <- vector()
for (sample in uniqueSampleIDs){
  sampleSizePerID <- c(sampleSizePerID, nrow(tableOfSNPs[ tableOfSNPs$ID == sample, ]))
}

# min and max number of variants per sample
min(sampleSizePerID)
max(sampleSizePerID)

length(unique(uniqueSNPID))/nrow(tableOfSNPs)

##################### $$$$$$$$$$$$$$$$$$$$$$ ######################
hist(tableOfSNPs[ tableOfSNPs$CHROM == "gi|383748882|gb|AJKG01000054.1|_consensus", ]$POS, breaks=100)
hist(tableOfNoCalls[ tableOfNoCalls$CHROM == "gi|383748882|gb|AJKG01000054.1|_consensus", ]$POS, breaks=100)

hist(tableOfSNPs[ tableOfSNPs$CHROM == "gi|383748912|gb|AJKG01000024.1|_consensus", ]$POS, breaks=100)
hist(tableOfNoCalls[ tableOfNoCalls$CHROM == "gi|383748912|gb|AJKG01000024.1|_consensus", ]$POS, breaks=100)

