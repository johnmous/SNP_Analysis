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
listOfFiles <- grep("vcf1$", dir(path), value=TRUE)

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
      
      # add aditonal variabel to filter out tri-allelic variants
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

oskarTableOfSNPs <- FetchFilesInTable(listOfFiles, path)
uniqueSNPID <- paste(oskarTableOfSNPs[,1], oskarTableOfSNPs[,2], oskarTableOfSNPs[,3], sep="_")
SNPIDsTableOskar <- data.frame(uniqueSNPID=uniqueSNPID, ID=oskarTableOfSNPs$ID) 

oskarTableOfSNPsStats <- extraStatsTVC(oskarTableOfSNPs)
oskarSNPFreq <- as.numeric(oskarTableOfSNPsStats$SNPFreq)
mean(oskarSNPFreq, na.rm = T)
hist(oskarSNPFreq, breaks = 100)
oskarSNPQual <- as.numeric(oskarTableOfSNPsStats$QUAL)
mean(oskarSNPQual)
hist(oskarSNPQual, breaks = 100)

IoannisTableOfSNPs <- FetchFilesInTable("Yanfang_Ioannis_V")
uniqueSNPID <- paste(IoannisTableOfSNPs[,1], IoannisTableOfSNPs[,2], IoannisTableOfSNPs[,3], sep="_")
SNPIDsTableIoannis <- data.frame(uniqueSNPID=uniqueSNPID, ID=IoannisTableOfSNPs$ID) 

IoannisTableOfSNPsStats <- extraStatsTVC(IoannisTableOfSNPs)
ioannisSNPFreq <- as.numeric(IoannisTableOfSNPsStats$SNPFreq)
mean(ioannisSNPFreq, na.rm = T)
hist(ioannisSNPFreq, breaks = 100)
ioannisSNPQual <- as.numeric(IoannisTableOfSNPsStats$QUAL)
mean(ioannisSNPQual)
hist(ioannisSNPQual, breaks = 100)

a <- SNPIDsTable[ SNPIDsTable$ID=="V7",][,1]
b <- SNPIDsTable[ SNPIDsTable$ID=="V8",][,1]
length(b)

pairWise <- function(setA, setB){
  boolean <- setA %in% setB  
  common <- setA[boolean]
  setAOnly <- setA[!boolean]
  setBOnly <- setB[! setB %in% setA ]
  overlap = list(common = common, setAOnly = setAOnly, setBOnly = setBOnly)
  return(overlap)
}

common <- as.character(SNPIDsTable[ SNPIDsTable$ID=="V1",][,1])
for (i in 2:5) {
  ID <- paste0("V",i)
  setB <- as.character(SNPIDsTable[ SNPIDsTable$ID == ID,][,1])
  common <- pairWise(common, setB)
  print(length(common))
}
  

  for (i in 2:5) {
    IDA <- paste0("V",i-1)
    IDB <- paste0("V",i)
    setA <- as.character(SNPIDsTable[ SNPIDsTable$ID == IDA,][,1])
    setB <- as.character(SNPIDsTable[ SNPIDsTable$ID == IDB,][,1])
    common <- pairWise(setA, setB)
    print(paste("common SNPs between version", i-1 , "and", i, "is", length(common), ". Version", i-1 , 
                "length:", length(setA), "Version",i , "length", length(setB)))
  }

############### $$$$$$$$$$$$$$$$$$$$$$$$$$ ####################
# compare SNP sets and decide on their qiality
setA <- as.character(SNPIDsTableOskar[ SNPIDsTableOskar$ID == "V1",][,1])
setB <- as.character(SNPIDsTableIoannis[ SNPIDsTableIoannis$ID == "V5",][,1])

common <- pairWise(setA, setB)

# count SNPs and indels occurences in the sets

setAVar <- common$setAOnly
setBVar <- common$setBOnly
common <- common$common

varTypeSetA <- substr(setAVar, nchar(setAVar)-2, nchar(setAVar))
varTypeSetB <- substr(setBVar, nchar(setBVar)-2, nchar(setBVar))
varTypeCommon <- substr(common, nchar(common)-2, nchar(common))

# Number of SNPs in set A 
length(which(varTypeSetA == "SNP"))
length(which(varTypeSetA == "SNP"))/length(varTypeSetA)
# number of Indels in set A
length(varTypeSetA) - length(which(varTypeSetA == "SNP"))
(length(varTypeSetA) - length(which(varTypeSetA == "SNP")))/length(varTypeSetA)

# Number of SNPs in set B 
length(which(varTypeSetB == "SNP"))
length(which(varTypeSetB == "SNP"))/length(varTypeSetB)
# Number of Indels in set B 
length(varTypeSetB) - length(which(varTypeSetB == "SNP"))
(length(varTypeSetB) - length(which(varTypeSetB == "SNP")))/length(varTypeSetB)

# Number of SNPs in overlap
length(which(varTypeCommon == "SNP"))
length(which(varTypeCommon == "SNP"))/length(varTypeCommon)
# Number of Indels in overlap 
length(varTypeCommon) - length(which(varTypeCommon == "SNP"))
(length(varTypeCommon) - length(which(varTypeCommon == "SNP")))/length(varTypeCommon)
