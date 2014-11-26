##########################################################################
##  		R version 3.0.1 (2013-05-16) 			##
##	Project 06.) Make overview of found variants			##
##	Data: 6 Samples of E.coli on IonProton (.bam)			##
##	Mapped: TMAP directly from IonReport				##
##	ref: ec_K12_MG1655.FASTA					##
##	Variants are called with Unified Genotyper(GATK)		##
##	Author script: Iris Kolder 
##  Modified by Ioannis Moustakas (i.moustakas@uva.nl)					##
##	last edited: 18-09-2014						##
##########################################################################

options(stringAsFactors = FALSE)

# The directory where vcf files are stored
path="/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1117-Ter_Kuile/MAD1117-P001-gDNA_sequencing/MAD1117-P001-E001_2014_32x_gDNA_Yanfang_svleeuw1/Results/"

# Load reference sequence
reference <- read.fasta(paste(path, "ps_LCT_PA102.fasta", sep=""), seqtype= "DNA")
##### reference <- unlist(reference[[1]])

# Name for annotation file 
# annotFile = "BSubltilis168_NCBI.gff"

# a function to produce few extra columns with variant statistics
extraStats <- function(tableOfSNPs) {
  statistics <- strsplit(as.character(tableOfSNPs[["2-1"]]), ":")
  allelicDepth <- unlist(lapply(statistics, function(row) { row[2] }))
  readDepth <- unlist(lapply(statistics, function(row) { row[3] }))
  referDepth <- unlist(lapply(strsplit(allelicDepth, ","), function(row) { row[1] } ))
  referDepth <- as.numeric(referDepth)
  altDepth <- unlist(lapply(strsplit(allelicDepth, ","), function(row) { row[2] } ))
  altDepth <- as.numeric(altDepth)
  totalDepth <- referDepth + altDepth
  altFrequency <- altDepth/totalDepth
  tableOfSNPs$SNPFreq <- altFrequency
  tableOfSNPs$Depth <- totalDepth
  return(tableOfSNPs)
}

# A function needed tp manipulate the annotation table
insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  existingDF
}

#  a function that adds annotation to the table of SNPs
annotateSNPTable <- function (path, annotFile, tableOfSNPs) {
    # read the annotation file & give columns more descriptive names
    annotation <- read.table(paste(path, annotFile, sep = ""), sep="\t", head=TRUE, stringsAsFactors=FALSE)
    names(annotation)[4] <- "Start"
    names(annotation)[5] <- "Stop"
    names(annotation)[9] <- "Description"
    
    # modify the annotation data.frame so a new line is added for every Non-annotated region. 
    # Each of them should have a unique name inthe "region" column
    previousElementStopCoordinates = 0 
    nonAnnotRegionCount = 0
  
    # save the last element's start coordinates to use it as a break trigger for the repeat loop
    lastElementStartCoord = tail(annotation, 1)[[4]]
  
    i = 1
    
    # Add an extra line in the annotation table to create an entry for non annotated region elements
    # Each elements bears a unique number to teel each other appart
    # Repeat the loop until the last element of the the annotation table is reached
    # Note the annotation table is mutatted as the loop proceeds by adding an extra lines
    repeat {    
      elementStartCoordinates <- as.numeric(annotation$Start[i])  
      #  break the loop if the last element in the gff file is reached
      if (elementStartCoordinates == lastElementStartCoord){
        break
      }
      gapBetweenElements <- elementStartCoordinates - previousElementStopCoordinates      
      if (gapBetweenElements > 1) {
        nonAnnotRegionCount <- nonAnnotRegionCount +1
        regionName <- paste("NonAnnotRegion", nonAnnotRegionCount, sep="_n_" ) 
        newRow <-  c( names(annotation)[1], names(annotation)[2], regionName, previousElementStopCoordinates+1, elementStartCoordinates-1, NA, NA, NA, NA ) 
        annotation <- insertRow(annotation, newRow, i)
      }
      previousElementStopCoordinates <- as.numeric(annotation$Stop[i])
      i = i+1
    }

    SNPPos <- as.numeric(tableOfSNPs[["POS"]])
    annotation$Start <- as.numeric(annotation$Start)
    annotation$Stop <- as.numeric(annotation$Stop)
    
    # get the type of element(s) for each SNP 
    elements <- sapply(SNPPos, function(x) {
                      positionOnTable <- which(annotation$Start<=x & annotation$Stop>=x)
                      elements <- annotation$region[positionOnTable]
                      # drop the NonAnnotRegion if other genomic elements are found along with it
                      if (length(elements) > 1 && length(grep("NonAnnotRegion", elements))>=1 ) {
                        index <- which(sapply("NonAnnotRegion", grepl, elements))
                        elements <- elements[-index]
                      }
                      return(elements)
                    })
    
    # replace empty vectors with "Non Annotated Region"
    elements[which(sapply(c(1:length(elements)), function (x) {length(elements[[x]])==0}))] <- "Non Annotated Region"
    
    # combine all elements of a vector into a single string, separated by comma
    elements <- sapply(c(1:length(elements)), function (x) {
                        as.vector(paste(elements[[x]], collapse=", "))
                      })
    
    tableOfSNPs$Elements <- elements

    
    # Get the gene name for each SNP. Names are in the descreption of the "gene" field of the annotation file 
    geneNames <- sapply(SNPPos, function(x) {
                        positionOnTable <- which(annotation$Start<=x & annotation$Stop>=x & annotation$region=="gene")
                        geneNames <- str_match(annotation$Description[positionOnTable], "Name=(.*?);")
                        return(geneNames)
                      })
    

    # go through the list and get the second element of the vector, that is the name of the gene
    geneNames <- sapply(c(1:length(geneNames)), function (x) {
                        geneNames[[x]] <- geneNames[[x]][2]
                      })
    
    # replace empty vectors with the description from "element" field 
    emptyGeneNameIndex <-  which(is.na(geneNames)) # which(sapply(c(1:length(geneNames)), function (x) {length(geneNames[[x]])==0 }))
    geneNames[emptyGeneNameIndex] <- elements[emptyGeneNameIndex]
    geneNames[is.na(geneNames)] <- "Non Annotated Region"
    
    tableOfSNPs$GeneName <-  geneNames
    
    # count the occurence of a gene and put it in a vector.
    count=0
    geneCount <-vector()
    previousGene = "anyGene"
    geneNames <- c(geneNames, "anyGene")
    for (gene in geneNames) {
      if (gene != previousGene) {
        geneCount = c(geneCount, rep(count, count))
        previousGene = gene
        count=1
      }
      else {
        count = count + 1
      }
    }
    
    tableOfSNPs$GeneCount <-  geneCount
    
    # Get the gene description for each SNP. Names are in the descreption of the "CDS" field of the annotation file 
    geneDescr <- sapply(SNPPos, function(x) {
                        positionOnTable <- which(annotation$Start<=x & annotation$Stop>=x & annotation$region=="CDS" )
                        geneDescr <- str_match(annotation$Description[positionOnTable], "Note=(.*?);")
                        return(geneDescr)
                      })
    
    # replace empty vectors with "Non Annotated Region"
    geneDescr[which(sapply(c(1:length(geneDescr)), function (x) {length(geneDescr[[x]])==0}))] <- "Non Annotated Region"
    
    # go through the list and get the second element of the vector, that is the description of the gene function
    geneDescr <- sapply(c(1:length(geneDescr)), function (x) {
                        geneDescr[[x]] <- geneDescr[[x]][2]
                      })  
    tableOfSNPs$Function <-  geneDescr
  
    return(tableOfSNPs)
}


# get the names of all vcf files in the directory
listOfFiles = grep("vcf$", dir(path), value=TRUE)
print 

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
    Sample$ID<-listOfFiles[i]
    
    # add aditonal variabel to filter out tri-allelic variants
    Sample$POS2<-paste(Sample$POS,Sample$REF,Sample$ALT,sep="_")
    
    # add variant type (SNP/DEL/INS)
    Sample$TYPE<-with(Sample,ifelse(Sample$REF!="A" & Sample$REF!="T" & Sample$REF!="G"& Sample$REF!="C",c("DEL"),
                                    ifelse(Sample$ALT=="A"|Sample$ALT=="T"|Sample$ALT=="C"|Sample$ALT=="G",c("SNP"),c("INS"))))
    
    # put the sample in the data.frame 
    tableOfSNPs <- rbind(tableOfSNPs, Sample, deparse.level = 1)
    }
}


# make barplot variant count
count <- table(tableOfSNPs$TYPE, tableOfSNPs$ID)
pdf(paste(path, "Barplot of variant counts.pdf", sep=""))
barplot(count,legend=rownames(count),main="Barplot of variant counts in all the samples",ylab=c("Variant counts"),xlab=c("Sample names"),ylim=c(0,600), las = 2)
dev.off()

# print out table
write.table(count,file="Variant_Counts_table.txt", sep="\t", row.names=T,col.names=T)

# Count overlapping variants
temp <- data.frame(table(tableOfSNPs$POS2,tableOfSNPs$TYPE))
overlap_counts<-table(temp$Var2,temp$Freq)[,-1]
pdf(paste(path, "Overlap of variants foun in samples.pdf", sep=""))
par(mfrow=c(3,1))
barplot(overlap_counts[1,],beside=T,main="Overlap of deletions found in the samples",ylab=c("Variant counts"),xlab=c("Number of samples variant was found in"))
barplot(overlap_counts[2,],beside=T,main="Overlap of insertions found in the samples",ylab=c("Variant counts"),xlab=c("Number of samples variant was found in"))
barplot(overlap_counts[3,],beside=T,main="Overlap of SNPs found in the samples",ylab=c("Variant counts"),xlab=c("Number of samples variant was found in"))
dev.off()

# print out table
write.table(overlap_counts, file=paste(path, "Overlap_Counts_table.txt", sep=""), sep="\t", row.names=T,col.names=T)

# Function to add a few more columns with extra statistics
 tableOfSNPs <- extraStats(tableOfSNPs)

# print variance table
# Remove useless columns from table of SNPs
drops <- c("INFO", "FORMAT", "2-1", "POS2")
tableOfSNPs <- tableOfSNPs[,!(names(tableOfSNPs) %in% drops)]

# Fill in the filter column, starting with low quality and going up
tableOfSNPs$FILTER <- "LowQual"
tableOfSNPs$QUAL <- as.numeric(tableOfSNPs$QUAL)
tableOfSNPs$FILTER[tableOfSNPs$QUAL > 30] <- "MediumQual"
tableOfSNPs$FILTER[tableOfSNPs$QUAL > 70] <- "HighQual"
tableOfSNPs$FILTER[tableOfSNPs$QUAL > 100] <- "VeryHighQual"

# run the annotation function
# tableOfSNPs <- annotateSNPTable(path, annotFile, tableOfSNPs)
#   A function to discover hotspots in the variants table and output their indexes + statistics table
#   User must imptut a vector with the positions of variants on the genome, the window size and the 
# number of variants
#   Algorithm will search inside the variant vector areas of length = window size that contain the
# number of variants specified by the user, or more

HotSpots <- function(variantsPos, windowSize, variants){
  stopifnot(length(variantsPos) > 3 )
  lastPos <- tail(variantsPos, 1) # last variant in vector
  # walk through the variant vector and output the variants inside Hotspots
  hotspotPosList <- sapply(seq(1, lastPos, by=2), function(walker){
    variantsInsideWindow <- variantsPos[(variantsPos > walker) & (variantsPos < walker+windowSize)]
    numberOfVariantsInsideWindow <- length(variantsInsideWindow)
    if (numberOfVariantsInsideWindow >= variants){
      return(variantsInsideWindow)
    }
  }
  )
  
  hotspotPosList <- unique(unlist(hotspotPosList))
  hotspotIndexes <- match(hotspotPosList, variantsPos)
  
  # Hotspot table of statistics: compile a table with hotspots start and stop positions, 
  # length and number of variants on them
  variantsDistance <- diff(hotspotPosList)
  hotspotStopPos <- hotspotPosList[variantsDistance>windowSize]
  hotspotStopPosIndex <- match(hotspotStopPos, hotspotPosList)
  hotspotStartPosIndex <- hotspotStopPosIndex +1
  hotspotStopPos <- c(hotspotStopPos, tail(hotspotPosList, 1))
  hotspotStartPos <- hotspotPosList[1]
  hotspotStartPos <- c(hotspotStartPos, hotspotPosList[hotspotStartPosIndex])
  hotspotLength <- hotspotStopPos - hotspotStartPos
  numberOfHotspots <- length(hotspotStartPos)
  numberOfVariantsInHotspot <- sapply(1:numberOfHotspots, function(x) {
    length(which( hotspotPosList <= hotspotStopPos[x] & 
                    hotspotPosList >= hotspotStartPos[x]))
  }
  )
  
  # compile information in a data.frame
  # if there is a hotspot detected, proceed
  if(length(hotspotLength)>0) {
    hotspotStats <- data.frame(StartPosition = hotspotStartPos, 
                               StopPosition = hotspotStopPos, 
                               Length = hotspotLength,
                               NumberOfVariants = numberOfVariantsInHotspot)
    #return a list with the indexes (position on variants vector) and the statistics
    returnList <- list(Indexes = hotspotIndexes, Stats =  hotspotStats)
  }
  # if no hotspot is detected, return NA
  else {
    returnList <- list(Indexes = NA, Stats =  NA)
  }
  return(returnList)
}

#  a function to check for homopolymer regions close to variants
variantsInsideHomop <- function(posOfVariants, reference, radious, minHomopLength){
  # Extract a sequence of length = radious around each deletion from the reference sequence
  listOfRegions <- sapply(posOfVariants, function(pos){
    left = pos-radious+1
    right = pos+radious+1
    reference[left:right]
  }
  )
  
  # check for homopolymer regions in the list of sequences
  homopolymerLength <- apply(listOfRegions, 2, function(sequence) {
    length = length(sequence)
    windowSize = minHomopLength
    homopolymerPos <- sapply(windowSize:length, function(right){
      left = right-windowSize+1
      window <- sequence[left:right]
      uniqueBasesInWindow <- unique(window)
      numberOfUniqueBases <- length(uniqueBasesInWindow)
      if (numberOfUniqueBases == 1){
        return(T)
      }
      else{
        return(F)
      }
    }
    )
    homopolymerPos <- which(homopolymerPos)
    homopolymerLength <- length(homopolymerPos)
    if (homopolymerLength >= 1) {
      return(homopolymerLength+windowSize-1)
    }
    else {
      return(0)
    }
  }  
  )
  # hist(homopolymerLength)
  # print(homopolymerLength)
  variantsInsideHomopolBoolean <- homopolymerLength > 0
  return(variantsInsideHomopolBoolean) 
}
hotSpotRegion <- vector()
homopolRegion <- vector()
listOfIDs <- unique(tableOfSNPs$ID)
for (id in listOfIDs){
  sliceOfTableOfSNPs <- tableOfSNPs[ tableOfSNPs$ID == id, ]
  chromosomeNames <- levels(sliceOfTableOfSNPs$CHROM)
  hotSpotRowIndex <- vector()
  delInsideHomopolIndex <- vector()
  for (i in chromosomeNames) {
    # Hotspots
    varPos <- sliceOfTableOfSNPs[sliceOfTableOfSNPs$CHROM == i, ]$POS
    varRowIndex <- as.numeric(row.names(sliceOfTableOfSNPs[sliceOfTableOfSNPs$CHROM == i, ]))
    # in there are more than 3 variants in this chromosome, proceed 
    if (length(varPos)>3 ){
      hotSpots <- HotSpots(varPos, 30, 3)
      indexes <-hotSpots$Indexes
    }
    # else return an empty index vector
    else { indexes <- vector() }
    # if there are hotspots detected
    if(length(indexes) > 1){
      hotSpotRowIndex <- c(hotSpotRowIndex, varRowIndex[indexes])
    }
    
    # Homopolymers
    # check DELs only
    chromosomeSeq <- reference[[i]]
    delPos <- sliceOfTableOfSNPs[sliceOfTableOfSNPs$CHROM == i & sliceOfTableOfSNPs$TYPE == "DEL", ]$POS
    delRowIndex <- as.numeric(row.names(sliceOfTableOfSNPs[sliceOfTableOfSNPs$CHROM == i & sliceOfTableOfSNPs$TYPE == "DEL", ]))
    if (length(delPos)> 0 ){
      delInsideHomopolBool <- variantsInsideHomop(delPos, chromosomeSeq, 10, 4)
      # delInsideHomopolIndex <- delRowIndex[delInsideHomopolBool]
    }
    else {delInsideHomopolBool <- vector()}
    if(length(delInsideHomopolBool) > 0){
      delInsideHomopolIndex <- c(delInsideHomopolIndex, delRowIndex[delInsideHomopolBool])
    }
    
  }
  # a true /false vector to be implemented as a 
  rows <- row.names(sliceOfTableOfSNPs)
  hotSpotRegion <- c(hotSpotRegion, rows %in% hotSpotRowIndex)
  homopolRegion <- c(homopolRegion, rows %in% delInsideHomopolIndex)
}
  tableOfSNPs$HotSpotRegion <- hotSpotRegion
  tableOfSNPs$DelInHomopolymerRegion <- homopolRegion

# write the data to a tab separated file
 write.table(tableOfSNPs, paste(path, "variantsCompiledHotspotsHomopolymers.txt", sep=""), sep="\t", row.names=FALSE)

tableOfSNPs
