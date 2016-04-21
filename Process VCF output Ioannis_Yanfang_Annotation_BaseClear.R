##########################################################################
##    	R version 3.0.1 (2013-05-16) 			##
##	Project 06.) Make overview of found variants			##
##	Data: 6 Samples of E.coli on IonProton (.bam)			##
##	Mapped: TMAP directly from IonReport				##
##	ref: ec_K12_MG1655.FASTA					##
##	Variants are called with Unified Genotyper(GATK)		##
##	Author script: Iris Kolder 
##  Modified by Ioannis Moustakas (i.moustakas@uva.nl)					##
##	last edited: 18-09-2014						##
##########################################################################

############# **************** ################
# libraries to load
library(seqinr)


options(stringAsFactors = FALSE)

# The directory where vcf files are stored
path="/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1117-Ter_Kuile/MAD1117-P001-gDNA_sequencing/MAD1117-P001-E001_2014_32x_gDNA_Yanfang_svleeuw1/Results/outputBaseClearReference/TVC_Output/"

# Load reference sequence
reference <- read.fasta(paste(path, "Wild_Type_scaffold-sequences.fa", sep=""), seqtype= "DNA")
# process the names of the contigs in reference so they match the annotation
#contigNames <- sapply(strsplit(names(reference), split="|", fixed =T), function(string) string[4] )
#contigNames <- sapply(strsplit(contigNames, split=".", fixed =T), function(string) string[1]) 
#names(reference) <- contigNames
##### reference <- unlist(reference[[1]])

#################### $$$$$$$$$$$$$$$$$$$$ #####################
#### If there are more than two altermative variants reported in the same POSition in the VCF file,
#### break it into separate lines and return the new table 
SplitMultipleVariants <- function(tableOfSNPS) {
  tableLinesSplit <- data.frame()
  for (i in 1:nrow(tableOfSNPs)) {
    row<- tableOfSNPs[i, ]
    alternativeBases <- unlist(strsplit(row$ALT, split=',', fixed=T))
    alternativeFreq <-  unlist(strsplit(row$SNPFreq, split=',', fixed=T))
    
    ### Save the info for Alternative Variants in a separate column
    row$MultAltVariants <- F
    if (length(alternativeBases)>1) {
      row$MultAltVariants <- T  
    }
    
    # now iterate over the alternative bases and freq, alter the row and save it
    for (j in 1:length(alternativeBases)){   
      row$ALT <- alternativeBases[j]
      row$SNPFreq <- alternativeFreq[j]   
      ### Now decide on the type of variant (SNP, INS, DEL)
      if (nchar(row$REF) == nchar(row$ALT)) {
        row$TYPE <- "SNP"
      } else if (nchar(row$REF) > nchar(row$ALT)) {
        row$TYPE <- "DEL"
      } else {
        row$TYPE <- "INS"
      }
      
      # Build a Variant identifier string (CHROM + POS + ALT ) and save it in a separate column
      row$VariantID <- paste(row$CHROM, row$POS, row$ALT, sep='_')  
      
      tableLinesSplit <- rbind(tableLinesSplit, row) 
    } 
  }
  return(tableLinesSplit)
}

# Name for annotation file 
annotFile = "Wild_Type.gff"

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

# a function to produce few extra columns with variant statistics with TVC input
extraStatsTVC <- function(tableOfSNPs) {
  statistics <- strsplit(as.character(tableOfSNPs[["INFO"]]), ";")
  allelicFreq <- sapply(statistics, function(row) { strsplit(row[1], "=")[[1]][2]})
  readFreq <- sapply(statistics, function(row) { strsplit(row[3], "=")[[1]][2]})
  tableOfSNPs$SNPFreq <- allelicFreq
  tableOfSNPs$Depth <- readFreq
  return(tableOfSNPs)
}


# A function needed tp manipulate the annotation table
insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  existingDF
}

#################### $$$$$$$$$$$$$$$$$$$ ##################### 
# A function that adds annotation to the table of SNPs
annotateSNPTable <- function (path, annotFile, tableOfSNPs) {
  # read the annotation file & give columns more descriptive names
  annotation <- read.table(paste0(path, annotFile), sep="\t", head=F, quote="", stringsAsFactors=FALSE)
  # rename the columns of annotation table
  names(annotation)[1] <- "seqname"
  names(annotation)[2] <- "RefSeq"
  names(annotation)[3] <- "feature"
  names(annotation)[4] <- "start"
  names(annotation)[5] <- "end"
  names(annotation)[7] <- "strand"
  names(annotation)[9] <- "attribute"
  
  # A function to extract the gene name from the attribute field of gff file
  extractGeneName <- function(attribute) {
    if (length(attribute)==0){
      return("Nothing")
    } else { 
      nameField <- strsplit(attribute, ";", fixed=T)[[1]][1]
      strsplit(nameField, "=", fixed=T)[[1]][2]
    }
  }
  
  # A function to extract the gene description from the attribute field of gff file
  extractGeneAttribute <- function(attribute) {
    if (length(attribute)==0){
      return("Nothing")
    } else { 
      nameField <- strsplit(attribute, ";", fixed=T)[[1]][4]
      strsplit(nameField, "=", fixed=T)[[1]][2]
    }
  }
  
  contigsInAnnot <- unique(annotation$seqname)
  annotationEnhanced <- data.frame()
  #names(annotationEnhanced) <- names(annotation)
  # count the non coding regions found in anotation table
  nonCodingCount = 1
  
  for (contigID in contigsInAnnot) {
    annotationSlice <- annotation[ annotation$seqname == contigID, ]
    row.names(annotationSlice) <- row.names(annotation)
    
    # select only the gene sections
    ###  annotationSlice <- annotationSlice [ annotationSlice$feature=="gene", ]
    
    # modify the annotation data.frame so a new line is added for every Non-annotated region. 
    # Each of them should have a unique name in the "region" column
    previousElementStopCoordinates = 0 
    nonAnnotRegionCount = 0
    previousElementStrand = "+"
    
    # save the last element's start coordinates to use it as a break trigger for the repeat loop
    lastElementStartCoord = tail( annotationSlice$start, 1)
    # this index is used to run over annotation table rows
    i = 1
    geneBefore <- "Nothing"
    
    # Add an extra line in the annotation table to create an entry for non annotated region elements
    # Each element bears a unique number to tell each other appart
    # Repeat the loop until the last element of the the annotation table is reached
    # Note the annotation table is mutatted as the loop proceeds adding extra lines
    repeat {
      elementStartCoordinates <- as.numeric(annotationSlice$start[i]) 
      elementStopCoordinates <- as.numeric(annotationSlice$end[i])   
      #  break the loop if the last element in the gff file is reached
      if (elementStartCoordinates == lastElementStartCoord){
        break
      }
      strand = "+"
      gapBetweenElements <- elementStartCoordinates - previousElementStopCoordinates    
      elementStrand <- annotationSlice$strand[i]
      element <- annotationSlice$feature[i]
      # extract the gene name if element is gene
#       if (element == "gene"){
#         geneAfter <- extractGeneName( annotationSlice$attribute[i] )
#         # print("GeneAfter :") 
#         # print(geneAfter)
#       }
      
      # if there is a gap between genes => non coding region
      if ( gapBetweenElements > 1 ) { 
        dummyRow <-  c( contigID, "MADInHouseScript", "NonCoding", previousElementStopCoordinates+1, 
                        elementStartCoordinates-1, "0", strand, "0", "attributes" ) 
        annotationInsert <- annotationSlice[i, ]
        annotationInsert <- insertRow(annotationInsert, dummyRow, 1)
        annotationEnhanced <- rbind(annotationEnhanced, annotationInsert)
      } else {
        annotationEnhanced <- rbind(annotationEnhanced, annotationSlice[i, ])
      }
    i = i+1
    previousElementStopCoordinates <- elementStopCoordinates
    }
      
####################### $$$$$$$$$$$$$$$$$$$$ ########################

      # Go through the enhansed annotation table and rpoperly name the attributes field
      previousElementStopCoordinates = 0 
      nonAnnotRegionCount = 0
      previousElementStrand = "+"
      
      # save the last element's start coordinates to use it as a break trigger for the repeat loop
      lastElementStartCoord = tail( annotationEnhanced$start, 1)
      # this index is used to run over annotation table rows
      i = 1
      geneBefore <- "Nothing"
      nonCodingPositions <- annotationEnhanced$feature == "NonCoding"

      repeat {
        
        elementStartCoordinates <- as.numeric(annotationEnhanced$start[i]) 
        elementStopCoordinates <- as.numeric(annotationEnhanced$end[i])   
        #  break the loop if the last element in the gff file is reached
        if (elementStartCoordinates == lastElementStartCoord){
          break
        }
        strand = "+"
        gapBetweenElements <- elementStartCoordinates - previousElementStopCoordinates    
        elementStrand <- annotationEnhanced$strand[i]
        element <- annotationEnhanced$feature[i]
        
      if (element ==  "gene") {  
        geneAfter <- extractGeneName( annotationEnhanced$attribute[i] )
        if ( elementStrand == "+" & previousElementStrand == "+") {
          regionName <- paste("Downstream (+) Of", geneBefore, "||| Upstream Of (+)", geneAfter, sep=" " )
          strand = "++"
        }
        else if ( elementStrand == "-" & previousElementStrand == "-") {
          strand ="--"
          regionName <- paste("Upstream (-) Of", geneBefore, "||| Downstream (-) Of", geneAfter, sep=" " ) 
        }
        else if ( elementStrand == "-" & previousElementStrand == "+") {
          strand = "+-"
          regionName <- paste("Downstream (+) Of", geneBefore, "||| Downstream (-) Of", geneAfter, sep=" " ) 
        }
        else if ( elementStrand == "+" & previousElementStrand == "-") {
          strand = "-+"
          regionName <- paste("Upstream (-) Of", geneBefore, "||| Upstream (+) Of", geneAfter, sep=" " ) 
          
        } else {
          print("Something wrong there!!!!")
        }
        
        previousElementStrand <- elementStrand
        nonCodingGeneID <- paste0("Name=NonCodingRegion_", nonCodingCount)
        nonCodingGeneAttribute <- paste0("Note=", regionName)
        newAttribute <- paste(nonCodingGeneID, nonCodingGeneAttribute, sep=";")
        # pick the non coding element(s) standing between previous and current gene
        elementPosition <- annotationEnhanced$start > previousElementStopCoordinates & annotationEnhanced$start < elementStartCoordinates
        # position of the noncoding element in the table
        noncodingIndex <- which(nonCodingPositions & elementPosition)
        # if there is a non coding element bwtween the two gene entries (noncodingIndex non empty)
        # substitute the strand and attribute
        if (length(noncodingIndex) >0 ){
          annotationEnhanced[ noncodingIndex, ]$strand <- strand
          annotationEnhanced[ noncodingIndex, ]$attribute <- newAttribute        
          nonCodingCount = nonCodingCount +1 
        }
        previousElementStopCoordinates <- elementStopCoordinates
        #annotationSlice <- insertRow(annotationSlice, newRow, dummyRowIndex)
        
        geneBefore <- geneAfter
      }
        i = i+1

    }
    # annotationEnhanced <- rbind(annotationEnhanced, annotationSlice)
  }
  
  chromIDsTableOfSNPs  <- unique(tableOfSNPs$CHROM)
  tableOfSNPsAnnotated <- data.frame()
  
  for (chromID in chromIDsTableOfSNPs) {
    
    sliceTableOfSNPs <- tableOfSNPs[tableOfSNPs$CHROM== chromID, ]
    SNPPos <- as.numeric(sliceTableOfSNPs$POS)
    sliceAnnotationEnhanced <- annotationEnhanced[ annotationEnhanced$seqname== chromID, ]
    
    sliceAnnotationEnhanced$start <- as.numeric(sliceAnnotationEnhanced$start)
    sliceAnnotationEnhanced$end <- as.numeric(sliceAnnotationEnhanced$end )
    
    # get the type of element(s) for each SNP 
    elements <- sapply(SNPPos, function(x) {
      positionOnTable <- which(sliceAnnotationEnhanced$start<=x & sliceAnnotationEnhanced$end>=x)
      elements <- sliceAnnotationEnhanced$feature[positionOnTable]
      # drop the NonAnnotRegion if other genomic elements are found along with it
      if (length(elements) > 1 && length(grep("NonCoding", elements))>=1 ) {
        index <- which(sapply("NonCoding", grepl, elements))
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
    
    sliceTableOfSNPs$Elements <- elements
    
    
    # Get the gene name for each SNP. Names are in the descreption of the "gene" field of the annotation file 
    geneNames <- sapply(SNPPos, function(x) {
      positionOnTable <- which(sliceAnnotationEnhanced$start<=x & sliceAnnotationEnhanced$end>=x )
      fractionTable <- sliceAnnotationEnhanced[positionOnTable, ]
      attributeString <- fractionTable[fractionTable$feature=="gene", ]$attribute
      geneNames <- extractGeneName(attributeString)
      # if the element is in non coding region, extract the non coding region Id form the 
      if (length(positionOnTable) <= 1) {
        if (fractionTable$feature== "NonCoding") {   
          attributeString <- fractionTable$attribute
          nameField <- strsplit(attributeString, ";", fixed=T)[[1]][1]
          geneNames <- strsplit(nameField, "=", fixed=T)[[1]][2]
        }
      }
      return(geneNames)
    })
    
    
    
    # go through the list and get the second element of the vector, that is the name of the gene
    #      geneNames <- sapply(c(1:length(geneNames)), function (x) {
    #                         geneNames[[x]] <- geneNames[[x]][2]
    #                       })
    
    # combine all elements of a vector into a single string, separated by comma
    geneNames <- sapply(c(1:length(geneNames)), function (x) {
      as.vector(paste(geneNames[[x]], collapse=", "))
    })
    
    # replace empty vectors with the description from "element" field 
    emptyGeneNameIndex <-  which(is.na(geneNames)) # which(sapply(c(1:length(geneNames)), function (x) {length(geneNames[[x]])==0 }))
    geneNames[emptyGeneNameIndex] <- elements[emptyGeneNameIndex]
    geneNames[is.na(geneNames)] <- "Non Annotated Region"
    geneNames <- unlist(geneNames)   
    
    sliceTableOfSNPs$GeneName <-  geneNames
    
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
    
    sliceTableOfSNPs$GeneCount <-  geneCount
    
    geneAttributes <- sapply(SNPPos, function(x) {
      positionOnTable <- which(sliceAnnotationEnhanced$start<=x & sliceAnnotationEnhanced$end>=x )
      fractionTable <- sliceAnnotationEnhanced[positionOnTable, ]
      attributeString <- fractionTable[fractionTable$feature=="CDS", ]$attribute
      geneAttributes <- extractGeneAttribute(attributeString)
      
      if (length(positionOnTable) <=1) {
        if ( fractionTable$feature == "NonCoding") {   
          start <- sliceAnnotationEnhanced$start[positionOnTable]
          end <- sliceAnnotationEnhanced$end[positionOnTable]
          distanceFromPrevious <- abs(x-start)
          distanceFromNext <- abs(x-end)
          attributeString <- fractionTable$attribute
          nameField <- strsplit(attributeString, ";", fixed=T)[[1]][2]
          geneAttributes <- strsplit(nameField, "=", fixed=T)[[1]][2]
          genes <- strsplit(geneAttributes, split= "|||", fixed=T)
          geneAttributes <- paste( distanceFromPrevious, "bp", genes[[1]][1], "AND", distanceFromNext, "bp", genes[[1]][2], sep=" ")
        }
      }
      return(geneAttributes)
    })
    
    sliceTableOfSNPs$GeneFunction <-  geneAttributes
    
    
    tableOfSNPsAnnotated <- rbind (tableOfSNPsAnnotated, sliceTableOfSNPs)
  }
  
  return(tableOfSNPsAnnotated)
}

# read in and process all samples (vcf files) 
ReadVCF <- function(listOfFiles) {
  # Make a data.frame variable to store all samples
  tableOfSNPs <- data.frame()
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
  # change the first column to keep only the contig ID in, so it matches with the annotation table
  contigs <- as.character(tableOfSNPs$CHROM)
  
#   contigIDs <- sapply(contigs, function(row) {
#     contig <- strsplit(row, split="|", fixed=T)[[1]][4]
#     strsplit(contig, split=".", fixed=T)[[1]][1]       
#   })
  
  tableOfSNPs$CHROM <- contigs
  return (tableOfSNPs)
}

# get the names of all vcf files in the directory and run the function that resds them 
# and then stores them in a data.frame
listOfFiles = grep("vcf$", dir(path), value=TRUE)
tableOfSNPs <- ReadVCF(listOfFiles)
listOfFilesSmallVariants <- grep("vcf1$", dir(path), value=TRUE)
tableOfSNPsSmallVariants <- ReadVCF(listOfFilesSmallVariants)
# remove the columns that are unnecessary
tableOfSNPsSmallVariants <- tableOfSNPsSmallVariants[ ,c("CHROM", "POS", "TYPE", "ALT", "ID")]
# rename the ID column of the small variants table to be identical to the other One
uniqueIDs <- sort(unique(tableOfSNPs$ID))
uniqueSmallVariantIDs <- sort(unique(tableOfSNPsSmallVariants$ID))
listOfIDs <- sapply(uniqueSmallVariantIDs, function(ID){  
  nOfRows<- nrow(tableOfSNPsSmallVariants[ tableOfSNPsSmallVariants$ID == ID, ])
  index <- match(ID, uniqueSmallVariantIDs)
  return(rep(uniqueIDs[index], nOfRows))
  
})
vectorOfIDs <- unlist(listOfIDs)
tableOfSNPsSmallVariants$ID <- vectorOfIDs

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
tableOfSNPs <- extraStatsTVC(tableOfSNPs)

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

### split multiple variants
tableOfSNPs <- SplitMultipleVariants(tableOfSNPs)

# run the annotation function and sort the table on the sample ID
tableOfSNPs <- annotateSNPTable(path, annotFile, tableOfSNPs)
tableOfSNPs <- tableOfSNPs[with(tableOfSNPs, order(ID)), ]
row.names(tableOfSNPs) <- NULL


hotSpotRegion <- vector()
homopolRegion <- vector()
listOfIDs <- unique(tableOfSNPs$ID)
for (id in listOfIDs){
  sliceOfTableOfSNPs <- tableOfSNPs[ tableOfSNPs$ID == id, ]
  chromosomeNames <- unique(sliceOfTableOfSNPs$CHROM)
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
    #     chromosomeSeq <- reference[[i]]
    #     delPos <- sliceOfTableOfSNPs[sliceOfTableOfSNPs$CHROM == i & sliceOfTableOfSNPs$TYPE == "DEL", ]$POS
    #     delRowIndex <- as.numeric(row.names(sliceOfTableOfSNPs[sliceOfTableOfSNPs$CHROM == i & sliceOfTableOfSNPs$TYPE == "DEL", ]))
    #     if (length(delPos)> 0 ){
    #       delInsideHomopolBool <- variantsInsideHomop(delPos, chromosomeSeq, 10, 4)
    #       # delInsideHomopolIndex <- delRowIndex[delInsideHomopolBool]
    #     }
    #     else {delInsideHomopolBool <- vector()}
    #     if(length(delInsideHomopolBool) > 0){
    #       delInsideHomopolIndex <- c(delInsideHomopolIndex, delRowIndex[delInsideHomopolBool])
    #     }
    
  }
  # a true /false vector to be implemented as a 
  rows <- row.names(sliceOfTableOfSNPs)
  hotSpotRegion <- c(hotSpotRegion, rows %in% hotSpotRowIndex)
  #homopolRegion <- c(homopolRegion, rows %in% delInsideHomopolIndex)
}

tableOfSNPs$HotSpotRegion <- hotSpotRegion
tableOfSNPs$DelInHomopolymerRegion <- homopolRegion

# A function to build a combined table with all samples from one experiment
# a list with all the samples Ids to combine into one table

CompileExperiment <- function(sampleIDs, tableOfSNPs, tableOfSNPsSmallVariants, experimentName){
  # build the table with the desired samples only
  tableForComparisson <- data.frame()
  tableForComparissonSmallVariants <- data.frame()
  for (sample in sampleIDs) {
    tableForComparisson <- rbind(tableForComparisson , tableOfSNPs[tableOfSNPs$ID ==  sample, ]) 
    tableForComparissonSmallVariants <- rbind(tableForComparissonSmallVariants , 
                                              tableOfSNPsSmallVariants[tableOfSNPsSmallVariants$ID ==  sample, ]) 
  }
  
  # create an extra colunt concisting of the Contig+POS+Type of SNP. 
  # If there are two different var types in the same pos, they are both included
  contigPosAlt <- paste(tableForComparisson$CHROM, tableForComparisson$POS, tableForComparisson$ALT, sep="_")
  tableForComparisson$contigPosAlt <- contigPosAlt
  contigPosAltSmallVar <- paste(tableForComparissonSmallVariants$CHROM, tableForComparissonSmallVariants$POS, 
                                 tableForComparissonSmallVariants$ALT, sep="_")
  tableForComparissonSmallVariants$contigPosAlt <- contigPosAltSmallVar
  
  uniqueVars <- unique(tableForComparisson$contigPosAlt)
  # sampleFreqList <- setNames(vector("data.frame", length( names(sampleIDs))), names(sampleIDs))
  # sampleDepthList <- setNames(vector("list", length( names(sampleIDs))), names(sampleIDs))
  # sampleQualList <- setNames(vector("list", length( names(sampleIDs))), names(sampleIDs))
  sampleFreqMatrix <-     data.frame(matrix(NA, nrow = length(uniqueVars), ncol = length(sampleIDs) ,
                                            dimnames = list(1:length(uniqueVars), names(sampleIDs))))
  sampleQualMatrix <-     sampleFreqMatrix
  sampleDepthMatrix <-    sampleFreqMatrix
  sampleHotSpotMatrix <-  sampleFreqMatrix
  
  # Go through the table and compile a vector for each sample for Quality, Frequency, Depth and hotspot. 
  # In the end, put the vector in a Data Frame
  i = 1
  for (sample in sampleIDs){
    sliceOfTableForComparisson <- tableForComparisson[ tableForComparisson$ID == sample, ]
    sliceOfTableOfSmallVarForComparisson <- tableForComparissonSmallVariants[ tableForComparissonSmallVariants$ID == sample, ]
    SNPFreq <- sapply(uniqueVars, function(variant){                        
      SNPFreq <- sliceOfTableForComparisson[ sliceOfTableForComparisson$contigPosAlt  == variant, ]$SNPFreq
      if (length(SNPFreq) != 0) {
        return(SNPFreq)
      } else if (variant %in% sliceOfTableOfSmallVarForComparisson$contigPosAlt) { 
        return (-1111)
      } else {
        return(NA)
      }
    })
    sampleFreqMatrix[ ,i] <- SNPFreq
    
    SNPDepth <- sapply(uniqueVars, function(variant) {                        
      SNPDepth <- sliceOfTableForComparisson[ sliceOfTableForComparisson$contigPosAlt  == variant, ]$Depth
      if (length(SNPDepth) != 0) {
        return(SNPDepth)
      } else {
        return(NA)
      }
    })
    sampleDepthMatrix[,i] <- SNPDepth
    
    SNPQual <- sapply(uniqueVars, function(variant) {                        
      SNPQual <- sliceOfTableForComparisson[ sliceOfTableForComparisson$contigPosAlt  == variant, ]$QUAL
      if (length(SNPQual) != 0) {
        return(SNPQual)
      } else {
        return(NA)
      }
    })
    sampleQualMatrix[,i] <- SNPQual
    
    hotspots <- sapply(uniqueVars, function(variant) {                        
      hotSpot<- sliceOfTableForComparisson[ sliceOfTableForComparisson$contigPosAlt  == variant, ]$HotSpotRegion
      if (length(hotSpot) != 0) {
        return(hotSpot)
      }else{ 
        return(NA)}
    })
    sampleHotSpotMatrix[,i] <- hotspots
    
    i = i+1
  }
  
  # go throught the hotspot list and call hotspot consensus
  hotSpotsConsensus <- apply(sampleHotSpotMatrix, 1, function(row){
    rowNARemoved <- row[!is.na(row)]
    uniqueRowNARemoved <- unique(rowNARemoved)
    if (length(uniqueRowNARemoved) ==1 ) {
      rowNARemoved[1]
    }else{
      "NoConsensus"
    }
  })
  
  compiledExperimentTable <- data.frame()
  for (var in uniqueVars) {
    compiledExperimentTable <- rbind(compiledExperimentTable, tableForComparisson[ tableForComparisson$contigPosAlt == var, ] )
  }
  uniqueVarIndex <- match(uniqueVars, tableForComparisson$contigPosAlt)
  experimentTableUniqueVar <- tableForComparisson[uniqueVarIndex, ]
  compiledExperimentTable <- experimentTableUniqueVar[ c("CHROM", "POS", "TYPE", "REF", "ALT")] #, "Elements", "GeneName") ]
  
  # Append Quality, Depth and frequency information for all samples. Name the columns accordingly 
  compiledExperimentTable <- cbind(compiledExperimentTable, sampleQualMatrix)
  qualityColIndex <- tail(1:ncol(compiledExperimentTable), ncol(sampleQualMatrix))
  colnames(compiledExperimentTable)[qualityColIndex] <- paste("Qual", names(sampleQualMatrix), sep="_")
  
  compiledExperimentTable <- cbind(compiledExperimentTable, sampleDepthMatrix)
  depthColIndex <- tail(1:ncol(compiledExperimentTable), ncol(sampleDepthMatrix))
  colnames(compiledExperimentTable)[depthColIndex] <- paste("Depth", names(sampleDepthMatrix), sep="_")
  
  compiledExperimentTable <- cbind(compiledExperimentTable, sampleFreqMatrix)
  freqColIndex <- tail(1:ncol(compiledExperimentTable), ncol(sampleFreqMatrix))
  colnames(compiledExperimentTable)[freqColIndex]  <- paste("Freq", names(sampleFreqMatrix), sep="_")
  
  
  # Go through the depth list and calculate the average depth per SNP. Remove non-numeric ("Absent") elements
  averageRowDepth <- apply(sampleDepthMatrix, 1, function(row) {
    numericRow <- as.numeric(row)
    numericRow <- numericRow[!is.na(numericRow)]
    mean(numericRow)
  })
  
  # Go through the depth list and calculate the average depth per SNP. Remove non-numeric ("Absent") elements
  averageRowQual <- apply(sampleQualMatrix, 1, function(row) {
    numericRow <- as.numeric(row)
    numericRow <- numericRow[!is.na(numericRow)]
    mean(numericRow)
  })
  
  freqAndAverageDepthMatrix <- sampleFreqMatrix
  # freqAndAverageDepthMatrix <- freqAndAverageDepthMatrix[ ,c(1:6) ]
  freqAndAverageDepthMatrix$AverageRowDepth  <- averageRowDepth
  
  # Append, averageRowDepth,  averageRowQual, Elements, GeneName, and Hotspot consensus  columns
  compiledExperimentTable$averageRowDepth <- averageRowDepth
  compiledExperimentTable$averageRowQual <- averageRowQual
  compiledExperimentTable$Element <- experimentTableUniqueVar$Elements
  compiledExperimentTable$GeneName <- experimentTableUniqueVar$GeneName
  compiledExperimentTable$HotSpotConsensus <- hotSpotsConsensus
  
  timoFilterAndSlope <- apply(freqAndAverageDepthMatrix, 1, function(row) {
    averageRowDepth <- as.numeric(tail(row, 1))
    # In the case multiple alernative alleles are reported, use the highest frequency
    rowFreq <- head(row, -1)
    rowFreq <- sapply(rowFreq, function(string){
      splitString <- as.numeric(unlist(strsplit(string, split=",", fixed=T)))
      if (splitString < 0 || is.na(splitString)) splitString <- NA
      if(all(is.na(splitString))){
        splitString
      } else{
        max(splitString)
      }
    })
    
    
    
    maxRowFreq <- max(rowFreq, na.rm = T)
    
    # Timo s filter                          
    diff <- round((max(rowFreq) - min(rowFreq)), digits=4)
    timoFilter <- "PASS"
    if (averageRowDepth < 100){
      if ( sum(is.na(rowFreq)) == 0 & diff < 0.1) {
        timoFilter <- paste("FAIL, Min-Max diff:", diff, "(is bellow 0.1)", sep=" ")
      }
    } else {
      if ( sum(is.na(rowFreq)) == 0 & diff < 0.1) {
        timoFilter <- paste("FAIL, Min-Max diff:", diff, "(is bellow 0.1)", sep=" ")
      }
    }
    
    # Minimium Maximum frequency filter
    if (maxRowFreq < 0.2) {
      maxFreqInRow <- "FAIL"
    } else {
      maxFreqInRow <- "PASS"
    }
    
    # Coverage filter
    coverage <- "WithinLimits"
    if ( averageRowDepth < 100) {
      coverage <- "LowCoverage"
    } else if ( averageRowDepth > 500 ) {
      coverage <- "HighCoverage"
    }
    
    freqWO <- tail(rowFreq, 1) # sample WithOut Antibiotic
    freqWO[is.na(freqWO)] <- 0 # absent SNPs should be substituted with 0
    rowFreq <- head(rowFreq, -1)
    # rowFreq <- rowFreq[which(!is.na(rowFreq))]
    underSelection <- rowFreq[2:length(rowFreq)] # Samples under selection
    pattern <- NA
    slope <- NA
    
    # Slope
    # if there are 2 or more genes 
    if (sum(!is.na(underSelection)) >= 2) {
      #print(paste("average Ant: ", averageAntibTreated, sep=""))
      averageAntibTreated <- mean(underSelection[2:length(underSelection)], na.rm = T)
      slope <- round(lm(underSelection ~ c(1:length(underSelection)))$coefficients[[2]], digits=4)
      freqAverTreatedMinusFreqWO <- averageAntibTreated - freqWO
      # Check if data follows pattern
      if (slope >= 0 & freqAverTreatedMinusFreqWO < 0){
        pattern <- "GoesUpStaysUp"                                      
      } else if (slope >= 0 & freqAverTreatedMinusFreqWO > 0){
        pattern <- "GoesUpThenDown"  
      } else if (slope < 0 & freqAverTreatedMinusFreqWO > 0) {
        pattern <- "GoesDownStaysDown"
      } else if (slope < 0 & freqAverTreatedMinusFreqWO < 0) {
        pattern <- "GoesDownThenUp"
      } else {
        pattern <- "NoPattern"
      }
      
    } else if(!is.na(tail(underSelection, 1))) {
      timoFilter <- "SNP in Last Selection Day Only"
    } else {
      timoFilter <- "SNP in Single Sample"
    }
    
    
    # Essential Mutation
    dayZero <- underSelection[1]
    if (is.na(dayZero)) { dayZero <- 0 }
    underAntibioticAndWO <- c(underSelection[2:length(underSelection)], freqWO)                           
    allAboveThershold <- all(underAntibioticAndWO > 0.8, na.rm = T)
    if ( !all(!is.na(underAntibioticAndWO) )) {allAboveThershold <- F}
    if (dayZero < 0.2 & allAboveThershold ) {
      pattern <- "EssentialMutation"
    }
    
    list(coverage, timoFilter, maxFreqInRow, slope, pattern)
  })
  
  timoFilterAndSlopeMatrix <- as.data.frame(matrix(unlist(timoFilterAndSlope), nrow=length(timoFilterAndSlope), byrow=T))
  colnames(timoFilterAndSlopeMatrix ) <- c("Coverage", "Timo's Filter", "MaxRowFreq>0.2", "Slope, Antibiotic Treated only", "Pattern")
  # freqAndAverageDepthMatrix$TimoFilter <- timoFilter
  compiledExperimentTable <- cbind(compiledExperimentTable, timoFilterAndSlopeMatrix)
  
  # print the table
  write.table(compiledExperimentTable, paste(path, "TableFor", experimentName, ".txt", sep=""), sep="\t", row.names=FALSE) 
}

######### ******* ######## 
# compile the tables for all experiments
# CAZ256
sampleIDsForCAZ256 = list(WT_Zero = listOfIDs[1],
                          Broth_WT_D30 = listOfIDs[2],
                          Broth_CAZ256_D6 = listOfIDs[4],
                          Broth_CAZ256_D9 = listOfIDs[5],
                          Broth_CAZ256_D12 = listOfIDs[6],
                          Broth_CAZ256_D21 = listOfIDs[7],
                          Broth_CAZ256_WO_D15 = listOfIDs[8]
)
CompileExperiment(sampleIDsForCAZ256, tableOfSNPs, tableOfSNPsSmallVariants, "CAZ256_test")

# MPN256
sampleIDsForMPN256 = list(WT_Zero = listOfIDs[1],
                          Broth_WT_D30 = listOfIDs[2],
                          Broth_MPN256_D5 = listOfIDs[9],
                          Broth_MPN256_D10 = listOfIDs[10],
                          Broth_MPN256_D25 = listOfIDs[11],
                          Broth_MPN256_WO_D15 = listOfIDs[12]
)
CompileExperiment(sampleIDsForMPN256, tableOfSNPs, tableOfSNPsSmallVariants, "MPN256_test")

# PT512
sampleIDsForPT512 = list(WT_Zero = listOfIDs[1],
                         Broth_WT_D30 = listOfIDs[2],
                         Broth_PT512_D5 = listOfIDs[14],
                         Broth_PT512_D12 = listOfIDs[15],
                         Broth_PT512_D22 = listOfIDs[16],
                         Broth_PT512_WO_D15 = listOfIDs[17]
)
CompileExperiment(sampleIDsForPT512, tableOfSNPs, tableOfSNPsSmallVariants, "PT512_new")

# TBM32
sampleIDsForTBM32 = list(WT_Zero = listOfIDs[1],
                         Broth_WT_D30 = listOfIDs[2],
                         Broth_TBM32_D6 = listOfIDs[18],
                         Broth_TBM32_D13 = listOfIDs[19],
                         Broth_TBM32_D30 = listOfIDs[20],
                         Broth_TBM32_WO_D15 = listOfIDs[21]
)
CompileExperiment(sampleIDsForTBM32, tableOfSNPs, tableOfSNPsSmallVariants, "TBM32_new")

# CIP128
sampleIDsForCIP128 = list(WT_Zero = listOfIDs[1],
                          Broth_WT_D30 = listOfIDs[2],
                          Broth_CIP128_D6 = listOfIDs[22],
                          Broth_CIP128_D12 = listOfIDs[23],
                          Broth_CIP128_D21 = listOfIDs[24],
                          Broth_CIP128_WO_D15 = listOfIDs[25]
)
CompileExperiment(sampleIDsForCIP128, tableOfSNPs, tableOfSNPsSmallVariants, "CIP128_new")

# TBM1024
sampleIDsForTBM1024 = list(WT_Zero = listOfIDs[1],
                           Evans_WT_D30 = listOfIDs[3],
                           Evans_TBM1024_D14 = listOfIDs[26],
                           Evans_TBM1024_D30 = listOfIDs[27],
                           Evans_TBM1024_WO_D15 = listOfIDs[28]
)
CompileExperiment(sampleIDsForTBM1024, tableOfSNPs, tableOfSNPsSmallVariants, "TBM1024_new")

# CIP256
sampleIDsForCIP256 = list(WT_Zero = listOfIDs[1],
                          Evans_WT_D30 = listOfIDs[3],
                          Evans_TCIP256_D11 = listOfIDs[29],
                          Evans_CIP256_D19 = listOfIDs[30]
)
CompileExperiment(sampleIDsForCIP256, tableOfSNPs, tableOfSNPsSmallVariants, "CIP256_new")

### Draw some lines 
# experimentPoints <- names(compiledExperimentTable)[c(8, 10, 12, 14, 16)]
slopes <- vector()
for (i in 1:nrow(compiledExperimentTable)){
  frequences <- as.numeric(compiledExperimentTable[,c(8, 10, 12, 14, 16)][i,])
  frequences <- frequences[!is.na(frequences)]
  xaxis <- c(1:length(frequences))
  if (length(frequences)>2){
    plot(xaxis, frequences, ylim=c(0, 1), main = i)
    abline(lm(frequences~xaxis ), col="red")
    slopes <- c(slopes, lm(frequences~xaxis)$coefficients[[2]])
  }
}

hist(slopes, breaks = 200)




# write the data to a tab separated file
# put output in separate files
sampleIDs = list(Control =  listOfIDs[1:3], #c("BS_variants_01.vcf", "BS_variants_02.vcf", "BS_variants_03.vcf" ),
                 CAZ256 =   listOfIDs[4:8],   #c("BS_variants_04.vcf", "BS_variants_05.vcf", "BS_variants_06.vcf", "BS_variants_07.vcf", "BS_variants_08.vcf" ),
                 MPN256 =   listOfIDs[9:13],  # c("BS_variants_09.vcf", "BS_variants_10.vcf", "BS_variants_11.vcf", "BS_variants_12.vcf", "BS_variants_13.vcf" ),
                 PT512 =    listOfIDs[14:17], # c("BS_variants_14.vcf", "BS_variants_15.vcf", "BS_variants_16.vcf", "BS_variants_17.vcf"),
                 TBM32 =    listOfIDs[18:21], # c("BS_variants_18.vcf", "BS_variants_19.vcf", "BS_variants_20.vcf", "BS_variants_21.vcf"),
                 CIP128 =   listOfIDs[22:25], # c("BS_variants_22.vcf", "BS_variants_23.vcf", "BS_variants_24.vcf", "BS_variants_25.vcf"),
                 TBM1024 =  listOfIDs[26:28], # c("BS_variants_26.vcf", "BS_variants_27.vcf", "BS_variants_28.vcf" ),
                 CIP256 =   listOfIDs[29:30]  #c("BS_variants_29.vcf", "BS_variants_30.vcf")
)

experimentNames <- names(sampleIDs)
for (experiment in experimentNames){
  IDs <- unlist(sampleIDs[experiment])
  sliceOfTableOfSNPs <- data.frame()
  for (id in IDs){
    sliceOfTableOfSNPs <- rbind(sliceOfTableOfSNPs, tableOfSNPs[tableOfSNPs$ID==id, ])
  }
  write.table(sliceOfTableOfSNPs, paste(path, "variantsFor_", experiment, "_Samples.txt", sep=""), sep="\t", row.names=FALSE)
  
}

############### *************** ################
# Dead Code

# else if (averageRowDepth < 300) {
#   if(diff < 0.02) {
#     timoFilter <- paste("FP, diff too small:", diff, sep=" ")
#   } else{
#     timoFilter <- "TrustWorthy"
#     slope <- lm(underSelection ~ c(1:length(underSelection)))$coefficients[[2]]
#     if (slope > 0 & freqAverTreatedMinusUntreated > 0){
#       pattern <- "FollowsPattern"                                      
#     } else if (slope < 0 & freqAverTreatedMinusUntreated < 0){
#       pattern <- "FollowsPattern"  
#     } else {
#       pattern <- "NoPattern"
#     }
#   }
# }

# # convert SNP depth list into a data frame
# sampleDepthMatrix <- data.frame(matrix(unlist(sampleDepthList), nrow= length(sampleDepthList[[1]]), byrow=F))
# # convert SNP frequency list into a data frame
# sampleFreqMatrix <- data.frame(matrix(unlist(sampleFreqList), nrow= length(sampleFreqList[[1]]), byrow=F))

# Go through all samples present in sub-experiment
# columnIndex = 5
# for (i in 1:length(sampleFreqList)) {
#   columnIndex = columnIndex + 1
#   compiledExperimentTable[ ,columnIndex] <- sampleFreqList[[i]]
#   colnames(compiledExperimentTable)[columnIndex] <- paste("Freq_", names(sampleFreqList)[i], sep="")
#   columnIndex = columnIndex + 1
#   compiledExperimentTable[ ,columnIndex] <- sampleDepthList[[i]]
#   colnames(compiledExperimentTable)[columnIndex] <- paste("Depth_", names(sampleDepthList)[i], sep="")
#   columnIndex = columnIndex + 1
#   compiledExperimentTable[ ,columnIndex] <- sampleQualList[[i]]
#   colnames(compiledExperimentTable)[columnIndex] <- paste("QualScore_", names(sampleDepthList)[i], " |||", sep="")
# }
# write.table(tableOfSNPs, paste(path, "variantsCompiledHotspotsHomopolymersAnnotation.txt", sep=""), sep="\t", row.names=FALSE)

# if (averageRowDepth < 20) {
#   timoFilter<- "LowCoverage"
# } else if (averageRowDepth < 100) {
#   if(diff > 0.03 & SNPsInAllSamplesBoolean) {
#     timoFilter <- "TrustWorthy"
#     slope <- round(lm(underSelection ~ c(1:length(underSelection)))$coefficients[[2]], digits=4)
#     # print(paste("slope", slope, sep=" "))
#     # print(paste("FreqAverTreatedMinusUntreaded ", freqAverTreatedMinusUntreated, sep=""))
#     if (slope > 0 & freqAverTreatedMinusUntreated > 0){
#       pattern <- "FollowsPattern"                                      
#     } else if (slope < 0 & freqAverTreatedMinusUntreated < 0){
#       pattern <- "FollowsPattern"  
#     } else {
#       pattern <- "NoPattern"
#     }
#   } else {
#     timoFilter <- paste("Max-min frequency diff:", diff, "(is bellow 0.03)", sep=" ")
#   }
# }  else if (averageRowDepth < 500) {
#   if(diff > 0.02 & SNPsInAllSamplesBoolean) {
#     timoFilter <- "TrustWorthy"
#     slope <- round(lm(underSelection ~ c(1:length(underSelection)))$coefficients[[2]], digits=4)
#     print(paste("slope", slope, sep=" "))
#     print(paste("FreqAverTreatedMinusUntreaded ", freqAverTreatedMinusUntreated, sep=""))
#     if (slope > 0 & freqAverTreatedMinusUntreated > 0){
#       pattern <- "FollowsPattern"                                      
#     } else if (slope < 0 & freqAverTreatedMinusUntreated < 0){
#       pattern <- "FollowsPattern"  
#     } else {
#       pattern <- "NoPattern"
#     }
#   } else{
#     timoFilter <- paste("Max-min frequency diff:", diff, "(is bellow 0.02)", sep=" ")
#   }
# } else {
#   timoFilter <- "HighCoverage"
# }