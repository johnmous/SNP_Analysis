# Test reading BAM files with Pysam
# Author: Ioannis Moustakas

import pysam
import matplotlib.pyplot as plt
import numpy as np
import pylab
import re
import csv


# Save file locations
BAMFile = "/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1095-Nadine_Handel/MAD1095-P002-Seq_E_Coli/MAD1095-P002-E001_2014_14x_gDNA_Seq_svleeuw1/Results/CLC/Data/BC_01/BC_01_aligned.sorted.bam.bam"
textFile = "../variant_calling/assembled_WT_ref/tmap/BC01_aligned_sort.mpileup.txt"
reference = "/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1095-Nadine_Handel/MAD1095-P002-Seq_E_Coli/MAD1095-P002-E001_2014_14x_gDNA_Seq_svleeuw1/Scratch/ec_K12_MG1655.fasta"

# Calculate the pileup from the Bam file
#pileup = pysam.mpileup("-f", 
#                       "/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1095-Nadine_Handel/MAD1095-P002-Seq_E_Coli/MAD1095-P002-E001_2014_14x_gDNA_Seq_svleeuw1/Scratch/ec_K12_MG1655.fasta", 
#                       BAMFile)
#                       
# A function to plot 
def Plot(frequencyList):
    count = 0
    threshold = 0 
    thresholdList = []
    percentageList= []
    while threshold < 0.10:
        valuesAboveThreshold = []
        for freq in frequencyList:
            if freq > threshold:
                valuesAboveThreshold.append(freq)
                count +=1 
        frequencyList = valuesAboveThreshold
        thresholdList.append(threshold)
        percentageList.append(len(valuesAboveThreshold))
        threshold += 0.001
    
       
    plt.plot(thresholdList, percentageList)    
    
# A function with a file name as input and relevant statistics as output
def BamFileStats (fileName):
                       
    # read it the text file with the saved pileup
    with open(fileName, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter="\t")
        pileup = [[e for e in r] for r in reader]
    
    frequencyOfIndels = []
    frequencyOfSNPs = []
    numberOfIndels = 0
    numberOfSNPs = 0
    numberOfLines = len(pileup)
    for line in pileup:
        elements = line
        #elements = line.split()
        # if the line has less than 6 elements, skip this iteration
        if len(elements) != 6:
            continue
        #elements = line
        #print(elements[1])
        depth = int(elements[3])
        mpileup = elements[4]
        #print(mpileup)
        # A regex to fish out the integers in the mpileup string 
        integers = re.findall('[\+-][0-9]+', mpileup)
        sumIndelLength = 0
        # If integers a non empty list, calculate the sum of their absulute values        
        if integers:
            #print(integers)
            for i in integers:
                sumIndelLength += abs(int(i))
            # Indels frequency, as a fraction of bases that differ from reference to total depth 
            percentIndel = sumIndelLength/depth
            frequencyOfIndels.append(percentIndel)
            numberOfIndels += 1
        
        letters = re.findall('[\.,\^\$]([ATGCNatgcn]+)\.', mpileup)  
        #print(letters)
        sumSNPLength = 0
        if letters:
            for letter in letters: 
                sumSNPLength += len(letter)            
                
            SNPFrequency = (sumSNPLength-sumIndelLength)/depth
            frequencyOfSNPs.append(SNPFrequency)
            if SNPFrequency>0:
                numberOfSNPs += 1
    return(numberOfIndels, numberOfSNPs, numberOfLines, frequencyOfIndels, frequencyOfSNPs)

numberOfIndels, numberOfSNPs, numberOfLines, frequencyOfIndelsYan, frequencyOfSNPsYan = BamFileStats("../variant_calling/assembled_WT_ref/tmap/BC01_aligned_sort.mpileup.txt")

#print(percentOfIndels[1:100])    
# number and percentage of positions where an indel or SNP is found    
print("Number of Indels: ", numberOfIndels)
print("Percentage of indels: ", numberOfIndels/numberOfLines)
print("Number of SNPs: ", numberOfSNPs)
print("Percentage ofSNPs: ", numberOfSNPs/numberOfLines)

  
numberOfIndels, numberOfSNPs, numberOfLines, frequencyOfIndelsNad, frequencyOfSNPsNad = BamFileStats("/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1095-Nadine_Handel/MAD1095-P002-Seq_E_Coli/MAD1095-P002-E001_2014_14x_gDNA_Seq_svleeuw1/Results/CLC/Data/BC_01/NadineBC01Mpileup.txt")

#print(percentOfIndels[1:100])    
# number and percentage of positions where an indel or SNP is found    
print("Number of Indels: ", numberOfIndels)
print("Percentage of indels: ", numberOfIndels/numberOfLines)
print("Number of SNPs: ", numberOfSNPs)
print("Percentage ofSNPs: ", numberOfSNPs/numberOfLines)  

Plot(frequencyOfIndelsYan)
Plot(frequencyOfIndelsNad)
plt.legend(["Yangang", "Nadine"])
plt.title("Indels")
plt.show()

Plot(frequencyOfSNPsYan)
Plot(frequencyOfSNPsNad)
plt.legend(["Yangang", "Nadine"])
plt.title("SNPs")
plt.show()

#pylab.savefig("NaadineIndelFreqGraph.png")         
#print("Number of SNPs with freq above 0.035: ", count)
  

    
#"-f", " ../raw_data/Wild-type_assembly.fasta",

#def interogateBAM(BAMFile):
#    samfile = pysam.AlignmentFile(BAMFile, "rb")
#    perfectAlignment = 0 
#    almostPerfect = 0
#    
#    reads = 0
#    insertions = []
#    deletions = []
#    SNPs = []
#    allVariants = []
#    for read in samfile.fetch():
#        readLength = read.infer_query_length()
#        # if NM tag is 0, then read is a perfect match to 
#        if read.opt("NM") == 0:
#            perfectAlignment +=1
#        # number of rts perfectly (1 or less mismatch)    
#        if read.opt("NM") <= 1:
#            almostPerfect +=1
#        else:
#            ins=0
#            dels=0
#            for tuple in read.cigartuples:
#                if tuple[0] == 1:
#                    ins += tuple[1]
#                if tuple[0] == 2:
#                    dels += tuple[1]
#            insertions.append(ins/readLength)
#            deletions.append(dels/readLength)  
#            SNPs.append(read.opt("NM")-ins-dels)
#            allVariants.append(read.opt("NM"))    
#        reads += 1
#    
#    samfile.close()
#    return (insertions, deletions, SNPs, allVariants, reads, perfectAlignment, almostPerfect)
#    
#insertions, deletions, SNPs, allVariants, reads, perfectAlignment, almostPerfect = interogateBAM("readsOnATtranscriptTAIR10_No_rRNA.sorted.bam")
#
#print("Perfect Aligned Reads, % of all: ", (perfectAlignment/reads)*100)
#print ("Reads with 1 or 0 mismatches, % of all:", (almostPerfect/reads)*100)
#print("Reads: ", reads)
#
#hist, bins = np.histogram(insertions, bins=100)
#width = 0.7 * (bins[1] - bins[0])
#center = (bins[:-1] + bins[1:]) / 2
#plt.bar(center, hist, align='center', width=width, log=True)
#pylab.savefig("InsertionsHistNadine.png")
#
#hist, bins = np.histogram(deletions, bins=100)
#width = 0.7 * (bins[1] - bins[0])
#center = (bins[:-1] + bins[1:]) / 2
#plt.bar(center, hist, align='center', width=width, log=True)
#pylab.savefig("DeletionsHistNadine.png")
#
#hist, bins = np.histogram(SNPs, bins=100)
#width = 0.7 * (bins[1] - bins[0])
#center = (bins[:-1] + bins[1:]) / 2
#plt.bar(center, hist, align='center', width=width, log=True)
#pylab.savefig("SNPsHistNadine.png")
# 
#maxMismatches = max(mismatchesInRead)
#hist, bins = np.histogram(mismatchesInRead, bins = maxMismatches)
#width = 0.7 * (bins[1] - bins[0])
#center = (bins[:-1] + bins[1:]) / 2
#plt.bar(center, hist, align='center', width=width, log=True)
#pylab.savefig("AllVariantsHistNadine.png")