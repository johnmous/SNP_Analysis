# Title: Use the mpileup ofa BAM file to draw graphs
# Author: Ioannis Moustakas, i.moustakas@uva.nl

# read the text file where the samtools mpileup is saved (as tab delimited)
mpileupFile = "/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1117-Ter_Kuile/MAD1117-P001-gDNA_sequencing/MAD1117-P001-E001_2014_32x_gDNA_Yanfang_svleeuw1/Scratch/Oskar/variant_calling/assembled_WT_ref/tmap/BC01_aligned_sort.mpileup.txt"

# mind the # and " characters in the file that mess up the table read
mpileup <- read.delim(mpileupFile, sep="\t", fill=T, comment.char="", stringsAsFactors=FALSE, quote="", header=F)

depth = mpileup[,4]
mpileupStr = mpileup[,5]
