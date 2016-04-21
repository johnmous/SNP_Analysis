library(ShortRead)
fq <- readFastq(file)
fq <- fq[!duplicated(sread(fq)]
writeFastq
writeXStringset(fq,new_filename,format="fastq")