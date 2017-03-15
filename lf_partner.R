#! /bin/bash/Rscript

### lf_partner.R by Clément Goubert ###
### Part of the High-Troughput Sequencing Transposon Display Analysis Suite ###

### Dependency: locus_factory.sh

### Description: This scripts is a dependency of the main script locus_factory.sh


Args = commandArgs()
rloc = Args[6]

reads=read.table(rloc)
countreads=as.data.frame(tapply(reads$V1, reads$V1, length))

write.table(countreads, "counts.tmp", quote=F)