#! /bin/bash/Rscript

### count_blat.R by Clément Goubert ###
### Part of the High-Troughput Sequencing Transposon Display Analysis Suite ###

### Dependency: Retreive_R2_reads.sh

### Description: This scripts is a dependency of the main script Retreive_R2_reads.sh

Args = commandArgs()
blat = Args[6]

tblat=read.table(blat)
countblat=as.data.frame(tapply(tblat$V1, tblat$V1, length))

write.table(countblat, "countblat.tmp", quote=F)