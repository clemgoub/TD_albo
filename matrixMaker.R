#! /usr/bin/Rscript

### matrixMaker.R by Clément Goubert ###
### Part of the High-Troughput Sequencing Transposon Display Analysis Suite ###

### Description: # This R function is to copy in paste into a R console (or to source into R). It allows to filter the raw coverage matrix produces by Retreive_from_blat.sh

###	Input:		 Raw coverage matrix for one TE family
###	Output:		 Filtered 1/0 matrix for one TE family, in a R object.

###PARAMETERS:
# path - path to the matrix file
# ind - min number of individual that must share a loci to keep it
# mcov - minimum quantile of coverage to keep a loci (e.g. 0.05 will discard the 5% of loci the less covered in average)
# sd_seuil - maximum quantile of standard deviation of coverage (e.g. 0.05 will discard loci in the extreme 5% of s.d. of the mean coverage between individuals)

matrixMaker=function(path, ind, mcov, sd_seuil){

mat=read.table(path, h=T)
mat3=as.matrix(mat)

use_i=(rowSums((mat3>0)+0) >= ind) ### select loci to nb ind

mat3=mat3[use_i,]

#filtering on min coverage or O
mat3[mat3 == 0] <- NA
mean_cov=quantile(rowMeans(mat3, na.rm=T), mcov, na.rm=T)
std.dev=apply(mat3, 1, function(x) sqrt(var(na.omit(x))))
sd=( std.dev <= quantile(std.dev, sd_seuil) )
use_cov=(rowMeans(mat3, na.rm=T) >= mean_cov )

print(quantile(rowMeans(mat3, na.rm=T),c(0.5,0.75,0.9,0.750, mcov), na.rm=T))
print(quantile(std.dev, c(0.5,0.750,0.975, sd_seuil)))
print(table(sd))
print(table(use_cov))

mat3[is.na(mat3)] <- 0
mat4=mat3[use_cov & sd,]

mat4=mat4[!is.na(rowSums(mat4)),]


#turning into binary and filtering for individual number
binM2=((mat4>0)+0)
# head(binM2)

# use_i=(rowSums((mat4>0)+0) >= ind)

f_matrix=binM2

unuse_fix=(rowSums(f_matrix)==length(f_matrix[1,]))

f_matrix[!unuse_fix,]

}