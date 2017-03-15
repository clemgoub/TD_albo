#! /bin/bash

### Retreive_R2_reads.sh by Clément Goubert ###
### Part of the High-Troughput Sequencing Transposon Display Analysis Suite ###

### Description: This script re-map the raw R2 reads over the high-quality loci produced by locus_factory.sh

###	Input:		 Reference loci produced by locus_factory. To produce this file, concatenate each fasta file from locus_factory (one per TE family, "All.*.clustering" files) containing the reference sequence of each cluster.
###				 Before concatenation, be careful to add the TE name in the sequence header to allow sorting in the next step.
###	Output:		 .out table containing the best hit per read on one locus (one file per individual)


### Directions: Execute in the sampled raw R2 reads file

###PARAMETERS:
### $1 = reference loci file (concatednaed output of locus_factory.sh [All.*.clustering files])     ###
### $2 = outdir                ###



FILE_LIST=""

mkdir -p $2

for file in ${FILE_LIST}
do

rm retreive""$file"".pbs

echo "#! /bin/bash" >> retreive""$file"".pbs

echo "#PBS -N ret_$file" >> retreive""$file"".pbs
echo "#PBS -q q1hour" >> retreive""$file"".pbs
echo "#PBS -l nodes=1:ppn=8" >> retreive""$file"".pbs
echo "#PBS -l mem=5gb" >> retreive""$file"".pbs
echo "#PBS -e /pandata/goubert/TD_albo/ret""$file"".err" >> retreive""$file"".pbs
echo "#PBS -o /pandata/goubert/TD_albo/ret""$file"".out" >> retreive""$file"".pbs

#run blat in parallel on one file
echo "cp /pandata/goubert/TD_albo/Raw_data/done/$file ." >> retreive""$file"".pbs
echo "/panhome/goubert/blat $1 $file -out=blast8 -noHead $file""_blat.out" >> retreive""$file"".pbs

# sort le blat.out keeping one best hit per read
echo "sort -k1,1 -k12,12nr $file""_blat.out | sort -u -k1,1 > $2/$file""_besthits_reads.blat.out" >> retreive""$file"".pbs

qsub retreive""$file"".pbs

done