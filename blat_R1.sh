#! /bin/bash

### Blat_R1.sh by Clément Goubert ###
### Part of the High-Troughput Sequencing Transposon Display Analysis Suite ###

### Description: This script filters the R1 reads for presence of the right TE family using a parallelized implementation of BLAT.

###	Input:		 Cleaned R1 reads (1 file / individual in 1 folder), TE consensuses (all families)
###	Output:		 .psl (Blat) file (1/individual) to use then with the script Sorting_blat.sh

###PARAMETERS:
#$1 input folder in fastq (R1)
#$2 outdir
#$3 TE database
#$4 minIdentity Blat



###CONFIG###
fastq_to_fasta="fastq_to_fasta"
blat="/home/clement/Downloads/blat"
parallel="~/pbil-panhome/parallel-20140622/src/parallel"
############

mkdir -p $2
mkdir -p $2/Fasta_R1
mkdir -p $2/Blat_out

cd $1

FILE_LIST="*R1*.fastq"

i=$((0))

for file in ${FILE_LIST}
do


i=$((1+$i))

echo ""
echo "processing file #$i/140 : $file"


$fastq_to_fasta -i $file -o $2/Fasta_R1/$file.fasta -Q33

echo "blat of $file..."

#parallel blat and send stdout to a blackhole
echo "cat  $2/Fasta_R1/$file.fasta | $parallel --round-robin --pipe --recstart '>' '$blat -noHead -minIdentity=$4 $3 stdin >(cat) >&2' > $2/Blat_out/$file""_blatout_""$4.psl #2>-" > blat.sh
chmod +x blat.sh

./blat.sh

echo "$file done"

done
