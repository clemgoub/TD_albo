#! /bin/bash

### Blat_R1.sh by Clément Goubert ###
### Part of the High-Troughput Sequencing Transposon Display Analysis Suite ###

### Description: This script sort the blat out from R1 read and go picking the cleaned (UrQt) R2 mate using the threshold described in the main manuscript. 
###				 It then put each R2.fasta file in one folder per TE.

###	Input:		 .psl (blat out) of the R1 reads (from Blat_R1.sh script), cleaned and trimmed R2 reads. 
###	Output:		 R2 reads for clustering (.fasta), 1 folder / TE, 1 file / individual / TE

###PARAMETERS:
# $1 = R2 files folder. 1 file per individual, fasta format

# EXECUTE in BLAT_OUT folder (where the .psl are)!!! 

cd ./

mkdir -p R2_mates_blat_OK
mkdir -p R1_headers_blat_OK


###################
### MAIN LOOP #####
###################

#compt to zero
i=$((0))

FILE_LIST="*.psl" #

for file in ${FILE_LIST}

do

#increment compt
i=$((1+$i))

echo ""
echo "processing file #$i/140 : $file"

#define the names of the file
file_name="$( ls $file | sed 's/_/ /g' | awk '{print $1}')"

echo "$file_name"
#sort 1 hit per read, and select according to threshold. Outputs only header for use in the read selection script
sort -u -k10,10 $file | grep -v 'RTE3' | awk '{if ($11 > 30 && (($13)-($12))/($11) >= 0.90 && $18 < 2) print $10 > "'R1_headers_blat_OK/$file_name'_"$14"_blat_OK_R1_headers"}'


echo ""
echo "picking R2 mate..."

#pick and format R2.fasta mate file for read picking
filetemp="$( find $1/*$file_name*.fasta)" 
awk '{print $1}' $filetemp > tmp

################
# second loop: # creates a fasta file per TE for this individual. Reads are put in R2_mates_blat_OK.fasta folder
################

TE="Lian1 IL1 L2B RTE5 RTE4"

for subfile in ${TE}
do

mkdir -p R2_mates_blat_OK/$subfile

echo "$subfile"

perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' R1_headers_blat_OK/$file_name""_""$subfile""_blat_OK_R1_headers tmp > R2_mates_blat_OK/$subfile/$file_name""_""$subfile""_R2_blat_OK.fasta

done

#cleaning
rm tmp

echo ""
echo "$file done"

done