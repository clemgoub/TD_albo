### Retreived reads counts report ###

### Retreive_from_blat.sh by Clément Goubert ###
### Part of the High-Troughput Sequencing Transposon Display Analysis Suite ###

### Description: This script build the raw matrices for each TE after the remaping of the raw R2 reads in the previous step (Retreive_R2_reads.sh)

###	Input:		 Reference loci produced by locus_factory.sh of one TE family (e.g. "All.IL1.clustering") containing the reference sequence of each cluster.
###				 sorted output from Retreive_R2_reads.sh
###	Output:		 Raw matrix of reads per loci per indivifual for one TE family

### Dependency: count_blat.R

### Directions: Execute in the sampled raw R2 reads file

###PARAMETERS:
# $1 = ref (full path) all_ET_clustering.bak.clstr
# $2 = ET

echo ""
echo "###################"
echo "# Retreive  readS #"
echo "# COVERAGE REPORT #"
echo "###################"
echo ""

grep '*' $1 | awk '{print $1" "$3}' | sed 's/>//g' | sed 's/\.\.\.//g' | sort -k2,2 > refs # sort the reference loci file and holds the headers

FILE_LIST=*besthits_reads.blat.out

for file in ${FILE_LIST}
do

grep "$2" $file > use_file # keeps only the reads mapping on the TE on interest

echo "$file..."

echo -ne '[#####------------------]   (33%)\r'
file_name="$( ls $file | sed 's/_/ /g' | awk '{print $1}')"

awk '{print $2}' use_file > file.R.tmp # create the file for the counts of reads per loci


Rscript ../../scripts/count_blat.R file.R.tmp # counts the reads per loci with R


echo -ne '[#############----------]   (66%)\r'
grep -v 'tapply' countblat.tmp | sort -k1,1 | sed "s/_$2//g" > $file_name""_counts # fetch the output
rm file.R.tmp

join -a1 -12 -21 -e 0 -o 2.2 refs $file_name""_counts > $file_name""_count_final # add the loci name to the counts
rm $file_name""_counts


echo "$file_name" > ind
cat ind $file_name""_count_final > $file_name""_matrix # creates the individual column

rm ind

echo -ne '[#######################]   (100%)\r'
echo -ne '\n'

done

echo "building matrix..." # assemble the final matrix

paste *_matrix > $2""_matrice_yeah 
cat <(echo "locus read") refs > refs2
paste refs2 $2""_matrice_yeah > $2""_matrix_R

rm refs
