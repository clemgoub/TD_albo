#! /bin/bash

#arguments $1 is the TE family

lof=$(ls /pandata/goubert/TD_albo/R2_t10-OK/RTE4_OK/*_$1""_R2_blat_OK.fasta_no_overlap_good_adapt_trimmed.fasta)

for file in $lof
do

ind=$(echo "$lof" | sed 's/OK\//\t/g;s/_/\t/g' | cut -f 6)

bwa aln reference_genome_index/Aedes-albopictus-Foshan_SCAFFOLDS_AaloF1.fa $file > $ind""_$1""_vs_AaloF1.sai
bwa samse reference_genome_index/Aedes-albopictus-Foshan_SCAFFOLDS_AaloF1.fa $ind""_$1""_vs_AaloF1.sai $file > $ind""_$1""_vs_AaloF1.sam
rm $ind""_$1""_vs_AaloF1.sai
samtools flagstat $ind""_$1""_vs_AaloF1.sam 

done
