### locus_factory.sh by Clément Goubert ###
### Part of the High-Troughput Sequencing Transposon Display Analysis Suite ###

### Dependency: lf_partner.R

### Description: This scripts uses the filtered R2 reads to assemble the insertion loci. 
###				 It realizes successively two clustering steps using the program CD-HIT-EST. 
###				 The first makes a intra-individual clustering of the R2 reads to retreive individual loci. 
###				 The second step clusterize individual loci to estimate the insertion polymorphism between individual
###				 This scrit is to run independantly for each TE famuly
###				 The last step filter out the loci which overlapp transposable elements found in Goubert et al. 2015.

### Directions:  Execute this script in each subfolder containing the individual R2 reads of each TE family

### Input:		 The cleaned and trimmed R2 reads (insertion loci) of a TE family, one file per individual, in fasta format
###	Outpout:	 It outputs a matrix of reads per loci per individual, with loci in lines and individuals in columns

### PARAMETERS (see CD-HIT-EST manual for details)
# -c = $1 (intra-individual indentity)
# -aS = $2 (intra-individual shorter sequence size)
# -l =  $3 (intra-individual read size cutoff in bp)
# ET name = $4 (e.g. IL1, RTE4,...)
# -c all = $5
# -aS all = $6 


echo ""
echo "##########################"
echo "#        WELCOME         #"
echo "#        to the          #"
echo "#     LOCUS FACTORY      #"
echo "##########################"

IND_LIST="*.fasta"

echo ""
echo "##########################"
echo "# INDIVIDUAL CLUSTERING  #"
echo "##########################"

for ind in ${IND_LIST}
do

file_name="$( ls $ind | sed 's/_/ /g' | awk '{print $1}')"
# run cd-hit individually

echo ""
echo "processing $file_name..."

cd-hit-est -i $ind -o $file_name""_clustering -G 1 -g 1 -c $1 -aS $2 -d 0 -l $3 -T 4

grep '*' $file_name""_clustering.bak.clstr | sed 's/>//g' | sed 's/\.\.\.//g' | awk '{print "CL-"$1 "\t" $3}' | sort -k2,2 > tmp

awk '{print $2}' tmp > tmp2

#part of the code sorting the fasta per header name is from Devon Ryan on SeqAnswer forum
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' tmp2 $ind | awk '{idx=index($0, " "); ID=$1; NAME=substr($0, idx+1); getline SEQ; printf("%s\t%s\t%s\n", NAME, ID, SEQ);}' | sort -k1d,1 | awk -F"\t" '{printf("%s %s\n%s\n", $2,$1,$3)}' | awk '{print $1}' | sed 's/>HWI/HWI/g' | awk '{getline second; print $0 "\t" second}' >  $ind""_refseq_tmp.fasta


# for line in $(grep -n '^\[.*\]$' $ind""_refseq_tmp.fasta | sort -k2 -t: |cut -f1 -d:);
# do
# tail -n +$line sections.txt | head -n 2 >> 
# done

join -a2 -12 -21 tmp $ind""_refseq_tmp.fasta  -o 2.1,1.1,2.2 | awk '{print $1"_"$2 ; print $3}' | sed 's/HWI/>HWI_'$file_name'_/g' > $ind""_refseq.fasta

rm *tmp*

echo ""
echo "$file_name done !!!"

done

echo ""
echo "##################"
echo "# CLUSTERING ALL #"
echo "##################"


cat *_refseq.fasta > all_ind_""$4"".fasta

echo "processing all individual clustering..."

cd-hit-est -i all_ind_""$4"".fasta -o all_""$4""_clustering -G 1 -g 1 -c $5 -aS $6 -d 0 -T 8

echo "done !!!"


echo ""
echo "###################"
echo "# COVERAGE REPORT #"
echo "###################"


FILE_LIST=*.bak.clstr

for file in ${FILE_LIST}
do

file_name="$( ls $file | sed 's/_/ /g' | awk '{print $1}')"


echo "$file_name..."

# sort -k1,1 $file | awk 'BEGIN {compt=0; cluster=$1} {if ($1==cluster) {cluster=$1; compt=compt+1; print $0 "\t" compt "\tREAD"} else {cluster=$1; compt=1; print $0 "\t" compt "\tNEW"}}' | grep -B 1 'NEW' | grep -v '^-' | sort -k1,1 -k5,5nr | sort -u -k1,1 | awk -v x=$file_name '/at/ {print x"_"$1 "\t" $6; next} {print x"_"$1 "\t" $5}' > $file_name""_locus_and_coverage.lf
awk '{print $1}' $file > R.tmp
Rscript ~/pbil-pandata/TD_albo/scripts/lf_partner.R R.tmp
awk -v x=$file_name '{print x"_"$0} ' counts.tmp > $file_name""_locus_and_coverage.lf

done

# add all individual loci and coverage and sort by name
cat *.lf | sort -k1,1 > tmp_lf

# display all individual loci name with individual loci name
awk '{print "All-"$1 "\t" $3}' all_""$4""_clustering.bak.clstr | sed 's/_/\t/g' | sed 's/-/_/g' | sed 's/\.\.\.//g' | awk '{print $1 "\t" $3 "_" $5}' | sed 's/CL_//g' | sort -k2,2 > tmp_all_locus

join -11 -22 tmp_lf tmp_all_locus -o 2.1,1.1,1.2 | sort -k1,1 -k2,2 > All_""$4""loci_and_coverage.lf

rm *tmp*


echo ""
echo "##################"
echo "#  MATRIX MAKER  #"
echo "##################"


awk '{print $1}' All*.lf | sort | uniq > locus_list
# tour=$((0))
indiv=$(sed 's/_/ /g' All*.lf | awk '{print $3}'| sort | uniq)


for ind in $indiv

do

echo "c'est le tour de $ind"
# tour=(($tour+1))



grep ' '$ind'_' All*.lf | sort -k1,1 > tmpR2
Rscript ~/pbil-pandata/TD_albo/scripts/lf_partner2.R tmpR2
awk 'NR > 1 { print }' Loc.tmp > Loc_$ind

join -a1 -11 -21 locus_list Loc_$ind -o 1.1,2.2 | awk '{if ($2=="") {print $1"\t 0"} else {print $O}}' | awk '{print $2}' > tmp


echo $ind > ind
cat ind tmp > Locus_$ind

rm *tmp* #IND_$locus

done


echo "locus" > loctext
cat loctext locus_list > loc2
paste loc2 Locus_* > Matrice_finale


# rm Loc_* loctext loc2

echo ""
echo "##################"
echo "# REMOVE TEs !!! #"
echo "##################"
echo ""

### This step filter out the loci that overlapp transposable element sequence, using here the file from the dnaPipeTE analysis (Goubert et al 2015, GBE) file and BLAT.

mkdir -p blat
#blating the locus on albo TEs froms dnaPipeTE
echo "blat on TE database..."
~/Downloads/blat ~/pbil-pandata/Papier_ET_albo/albo1.19_2x10X/Trinity.fasta all_""$4""_clustering -minIdentity=0.8 -out=blast8 -noHead -q=dna -t=dna ./blat/locus_vs_dnaPipeTE.out 

echo "sorting..."
#sort blatout file
sort -u -k1,1 ./blat/locus_vs_dnaPipeTE.out | awk '{print $1}' > ./blat/sorted_locus_vs_dnaPipeTE.out

echo "print new refs without TEs..."
#making reflocus list
perl -ne 'if(/^>(\S+)/){$c=!$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' ./blat/sorted_locus_vs_dnaPipeTE.out all_""$4""_clustering > noET_all""$4""_clustering.fasta


echo ""
echo "##############################"
echo "# LOCUS FACTORY FINISHED !!! #"
echo "##############################"