#! /bin/bash

### files_sampler.sh by Clément Goubert ###
### Part of the High-Troughput Sequencing Transposon Display Analysis Suite ###

### Description: # Sample read files individual per individual to a fixed value

###	Input:		 read files (fasta)
###	Output:		 sampled read files (fasta)

###PARAMETERS:
#$1 = n to sample

FILE_LIST=$"*.fasta"
mkdir -p samples_1_""$1
rm samples_1_$1/discarted

for file in ${FILE_LIST}

do

reads="$( grep -c '>' $file)"
file_name="$( ls $file | sed 's/_/ /g' | awk '{print $1}')"



	if (( $reads >= $1 ))

		then echo "sampling individual $file_name $reads to $1" ; perl /pandata/goubert/TD_albo/scripts/random_sequence_sample.pl -i $file -o samples_1_$1/$file_name"_"$1""_sample.fasta -n $1

		else echo "individual $file_name discarted" ; echo "$file_name $reads" >> samples_1_$1/discarted ### individuals with not enough reads

	fi

done




FILE_LIST=$"*.fasta"
mkdir -p samples_2_""$1
rm samples_2_$1/discarted

for file in ${FILE_LIST}

do

reads="$( grep -c '>' $file)"
file_name="$( ls $file | sed 's/_/ /g' | awk '{print $1}')"



	if (( $reads >= $1 ))

		then echo "sampling individual $file_name $reads to $1" ; perl /pandata/goubert/TD_albo/scripts/random_sequence_sample.pl -i $file -o samples_2_$1/$file_name"_"$1""_sample.fasta -n $1

		else echo "individual $file_name discarted" ; echo "$file_name $reads" >> samples_2_$1/discarted

	fi

done

FILE_LIST=$"*.fasta"
mkdir -p samples_3_""$1
rm samples_3_$1/discarted

for file in ${FILE_LIST}

do

reads="$( grep -c '>' $file)"
file_name="$( ls $file | sed 's/_/ /g' | awk '{print $1}')"



	if (( $reads >= $1 ))

		then echo "sampling individual $file_name $reads to $1" ; perl /pandata/goubert/TD_albo/scripts/random_sequence_sample.pl -i $file -o samples_3_$1/$file_name"_"$1""_sample.fasta -n $1

		else echo "individual $file_name discarted" ; echo "$file_name $reads" >> samples_3_$1/discarted

	fi

done

