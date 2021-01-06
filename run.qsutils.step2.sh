#!/bin/bash

#source activate qap
#export PATH=/gpfs1/scratch/30days/uqkchew2/03.KAT_ASTHMA/qap/:$PATH
threads=24

ulimit -v 240000000

#array1=$(cut -f1 DNA_diversity.csv| sed 's/"//g')
#array2=$(cut -f2 DNA_diversity.csv| sed 's/"//g')

cut -f1 DNA_diversity.csv | sed 's/"//g' > file1
cut -f2 DNA_diversity.csv | sed 's/"//g' > file2
while read X; do C[++i]=$X;done < file1
while read Y; do D[++i]=$Y;done < file2

array1=$(echo ${C[@]})
array2=$(echo ${D[@]}) 


for i in `cat text`; do base=`basename $i "_R1_001.fastq.gz"` 

for gene in A_HA_H1 A_NA_N1 A_PB1 A_PB2 A_NS A_MP A_NP A_PA 
do 

grep $gene DNA_diversity.csv | cut -f3 > ${gene}_size.txt

mybase=${base}
mygene=${gene}
mybam=${mybase}-V2/${gene}.bam
myfasta=${mybase}-V2/${gene}.fasta
mysize=${gene}_size.txt
value1=${base}
value2=${gene}

echo $base

if [[ " ${array1[@]} "  =~ " ${value1} " && " ${array2[@]} "  =~ " ${value2} " ]]
then
	echo "${value2} for ${value1} is already calculated"
	#echo "The ${value2} fro ${value1} is not present, calculating"
        #echo "Haplotypes from bam files : analysing ${mybam} and writing results to a file"
        Rscript --vanilla test2.R ${mybase} ${mygene} ${mybam} ${myfasta} ${mysize} > /dev/null	
else
	echo "The ${value2} in ${value1} is not present, calculating..."
	echo "Haplotypes from bam files : analysing ${mybam} and writing results to a file"
	Rscript --vanilla test1.R ${mybase} ${mygene} ${mybam} ${myfasta} > /dev/null
fi	

done 
done
