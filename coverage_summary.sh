#!/bin/bash

#$1 = directory containing statistic files

#when no input directory was supplied, assume "."
if [[ $# -eq 0 ]] ;
then
        dir='.'
else
        dir=$1
fi


#insert header
echo 'File name,# reads on target,Average coverage,No coverage,1x Coverage,10x Coverage,15x Coverage,50x Coverage,120x Coverage,200x Coverage,500x Coverage,1000x Coverage,1500x Coverage,2000x Coverage,2500x Coverage' > ${dir}/coverage_summary.csv

#extract coverage information from all BED.txt-files in $dir 
for file in ${dir}/*BED.txt
do
	echo -n "$(basename $file)," >> ${dir}/coverage_summary.csv
	sed 's/.*(\(.*\))/\1/g' $file | sed 's/on average \(.*\) reads.*/\1/g' | tr '\n' ',' | sed 's/,$/\n/g' >> ${dir}/coverage_summary.csv 
done
