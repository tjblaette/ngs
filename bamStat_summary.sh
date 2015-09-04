#!/bin/bash

#$1 = directory containing Picards "*_metrics.txt" statistic files

#when no input directory was supplied, assume "."
if [[ $# -eq 0 ]] ;
then
        dir='.'
else
        dir=$1
fi


#insert header
echo 'File name,# reads unmapped total,# reads mapped total,# paired reads mapped,# duplicate reads total,# paired duplicate reads,# paired optical duplicate reads,percent duplicate reads total' > ${dir}/bamStat_summary.csv

#extract mapping and duplicate statistics from all metrics.txt-files in $dir 
for file in ${dir}/*metrics.txt
do
        echo -n "$(basename $file)," >> ${dir}/bamStat_summary.csv
	head  -9 $file | tail -1 | sed 's/[a-zA-Z]//g' | awk '{print $3,$1+2*$2,$2*2,$4+$5*2,$5*2,$6*2,$7}' >> ${dir}/bamStat_summary.csv
done

