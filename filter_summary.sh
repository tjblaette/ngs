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
echo 'File name,Total calls,Exonic/Splicing,Nonsynonymous,Passed dbSNPi/COSMIC-filter,AML candidate genes' > ${dir}/filter_summary.txt

#extract variant filtering statistics from all *stat*-files in $dir
for file in ${dir}/*_filter_statistic.txt
do
        echo -n "$(basename ${file%_merged_filter_statistic.txt})," >> ${dir}/filter_summary.txt #file name
	grep 'total' $file | sed 's/.*total of \(.*\) calls.*/\1/g' | tr '\n' ',' >> ${dir}/filter_summary.txt #total calls
	grep 'total' $file | sed 's/\(.*\) out of.*/\1/g' | tr '\n' ',' >> ${dir}/filter_summary.txt #exonic/splicing variants in total
	grep 'remain' $file | sed 's/\(.*\) calls.*/\1/g' | tr '\n' ',' >> ${dir}/filter_summary.txt #calls remaining after removal of synonymous variants and those failing dbSNP/COSMIC-filter
	grep 'candidate' $file | sed 's/\(.*\) calls.*/\1/g' >> ${dir}/filter_summary.txt #calls affecting AML candidate genes
done
