#!/bin/bash


# $1 = input file for DESeq2
# $2 = biomart DB for gene id conversion
# $3 = alpha for DEG statistic (default:0.1)
# $4 = min absolute LFC tested for (default:0)

Rscript deseq2.R $1 $3 $4

head -1 ${1}_DESeq2results.txt > ${1}_DESeq2results_sorted.txt
tail -n +2 ${1}_DESeq2results.txt | sort >> ${1}_DESeq2results_sorted.txt
sed -i 's/"//g' ${1}_DESeq2results_sorted.txt

join --header ${1}_DESeq2results_sorted.txt $2 > ${1}_DESeq2results_annotated.txt

sort -g -k 7,7 ${1}_DESeq2results_annotated.txt > ${1}_DESeq2results_sorted.txt
grep -wv 'NA' ${1}_DESeq2results_sorted.txt >  ${1}_DESeq2results_sorted_woutNA.txt

rm ${1}_DESeq2results_annotated.txt
