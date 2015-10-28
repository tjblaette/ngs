#!/bin/bash


# $1 = input file for DESeq2
# $2 = alpha for DEG statistic (default:0.1)
# $3 = biomart DB for gene id conversion
# $4 = min absolute LFC tested for (default:0) -> does not take any effect because R function in deseq does not accept variable for lfc parameter

IN=$1
ALPHA=${2:-0.1}
BIOMART=${3:-'/NGS/known_sites/human_ensembl_biomart_gene_ID_to_symbol/mart_export_sorted_woutLRG.txt'}
LFC=${4:-0} # does not take any effect!

Rscript /NGS/myscripts/deseq2.R $IN $ALPHA $LFC

head -1 ${IN}_DESeq2results.txt > ${IN}_DESeq2results_sorted.txt
tail -n +2 ${IN}_DESeq2results.txt | sort >> ${IN}_DESeq2results_sorted.txt
sed -i 's/"//g' ${IN}_DESeq2results_sorted.txt

join --header ${IN}_DESeq2results_sorted.txt $BIOMART > ${IN}_DESeq2results_annotated.txt

sort -g -k 7,7 ${IN}_DESeq2results_annotated.txt > ${IN}_DESeq2results_sorted.txt
grep -wv 'NA' ${IN}_DESeq2results_sorted.txt >  ${IN}_DESeq2results_sorted_woutNA.txt

rm ${IN}_DESeq2results_annotated.txt
