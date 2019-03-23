#!/bin/bash

####
# T.J.Blätte
# 2015
####
#
# Collects all of the filtering stat files
#       created by bpipe filter scripts
#       from a specified folder, containing
#       alignment and duplicate statistics,
#       and saves them to a single TSV file.
#
# Args:
#   DIR: Name of / Path to the folder to
#       collect filtering stat files from.
#
# Output:
#   $DIR/filter_summary.tsv: containing all
#       of the collected and merged statistics.
#
####


DIR=${1:-'.'}
OUT="${DIR}/filter_summary.tsv"


# insert header
echo -e "File name\tTotal calls\tExonic/Splicing\tNonsynonymous\tPassed dbSNP/COSMIC-filter\tAML candidate genes" > "$OUT"

# extract variant filtering statistics from all *stat*-files in $DIR
for file in ${DIR}/*_filter_statistic.txt
do
        echo -ne "$(basename ${file%_merged_filter_statistic.txt})\t" >> "$OUT" #file name
        grep 'total' "$file" | sed 's/.*total of \(.*\) calls.*/\1/g' | tr '\n' '\t' >> "$OUT" #total calls
        grep 'total' "$file" | sed 's/\(.*\) out of.*/\1/g' | tr '\n' '\t' >> "$OUT" #exonic/splicing variants in total
        grep 'remain' "$file" | sed 's/\(.*\) calls.*/\1/g' | tr '\n' '\t' >> "$OUT" #calls remaining after removal of synonymous variants and those failing dbSNP/COSMIC-filter
        grep 'candidate' "$file" | sed 's/\(.*\) calls.*/\1/g' >> "$OUT" #calls affecting AML candidate genes
done
