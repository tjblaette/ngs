#!/bin/bash

####
# T.J.Blätte
# 2015
####
#
# Collects all of the coverage stat summary
#       files created using BEDTools
#       from a specified folder and saves
#       them to a single TSV file.
#
# Args:
#   DIR: Name of / Path to the folder to
#       collect coverage stat files from.
#
# Output:
#   $DIR/coverage_summary.tsv: containing all
#       of the collected and merged statistics.
#
####


DIR=${1:-'.'}
OUT="${DIR}/coverage_summary.tsv"


# insert header
echo -e "File name\tBp on target\tAverage coverage\tNo coverage\tMin 1x coverage\tMin 10x coverage\tMin 15x coverage\tMin 50x coverage\tMin 100x coverage\tMin 120x coverage\tMin 200x coverage\tMin 500x coverage\tMin 1000x coverage\tMin 1500x coverage\tMin 2000x coverage\tMin 2500x coverage" > "$OUT"

# extract coverage information from all BEDTools coverage stat files in $DIR
for file in ${DIR}/*BED.summary
do
        echo -ne "$(basename "$file")\t" >> "$OUT"
        sed 's/.*(\(.*\))/\1/g' $file | sed 's/on average \(.*\) reads.*/\1/g' | tr '\n' '\t' | sed 's/\t$/\n/g' >> "$OUT"
done

