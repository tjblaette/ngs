#!/bin/bash

# for -sorted option of BEDTools, input BAM and BED have to be sorted the same way
# BAM position sorting is required by GATK for example
# BED sorting has to be adjusted accordingly!
# -> -sorted reduces RAM consumption of coverageBed command
# contig order of BAM is defined by REF.fai = REFFAI in bpipe config files

# in a second step, 5bp padding is added to intervals to retain splicing variants
# when restricting pileup file to target BED

# $1 = BED file to be sorted
# $2 = REFFAI that defines sorting order

BED=$1
REFFAI=$2

## sort BED file
# create genome file from REFFAI
cut -f1-2 /NGS/refgenome/GATK/ucsc.hg19.fasta.fai > ${REFFAI%.fai}.genomeFile

rm -f ${BED%.bed}_sorted.bed

cut -f1 $REFFAI |
while read CONTIG
do
	echo "$CONTIG"
	grep -w "$CONTIG" $BED | sort -k1,1 -k2,2n | uniq >> ${BED%.bed}_sorted.bed
done 

echo "Done sorting!"

## add 5bp padding to each interval in the BED file so that we can safely restrict variant calling to that area
awk -v OFS='\t' '!/browser|track/ {$2=$2-5; $3=$3+5; print $0} /browser|track/' ${BED%.bed}_sorted.bed > ${BED%.bed}_sorted.bed_padded

echo "Done padding!"
