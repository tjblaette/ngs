#!/bin/bash

# for -sorted option of BEDTools, input BAM and BED have to be sorted the same way
# BAM position sorting is required by GATK for example
# BED sorting has to be adjusted accordingly!
# -> -sorted reduces RAM consumption of coverageBed command


# contig order of BAM is defined by REF.fai = REFFAI in bpipe config files

# $1 = BED file to be sorted
# $2 = REFFAI that defines sorting order

BED=$1
REFFAI=$2

# GENOME_FILE FROM REFFAI
cut -f1-2 /NGS/refgenome/GATK/ucsc.hg19.fasta.fai > ${REFFAI%.fai}.genomeFile

rm -f ${BED%.bed}_sorted.bed

cut -f1 $REFFAI |
while read CONTIG
do
	#echo "$CONTIG #############################################" >> ${BED%.bed}_sorted.bed
	grep -w "$CONTIG" $BED | sort -k1,1 -k2,2n >> ${BED%.bed}_sorted.bed
done 
