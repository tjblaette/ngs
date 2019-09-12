#!/bin/bash

####
# T.J.Blätte
# 2015
####
#
# Takes a BED file as input and creates a sorted
#       as well as a sorted and padded version.
#       The padded version's intervals are extended
#       by 5 bp in both directions compared to the original.
#       The BED file passed to the mutation calling
#       pipelines for coverage calculation (via BEDTools'
#       coverageBedd command) must be sorted to reduce RAM usage.
#       Since variant calling is restricted to this BED file as well,
#       the padded version is used for targeted enrichment protocols
#       to profit from enriched sequences flanking the actual targets,
#       which frequently support additional splice site mutations.
#
# Args:
#   BED: BED file to be sorted / padded
#   REFFAI: Reference genome index (.fai) with the order of its
#       contigs, as defined in the respective bpipe config file.
#
# Output:
#   ${BED%.bed}_sorted.bed: Sorted BED file.
#   ${BED%.bed}_sorted.bed_padded: Sorted & padded BED file.
#
####


BED=$1
REFFAI=$2


## sort BED file
rm -f ${BED%.bed}_sorted.bed

cut -f1 $REFFAI |
while read CONTIG
do
	echo "$CONTIG"
	grep -w "$CONTIG" $BED | sort -k1,1 -k2,2n | uniq >> ${BED%.bed}_sorted.bed
done
echo "Done sorting!"


## add 5bp padding to each interval in the BED file
#  so that we can safely restrict variant calling to that area
awk -v OFS='\t' '!/browser|track/ {$2=$2-5; $3=$3+5; print $0} /browser|track/' ${BED%.bed}_sorted.bed > ${BED%.bed}_sorted.bed_padded
echo "Done padding!"
