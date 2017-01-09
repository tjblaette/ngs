#!/bin/bash

# Script takes as input a "*ReadsPerGene*" counts file from STAR and converts to the exact format of HTSeq output
# Default column to extract from STAR is 2 which corresponds to read counts for an unstranded protocol

 INPUT=$1
 COLUMN=${2:-2} #for unstranded

 OUTPUT=$(echo ${INPUT%_R1ReadsPerGene.out.tab}.counts)
 
 NOFEAT=$(grep 'noFeature' $INPUT | cut -f2) 
 AMBIG=$(grep 'ambiguous' $INPUT | cut -f2) 
 MULTIMAPPED=$(grep 'multimap' $INPUT | cut -f2) 
 UNMAPPED=$(grep 'unmapped' $INPUT | cut -f2)

cat <(grep -v '^N_' $INPUT | sort | cut -f1-$COLUMN) <(echo "__no_feature	$NOFEAT") <(echo "__ambiguous	$AMBIG") <(echo "__too_low_aQual	0") <(echo "__not_aligned	$UNMAPPED") <(echo "__alignment_not_unique	$MULTIMAPPED") > $OUTPUT
