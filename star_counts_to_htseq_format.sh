#!/bin/bash

####
# T.J.Blätte
# 2016
####
#
# Converts read counts output by STAR into HTSeq-based
#       file format, to prepare them as input to DESeq2.
#       The main difference is that STAR's output files
#       contain four columns: The gene / feature ID that
#       was counted, the counts for unstranded libraries
#       and the complimentary counts for stranded protocols.
#       HTSeq outputs two columns, where the second column
#       depends on the library preparation protocol pursued.
#       (STAR outputs all three variants, HTSeq requires the
#       correct protocol be provided prior to counting.
#
# Args:
#   INPUT: Read counts output by STAR. Should be
#       "*ReadsPerGene.out.tab" when generated by our pipeline.
#   COLUMN: Column number to extract, defaults to 2 when none is given.
#       -> 2 for unstranded protocols
#       -> 3/4 for respective stranded protocl
#
# Output:
#   $(echo ${INPUT%_R1ReadsPerGene.out.tab}.counts): HTSeq-formatted
#       read counts.
#
####


 INPUT=$1
 COLUMN=${2:-2} #for unstranded


 OUTPUT=$(echo ${INPUT%_R1ReadsPerGene.out.tab}.counts)

 NOFEAT=$(grep 'noFeature' $INPUT | cut -f2)
 AMBIG=$(grep 'ambiguous' $INPUT | cut -f2)
 MULTIMAPPED=$(grep 'multimap' $INPUT | cut -f2)
 UNMAPPED=$(grep 'unmapped' $INPUT | cut -f2)

cat <(grep -v '^N_' $INPUT | sort | cut -f1,$COLUMN) <(echo "__no_feature	$NOFEAT") <(echo "__ambiguous	$AMBIG") <(echo "__too_low_aQual	0") <(echo "__not_aligned	$UNMAPPED") <(echo "__alignment_not_unique	$MULTIMAPPED") > $OUTPUT
