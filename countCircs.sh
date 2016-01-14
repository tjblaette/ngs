#!/bin/bash

PREFIX=$1
#SAMPLES= ALL OTHER INPUTS

KEY="${PREFIX}_key.tsv"
LIST="${PREFIX}_list.tsv"

mkdir -p intermediate_files

tail -q -n +2 "${@:2}" | cut -f1-3,8 | sed -e 's/\t/_/' -e 's/\t/_/' | sort | uniq > $KEY
awk -v OFS='\t' '{print $1,0}' $KEY > $LIST

for FILE in "${@:2}"
do
   echo "$FILE"
   # prepare sample
   SAMPLE="$(tail -n +2 "$FILE" | cut -f1-4 | sed -e 's/\t/_/' -e 's/\t/_/' | sort | uniq)"
   SAMPLE_LONG="$(cat <(echo "$SAMPLE") "$LIST" | sort)"

   UNIFORM="intermediate_files/$(basename ${FILE%.txt}_uniformIDs.txt)"
   ENSEMBL="intermediate_files/$(basename ${FILE%.txt}_uniformEnsemblIDs.txt)"

   echo "$SAMPLE_LONG" | awk -v OFS='\t' '{if(PREV==0){PREV=$1;} ID=$1; if(ID==PREV){SUM=SUM+$2}else{print(PREV,SUM); SUM=$2;}; PREV=ID;} END {print(ID,SUM);}' > $UNIFORM

   # create the same output but with ENSEMBL gene IDs as identifiers instead of custom junction coords
   join -t '	' $KEY $UNIFORM  | cut -f2,3 | sed 's/.*"\([^"]\+\)"/\1/' | sort | awk -v OFS='\t' 'NR==1 {PREV=$1; SUM=0} {ID=$1; if(ID==PREV){SUM=SUM+$2}else{print(PREV,SUM); SUM=$2;}; PREV=ID;} END {print(ID,SUM);}' > $ENSEMBL
   
   echo '__no_feature	0
__ambiguous	0
__too_low_aQual	0
__not_aligned	0
__alignment_not_unique	0' >> $UNIFORM

   echo '__no_feature	0
__ambiguous	0
__too_low_aQual	0
__not_aligned	0
__alignment_not_unique	0' >> $ENSEMBL
done

# create DESeq2 input files -> one with counts per circRNA junction, one with counts per ENSEMBL ID
echo 'sampleName      fileName        condition' > ${PREFIX}_circs.tsv
ls -C intermediate_files/*uniformIDs.txt | sed -e 's/.*\///' -e 's/\(.*\)/\1\t\1/' >> ${PREFIX}_circs.tsv
# add to last column: Group information of each sample


echo 'sampleName      fileName        condition' > ${PREFIX}_genes.tsv
ls -C intermediate_files/*uniformEnsemblIDs.txt | sed -e 's/.*\///' -e 's/\(.*\)/\1\t\1/' >> ${PREFIX}_genes.tsv
# add to last column: Group information of each sample
