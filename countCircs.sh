#!/bin/bash

LIST=$1
#SAMPLES= ALL OTHER INPUTS

mkdir -p intermediate_files

for FILE in "${@:2}"
do
   echo "$FILE"
   # prepare sample
   SAMPLE="$(cat "$FILE" | cut -f1-4 | sed -e 's/\t/_/' -e 's/\t/_/' | sort | uniq)"
   SAMPLE_LONG="$(cat <(echo "$SAMPLE") "$LIST" | sort)"

   echo "$SAMPLE_LONG" | awk -v OFS='\t' '{if(PREV==0){PREV=$1;} ID=$1; if(ID==PREV){SUM=SUM+$2}else{print(PREV,SUM); SUM=0;}; PREV=ID;} END {print(ID,SUM);}' > intermediate_files/$(basename ${FILE%.txt}_uniformIDs.txt)


   echo '__no_feature	0
__ambiguous	0
__too_low_aQual	0
__not_aligned	0
__alignment_not_unique	0' >> intermediate_files/$(basename ${FILE%.txt}_uniformIDs.txt)
done

# create DESeq2 input file
echo 'sampleName      fileName        condition' > forDESeq2.tsv
ls -C intermediate_files/*unif* | sed -e 's/.*\///' -e 's/\(.*\)/\1\t\1/' >> forDESeq2.tsv
# add to last column: Group information of each sample
