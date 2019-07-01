#!/bin/bash

SAMPLE="$1"

# make sure the most 5' filename prefix is used:
SAMPLE=$(echo "$SAMPLE" | cut -f1 -d'R')

# add columns to variant CSVs to collect information of manual inspection
if [ -e results_csv/ ]
then
  for file in results_csv/${SAMPLE}*csv
  do
    if [ -e "$file" ]
    then
        sed -e 's/^"Chr"/"Author","Project","BioID","Study","StudyID","PB-KM_ID","Time","TypeOfSequencing","Decision","Comment","Chr"/' -e 's/^"chr/"","","","","","","","","","","chr/' "$file" > "${file}_tmp" && mv "${file}_tmp" "$file"
    fi
  done
fi



#delete temporary files
rm -f ${SAMPLE}*.subsample
rm -f ${SAMPLE}*.success
#rm -f ${SAMPLE}*.pileup
rm -f intermediate_files/${SAMPLE}*.sam
rm -f intermediate_files/${SAMPLE}*.alignMEM.sortPIC.ba*
rm -f intermediate_files/${SAMPLE}*dedupPIC.ba*
rm -f intermediate_files/${SAMPLE}*alignMEM*dedupOptPIC.ba* # final BAM for RNA pipeline, delete for DNA only (which ends on GATK realignment)
rm -f intermediate_files/${SAMPLE}*.coverBED_exon.txt
rm -f intermediate_files/${SAMPLE}*otherinfo*
rm -f results_varscan/${SAMPLE}*fixFormat*
rm -f results_csv/${SAMPLE}*_PREcandidates.csv
#rm -rf tmp/
#rm -rf JAVA_TMP/
