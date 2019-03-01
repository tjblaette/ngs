#!/bin/bash

SAMPLE="$1"

# add columns to variant CSVs to collect information of manual inspection
if [ -e results_csv/ ]
then 
  for file in results_csv/${SAMPLE}*csv
  do
    sed -i -e 's/^"Chr"/"Author","Project","BioID","Study","StudyID","PB-KM_ID","Time","TypeOfSequencing","Decision","Comment","Chr"/' -e 's/^"chr/"","","","","","","","","","","chr/' "$file"
  done
fi


#delete temporary files
rm -f ${SAMPLE}*.subsample
rm -f ${SAMPLE}*.success
rm -f ${SAMPLE}*.pileup
rm -f intermediate_files/${SAMPLE}*.sam
rm -f intermediate_files/${SAMPLE}*somVARSC*
rm -f intermediate_files/${SAMPLE}*.alignMEM.sortPIC.ba*
rm -f intermediate_files/${SAMPLE}*dedupPIC.ba*
rm -f intermediate_files/${SAMPLE}*.coverBED_exon.txt
rm -f intermediate_files/${SAMPLE}*otherinfo*
rm -f results_varscan/${SAMPLE}*dummy
rm -f results_csv/${SAMPLE}*_PREcandidates.csv
rm -rf tmp/
rm -rf JAVA_TMP/
