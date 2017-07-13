
# add columns to variant CSVs to collect information of manual inspection
if [ -e results_csv/ ]
then 
  for file in results_csv/*csv
  do
    sed -i -e 's/^"Chr"/"Author","Project","BioID","Study","StudyID","PB-KM_ID","Time","TypeOfSequencing","Decision","Comment","Chr"/' -e 's/^"chr/"","","","","","","","","","","chr/' "$file"
  done
fi


#delete temporary files
rm -f *.subsample
rm -f *.success
rm -f *.pileup
rm -f intermediate_files/*.sam
rm -f intermediate_files/*somVARSC*
rm -f intermediate_files/*.alignMEM.sortPIC.ba*
rm -f intermediate_files/*dedupPIC.ba*
rm -f intermediate_files/*.coverBED_exon.txt
rm -f intermediate_files/*otherinfo*
rm -f results_varscan/*dummy
rm -f results_csv/*_PREcandidates.csv
rm -rf tmp/
rm -rf JAVA_TMP/
