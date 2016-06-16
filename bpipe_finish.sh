
#delete temporary files
rm -f *.subsample
rm -f *.success
rm -f *.pileup
rm -f intermediate_files/*.sam
rm -f intermediate_files/*somVARSC.*
rm -f intermediate_files/*.alignMEM.sortPIC.ba*
rm -f intermediate_files/*dedupPIC.ba*
rm -f intermediate_files/*.coverBED_exon.txt
rm -f intermediate_files/*otherinfo*
rm -f results_varscan/*dummy
rm -f results_csv/*_PREcandidates.csv
rm -rf tmp/
rm -rf JAVA_TMP/
