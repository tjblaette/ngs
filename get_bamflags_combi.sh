#!/bin/sh

# Description: This script extracts all bamflags from all files given as arguments (so more than one can be supplied). Flags of each read are saved in the file ending with "_flags_final_combi.txt". The combination of flags of each read and their order is preserved
# Duration: ~1h

for var in "$@"
do

cut -f2 "$var" | uniq -c | sed 's/^ *//' > ${var%.sam}_flagsCombi.txt #sed to remove leading whitespace
#samtools view $var | cut -f2 | uniq -c | sed 's/^ *//' > ${var}_flags_combi.txt #sed to remove leading whitespace
Rscript /home/tabl/scripts/get_flags_combi.R  ${var%.sam}_flagsCombi.txt ${var%.sam}_flags.txt  

done

rm -f ${var%.sam}_flagsCombi.txt
