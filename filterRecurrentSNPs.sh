#!/bin/bash

# script requested by Montse to filter (recurrent) variant calls (CSVs) 
# for a provided list of specific mutations

# Input is first a file defining certain mutations
# which must be in the following format with one variant per line
# exactly as they are in our variant calling output files (CSVs):

# "chr17","7579472","7579472","G","C"
# "chr1","120611960","120611960","C","T"

# Output are two files, one containing calls of these mutations, the other 
# containing all remaining variants 

# Ex.: filterRecurrentSNPs.sh recurrentSNPs.txt *final.csv 


RECURRENCES=$1

 cut -f1-5 "$RECURRENCES" | sed -e 's/^\([^"]\)/"\1/' -e 's/\(^"]\)$/\1"/' -e 's/\t/","/g' > tmp && mv tmp $1;


 for CSV in "${@:2}"
 do
   ls $CSV
   grep -f  "$RECURRENCES"  "$CSV" > ${CSV%.csv}_recurrent.csv
   grep -vf  "$RECURRENCES"  "$CSV" > ${CSV%.csv}_nonRecurrent.csv
 done


 
