#!/bin/bash

####
# T.J.Blätte
# 2017
####
#
# Script requested by Montse to filter a specified list of variants
#       from our mutation calling pipelines' outputs.
#
# Args:
#   RECURRENCES: File defining the specific mutations to filter.
#       These must be in the same format as they would be in the
#       files to be filtered, for example:
#
#       "chr17","7579472","7579472","G","C"
#       "chr1","120611960","120611960","C","T"
#
#   [...]: One or more CSV files to be filtered, as created
#       by our bpipe mutation calling pipelines.
#
# Output:
#   ${CSV%.csv}_recurrent.csv: File containing those mutations
#       matching the filter. One of these is written for each input CSV.
#   ${CSV%.csv}_non-recurrent.csv: File containing all mutations
#       *not* matching the filter. Again, one of these is written
#       for each input CSV.
#
####


RECURRENCES=$1

 cut -f1-5 "$RECURRENCES" | sed -e 's/^\([^"]\)/"\1/' -e 's/\(^"]\)$/\1"/' -e 's/\t/","/g' > tmp && mv tmp $1;


 for CSV in "${@:2}"
 do
   ls $CSV
   grep -f  "$RECURRENCES"  "$CSV" > ${CSV%.csv}_recurrent.csv
   grep -vf  "$RECURRENCES"  "$CSV" > ${CSV%.csv}_non-recurrent.csv
 done
