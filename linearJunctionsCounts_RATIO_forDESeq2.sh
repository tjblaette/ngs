# 12.12.2018
# Susanne Lux proudly presents
# commandline  ./linearJunctionsCounts_RATIO_forDEseq.sh , works for *linearJunctionsCounts.tsv files in same folder
# script for usage with output from linear junction read script
# !/bin/bash
# prepare as input for deseq: uniformIDs for uniform circRNA followed by ratio mean linear junction cout/circJunction count
#extract circRNA coordinates and linear junction counts from all linearJunctionsCounts.tsv-files in $dir
for file in ./*linearJunctionsCounts.tsv
do
  echo "$file"
  UNIFORM="$(basename ${file%.tsv}_RATIOuniformIDs.tsv)"
  echo "$file"
  paste <(cut -f 1-3 $file | sed 's/\t/_/g' ) <(cut -f12 $file) > $UNIFORM

done
