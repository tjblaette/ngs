#!/bin/sh

#columns in the csv that is to be sorted minus 1
somatic_status=25 #26-1
genomicSuperDups=10 #11-1
cosmic=14 #15-1
dbsnp=13 #14-1
tumor_reads2_plus=30 #31-1
tumor_reads2_minus=31 #32-1
normal_var_freq=19 #20-1

#Input
#$1 = csv-file that is to be sorted
#$2 = Prefix for output-files
#$3 = file containing candidate gene list (Format: one gene name in quotation marks per line)
#$4 = Path to reference genome
#$5 = number of flanking bases to print

#Variable to save line-number for statistic calculation
lines=0

echo "Filter statistics for $1" > ${2}_filter_statistic.txt
echo "" >> ${2}_filter_statistic.txt

#filter exonic, splicing and exonic;splicing (but not ncRNA_splicing, _exonic,...)
head -n 1 $1 > ${2}_filtered.csv
grep -e 'exonic' -e 'splicing' $1 | grep -v '_splicing' | grep -v '_exonic' >> ${2}_filtered.csv
echo "$(( $(wc -l ${2}_filtered.csv | cut -f1 -d' ') -1)) out of a total of $(wc -l $1 | cut -f1 -d' ') calls were classified as exonic, splicing or exonic;splicing and saved to ${2}_filtered.csv" >> ${2}_filter_statistic.txt
echo "" >> ${2}_filter_statistic.txt

#remove synonymous variants
head -n 1 $1 > ${2}_filtered_final.csv
grep 'chr' ${2}_filtered.csv | grep -v '"synonymous' >> ${2}_filtered_final.csv
echo "$(grep -c '"synonymous' ${2}_filtered.csv) of these were synonymous and dropped" >> ${2}_filter_statistic.txt
lines="$(( $(wc -l ${2}_filtered_final.csv | cut -f1 -d' ') -1))"
echo "$lines calls remain after filtering" >> ${2}_filter_statistic.txt

#take out calls that have an entry for dbsnp and not at least two for cosmic
head -n 1 $1 > ${2}_dbsnp_and_one_cosmic.csv
sed -i -e '/^\([^,]*,\)\{13\}"[^\(\.\)"][^,]*,"\."/d' ${2}_filtered_final.csv #delete variants with entry for dbsnp but not cosmic
sed -n -e '/^\([^,]*,\)\{13\}"[^\(\.\)"][^,]*,[^,]*OCCURENCE=1([^)]*)",/p' ${2}_filtered_final.csv >> ${2}_dbsnp_and_one_cosmic.csv #save variants with dbsnp entry but only one for cosmic in a separate file before discarding them
sed -i -e '/^\([^,]*,\)\{13\}"[^\(\.\)"][^,]*,[^,]*OCCURENCE=1([^)]*)",/d' ${2}_filtered_final.csv #delete variants with an entry for dbsnp but only one for cosmic
echo $(($lines - $(wc -l ${2}_filtered_final.csv | cut -f1 -d' ') +1)) " of these calls had an entry for SNP-DB but less than two for cosmic-DB and were dropped" >> ${2}_filter_statistic.txt
echo $(( $(wc -l ${2}_dbsnp_and_one_cosmic.csv | cut -f1 -d' ') -1)) "of these calls had an entry for SNP-DB and exactly one for cosmic-DB and were saved to ${2}_dbsnp_and_one_cosmic.csv" >> ${2}_filter_statistic.txt
lines="$(( $(wc -l ${2}_filtered_final.csv | cut -f1 -d' ') -1))"
echo "$lines calls remain after filtering" >> ${2}_filter_statistic.txt


#######################################################################################################################################################################################################

#add annotation for flanking sequences
get_flanking_sequence.sh ${2}_filtered_final.csv $4 $5


#######################################################################################################################################################################################################

#search for known AML candidate genes in germline calls, somatic_filtered_check and somatic_filtered_normalVarFreq30plus
head -n 1 $1 > ${2}_candidates.csv
grep -i -f $3 ${2}_filtered_final.csv >> ${2}_candidates.csv
echo "" >> ${2}_filter_statistic.txt
echo "$(( $(wc -l ${2}_candidates.csv | cut -f1 -d' ') -1)) calls in ${2}_filtered_final.csv were matched to an AML candidate gene in $3 and saved to ${2}_candidates.csv" >> ${2}_filter_statistic.txt
