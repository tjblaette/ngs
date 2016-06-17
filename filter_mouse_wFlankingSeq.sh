#!/bin/sh

#columns in the csv that is to be sorted minus 1
somatic_status=25 #26-1
genomicSuperDups=10 #11-1
#cosmic=14 #15-1
dbsnp=11 #12-1
tumor_reads2=17 #18-1
tumor_reads2_plus=25 #26-1 -> do not forget to subract 5 columns for snv selection below!
tumor_reads2_minus=26 #27-1 -> same here
normal_var_freq=14 #15-1

#Input
#$1 = input csv-file that is to be sorted, filtered and annotated with flanking sequences
#$2 = Prefix for output-files
#$3 = file containing candidate gene list (Format: one gene name in quotation marks per line)
#$4 = reference genome -> same as used for variant calling!
#$5 = length of flanking sequence to print

#Variable to save line-number for statistic calculation
lines=0

echo "Filter statistics for $1" > ${2}_filter_statistic.txt
echo "" >> ${2}_filter_statistic.txt

#filter exonic, splicing and exonic;splicing (but not ncRNA_splicing, _exonic,...)
head -n 1 $1 > ${2}_filtered.csv
grep -e 'exonic' -e 'splicing' $1 | grep -v '_splicing' | grep -v '_exonic' >> ${2}_filtered.csv
echo "$(( $(wc -l ${2}_filtered.csv | cut -f1 -d' ') -1)) out of a total of $(wc -l $1 | cut -f1 -d' ') calls were classified as exonic, splicing or exonic;splicing and saved to ${2}_filtered.csv" >> ${2}_filter_statistic.txt
echo "" >> ${2}_filter_statistic.txt

#filter Somatic calls and remove synonymous variants
head -n 1 $1 > ${2}_filtered_somatic.csv
grep 'Somatic' ${2}_filtered.csv | grep -v '"synonymous' >> ${2}_filtered_somatic.csv
echo "$(grep -c 'Somatic' ${2}_filtered.csv) of these calls were classified as Somatic" >> ${2}_filter_statistic.txt
echo "$(grep 'Somatic' ${2}_filtered.csv | grep -c '"synonymous') of these were synonymous and dropped" >> ${2}_filter_statistic.txt
lines="$(( $(wc -l ${2}_filtered_somatic.csv | cut -f1 -d' ') -1))"
echo "$lines calls remain after filtering" >> ${2}_filter_statistic.txt

#take out reads with tumor_reads2 < 4
#sed '/^\([^,]*,\)\{17\}[0-3],/d' ${2}_filtered_somatic.csv > tumorreads2min4
sed -i -e '/^\([^,]*,\)\{17\}"[0-3]",/d' ${2}_filtered_somatic.csv 
echo $(($lines - $(wc -l ${2}_filtered_somatic.csv | cut -f1 -d' ') +1)) " of these calls had tumor_reads2 < 4 and were dropped" >> ${2}_filter_statistic.txt
lines=$(( $(wc -l ${2}_filtered_somatic.csv | cut -f1 -d' ') -1))
echo "$lines calls remain after filtering" >> ${2}_filter_statistic.txt

#take out calls that have an entry for GenomicSuperDups-DB
#sed '/^\([^,]*,\)\{10\}"[^\.]/d' tumorreads2min4 > noSuperDups
sed -i -e '/^\([^,]*,\)\{10\}"[^\.]/d' ${2}_filtered_somatic.csv
echo $(($lines - $(wc -l ${2}_filtered_somatic.csv | cut -f1 -d' ') +1)) " of these calls had an entry for GenomicSuperDups-DB and were dropped" >> ${2}_filter_statistic.txt
lines="$(( $(wc -l ${2}_filtered_somatic.csv | cut -f1 -d' ') -1))"
echo "$lines calls remain after filtering" >> ${2}_filter_statistic.txt

##take out calls that have an entry for dbsnp and not for cosmic
##sed '/^\([^,]*,\)\{14\}"[^\(\.\)"][^,]*,"\."/d' noSuperDups > nodbsnpbutcosmic
#sed -i -e '/^\([^,]*,\)\{13\}"[^\(\.\)"][^,]*,"\."/d' ${2}_filtered_somatic.csv
#echo $(($lines - $(wc -l ${2}_filtered_somatic.csv | cut -f1 -d' ') +1)) " of these calls had an entry for SNP-DB but not for cosmic-DB and were dropped" >> ${2}_filter_statistic.txt
#lines="$(( $(wc -l ${2}_filtered_somatic.csv | cut -f1 -d' ') -1))"
#echo "$lines calls remain after filtering" >> ${2}_filter_statistic.txt

#take out snps with 0 tumor_reads2_plus and 0 tumor_reads2_minus
#sed '/^\([^,]*,\)\{3\}[A-Z],[A-Z],\([^,]*,\)\{20\}0/d' noSuperDups > noplus0
sed -i -e '/^\([^,]*,\)\{3\}"[A-Z]","[A-Z]",\([^,]*,\)\{20\}"0"/d' ${2}_filtered_somatic.csv
echo $(($lines - $(wc -l ${2}_filtered_somatic.csv | cut -f1 -d' ') +1)) " of these calls had tumor_reads2_plus = 0 and were dropped" >> ${2}_filter_statistic.txt
lines="$(( $(wc -l ${2}_filtered_somatic.csv | cut -f1 -d' ') -1))"
echo "$lines calls remain after filtering" >> ${2}_filter_statistic.txt

#sed '/^\([^,]*,\)\{3\}[A-Z],[A-Z],\([^,]*,\)\{21\}0/d' noplus0 > noplusminus0
sed -i -e '/^\([^,]*,\)\{3\}"[A-Z]","[A-Z]",\([^,]*,\)\{21\}"0"/d' ${2}_filtered_somatic.csv
echo $(($lines - $(wc -l ${2}_filtered_somatic.csv | cut -f1 -d' ') +1)) " of these calls had tumor_reads2_minus = 0 and were dropped" >> ${2}_filter_statistic.txt
lines="$(( $(wc -l ${2}_filtered_somatic.csv | cut -f1 -d' ') -1))"
echo "$lines calls remain after filtering and are saved to ${2}_filtered_somatic.csv" >> ${2}_filter_statistic.txt

#add annotation for flanking sequences
get_flanking_sequence.sh ${2}_filtered_somatic.csv $4 $5 

#split files according to normal_var_freq
head -n 1 $1 > ${2}_filtered_somatic_true.csv
sed -n '/^\([^,]*,\)\{14\}"[0-7]\(\.[0-9]\+\)\?%",/p' ${2}_filtered_somatic.csv >> ${2}_filtered_somatic_true.csv
echo "$(( $(wc -l ${2}_filtered_somatic_true.csv | cut -f1 -d' ') -1)) of these calls have normal_variant_frequency < 8 and were saved to ${2}_filtered_somatic_true.csv" >> ${2}_filter_statistic.txt

head -n 1 $1 > ${2}_filtered_somatic_check.csv
sed -n '/^\([^,]*,\)\{14\}"[8-9]\(\.[0-9]\+\)\?%",/p' ${2}_filtered_somatic.csv >> ${2}_filtered_somatic_check.csv
sed -n '/^\([^,]*,\)\{14\}"[1-2][0-9]\(\.[0-9]\+\)\?%",/p' ${2}_filtered_somatic.csv >> ${2}_filtered_somatic_check.csv
echo "$(( $(wc -l ${2}_filtered_somatic_check.csv | cut -f1 -d' ') -1)) of these calls have 7 < normal_variant_frequency < 30 and were saved to ${2}_filtered_somatic_check.csv" >> ${2}_filter_statistic.txt

head -n 1 $1 > ${2}_filtered_somatic_normalVarFreq30plus.csv
sed -n '/^\([^,]*,\)\{14\}"[3-9][0-9]\(\.[0-9]\+\)\?%"/p' ${2}_filtered_somatic.csv >> ${2}_filtered_somatic_normalVarFreq30plus.csv
sed -n '/^\([^,]*,\)\{14\}"100%"/p' ${2}_filtered_somatic.csv >> ${2}_filtered_somatic_normalVarFreq30plus.csv
echo "$(( $(wc -l ${2}_filtered_somatic_normalVarFreq30plus.csv | cut -f1 -d' ') -1)) of these calls have normal_variant_frequency > 29 and were saved to ${2}_filtered_somatic_normalVarFreq30plus.csv" >> ${2}_filter_statistic.txt





#same for Germline and LOH (except normal_var_freq-Splitting)
#filter Somatic calls and remove synonymous variants
head -n 1 $1 > ${2}_filtered_germline.csv
grep 'Germline' ${2}_filtered.csv | grep -v '"synonymous' >> ${2}_filtered_germline.csv
echo "" >> ${2}_filter_statistic.txt
echo "$(grep -c 'Germline' ${2}_filtered.csv) of these calls were classified as Germline" >> ${2}_filter_statistic.txt
echo "$(grep 'Germline' ${2}_filtered.csv | grep -c '"synonymous') of these were synonymous and dropped" >> ${2}_filter_statistic.txt
lines="$(( $(wc -l ${2}_filtered_germline.csv | cut -f1 -d' ') -1))"
echo "$lines calls remain after filtering" >> ${2}_filter_statistic.txt

#take out reads with tumor_reads2 < 4
sed -i -e '/^\([^,]*,\)\{17\}"[0-3]",/d' ${2}_filtered_germline.csv
echo $(($lines - $(wc -l ${2}_filtered_germline.csv | cut -f1 -d' ') +1)) " of these calls had tumor_reads2 < 4 and were dropped" >> ${2}_filter_statistic.txt
lines="$(( $(wc -l ${2}_filtered_germline.csv | cut -f1 -d' ') -1))"
echo "$lines calls remain after filtering" >> ${2}_filter_statistic.txt

#take out calls that have an entry for GenomicSuperDups-DB
sed -i -e '/^\([^,]*,\)\{10\}"[^\.]/d' ${2}_filtered_germline.csv
echo $(($lines - $(wc -l ${2}_filtered_germline.csv | cut -f1 -d' ') +1)) " of these calls had an entry for GenomicSuperDups-DB and were dropped" >> ${2}_filter_statistic.txt
lines="$(( $(wc -l ${2}_filtered_germline.csv | cut -f1 -d' ') -1))"
echo "$lines calls remain after filtering" >> ${2}_filter_statistic.txt

##take out calls that have an entry for dbsnp and not for cosmic
#sed -i -e '/^\([^,]*,\)\{13\}"[^\(\.\)"][^,]*,"\."/d' ${2}_filtered_germline.csv
#echo $(($lines - $(wc -l ${2}_filtered_germline.csv | cut -f1 -d' ') +1)) " of these calls had an entry for SNP-DB but not for cosmic-DB and were dropped" >> ${2}_filter_statistic.txt
#lines="$(( $(wc -l ${2}_filtered_germline.csv | cut -f1 -d' ') -1))"
#echo "$lines calls remain after filtering" >> ${2}_filter_statistic.txt

#take out snps with 0 tumor_reads2_plus and 0 tumor_reads2_minus
sed -i -e '/^\([^,]*,\)\{3\}"[A-Z]","[A-Z]",\([^,]*,\)\{20\}"0"/d' ${2}_filtered_germline.csv
echo $(($lines - $(wc -l ${2}_filtered_germline.csv | cut -f1 -d' ') +1)) " of these calls had tumor_reads2_plus = 0 and were dropped" >> ${2}_filter_statistic.txt
lines="$(( $(wc -l ${2}_filtered_germline.csv | cut -f1 -d' ') -1))"
echo "$lines calls remain after filtering" >> ${2}_filter_statistic.txt

sed -i -e '/^\([^,]*,\)\{3\}"[A-Z]","[A-Z]",\([^,]*,\)\{21\}"0"/d' ${2}_filtered_germline.csv
echo $(($lines - $(wc -l ${2}_filtered_germline.csv | cut -f1 -d' ') +1)) " of these calls had tumor_reads2_minus = 0 and were dropped" >> ${2}_filter_statistic.txt
lines="$(( $(wc -l ${2}_filtered_germline.csv | cut -f1 -d' ') -1))"
echo "$lines calls remain after filtering and are saved to ${2}_filtered_germline.csv" >> ${2}_filter_statistic.txt

#add annotation for flanking sequences
get_flanking_sequence.sh ${2}_filtered_germline.csv $4 $5

#add filtered germline calls with normal_variant_frequency <= 30 to _filtered_somatic_check.csv
prev="$(wc -l ${2}_filtered_somatic_check.csv | cut -f1 -d' ')"
sed -n '/^\([^,]*,\)\{14\}"[0-9]\(\.[0-9]\+\)\?%",/p' ${2}_filtered_germline.csv >> ${2}_filtered_somatic_check.csv
sed -n '/^\([^,]*,\)\{14\}"1[0-9]\(\.[0-9]\+\)\?%",/p' ${2}_filtered_germline.csv >> ${2}_filtered_somatic_check.csv
sed -n '/^\([^,]*,\)\{14\}"2[0-9]\(\.[0-9]\+\)\?%/p' ${2}_filtered_germline.csv >> ${2}_filtered_somatic_check.csv
sed -n '/^\([^,]*,\)\{14\}"30%",/p' ${2}_filtered_germline.csv >> ${2}_filtered_somatic_check.csv
echo "$(( $(wc -l ${2}_filtered_somatic_check.csv | cut -f1 -d' ') - $prev)) filtered germline calls had a normal_variant_frequency up to 30% and were also added to ${2}_filtered_somatic_check.csv" >> ${2}_filter_statistic.txt



#filter LOH calls and remove synonymous variants
head -n 1 $1 > ${2}_filtered_LOH.csv
grep 'LOH' ${2}_filtered.csv | grep -v '"synonymous' >> ${2}_filtered_LOH.csv
echo "" >> ${2}_filter_statistic.txt
echo "$(grep -c 'LOH' ${2}_filtered.csv) of these calls were classified as LOH" >> ${2}_filter_statistic.txt
echo "$(grep 'LOH' ${2}_filtered.csv | grep -c '"synonymous') of these were synonymous and dropped" >> ${2}_filter_statistic.txt
lines="$(( $(wc -l ${2}_filtered_LOH.csv | cut -f1 -d' ') -1))" 
echo "$lines calls remain after filtering" >> ${2}_filter_statistic.txt

#take out reads with tumor_reads2 < 4
sed -i -e '/^\([^,]*,\)\{17\}"[0-3]",/d' ${2}_filtered_LOH.csv
echo $(($lines - $(wc -l ${2}_filtered_LOH.csv | cut -f1 -d' ') +1)) " of these calls had tumor_reads2 < 4 and were dropped" >> ${2}_filter_statistic.txt
lines="$(( $(wc -l ${2}_filtered_LOH.csv | cut -f1 -d' ') -1))" 
echo "$lines calls remain after filtering" >> ${2}_filter_statistic.txt

#take out calls that have an entry for GenomicSuperDups-DB
sed -i -e '/^\([^,]*,\)\{10\}"[^\.]/d' ${2}_filtered_LOH.csv
echo $(($lines - $(wc -l ${2}_filtered_LOH.csv | cut -f1 -d' ') +1)) " of these calls had an entry for GenomicSuperDups-DB and were dropped" >> ${2}_filter_statistic.txt
lines="$(( $(wc -l ${2}_filtered_LOH.csv | cut -f1 -d' ') -1))" 
echo "$lines calls remain after filtering" >> ${2}_filter_statistic.txt

##take out calls that have an entry for dbsnp and not for cosmic
#sed -i -e '/^\([^,]*,\)\{13\}"[^\(\.\)"][^,]*,"\."/d' ${2}_filtered_LOH.csv
#echo $(($lines - $(wc -l ${2}_filtered_LOH.csv | cut -f1 -d' ') +1)) " of these calls had an entry for SNP-DB but not for cosmic-DB and were dropped" >> ${2}_filter_statistic.txt
#lines="$(( $(wc -l ${2}_filtered_LOH.csv | cut -f1 -d' ') -1))" 
#echo "$lines calls remain after filtering" >> ${2}_filter_statistic.txt

#take out snps with 0 tumor_reads2_plus and 0 tumor_reads2_minus
sed -i -e '/^\([^,]*,\)\{3\}"[A-Z]","[A-Z]",\([^,]*,\)\{20\}"0"/d' ${2}_filtered_LOH.csv
echo $(($lines - $(wc -l ${2}_filtered_LOH.csv | cut -f1 -d' ') +1)) " of these calls had tumor_reads2_plus = 0 and were dropped" >> ${2}_filter_statistic.txt
lines="$(( $(wc -l ${2}_filtered_LOH.csv | cut -f1 -d' ') -1))" 
echo "$lines calls remain after filtering" >> ${2}_filter_statistic.txt

sed -i -e '/^\([^,]*,\)\{3\}"[A-Z]","[A-Z]",\([^,]*,\)\{21\}"0"/d' ${2}_filtered_LOH.csv
echo $(($lines - $(wc -l ${2}_filtered_LOH.csv | cut -f1 -d' ') +1)) " of these calls had tumor_reads2_minus = 0 and were dropped" >> ${2}_filter_statistic.txt
lines="$(( $(wc -l ${2}_filtered_LOH.csv | cut -f1 -d' ') -1))" 
echo "$lines calls remain after filtering and are saved to ${2}_filtered_LOH.csv" >> ${2}_filter_statistic.txt

#add annotation for flanking sequences
get_flanking_sequence.sh ${2}_filtered_LOH.csv $4 $5


#search for known AML candidate genes in germline calls, somatic_filtered_check and somatic_filtered_normalVarFreq30plus
head -n 1 $1 > ${2}_candidates.csv
grep -i -f $3 ${2}_filtered_germline.csv > ${2}_PREcandidates.csv
grep -i -f $3 ${2}_filtered_somatic_check.csv >> ${2}_PREcandidates.csv
grep -i -f $3 ${2}_filtered_somatic_normalVarFreq30plus.csv >> ${2}_PREcandidates.csv
#remove duplicates (contained in both germline_filtered and because of normal_var_freq in somatic_check)
sort ${2}_PREcandidates.csv | uniq >> ${2}_candidates.csv
echo "" >> ${2}_filter_statistic.txt
echo "$(( $(wc -l ${2}_candidates.csv | cut -f1 -d' ') -1)) calls in ${2}_filtered_germline.csv, ${2}_filtered_somatic_check and ${2}_filtered_somatic_normalVarFreq30plus were matched to an AML candidate gene in $3 and saved to ${2}_candidates.csv" >> ${2}_filter_statistic.txt


#######################################################################################################################################################################################################
#######################################################################################################################################################################################################

#delete temporary files
bpipe_finish.sh
