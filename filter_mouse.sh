#!/bin/sh

####
# T.J.Blätte
# 2015
####
#
# Filter the mutations output by our bpipe
#       pipelines and annotate them with
#       flanking reference sequence of a
#       given length.
#
# Args:
#   IN: Input CSV file to be filtered and
#       annotated, as output by our pipelines.
#   OUT: Output prefix for filtered and
#       annotated files.
#   CANDIDATES: List of candidate genes,
#       quoted, one per line.
#   REF: Reference genome used during alignment
#       and variant calling. This is only used
#       when variants are annotated offline.
#       (Currently variants are annotated online
#       but the parameter has to be passed for
#       compatibility reasons, since all parameters
#       are positional.
#   FLANKING_SEQ_LEN: Number of bases of reference
#       flanking sequence to annotate on either
#       side of the variants.
#
# Output:
#   $OUT_*filtered*csv: Filtered CSV files at various steps
#       of filtering.
#   $OUT_filter_statistics.txt: Stats on the number of
#       mutations filtered and retained at each step.
#       of the collected and merged statistics.
#   $OUT_candidates.txt: Separate file with variants
#       affecting the provided candidate genes of interest.
#
####


IN="$1"
OUT="$2"
CANDIDATES="$3"
REF="$4"
FLANKING_SEQ_LEN="$5"

#columns in the csv that is to be sorted minus 1
somatic_status=25 #26-1
genomicSuperDups=10 #11-1
dbsnp=11 #12-1
tumor_reads2=17 #18-1
tumor_reads2_plus=25 #26-1 -> do not forget to subract 5 columns for snv selection below!
tumor_reads2_minus=26 #27-1 -> same here
normal_var_freq=14 #15-1


#Variable to save line-number for statistic calculation
lines=0

echo "Filter statistics for $IN" > ${OUT}_filter_statistic.txt
echo "" >> ${OUT}_filter_statistic.txt

#filter exonic, splicing and exonic;splicing (but not ncRNA_splicing, _exonic,...)
head -n 1 $IN > ${OUT}_filtered.csv
grep -e 'exonic' -e 'splicing' $IN | grep -v '_splicing' | grep -v '_exonic' >> ${OUT}_filtered.csv
echo "$(( $(wc -l ${OUT}_filtered.csv | cut -f1 -d' ') -1)) out of a total of $(wc -l $IN | cut -f1 -d' ') calls were classified as exonic, splicing or exonic;splicing and saved to ${OUT}_filtered.csv" >> ${OUT}_filter_statistic.txt
echo "" >> ${OUT}_filter_statistic.txt

#filter Somatic calls and remove synonymous variants
head -n 1 $IN > ${OUT}_filtered_somatic.csv
grep 'Somatic' ${OUT}_filtered.csv | grep -v '"synonymous' >> ${OUT}_filtered_somatic.csv
echo "$(grep -c 'Somatic' ${OUT}_filtered.csv) of these calls were classified as Somatic" >> ${OUT}_filter_statistic.txt
echo "$(grep 'Somatic' ${OUT}_filtered.csv | grep -c '"synonymous') of these were synonymous and dropped" >> ${OUT}_filter_statistic.txt
lines="$(( $(wc -l ${OUT}_filtered_somatic.csv | cut -f1 -d' ') -1))"
echo "$lines calls remain after filtering" >> ${OUT}_filter_statistic.txt

#take out reads with tumor_reads2 < 4
#sed '/^\([^,]*,\)\{17\}[0-3],/d' ${OUT}_filtered_somatic.csv > tumorreads2min4
sed -i -e '/^\([^,]*,\)\{17\}"[0-3]",/d' ${OUT}_filtered_somatic.csv
echo $(($lines - $(wc -l ${OUT}_filtered_somatic.csv | cut -f1 -d' ') +1)) " of these calls had tumor_reads2 < 4 and were dropped" >> ${OUT}_filter_statistic.txt
lines=$(( $(wc -l ${OUT}_filtered_somatic.csv | cut -f1 -d' ') -1))
echo "$lines calls remain after filtering" >> ${OUT}_filter_statistic.txt

#take out calls that have an entry for GenomicSuperDups-DB
#sed '/^\([^,]*,\)\{10\}"[^\.]/d' tumorreads2min4 > noSuperDups
sed -i -e '/^\([^,]*,\)\{10\}"[^\.]/d' ${OUT}_filtered_somatic.csv
echo $(($lines - $(wc -l ${OUT}_filtered_somatic.csv | cut -f1 -d' ') +1)) " of these calls had an entry for GenomicSuperDups-DB and were dropped" >> ${OUT}_filter_statistic.txt
lines="$(( $(wc -l ${OUT}_filtered_somatic.csv | cut -f1 -d' ') -1))"
echo "$lines calls remain after filtering" >> ${OUT}_filter_statistic.txt

##take out calls that have an entry for dbsnp and not for cosmic
##sed '/^\([^,]*,\)\{14\}"[^\(\.\)"][^,]*,"\."/d' noSuperDups > nodbsnpbutcosmic
#sed -i -e '/^\([^,]*,\)\{13\}"[^\(\.\)"][^,]*,"\."/d' ${OUT}_filtered_somatic.csv
#echo $(($lines - $(wc -l ${OUT}_filtered_somatic.csv | cut -f1 -d' ') +1)) " of these calls had an entry for SNP-DB but not for cosmic-DB and were dropped" >> ${OUT}_filter_statistic.txt
#lines="$(( $(wc -l ${OUT}_filtered_somatic.csv | cut -f1 -d' ') -1))"
#echo "$lines calls remain after filtering" >> ${OUT}_filter_statistic.txt

#take out snps with 0 tumor_reads2_plus and 0 tumor_reads2_minus
#sed '/^\([^,]*,\)\{3\}[A-Z],[A-Z],\([^,]*,\)\{20\}0/d' noSuperDups > noplus0
sed -i -e '/^\([^,]*,\)\{3\}"[A-Z]","[A-Z]",\([^,]*,\)\{20\}"0"/d' ${OUT}_filtered_somatic.csv
echo $(($lines - $(wc -l ${OUT}_filtered_somatic.csv | cut -f1 -d' ') +1)) " of these calls had tumor_reads2_plus = 0 and were dropped" >> ${OUT}_filter_statistic.txt
lines="$(( $(wc -l ${OUT}_filtered_somatic.csv | cut -f1 -d' ') -1))"
echo "$lines calls remain after filtering" >> ${OUT}_filter_statistic.txt

#sed '/^\([^,]*,\)\{3\}[A-Z],[A-Z],\([^,]*,\)\{21\}0/d' noplus0 > noplusminus0
sed -i -e '/^\([^,]*,\)\{3\}"[A-Z]","[A-Z]",\([^,]*,\)\{21\}"0"/d' ${OUT}_filtered_somatic.csv
echo $(($lines - $(wc -l ${OUT}_filtered_somatic.csv | cut -f1 -d' ') +1)) " of these calls had tumor_reads2_minus = 0 and were dropped" >> ${OUT}_filter_statistic.txt
lines="$(( $(wc -l ${OUT}_filtered_somatic.csv | cut -f1 -d' ') -1))"
echo "$lines calls remain after filtering and are saved to ${OUT}_filtered_somatic.csv" >> ${OUT}_filter_statistic.txt

#add annotation for flanking sequences
get_flanking_sequence.sh ${OUT}_filtered_somatic.csv $REF $FLANKING_SEQ_LEN

#split files according to normal_var_freq
head -n 1 $IN > ${OUT}_filtered_somatic_true.csv
sed -n '/^\([^,]*,\)\{14\}"[0-7]\(\.[0-9]\+\)\?%",/p' ${OUT}_filtered_somatic.csv >> ${OUT}_filtered_somatic_true.csv
echo "$(( $(wc -l ${OUT}_filtered_somatic_true.csv | cut -f1 -d' ') -1)) of these calls have normal_variant_frequency < 8 and were saved to ${OUT}_filtered_somatic_true.csv" >> ${OUT}_filter_statistic.txt

head -n 1 $IN > ${OUT}_filtered_somatic_check.csv
sed -n '/^\([^,]*,\)\{14\}"[8-9]\(\.[0-9]\+\)\?%",/p' ${OUT}_filtered_somatic.csv >> ${OUT}_filtered_somatic_check.csv
sed -n '/^\([^,]*,\)\{14\}"[1-2][0-9]\(\.[0-9]\+\)\?%",/p' ${OUT}_filtered_somatic.csv >> ${OUT}_filtered_somatic_check.csv
echo "$(( $(wc -l ${OUT}_filtered_somatic_check.csv | cut -f1 -d' ') -1)) of these calls have 7 < normal_variant_frequency < 30 and were saved to ${OUT}_filtered_somatic_check.csv" >> ${OUT}_filter_statistic.txt

head -n 1 $IN > ${OUT}_filtered_somatic_normalVarFreq30plus.csv
sed -n '/^\([^,]*,\)\{14\}"[3-9][0-9]\(\.[0-9]\+\)\?%"/p' ${OUT}_filtered_somatic.csv >> ${OUT}_filtered_somatic_normalVarFreq30plus.csv
sed -n '/^\([^,]*,\)\{14\}"100%"/p' ${OUT}_filtered_somatic.csv >> ${OUT}_filtered_somatic_normalVarFreq30plus.csv
echo "$(( $(wc -l ${OUT}_filtered_somatic_normalVarFreq30plus.csv | cut -f1 -d' ') -1)) of these calls have normal_variant_frequency > 29 and were saved to ${OUT}_filtered_somatic_normalVarFreq30plus.csv" >> ${OUT}_filter_statistic.txt





#same for Germline and LOH (except normal_var_freq-Splitting)
#filter Somatic calls and remove synonymous variants
head -n 1 $IN > ${OUT}_filtered_germline.csv
grep 'Germline' ${OUT}_filtered.csv | grep -v '"synonymous' >> ${OUT}_filtered_germline.csv
echo "" >> ${OUT}_filter_statistic.txt
echo "$(grep -c 'Germline' ${OUT}_filtered.csv) of these calls were classified as Germline" >> ${OUT}_filter_statistic.txt
echo "$(grep 'Germline' ${OUT}_filtered.csv | grep -c '"synonymous') of these were synonymous and dropped" >> ${OUT}_filter_statistic.txt
lines="$(( $(wc -l ${OUT}_filtered_germline.csv | cut -f1 -d' ') -1))"
echo "$lines calls remain after filtering" >> ${OUT}_filter_statistic.txt

#take out reads with tumor_reads2 < 4
sed -i -e '/^\([^,]*,\)\{17\}"[0-3]",/d' ${OUT}_filtered_germline.csv
echo $(($lines - $(wc -l ${OUT}_filtered_germline.csv | cut -f1 -d' ') +1)) " of these calls had tumor_reads2 < 4 and were dropped" >> ${OUT}_filter_statistic.txt
lines="$(( $(wc -l ${OUT}_filtered_germline.csv | cut -f1 -d' ') -1))"
echo "$lines calls remain after filtering" >> ${OUT}_filter_statistic.txt

#take out calls that have an entry for GenomicSuperDups-DB
sed -i -e '/^\([^,]*,\)\{10\}"[^\.]/d' ${OUT}_filtered_germline.csv
echo $(($lines - $(wc -l ${OUT}_filtered_germline.csv | cut -f1 -d' ') +1)) " of these calls had an entry for GenomicSuperDups-DB and were dropped" >> ${OUT}_filter_statistic.txt
lines="$(( $(wc -l ${OUT}_filtered_germline.csv | cut -f1 -d' ') -1))"
echo "$lines calls remain after filtering" >> ${OUT}_filter_statistic.txt

##take out calls that have an entry for dbsnp and not for cosmic
#sed -i -e '/^\([^,]*,\)\{13\}"[^\(\.\)"][^,]*,"\."/d' ${OUT}_filtered_germline.csv
#echo $(($lines - $(wc -l ${OUT}_filtered_germline.csv | cut -f1 -d' ') +1)) " of these calls had an entry for SNP-DB but not for cosmic-DB and were dropped" >> ${OUT}_filter_statistic.txt
#lines="$(( $(wc -l ${OUT}_filtered_germline.csv | cut -f1 -d' ') -1))"
#echo "$lines calls remain after filtering" >> ${OUT}_filter_statistic.txt

#take out snps with 0 tumor_reads2_plus and 0 tumor_reads2_minus
sed -i -e '/^\([^,]*,\)\{3\}"[A-Z]","[A-Z]",\([^,]*,\)\{20\}"0"/d' ${OUT}_filtered_germline.csv
echo $(($lines - $(wc -l ${OUT}_filtered_germline.csv | cut -f1 -d' ') +1)) " of these calls had tumor_reads2_plus = 0 and were dropped" >> ${OUT}_filter_statistic.txt
lines="$(( $(wc -l ${OUT}_filtered_germline.csv | cut -f1 -d' ') -1))"
echo "$lines calls remain after filtering" >> ${OUT}_filter_statistic.txt

sed -i -e '/^\([^,]*,\)\{3\}"[A-Z]","[A-Z]",\([^,]*,\)\{21\}"0"/d' ${OUT}_filtered_germline.csv
echo $(($lines - $(wc -l ${OUT}_filtered_germline.csv | cut -f1 -d' ') +1)) " of these calls had tumor_reads2_minus = 0 and were dropped" >> ${OUT}_filter_statistic.txt
lines="$(( $(wc -l ${OUT}_filtered_germline.csv | cut -f1 -d' ') -1))"
echo "$lines calls remain after filtering and are saved to ${OUT}_filtered_germline.csv" >> ${OUT}_filter_statistic.txt

#add annotation for flanking sequences
get_flanking_sequence.sh ${OUT}_filtered_germline.csv $4 $5

#add filtered germline calls with normal_variant_frequency <= 30 to _filtered_somatic_check.csv
prev="$(wc -l ${OUT}_filtered_somatic_check.csv | cut -f1 -d' ')"
sed -n '/^\([^,]*,\)\{14\}"[0-9]\(\.[0-9]\+\)\?%",/p' ${OUT}_filtered_germline.csv >> ${OUT}_filtered_somatic_check.csv
sed -n '/^\([^,]*,\)\{14\}"1[0-9]\(\.[0-9]\+\)\?%",/p' ${OUT}_filtered_germline.csv >> ${OUT}_filtered_somatic_check.csv
sed -n '/^\([^,]*,\)\{14\}"2[0-9]\(\.[0-9]\+\)\?%/p' ${OUT}_filtered_germline.csv >> ${OUT}_filtered_somatic_check.csv
sed -n '/^\([^,]*,\)\{14\}"30%",/p' ${OUT}_filtered_germline.csv >> ${OUT}_filtered_somatic_check.csv
echo "$(( $(wc -l ${OUT}_filtered_somatic_check.csv | cut -f1 -d' ') - $prev)) filtered germline calls had a normal_variant_frequency up to 30% and were also added to ${OUT}_filtered_somatic_check.csv" >> ${OUT}_filter_statistic.txt



#filter LOH calls and remove synonymous variants
head -n 1 $IN > ${OUT}_filtered_LOH.csv
grep 'LOH' ${OUT}_filtered.csv | grep -v '"synonymous' >> ${OUT}_filtered_LOH.csv
echo "" >> ${OUT}_filter_statistic.txt
echo "$(grep -c 'LOH' ${OUT}_filtered.csv) of these calls were classified as LOH" >> ${OUT}_filter_statistic.txt
echo "$(grep 'LOH' ${OUT}_filtered.csv | grep -c '"synonymous') of these were synonymous and dropped" >> ${OUT}_filter_statistic.txt
lines="$(( $(wc -l ${OUT}_filtered_LOH.csv | cut -f1 -d' ') -1))"
echo "$lines calls remain after filtering" >> ${OUT}_filter_statistic.txt

#take out reads with tumor_reads2 < 4
sed -i -e '/^\([^,]*,\)\{17\}"[0-3]",/d' ${OUT}_filtered_LOH.csv
echo $(($lines - $(wc -l ${OUT}_filtered_LOH.csv | cut -f1 -d' ') +1)) " of these calls had tumor_reads2 < 4 and were dropped" >> ${OUT}_filter_statistic.txt
lines="$(( $(wc -l ${OUT}_filtered_LOH.csv | cut -f1 -d' ') -1))"
echo "$lines calls remain after filtering" >> ${OUT}_filter_statistic.txt

#take out calls that have an entry for GenomicSuperDups-DB
sed -i -e '/^\([^,]*,\)\{10\}"[^\.]/d' ${OUT}_filtered_LOH.csv
echo $(($lines - $(wc -l ${OUT}_filtered_LOH.csv | cut -f1 -d' ') +1)) " of these calls had an entry for GenomicSuperDups-DB and were dropped" >> ${OUT}_filter_statistic.txt
lines="$(( $(wc -l ${OUT}_filtered_LOH.csv | cut -f1 -d' ') -1))"
echo "$lines calls remain after filtering" >> ${OUT}_filter_statistic.txt

##take out calls that have an entry for dbsnp and not for cosmic
#sed -i -e '/^\([^,]*,\)\{13\}"[^\(\.\)"][^,]*,"\."/d' ${OUT}_filtered_LOH.csv
#echo $(($lines - $(wc -l ${OUT}_filtered_LOH.csv | cut -f1 -d' ') +1)) " of these calls had an entry for SNP-DB but not for cosmic-DB and were dropped" >> ${OUT}_filter_statistic.txt
#lines="$(( $(wc -l ${OUT}_filtered_LOH.csv | cut -f1 -d' ') -1))"
#echo "$lines calls remain after filtering" >> ${OUT}_filter_statistic.txt

#take out snps with 0 tumor_reads2_plus and 0 tumor_reads2_minus
sed -i -e '/^\([^,]*,\)\{3\}"[A-Z]","[A-Z]",\([^,]*,\)\{20\}"0"/d' ${OUT}_filtered_LOH.csv
echo $(($lines - $(wc -l ${OUT}_filtered_LOH.csv | cut -f1 -d' ') +1)) " of these calls had tumor_reads2_plus = 0 and were dropped" >> ${OUT}_filter_statistic.txt
lines="$(( $(wc -l ${OUT}_filtered_LOH.csv | cut -f1 -d' ') -1))"
echo "$lines calls remain after filtering" >> ${OUT}_filter_statistic.txt

sed -i -e '/^\([^,]*,\)\{3\}"[A-Z]","[A-Z]",\([^,]*,\)\{21\}"0"/d' ${OUT}_filtered_LOH.csv
echo $(($lines - $(wc -l ${OUT}_filtered_LOH.csv | cut -f1 -d' ') +1)) " of these calls had tumor_reads2_minus = 0 and were dropped" >> ${OUT}_filter_statistic.txt
lines="$(( $(wc -l ${OUT}_filtered_LOH.csv | cut -f1 -d' ') -1))"
echo "$lines calls remain after filtering and are saved to ${OUT}_filtered_LOH.csv" >> ${OUT}_filter_statistic.txt

#add annotation for flanking sequences
get_flanking_sequence.sh ${OUT}_filtered_LOH.csv $4 $5


#search for known AML candidate genes in germline calls, somatic_filtered_check and somatic_filtered_normalVarFreq30plus
head -n 1 $IN > ${OUT}_candidates.csv
grep -i -f $CANDIDATES ${OUT}_filtered_germline.csv > ${OUT}_PREcandidates.csv
grep -i -f $CANDIDATES ${OUT}_filtered_somatic_check.csv >> ${OUT}_PREcandidates.csv
grep -i -f $CANDIDATES ${OUT}_filtered_somatic_normalVarFreq30plus.csv >> ${OUT}_PREcandidates.csv
#remove duplicates (contained in both germline_filtered and because of normal_var_freq in somatic_check)
sort ${OUT}_PREcandidates.csv | uniq >> ${OUT}_candidates.csv
echo "" >> ${OUT}_filter_statistic.txt
echo "$(( $(wc -l ${OUT}_candidates.csv | cut -f1 -d' ') -1)) calls in ${OUT}_filtered_germline.csv, ${OUT}_filtered_somatic_check and ${OUT}_filtered_somatic_normalVarFreq30plus were matched to an AML candidate gene in $CANDIDATES and saved to ${OUT}_candidates.csv" >> ${OUT}_filter_statistic.txt
