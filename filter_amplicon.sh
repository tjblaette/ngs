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


# whether to annotate flanking sequence on- or offline
# --> pass REF as additional parameter to script call below
#     when annotating offline!
GET_FLANKING_SEQ="get_flanking_sequence_online.sh"
#GET_FLANKING_SEQ=get_flanking_sequence.sh


#columns in the csv that is to be sorted minus 1
somatic_status=25 #26-1
genomicSuperDups=10 #11-1
cosmic=14 #15-1
dbsnp=13 #14-1
tumor_reads2_plus=30 #31-1
tumor_reads2_minus=31 #32-1
normal_var_freq=19 #20-1


#Variable to save line-number for statistic calculation
lines=0

echo "Filter statistics for $IN" > ${OUT}_filter_statistic.txt
echo "" >> ${OUT}_filter_statistic.txt

#filter exonic, splicing and exonic;splicing (but not ncRNA_splicing, _exonic,...)
head -n 1 $IN > ${OUT}_filtered.csv
grep -e 'exonic' -e 'splicing' $IN | grep -v '_splicing' | grep -v '_exonic' >> ${OUT}_filtered.csv
echo "$(( $(wc -l ${OUT}_filtered.csv | cut -f1 -d' ') -1)) out of a total of $(wc -l $IN | cut -f1 -d' ') calls were classified as exonic, splicing or exonic;splicing and saved to ${OUT}_filtered.csv" >> ${OUT}_filter_statistic.txt
echo "" >> ${OUT}_filter_statistic.txt

#remove synonymous variants
head -n 1 $IN > ${OUT}_filtered_final.csv
grep 'chr' ${OUT}_filtered.csv | grep -v '"synonymous' >> ${OUT}_filtered_final.csv
echo "$(grep -c '"synonymous' ${OUT}_filtered.csv) of these were synonymous and dropped" >> ${OUT}_filter_statistic.txt
lines="$(( $(wc -l ${OUT}_filtered_final.csv | cut -f1 -d' ') -1))"
echo "$lines calls remain after filtering" >> ${OUT}_filter_statistic.txt

#take out calls that have an entry for dbsnp and not at least two for cosmic
head -n 1 $IN > ${OUT}_dbsnp_and_one_cosmic.csv
sed -e '/^\([^,]*,\)\{13\}"[^\(\.\)"][^,]*,"\."/d' ${OUT}_filtered_final.csv > ${OUT}_tmp && mv ${OUT}_tmp ${OUT}_filtered_final.csv #delete variants with entry for dbsnp but not cosmic
sed -n -e '/^\([^,]*,\)\{13\}"[^\(\.\)"][^,]*,[^,]*OCCURENCE=1([^)]*)",/p' ${OUT}_filtered_final.csv >> ${OUT}_dbsnp_and_one_cosmic.csv #save variants with dbsnp entry but only one for cosmic in a separate file before discarding them
sed -e '/^\([^,]*,\)\{13\}"[^\(\.\)"][^,]*,[^,]*OCCURENCE=1([^)]*)",/d' ${OUT}_filtered_final.csv > ${OUT}_tmp && mv ${OUT}_tmp ${OUT}_filtered_final.csv #delete variants with an entry for dbsnp but only one for cosmic
echo $(($lines - $(wc -l ${OUT}_filtered_final.csv | cut -f1 -d' ') +1)) " of these calls had an entry for SNP-DB but less than two for cosmic-DB and were dropped" >> ${OUT}_filter_statistic.txt
echo $(( $(wc -l ${OUT}_dbsnp_and_one_cosmic.csv | cut -f1 -d' ') -1)) "of these calls had an entry for SNP-DB and exactly one for cosmic-DB and were saved to ${OUT}_dbsnp_and_one_cosmic.csv" >> ${OUT}_filter_statistic.txt
lines="$(( $(wc -l ${OUT}_filtered_final.csv | cut -f1 -d' ') -1))"
echo "$lines calls remain after filtering" >> ${OUT}_filter_statistic.txt


#add annotation for flanking sequences
$GET_FLANKING_SEQ ${OUT}_filtered_final.csv $FLANKING_SEQ_LEN
#$GET_FLANKING_SEQ ${OUT}_filtered_final.csv $REF $FLANKING_SEQ_LEN


#search for known AML candidate genes in germline calls, somatic_filtered_check and somatic_filtered_normalVarFreq30plus
head -n 1 ${OUT}_filtered_final.csv > ${OUT}_candidates.csv
grep -i -f $CANDIDATES ${OUT}_filtered_final.csv >> ${OUT}_candidates.csv
echo "" >> ${OUT}_filter_statistic.txt
echo "$(( $(wc -l ${OUT}_candidates.csv | cut -f1 -d' ') -1)) calls in ${OUT}_filtered_final.csv were matched to an AML candidate gene in $CANDIDATES and saved to ${OUT}_candidates.csv" >> ${OUT}_filter_statistic.txt
