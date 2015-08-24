#!/bin/bash

#$1 = sample prefix of fastq files -> everything except terminal "_R*.fastq"

#extract scores from both fastqs per sample
sed -n '4~4p' $1_R1*.fastq > ${1}_bqs.txt
sed -n '4~4p' $1_R2*.fastq >> ${1}_bqs.txt
#to make sure all signs are included, add each once
echo '!"#$%&'\''()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~' >> ${1}_bqs.txt

#insert linebreak after each character
sed -i 's/\(.\)/\1\n/g' ${1}_bqs.txt 

#sort & uniq - LC_ALL=C required for sorting according to ASCII
LC_ALL=C sort ${1}_bqs.txt | sed -e '/^$/d' | sed -e '/^ $/d'> ${1}_bqs_sorted.txt
#uniq -c ${1}_bqs_sorted.txt | sed -e '/^$/d' | sed -e '/^ $/d' > ${1}_bqs_counted.txt

#statistics
# have to subtract 1 because grep finds 11th symbol for min10 (and so on) and then 1 for each min-bqs and 1 because 0-based
echo 'Minimum BQS	Percentage of Bases' > ${1}_bqs_statistics.txt
total=$(( $(wc -l ${1}_bqs_sorted.txt | cut -f1 -d' ') -94))
min10=$(( $(( $total - $(grep -n -m1 ','  ${1}_bqs_sorted.txt | cut -f1 -d':') +12 )) * 100 / $total ))
min20=$(( $(( $total - $(grep -n -m1 '6'  ${1}_bqs_sorted.txt | cut -f1 -d':') +22 )) * 100 / $total ))
min30=$(( $(( $total - $(grep -n -m1 '@'  ${1}_bqs_sorted.txt | cut -f1 -d':') +32 )) * 100 / $total ))
min40=$(( $(( $total - $(grep -n -m1 'J'  ${1}_bqs_sorted.txt | cut -f1 -d':') +42 )) * 100 / $total ))
min50=$(( $(( $total - $(grep -n -m1 'T'  ${1}_bqs_sorted.txt | cut -f1 -d':') +52 )) * 100 / $total ))
min60=$(( $(( $total - $(grep -n -m1 '\^' ${1}_bqs_sorted.txt | cut -f1 -d':') +62 )) * 100 / $total ))

echo "10	$min10" >> ${1}_bqs_statistics.txt
echo "20	$min20" >> ${1}_bqs_statistics.txt
echo "30	$min30" >> ${1}_bqs_statistics.txt
echo "40	$min40" >> ${1}_bqs_statistics.txt
echo "50	$min50" >> ${1}_bqs_statistics.txt
echo "60	$min60" >> ${1}_bqs_statistics.txt
