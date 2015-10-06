#!/bin/bash


IN="$1"


# extract circRNAs from STAR chimeric output
awk -v OFS='\t' -v INDEX=1 '$1==$4 && $3==$6 && $7>=0 && (($3=="-" && $5>$2 && $5-$2<1000000) || ($3=="+" && $2>$5 && $2-$5<1000000)) {print $0,INDEX; INDEX++;}' $IN | sort -k1,1 -k2,2n > ${IN}_circs.txt


######

# calculate start and stop positions of the chimeric donor segment, acceptor segment and their paired mate
FILE="$(cat "${IN}_circs.txt")"

# + strand, case1 (inter-mate gap within segment2):
# Apppend: SS1 (Start segment1 = Donor segment), LS1 (Length segment1), SS2a, LS2a, SS2b,LS2b, IMG (Inter-mate Gap)
P1="$(echo "$FILE" | awk -v OFS='\t' '$3 == "+" && $14 ~ /p/ {print $0}')"
P1_SS1="$(echo "$P1" | cut -f11)"
P1_LS1="$(echo "$P1" | cut -f12 | sed -e 's/[MDN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}')"
P1_SS2a="$(echo "$P1" | cut -f13)"
P1_LS2a="$(echo "$P1" | cut -f14 | sed -e 's/[0-9-]\+p.*$//' -e 's/[MDN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}')"
P1_IMG="$(echo "$P1" | cut -f14 | sed -e 's/.*[MIDNSH]\([0-9-]\+\)p.*$/\1/')"
# P1_LS2b = Mate2_Length
P1_LS2b="$(echo "$P1" | cut -f14 | sed -e 's/.*p\(.*\)$/\1/' -e 's/[MDN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}')"

# calculate donor start/stop, acceptor start/stop and second mate start/stop
# CAREFUL: start (inclusive) + segment length = end (exclusive)!!!
P1_POSCOORDS="$(paste <(echo "$P1") <(paste <(echo "$P1_SS1") <(echo "$P1_LS1") <(echo "$P1_SS2a") <(echo "$P1_LS2a") <(echo "$P1_IMG") <(echo "$P1_LS2b") | awk -v OFS='\t' '{print $1,$1+$2,$3,$3+$4,$3+$4+$5,$3+$4+$5+$6}'))"


# + strand, case2 (intermate-gap within segment1):
P2="$(echo "$FILE" | awk -v OFS='\t' '$3 == "+" && $12 ~ /p/ {print $0}')"
P2_SS1a="$(echo "$P2" | cut -f11)"
# P2_LS1a = Mate2_Length
P2_LS1a="$(echo "$P2" | cut -f12 | sed -e 's/[0-9-]\+p.*$//' -e 's/[MDN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}')"
P2_IMG="$(echo "$P2" | cut -f12 | sed -e 's/.*[MIDNSH]\([0-9-]\+\)p.*$/\1/')"
P2_LS1b="$(echo "$P2" | cut -f12 | sed -e 's/.*p\(.*\)$/\1/' -e 's/[MDN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}')"
P2_SS2="$(echo "$P2" | cut -f13)"
P2_LS2="$(echo "$P2" | cut -f14 | sed -e 's/[MDN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}')"

# calculate donor start/stop, acceptor start/stop and second mate start/stop
P2_POSCOORDS="$(paste <(echo "$P2") <(paste <(echo "$P2_SS1a") <(echo "$P2_LS1a") <(echo "$P2_IMG") <(echo "$P2_LS1b") <(echo "$P2_SS2") <(echo "$P2_LS2") | awk -v OFS='\t' '{print $1+$2+$3,$1+$2+$3+$4,$5,$5+$6,$1,$1+$2}'))"


# - strand, case1 (inter-mate gap within segment2):
N1="$(echo "$FILE" | awk -v OFS='\t' '$3 == "-" && $14 ~ /p/ {print $0}')"
N1_SS1="$(echo "$N1" | cut -f11)"
N1_LS1="$(echo "$N1" | cut -f12 | sed -e 's/[MDN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}')"
N1_SS2a="$(echo "$N1" | cut -f13)"
# N1_LS2a = Mate2_Length
N1_LS2a="$(echo "$N1" | cut -f14 | sed -e 's/[0-9-]\+p.*$//' -e 's/[MDN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}')"
N1_IMG="$(echo "$N1" | cut -f14 | sed -e 's/.*[MIDNSH]\([0-9-]\+\)p.*$/\1/')"
N1_LS2b="$(echo "$N1" | cut -f14 | sed -e 's/.*p\(.*\)$/\1/' -e 's/[MDN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}')"

# calculate donor start/stop, acceptor start/stop and second mate start/stop
N1_POSCOORDS="$(paste <(echo "$N1") <(paste <(echo "$N1_SS1") <(echo "$N1_LS1") <(echo "$N1_SS2a") <(echo "$N1_LS2a") <(echo "$N1_IMG") <(echo "$N1_LS2b") | awk -v OFS='\t' '{print $1,$1+$2,$3+$4+$5,$3+$4+$5+$6,$3,$3+$4}'))"


# - strand, case2 (intermate-gap within segment1):
N2="$(echo "$FILE" | awk -v OFS='\t' '$3 == "-" && $12 ~ /p/ {print $0}')"
N2_SS1a="$(echo "$N2" | cut -f11)"
N2_LS1a="$(echo "$N2" | cut -f12 | sed -e 's/[0-9-]\+p.*$//' -e 's/[MDN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}')"
N2_IMG="$(echo "$N2" | cut -f12 | sed -e 's/.*[MIDNSH]\([0-9-]\+\)p.*$/\1/')"
# N2_LS1b = Mate2_Length
N2_LS1b="$(echo "$N2" | cut -f12 | sed -e 's/.*p\(.*\)$/\1/' -e 's/[MDN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}')"
N2_SS2="$(echo "$N2" | cut -f13)"
N2_LS2="$(echo "$N2" | cut -f14 | sed -e 's/[MDN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}')"

# calculate donor start/stop, acceptor start/stop and second mate start/stop
N2_POSCOORDS="$(paste <(echo "$N2") <(paste <(echo "$N2_SS1a") <(echo "$N2_LS1a") <(echo "$N2_IMG") <(echo "$N2_LS1b") <(echo "$N2_SS2") <(echo "$N2_LS2") | awk -v OFS='\t' '{print $1,$1+$2,$5,$5+$6,$1+$2+$3,$1+$2+$3+$4}'))"

# save all of these coords into one file
cat <(echo "$P1_POSCOORDS") <(echo "$P2_POSCOORDS") <(echo "$N1_POSCOORDS") <(echo "$N2_POSCOORDS") > ${IN}_circsWcoords.txt


#######

# count supporting reads of each junction and output as vector of equal length as input
cut -f1-2,5 "${IN}_circsWcoords.txt" | uniq -c | awk '{ for (i=1; i<= $1; i++) print $1}' > ${IN}_circsSupportingReadsCount.txt

# combine read counts and junction file
paste "${IN}_circsSupportingReadsCount.txt" "${IN}_circsWcoords.txt" | sort -k 16b,16 > ${IN}_circsWcounts.txt

# create BED files for exon annotation (chr,exon base adjacent to backsplice,exon base two coordinates away from backplice (interval end coordinate is not included in BED analysis))
# consider splice donors and acceptors on different strands (assign complementary indixes (INDEX & -INDEX))
# add padding according to left and right flanking repeats ($9 and $10) and remember the shift applied to correct other coords later for beyond/within backsplice testing
# remember that BED is 0-based while STAR coordinates are 1-based (subtract all start coords by 1! (end is inclusive for STAR but not for BED))
awk -v OFS='\t' -v INDEX=0 '$4=="+" {for(i=$3-$9-2;i<$3+$10-1;i++){INDEX++; print $2,i,i+1,$4,i-$3+2,$16,INDEX}}' ${IN}_circsWcounts.txt  > ${IN}_donor+.bed
awk -v OFS='\t' -v INDEX=0 '$4=="-" {for(i=$3-$10;i<$3+$9+1;i++){INDEX++; print $2,i,i+1,$4,i-$3,$16,INDEX*-1}}' ${IN}_circsWcounts.txt  > ${IN}_donor-.bed
awk -v OFS='\t' -v INDEX=0 '$4=="+" {for(i=$6-$9;i<$6+$10+1;i++){INDEX++; print $5,i,i+1,$4,i-$6,$16,INDEX}}' ${IN}_circsWcounts.txt > ${IN}_acceptor+.bed
awk -v OFS='\t' -v INDEX=0 '$4=="-" {for(i=$6-$10-2;i<$6+$9-1;i++){INDEX++; print $5,i,i+1,$4,i-$6+2,$16,INDEX*-1}}' ${IN}_circsWcounts.txt  > ${IN}_acceptor-.bed

### ADDED COLUMN TO donor+/- & acceptor+/- BEDs -> adjust downstream column indices if necessary!


# annotate splice donors and acceptors with exon and gene information
/NGS/links/bedtools/intersectBed -wa -wb -a ${IN}_donor+.bed -b /NGS/known_sites/hg19/gencode.v19.exons.toUCSC.END.uniq.bed > ${IN}_exonicDonor.bed
/NGS/links/bedtools/intersectBed -wa -wb -a ${IN}_donor-.bed -b /NGS/known_sites/hg19/gencode.v19.exons.toUCSC.START.uniq.bed >> ${IN}_exonicDonor.bed
/NGS/links/bedtools/intersectBed -wa -wb -a ${IN}_acceptor+.bed -b /NGS/known_sites/hg19/gencode.v19.exons.toUCSC.START.uniq.bed > ${IN}_exonicAcceptor.bed
/NGS/links/bedtools/intersectBed -wa -wb -a ${IN}_acceptor-.bed -b /NGS/known_sites/hg19/gencode.v19.exons.toUCSC.END.uniq.bed >> ${IN}_exonicAcceptor.bed

# sort results as required by join
sort -k6b,6 ${IN}_exonicDonor.bed > ${IN}_exonicDonorSorted.bed
sort -k6b,6 ${IN}_exonicAcceptor.bed > ${IN}_exonicAcceptorSorted.bed

# merge corresponding donors and acceptors after annotation based on INDEX
join -j6 -t'	' ${IN}_exonicDonorSorted.bed ${IN}_exonicAcceptorSorted.bed > ${IN}_exonicJunctions.bed

# extract intra- & intergenic junctions
awk -v OFS='\t' '$10==$19' ${IN}_exonicJunctions.bed > ${IN}_exonicJunctionsIntragenic.bed
awk -v OFS='\t' '$10!=$19' ${IN}_exonicJunctions.bed > ${IN}_exonicJunctionsIntergenic.bed

# extract ambiguous and nonambiguous annotations
grep -wf <(cut -f6 ${IN}_exonicJunctionsIntragenic.bed | sort | uniq -c | awk '$1 > 1 {print $2}') ${IN}_exonicJunctionsIntragenic.bed > ${IN}_exonicJunctionsIntragenic_ambiguous.bed 
grep -wf <(cut -f6 ${IN}_exonicJunctionsIntragenic.bed | sort | uniq -c | awk '$1 == 1 {print $2}') ${IN}_exonicJunctionsIntragenic.bed | sort -k 6b,6 > ${IN}_exonicJunctionsIntragenic_nonambiguous.bed 
# join annotated exonic intragenic junctions with original circRNA junctions based on first index
join -1 16 -2 6 -t'	' ${IN}_circsWcoords.txt ${IN}_exonicJunctionsIntragenic_nonambiguous.bed > ${IN}_circsAnnotated.txt

# create nicely formatted final output file with header
#echo -e "supportingReads\tchr\tspliceDonor\tspliceAcceptor\tstrand\tspliceSignal\t5'shift\t3'shift\tsupportingReadID\tdonorSegmentStart\tdonorSegmentCIGAR\tacceptorSegmentStart\tacceptorSegmentCIGAR\tgeneSymbol;geneID" > ${IN}_circsAnnotatedFinal.txt
awk -v OFS='\t' '{print $2,$3,$4,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$25}' ${IN}_circsAnnotated.txt | sort -k1,1nr -k2,2V -k3,3n -k4,4n > ${IN}_circsAnnotatedFinal.txt

# filter out paired-end reads that span regions beyond the backsplice
# for that, calculate the reference length spanned by each segment
# also add an index to number unique junctions
echo -e "Index\tsupportingReads\tchr\tspliceDonor\tspliceAcceptor\tstrand\tspliceSignal\t5'shift\t3'shift\tsupportingReadID\tdonorSegmentStart\tdonorSegmentCIGAR\tacceptorSegmentStart\tacceptorSegmentCIGAR\tgeneSymbol;geneID\tdonorSegmentLength\tacceptorSegmentLength" > ${IN}_circsAnnotatedFinal_withinBS.txt
echo -e "supportingReads\tchr\tspliceDonor\tspliceAcceptor\tstrand\tspliceSignal\t5'shift\t3'shift\tsupportingReadID\tdonorSegmentStart\tdonorSegmentCIGAR\tacceptorSegmentStart\tacceptorSegmentCIGAR\tgeneSymbol;geneID\tdonorSegmentLength\tacceptorSegmentLength" > ${IN}_circsAnnotatedFinal_beyondBS.txt

paste <(cut -f2-4 ${IN}_circsAnnotatedFinal.txt | uniq -c | awk '{ for (i=1; i<= $1; i++) print NR}')  ${IN}_circsAnnotatedFinal.txt <(cut -f11 ${IN}_circsAnnotatedFinal.txt | sed -e 's/[MDpN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}') <(cut -f13 ${IN}_circsAnnotatedFinal.txt | sed -e 's/[MDpN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}') | awk '($6 == "+" && $11 >= $5 && $13 + $17 <=$4) || ($6 == "-" && $13 >= $4 && $11 + $16 <= $5)' >> ${IN}_circsAnnotatedFinal_withinBS.txt 

paste ${IN}_circsAnnotatedFinal.txt <(cut -f11 ${IN}_circsAnnotatedFinal.txt | sed -e 's/[MDpN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}') <(cut -f13 ${IN}_circsAnnotatedFinal.txt | sed -e 's/[MDpN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}') | awk '!(($5 == "+" && $10 >= $4 && $12 + $16 <=$3) || ($5 == "-" && $12 >= $3 && $10 + $15 <= $4))' >> ${IN}_circsAnnotatedFinal_beyondBS.txt


# extract splice donors and acceptors not overlapping any exon
/NGS/links/bedtools/intersectBed -wa -v -a ${IN}_donor+.bed -b /NGS/known_sites/hg19/gencode.v19.exons.toUCSC.END.uniq.bed > ${IN}_nonExonicDonor.bed
/NGS/links/bedtools/intersectBed -wa -v -a ${IN}_donor-.bed -b /NGS/known_sites/hg19/gencode.v19.exons.toUCSC.START.uniq.bed >> ${IN}_nonExonicDonor.bed
/NGS/links/bedtools/intersectBed -wa -v -a ${IN}_acceptor+.bed -b /NGS/known_sites/hg19/gencode.v19.exons.toUCSC.START.uniq.bed > ${IN}_nonExonicAcceptor.bed
/NGS/links/bedtools/intersectBed -wa -v -a ${IN}_acceptor-.bed -b /NGS/known_sites/hg19/gencode.v19.exons.toUCSC.END.uniq.bed >> ${IN}_nonExonicAcceptor.bed


# delete temporary files
rm -f ${IN}_*Count*
rm -f ${IN}_*count*
#rm -f ${IN}_*+*
#rm -f ${IN}_*-*
#rm -f ${IN}_*exonicDonor.bed
rm -f ${IN}_*exonicAcceptor.bed
rm -f ${IN}_*Sorted.bed
#rm -f ${IN}_exonicJunctions.bed
rm -f ${IN}_*genic.bed
rm -f ${IN}_exonicJunctionsIntragenic_nonambiguous.bed
rm -f ${IN}_circsAnnotated.txt
#rm -f ${IN}_circsAnnotatedFinal.txt
