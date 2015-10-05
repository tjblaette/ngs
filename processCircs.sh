#!/bin/bash


IN="$1"


# extract circRNAs from STAR chimeric output
awk -v OFS='\t' '$1==$4 && $3==$6 && $7>=0 && (($3=="-" && $5>$2 && $5-$2<1000000) || ($3=="+" && $2>$5 && $2-$5<1000000)) {print $0,NR}' $IN | sort -k1,1 -k2,2n > ${IN}_circs.txt

# count supporting reads of each junction and output as vector of equal length as input
cut -f1-2,5 ${IN}_circs.txt | uniq -c | awk '{ for (i=1; i<= $1; i++) print $1}' > ${IN}_circsSupportingReadsCount.txt
     
# combine read counts and junction file
paste ${IN}_circsSupportingReadsCount.txt ${IN}_circs.txt | sort -k 16b,16 > ${IN}_circsWcounts.txt

# create BED files for exon annotation (chr,exon base adjacent to backsplice,exon base two coordinates away from backplice (interval end coordinate is not included in BED analysis))
# consider splice donors and acceptors on different strands
# add padding according to left and right flanking repeats ($9 and $10)
# remember that BED is 0-based while STAR coordinates are 1-based (subtract all by 1!)
awk -v OFS='\t' -v INDEX=0 '$4=="+" {for(i=$3-$9-2;i<$3+$10-1;i++){INDEX++; print $2,i,i+1,$4,$16,INDEX}}' ${IN}_circsWcounts.txt  > ${IN}_donor+.bed
awk -v OFS='\t' -v INDEX=0 '$4=="-" {for(i=$3-$10;i<$3+$9+1;i++){INDEX++; print $2,i,i+1,$4,$16,INDEX*-1}}' ${IN}_circsWcounts.txt  > ${IN}_donor-.bed
awk -v OFS='\t' -v INDEX=0 '$4=="+" {for(i=$6-$9;i<$6+$10+1;i++){INDEX++; print $5,i,i+1,$4,$16,INDEX}}' ${IN}_circsWcounts.txt > ${IN}_acceptor+.bed
awk -v OFS='\t' -v INDEX=0 '$4=="-" {for(i=$6-$10-2;i<$6+$9-1;i++){INDEX++; print $5,i,i+1,$4,$16,INDEX*-1}}' ${IN}_circsWcounts.txt  > ${IN}_acceptor-.bed

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
join -1 16 -2 6 -t'	' ${IN}_circsWcounts.txt ${IN}_exonicJunctionsIntragenic_nonambiguous.bed > ${IN}_circsAnnotated.txt

# create nicely formatted final output file with header
#echo -e "supportingReads\tchr\tspliceDonor\tspliceAcceptor\tstrand\tspliceSignal\t5'shift\t3'shift\tsupportingReadID\tdonorSegmentStart\tdonorSegmentCIGAR\tacceptorSegmentStart\tacceptorSegmentCIGAR\tgeneSymbol;geneID" > ${IN}_circsAnnotatedFinal.txt
awk -v OFS='\t' '{print $2,$3,$4,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$25}' ${IN}_circsAnnotated.txt | sort -k1,1nr -k2,2V -k3,3n -k6,6n > ${IN}_circsAnnotatedFinal.txt

# filter out paired-end reads that span regions beyond the backsplice
# first calculate the reference length spanned by each segment
echo -e "supportingReads\tchr\tspliceDonor\tspliceAcceptor\tstrand\tspliceSignal\t5'shift\t3'shift\tsupportingReadID\tdonorSegmentStart\tdonorSegmentCIGAR\tacceptorSegmentStart\tacceptorSegmentCIGAR\tgeneSymbol;geneID\tdonorSegmentLength\tacceptorSegmentLength" > ${IN}_circsAnnotatedFinal_withinBS.txt
echo -e "supportingReads\tchr\tspliceDonor\tspliceAcceptor\tstrand\tspliceSignal\t5'shift\t3'shift\tsupportingReadID\tdonorSegmentStart\tdonorSegmentCIGAR\tacceptorSegmentStart\tacceptorSegmentCIGAR\tgeneSymbol;geneID\tdonorSegmentLength\tacceptorSegmentLength" > ${IN}_circsAnnotatedFinal_beyondBS.txt
paste ${IN}_circsAnnotatedFinal.txt <(cut -f11 ${IN}_circsAnnotatedFinal.txt | sed -e 's/[MDpN]/\t/g' -e 's/[1-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}') <(cut -f13 ${IN}_circsAnnotatedFinal.txt | sed -e 's/[MDpN]/\t/g' -e 's/[1-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}') | awk '($5 == "+" && $10 >= $4 && $12 + $16 <=$3) || ($5 == "-" && $12 >= $3 && $10 + $15 <= $4)' >> ${IN}_circsAnnotatedFinal_withinBS.txt 

paste ${IN}_circsAnnotatedFinal.txt <(cut -f11 ${IN}_circsAnnotatedFinal.txt | sed -e 's/[MDpN]/\t/g' -e 's/[1-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}') <(cut -f13 ${IN}_circsAnnotatedFinal.txt | sed -e 's/[MDpN]/\t/g' -e 's/[1-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}') | awk '!(($5 == "+" && $10 >= $4 && $12 + $16 <=$3) || ($5 == "-" && $12 >= $3 && $10 + $15 <= $4))' >> ${IN}_circsAnnotatedFinal_beyondBS.txt

# extract splice donors and acceptors not overlapping any exon
/NGS/links/bedtools/intersectBed -wa -v -a ${IN}_donor+.bed -b /NGS/known_sites/hg19/gencode.v19.exons.toUCSC.END.uniq.bed > ${IN}_nonExonicDonor.bed
/NGS/links/bedtools/intersectBed -wa -v -a ${IN}_donor-.bed -b /NGS/known_sites/hg19/gencode.v19.exons.toUCSC.START.uniq.bed >> ${IN}_nonExonicDonor.bed
/NGS/links/bedtools/intersectBed -wa -v -a ${IN}_acceptor+.bed -b /NGS/known_sites/hg19/gencode.v19.exons.toUCSC.START.uniq.bed > ${IN}_nonExonicAcceptor.bed
/NGS/links/bedtools/intersectBed -wa -v -a ${IN}_acceptor-.bed -b /NGS/known_sites/hg19/gencode.v19.exons.toUCSC.END.uniq.bed >> ${IN}_nonExonicAcceptor.bed


# delete temporary files
rm -f ${IN}_*Count*
rm -f ${IN}_*count*
rm -f ${IN}_*+*
rm -f ${IN}_*-*
rm -f ${IN}_*exonicDonor.bed
rm -f ${IN}_*exonicAcceptor.bed
rm -f ${IN}_*Sorted.bed
rm -f ${IN}_exonicJunctions.bed
rm -f ${IN}_*genic.bed
rm -f ${IN}_exonicJunctionsIntragenic_nonambiguous.bed
rm -f ${IN}_circsAnnotated.txt
rm -f ${IN}_circsAnnotatedFinal.txt
