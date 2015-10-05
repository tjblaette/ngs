#!/bin/bash

#cat *circs.txt | awk -v OFS='\t' '$14 ~ /p/ {print $0,$11}' | paste - <( awk '$14 ~ /p/ {print $12}' *circs.txt | sed -e 's/[MDpN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}') <(awk '$14 ~ /p/ {print $13}' *circs.txt) <( awk '$14 ~ /p/ {print $14}' *circs.txt | sed -e 's/[0-9-]\+p.*$//' -e 's/[MDN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}') | grep -w '+' | head


FILE=$1

# start with + strand, case1 (inter-mate gap within segment2):
# Apppend: SS1 (Start segment1), LS1 (Length segment1), SS2a, LS2a, SS2b,LS2b, IMG (Inter-mate Gap)

POS1_and_SS1="$(awk -v OFS='\t' '$14 ~ /p/ {print $0,$11}' "$FILE")"
LS1="$( awk '$14 ~ /p/ {print $12}' "$FILE" | sed -e 's/[MDpN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}')"
SS2a="$(awk '$14 ~ /p/ {print $13}' "$FILE")"
LS2a="$(awk '$14 ~ /p/ {print $14}' "$FILE" | sed -e 's/[0-9-]\+p.*$//' -e 's/[MDN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}')"
IMG="$(awk '$14 ~ /p/ {print $14}' "$FILE" | sed -e 's/.*[MIDNSH]\([0-9-]\+\)p.*$/\1/')"
LS2b="$(awk '$14 ~ /p/ {print $14}' "$FILE" | sed -e 's/.*p\(.*\)$/\1/' -e 's/[MDN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}')"


#echo "$POS1_and_SS1" | head
#echo "$LS1" | head 

paste <(echo "$POS1_and_SS1") <(echo "$LS1") <(echo "$SS2a") <(echo "$LS2a") <(echo "$IMG") <(echo "$LS2b")> coordsAdded.txt

head "coordsAdded.txt"
