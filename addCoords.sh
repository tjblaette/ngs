#!/bin/bash

#cat *circs.txt | awk -v OFS='\t' '$14 ~ /p/ {print $0,$11}' | paste - <( awk '$14 ~ /p/ {print $12}' *circs.txt | sed -e 's/[MDpN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}') <(awk '$14 ~ /p/ {print $13}' *circs.txt) <( awk '$14 ~ /p/ {print $14}' *circs.txt | sed -e 's/[0-9-]\+p.*$//' -e 's/[MDN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}') | grep -w '+' | head


FILE=$1

# start with + strand, case1 (inter-mate gap within segment2):
# Apppend: SS1 (Start segment1), LS1 (Length segment1), SS2a, LS2a, SS2b,LS2b, IMG (Inter-mate Gap)

P1="$(awk -v OFS='\t' '$14 ~ /p/ {print $0}' "$FILE")"
P1_SS1="$(awk -v OFS='\t' '$14 ~ /p/ {print $11}' "$FILE")"
P1_LS1="$( awk '$14 ~ /p/ {print $12}' "$FILE" | sed -e 's/[MDpN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}')"
P1_SS2a="$(awk '$14 ~ /p/ {print $13}' "$FILE")"
P1_LS2a="$(awk '$14 ~ /p/ {print $14}' "$FILE" | sed -e 's/[0-9-]\+p.*$//' -e 's/[MDN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}')"
P1_IMG="$(awk '$14 ~ /p/ {print $14}' "$FILE" | sed -e 's/.*[MIDNSH]\([0-9-]\+\)p.*$/\1/')"
P1_LS2b="$(awk '$14 ~ /p/ {print $14}' "$FILE" | sed -e 's/.*p\(.*\)$/\1/' -e 's/[MDN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}')"

P1_POSCOORDS="$(paste <(echo "$P1_SS1") <(echo "$P1_LS1") <(echo "$P1_SS2a") <(echo "$P1_LS2a") <(echo "$P1_IMG") <(echo "$P1_LS2b") | awk -v OFS='\t' '{print $1,$1+$2,$3,$3+$4,$3+$4+$5,$3+$4+$5+$6}')"

#echo "$POSCOORDS" | head
#head "coordsAdded.txt"



