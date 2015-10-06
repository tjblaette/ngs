#!/bin/bash

#cat *circs.txt | awk -v OFS='\t' '$14 ~ /p/ {print $0,$11}' | paste - <( awk '$14 ~ /p/ {print $12}' *circs.txt | sed -e 's/[MDpN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}') <(awk '$14 ~ /p/ {print $13}' *circs.txt) <( awk '$14 ~ /p/ {print $14}' *circs.txt | sed -e 's/[0-9-]\+p.*$//' -e 's/[MDN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}') | grep -w '+' | head


FILE=$1

# start with + strand, case1 (inter-mate gap within segment2):
# Apppend: SS1 (Start segment1), LS1 (Length segment1), SS2a, LS2a, SS2b,LS2b, IMG (Inter-mate Gap)

P1="$(awk -v OFS='\t' '$3 == "+" && $14 ~ /p/ {print $0}' "$FILE")"
P1_SS1="$(awk -v OFS='\t' '$3 == "+" && $14 ~ /p/ {print $11}' "$FILE")"
P1_LS1="$( awk '$3 == "+" && $14 ~ /p/ {print $12}' "$FILE" | sed -e 's/[MDN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}')"
P1_SS2a="$(awk '$3 == "+" && $14 ~ /p/ {print $13}' "$FILE")"
P1_LS2a="$(awk '$3 == "+" && $14 ~ /p/ {print $14}' "$FILE" | sed -e 's/[0-9-]\+p.*$//' -e 's/[MDN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}')"
P1_IMG="$(awk '$3 == "+" && $14 ~ /p/ {print $14}' "$FILE" | sed -e 's/.*[MIDNSH]\([0-9-]\+\)p.*$/\1/')"
# P1_LS2b = Mate2_Length
P1_LS2b="$(awk '$3 == "+" && $14 ~ /p/ {print $14}' "$FILE" | sed -e 's/.*p\(.*\)$/\1/' -e 's/[MDN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}')"

# calculate donor start/stop, acceptor start/stop and second mate start/stop
P1_POSCOORDS="$(paste <(echo "$P1_SS1") <(echo "$P1_LS1") <(echo "$P1_SS2a") <(echo "$P1_LS2a") <(echo "$P1_IMG") <(echo "$P1_LS2b") | awk -v OFS='\t' '{print $1,$1+$2,$3,$3+$4,$3+$4+$5,$3+$4+$5+$6}')"

#echo "$P1_POSCOORDS" | head
#head "coordsAdded.txt"

#################################################################################
#################################################################################

# continue with + strand, case2 (intermate-gap within segment1):
P2="$(awk -v OFS='\t' '$3 == "+" && $12 ~ /p/ {print $0}' "$FILE")"
P2_SS1a="$(awk -v OFS='\t' '$3 == "+" && $12 ~ /p/ {print $11}' "$FILE")"
# P2_LS1a = Mate2_Length
P2_LS1a="$(awk '$3 == "+" && $12 ~ /p/ {print $12}' "$FILE" | sed -e 's/[0-9-]\+p.*$//' -e 's/[MDN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}')"
P2_IMG="$(awk '$3 == "+" && $12 ~ /p/ {print $12}' "$FILE" | sed -e 's/.*[MIDNSH]\([0-9-]\+\)p.*$/\1/')"
P2_LS1b="$(awk '$3 == "+" && $12 ~ /p/ {print $12}' "$FILE" | sed -e 's/.*p\(.*\)$/\1/' -e 's/[MDN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}')"
P2_SS2="$(awk '$3 == "+" && $12 ~ /p/ {print $13}' "$FILE")"
P2_LS2="$( awk '$3 == "+" && $12 ~ /p/ {print $14}' "$FILE" | sed -e 's/[MDN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}')"

# calculate donor start/stop, acceptor start/stop and second mate start/stop
P2_POSCOORDS="$(paste <(echo "$P2_SS1a") <(echo "$P2_LS1a") <(echo "$P2_IMG") <(echo "$P2_LS1b") <(echo "$P2_SS2") <(echo "$P2_LS2") | awk -v OFS='\t' '{print $1+$2+$3,$1+$2+$3+$4,$5,$5+$6,$1,$1+$2}')"

#echo "$P2" | head
#echo "$P2_POSCOORDS" | head

#################################################################################
#################################################################################

# start with - strand, case1 (inter-mate gap within segment2):
N1="$(awk -v OFS='\t' '$3 == "-" && $14 ~ /p/ {print $0}' "$FILE")"
N1_SS1="$(awk -v OFS='\t' '$3 == "-" && $14 ~ /p/ {print $11}' "$FILE")"
N1_LS1="$( awk '$3 == "-" && $14 ~ /p/ {print $12}' "$FILE" | sed -e 's/[MDN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}')"
N1_SS2a="$(awk '$3 == "-" && $14 ~ /p/ {print $13}' "$FILE")"
# N1_LS2a = Mate2_Length
N1_LS2a="$(awk '$3 == "-" && $14 ~ /p/ {print $14}' "$FILE" | sed -e 's/[0-9-]\+p.*$//' -e 's/[MDN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}')"
N1_IMG="$(awk '$3 == "-" && $14 ~ /p/ {print $14}' "$FILE" | sed -e 's/.*[MIDNSH]\([0-9-]\+\)p.*$/\1/')"
N1_LS2b="$(awk '$3 == "-" && $14 ~ /p/ {print $14}' "$FILE" | sed -e 's/.*p\(.*\)$/\1/' -e 's/[MDN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}')"

# calculate donor start/stop, acceptor start/stop and second mate start/stop
N1_POSCOORDS="$(paste <(echo "$N1_SS1") <(echo "$N1_LS1") <(echo "$N1_SS2a") <(echo "$N1_LS2a") <(echo "$N1_IMG") <(echo "$N1_LS2b") | awk -v OFS='\t' '{print $1,$1+$2,$3+$4+$5,$3+$4+$5+$6,$3,$3+$4}')"

#echo "$N1" | head
#echo "$N1_POSCOORDS" | head

#################################################################################
#################################################################################

# continue with - strand, case2 (intermate-gap within segment1):
N2="$(awk -v OFS='\t' '$3 == "-" && $12 ~ /p/ {print $0}' "$FILE")"
N2_SS1a="$(awk -v OFS='\t' '$3 == "-" && $12 ~ /p/ {print $11}' "$FILE")"
N2_LS1a="$(awk '$3 == "-" && $12 ~ /p/ {print $12}' "$FILE" | sed -e 's/[0-9-]\+p.*$//' -e 's/[MDN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}')"
N2_IMG="$(awk '$3 == "-" && $12 ~ /p/ {print $12}' "$FILE" | sed -e 's/.*[MIDNSH]\([0-9-]\+\)p.*$/\1/')"
# N2_LS1b = Mate2_Length
N2_LS1b="$(awk '$3 == "-" && $12 ~ /p/ {print $12}' "$FILE" | sed -e 's/.*p\(.*\)$/\1/' -e 's/[MDN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}')"
N2_SS2="$(awk '$3 == "-" && $12 ~ /p/ {print $13}' "$FILE")"
N2_LS2="$( awk '$3 == "-" && $12 ~ /p/ {print $14}' "$FILE" | sed -e 's/[MDN]/\t/g' -e 's/[0-9]\+[SHI]//g' | awk -v OFS='\t'  '{SUM=$1; for(i=2;i<=NF;i++)SUM=SUM+$i; print SUM}')"

# calculate donor start/stop, acceptor start/stop and second mate start/stop
N2_POSCOORDS="$(paste <(echo "$N2_SS1a") <(echo "$N2_LS1a") <(echo "$N2_IMG") <(echo "$N2_LS1b") <(echo "$N2_SS2") <(echo "$N2_LS2") | awk -v OFS='\t' '{print $1,$1+$2,$5,$5+$6,$1+$2+$3,$1+$2+$3+$4}')"

echo "$N2" | head
echo "$N2_POSCOORDS" | head





