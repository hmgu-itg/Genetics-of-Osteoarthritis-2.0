#!/usr/bin/env bash

sumstats=$1
output=$2

get_colnum () {
    column=$1
    sumstats=$2
    zcat $sumstats | head -n 1 | tr '\t' '\n' | awk '{print NR" "$1}' | grep -w $column | awk -F' ' '{print $1}'
}

CPTID=$(get_colnum CPTID $sumstats)
EA=$(get_colnum EA $sumstats)
NEA=$(get_colnum NEA $sumstats)
EAF=$(get_colnum EAF $sumstats)
BETA=$(get_colnum BETA $sumstats)
SE=$(get_colnum SE $sumstats)
P=$(get_colnum P $sumstats)
N=$(get_colnum N $sumstats)

echo SNP A1 A2 freq b se p N > $output
awk -F$'\t' \
  -v CPTID="$CPTID" \
  -v EA="$EA" \
  -v NEA="$NEA" \
  -v EAF="$EAF" \
  -v BETA="$BETA" \
  -v SE="$SE" \
  -v P="$P" \
  -v N="$N" \
  '{if (NR!=1){print $CPTID,toupper($EA),toupper($NEA),$EAF,$BETA,$SE,$P,$N}}' <(zcat $sumstats) \
  >> $output