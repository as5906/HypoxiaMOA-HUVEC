#!/bin/bash

File=${1}
MOA1=${2}
MOA2=${3}

awk 'BEGIN{OFS="\t"} NR == FNR {a[$1] = 1; next} {if ($1 in a) $4 = "gain"; else $4 = "na"; print}' $MOA1 $File > out
awk 'BEGIN{OFS="\t"} NR == FNR {a[$1] = 1; next} {if ($1 in a) $5 = "loss"; else $5 = "na"; print}' $MOA2 out > out2

awk '{OFS="\t" ; print $1,$2,$3,$4","$5}' out2 > out4

sed -i 's/na,na/na/g' out4
sed -i 's/gain,na/gain/g' out4

sed -i 's/na,loss/loss/g' out4

rm out out2
mv out4 out
