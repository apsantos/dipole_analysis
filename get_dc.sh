#!/bin/bash

module load anaconda
# INPUT FOR JOB
header=$1
end='BHM'

nB=$2
NAVG=$3
half_nB=`python -c "print int(0.5*${nB})"`
onefive_nB=`python -c "print int(1.5*${nB})"`
ntot=`python -c "print int(2.5*${nB})"`
echo 0 0 0 0 0 0 0 > ${header}allH.BHM
for f in ${header}AH*.${end} ${header}BH*.${end} ${header}CH*.${end}; do
   t=`tail -n 1 $f | awk '{print $1}'`
   H=`tail -n 1 $f | awk '{print $2}'`
   B=`tail -n 1 $f | awk '{print $3}'`
   tail -n $NAVG $f > tmp.BHM
   M=`dataFunctions.py -i tmp.BHM -a -d -l 0 -c 3`
   echo $t $H $B $M >> ${header}allH.BHM
done
head -n1 $f > ${header}loop.BHM
sort_dc.py -f ${header}allH.BHM -l 1 -o tmp.BHM
cat tmp.BHM >> ${header}loop.BHM

#for (( c=1; c<=$ntot; c++ ))
#do
#  if [ "$c" -lt "${half_nB}" ]; then
#    cycle="A"
#    Bmag=`python -c "print ${c}*${Bscale}/${nB}"`
#  elif [ "$c" -lt "${onefive_nB}" ]; then
#    cycle="B"
#    Bmag=`python -c "print (${nB}-${c})*${Bscale}/${nB}"`
#  else
#    cycle="C"
#    Bmag=`python -c "print (${c}-2*${nB})*${Bscale}/${nB}"`
#   fi
#   if [ "$Bmag" = 0.0 ]; then
#     Bmag='0'
#   fi
#   tail -n 1 ${header}${cycle}'B'${Bmag}'.'${end}
#done
