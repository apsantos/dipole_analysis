#!/bin/bash

# INPUT FOR JOB
header='dnp155rl24_'
end='BHM'

Bscale=0.001
nB=10
half_nB=`python -c "print int(0.5*${nB})"`
onefive_nB=`python -c "print int(1.5*${nB})"`
ntot=`python -c "print int(2.5*${nB})"`
#echo "# TimeStep v_B v_H v_M v_Chi" > ${header}allB.BHM
echo 0 0 0 0 0 >> ${header}allB.BHM
files=(${header}AB*.${end})
for f in ${header}AB*.${end} ${header}BB*.${end} ${header}CB*.${end}; do
   tail -n 1 $f >> ${header}allB.BHM
done
python sort.py -f ${header}allB.BHM -l 1 -o ${header}loop.BHM

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
