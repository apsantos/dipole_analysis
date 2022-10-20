#!/bin/bash

module load anaconda
# INPUT FOR JOB
header=$1
end='BHM'

nB=$2
NAVG=$3
nsims=$4
half_nB=`python -c "print int(0.5*${nB})"`
onefive_nB=`python -c "print int(1.5*${nB})"`
ntot=`python -c "print int(2.5*${nB})"`
echo 0 0 0 0 0 0 0 > ${header}allH.BHM
for f in ${header}AH*_i1.${end} ${header}BH*_i1.${end} ${header}CH*_i1.${end}; do
   t=`tail -n 1 $f | awk '{print $1}'`
   H=`tail -n 1 $f | awk '{print $2}'`
   B=`tail -n 1 $f | awk '{print $3}'`
   rm ${header}nsims_allH.BHM
   dM=0.0
   for ((isim=1; isim<=${nsims}; isim++)); do
     tail -n $NAVG $f > tmp.BHM
     ave=`dataFunctions.py -i tmp.BHM -a -d -l 0 -c 3`
     M=`echo $ave | awk '{print $1}'`
     ddM=`echo $ave | awk '{print $2}'`
     dM=`python -c "print(float(${dM})+float(${ddM})**2.0)"`
     echo $t $H $B $M >> ${header}nsims_allH.BHM
   done
   M=`dataFunctions.py -i ${header}nsims_allH.BHM -a -l 0 -c 3`
   dM=`python -c "print(float(${dM})**0.5/float(${nsims}))"`
   echo $t $H $B $M $dM >> ${header}allH.BHM
done
head -n1 $f > ${header}loop.BHM
sed -i '/starting/d' ${header}allH.BHM
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

