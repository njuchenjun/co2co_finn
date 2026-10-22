#!/bin/sh


##    file="/home/fengyq/vdw-co-co2/3-nnfit/fit-08/finn.dat" #########

    for i in {-12..25}
do
     cat total-65330.dat | awk '$14>'$i'*100 && $14<='$i'*100+100 {sum+=1} END {printf "%.1f-%.1f = %d\n", '$i'*100,'$i'*100+100, sum}'
done










