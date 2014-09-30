#!/bin/sh
#PBS -l walltime=600:00:00
#PBS -l nodes=node18:ppn=3
#cd $PBS_O_WORKDIR
cd ./data
for dir in *;do
	for i in $dir/C*.tss.nr.bed;do 
		k=$(sed -e 's/.*\/C\(.*\)\.tss\.nr\.bed/\1/' <<< $i)
		awk '{print $1,$2,$3}' $i > ./$dir/$k.bed
		echo $k ;done|sort -n >./$dir/modules_list.txt
done
