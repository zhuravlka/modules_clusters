#!/bin/sh
#PBS -l walltime=600:00:00
#PBS -l nodes=node11:ppn=20
cd $PBS_O_WORKDIR
cd ~zhuravlka/modules_clusters/
if [ ! -d "results" ]; then
  mkdir results
fi

R CMD BATCH --no-save --no-restore "--args  modules_list.txt test_brain_10000_10000.RData brain best1 -u 10000 -d 10000 -g" count_Corr_27.08.R 
R CMD BATCH --no-save --no-restore "--args  modules_list.txt test_active_10000_10000.RData active best1 -u 10000 -d 10000 -g" count_Corr_27.08.R
R CMD BATCH --no-save --no-restore "--args  modules_list.txt test_all_10000_10000.RData all best1 -u 10000 -d 10000 -g" count_Corr_27.08.R