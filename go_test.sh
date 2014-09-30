#!/bin/sh
#PBS -l walltime=600:00:00
#PBS -l nodes=node15:ppn=20
cd $PBS_O_WORKDIR
cd ~zhuravlka/modules_clusters/
if [ ! -d "results" ]; then
  mkdir results
fi

R CMD BATCH --no-save --no-restore "--args  modules_list.txt test_brain.RData brain best1 -g" count_Corr_27.08.R testbrain.RLog
R CMD BATCH --no-save --no-restore "--args  modules_list.txt test_active.RData active best1 -g" count_Corr_27.08.R testactive.RLog
R CMD BATCH --no-save --no-restore "--args  modules_list.txt test_all.RData all best1 -g" count_Corr_27.08.R testall.RLog