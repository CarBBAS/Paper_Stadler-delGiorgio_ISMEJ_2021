#!/bin/bash
#SBATCH --account=mstadler
#SBATCH --time=0-24:00
#SBATCH --mem-per-cpu=100G
#SBATCH --cpus-per-task=15
#SBATCH --job-name="OTU_clustering_99"
#SBATCH --mail-user=m.stadler.jp.at@gmail.com
#SBATCH --mail-type=ALL
#creating a variable for the blast database

#to load the executables
module load StdEnv/2020 gcc/9.3.0 r/4.0.2 openmpi/4.0.3

#command to execute
Rscript /home/mstadler/Jobs/0_clusterOTU_slurm.R
#choose --save if you want #to #save your results on the R interface
