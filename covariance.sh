#!/bin/bash
# Job name:
#SBATCH --job-name=freya
#
# Account:
#SBATCH --account=fc_deans
#
# Partition:
#SBATCH --partition=savio2
#
# Wall clock limit:
#SBATCH --time=72:00:00
#
## Command(s) to run:
module load cmake
module load python/2.7
module load python/2.7/numpy
module load python/2.7/scipy
module load python/2.7/cython
python cluster.py covar e x
python cluster.py covar e c
python cluster.py covar e T
python cluster.py covar e d
python cluster.py covar x c
python cluster.py covar x T
python cluster.py covar x d
python cluster.py covar c T
python cluster.py covar c d
python cluster.py covar T d
sacct -j $SLURM_JOB_ID --format=JobID%16,Elapsed
