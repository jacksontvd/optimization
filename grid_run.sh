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
python cluster.py grid e x
python cluster.py grid e c
python cluster.py grid e T
python cluster.py grid e d
python cluster.py grid x c
python cluster.py grid x T
python cluster.py grid x d
python cluster.py grid c T
python cluster.py grid c d
python cluster.py grid T d
sacct -j $SLURM_JOB_ID --format=JobID%16,Elapsed
