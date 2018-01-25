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
python cluster.py var e
python cluster.py var x
python cluster.py var c
python cluster.py var T
python cluster.py var d
