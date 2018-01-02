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
#SBATCH --time=10:00:00
#
## Command(s) to run:
module load cmake
module load python/2.7.8
module load numpy
module load scipy
module load cython
python cluster.py
