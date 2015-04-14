#!/bin/bash
#
#SBATCH -J PMCAMxUF_OMP
#SBATCH -t 6-3:00:00
#SBATCH --mem=4000
#SBATCH -N 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jan.julin@aces.su.se
#SBATCH --constraint=vtune
#
# Run a single task in the foreground
./run_days 
#
# Script ends here
