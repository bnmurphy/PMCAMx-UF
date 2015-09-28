#!/bin/bash
#
#SBATCH -J PMCAMxUF_OMP
#SBATCH -t 3-12:00:00
#SBATCH --mem=4000
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jan.julin@aces.su.se
#
# Run a single task in the foreground
./run_days 
#
# Script ends here
