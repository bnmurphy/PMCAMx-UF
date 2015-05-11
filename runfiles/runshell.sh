#!/bin/bash
#
#SBATCH -J test_def
#SBATCH -t 2-20:00:00
#SBATCH --mem=4000
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=davidpat@chemeng.upatras.gr
#
# Run a single task in the foreground
./run_days 
#
# Script ends here
