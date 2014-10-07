#!/bin/bash
#
#SBATCH -J test_def
#SBATCH -t 5-20:00:00
#SBATCH --mem=4000
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elham.baranizadeh@uef.fi
#
# Run a single task in the foreground
./run_days 
#
# Script ends here
