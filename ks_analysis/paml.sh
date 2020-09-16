#!/bin/bash
#
#SBATCH -n 1                 # Number of cores
#SBATCH -N 1                 # Number of nodes for the cores
#SBATCH -t 0-20:05           # Runtime in D-HH:MM format
#SBATCH -p serial_requeue    # Partition to submit to
#SBATCH --mem=200            # Memory pool for all CPUs
#SBATCH -o lmcai.out      # File to which standard out will be written
#SBATCH -e lmcai.err      # File to which standard err will be written

source new-modules.sh
module load paml/4.8-fasrc01
codeml paml.clt
