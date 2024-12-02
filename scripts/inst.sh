#! /bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:30:00
#SBATCH --partition=short
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mark.pritchard@ndm.ox.ac.uk

export JULIA_NUM_THREADS=1
module load Julia/1.9.3-linux-x86_64

julia scripts/instantiate.jl

wait

