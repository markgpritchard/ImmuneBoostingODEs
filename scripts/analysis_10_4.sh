#! /bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --partition=medium
#SBATCH --clusters=all
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mark.pritchard@ndm.ox.ac.uk

export JULIA_NUM_THREADS=20
module load Julia/1.8.5-linux-x86_64

n_rounds=10
chainid=4

julia scripts/rsvanalysis.jl "$chainid" "$n_rounds"

wait


