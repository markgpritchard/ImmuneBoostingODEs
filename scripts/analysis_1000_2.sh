#! /bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=24:00:00
#SBATCH --partition=medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mark.pritchard@ndm.ox.ac.uk

export JULIA_NUM_THREADS=20
module load Julia/1.11.3-linux-x86_64

n_rounds=1000
omega=2.0

julia scripts/rsvanalysis.jl "$omega" "$n_rounds" & 

wait



