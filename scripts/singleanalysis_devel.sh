#! /bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:10:00
#SBATCH --partition=devel
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mark.pritchard@ndm.ox.ac.uk

export JULIA_NUM_THREADS=1
module load Julia/1.11.3-linux-x86_64

julia scripts/instantiate.jl  

wait



