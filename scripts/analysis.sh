#! /bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
#SBATCH --time=24:00:00
#SBATCH --partition=medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mark.pritchard@ndm.ox.ac.uk

export JULIA_NUM_THREADS=8
module load Julia/1.11.3-linux-x86_64

julia scripts/rsvanalysis.jl "2.0" & 
julia scripts/rsvanalysis.jl "1.0" & 
julia scripts/rsvanalysis.jl "0.5" &
julia scripts/rsvanalysis.jl "0.2" & 

wait



