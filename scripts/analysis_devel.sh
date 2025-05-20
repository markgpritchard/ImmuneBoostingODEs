#! /bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=00:10:00
#SBATCH --partition=devel
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mark.pritchard@ndm.ox.ac.uk

export JULIA_NUM_THREADS=8
module load Julia/1.11.3-linux-x86_64

julia scripts/rsvanalysis.jl "25" "2500" & 
julia scripts/rsvanalysis.jl "1000" "1000000" & 
julia scripts/rsvanalysis.jl "2500" "1000000" & 

wait



