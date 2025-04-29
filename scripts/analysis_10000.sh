#! /bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --time=24:00:00
#SBATCH --partition=medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mark.pritchard@ndm.ox.ac.uk

export JULIA_NUM_THREADS=4
module load Julia/1.11.3-linux-x86_64

n_rounds=10000

julia scripts/rsvanalysis.jl "0.1" "$n_rounds" & 
julia scripts/rsvanalysis.jl "0.2" "$n_rounds" & 
julia scripts/rsvanalysis.jl "0.4" "$n_rounds" & 
julia scripts/rsvanalysis.jl "1.0" "$n_rounds" & 
julia scripts/rsvanalysis.jl "2.0" "$n_rounds" & 
julia scripts/rsvanalysis.jl "4.0" "$n_rounds" & 
julia scripts/rsvanalysis.jl "6.0" "$n_rounds" &  

wait



