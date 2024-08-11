#! /bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --time=12:00:00
#SBATCH --partition=short
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mark.pritchard@ndm.ox.ac.uk

export JULIA_NUM_THREADS=24
module load Julia/1.8.5-linux-x86_64

omega=10.0
n_rounds=12

for n in {1..5}
do
	julia scripts/rsvanalysis.jl "$omega" "$n" "$n_rounds" &
done

wait



