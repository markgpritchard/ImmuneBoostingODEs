#! /bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
#SBATCH --time=24:00:00
#SBATCH --partition=medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mark.pritchard@ndm.ox.ac.uk

export JULIA_NUM_THREADS=36
module load Julia/1.9.3-linux-x86_64

omega=6.0
n_rounds=12

for n in {1..5}
do
	julia scripts/rsvanalysis.jl "$omega" "$n" "$n_rounds"
done

wait


