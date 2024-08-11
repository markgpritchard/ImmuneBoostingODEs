#! /bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --time=00:10:00
#SBATCH --partition=devel
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mark.pritchard@ndm.ox.ac.uk

export JULIA_NUM_THREADS=24
module load Julia/1.8.5-linux-x86_64

omega=0.1
n_rounds=8

for n in {1..5}
do
	julia scripts/rsvanalysis.jl "$omega" "$n" "$n_rounds" &
done

wait



