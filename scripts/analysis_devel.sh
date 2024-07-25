#! /bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=00:10:00
#SBATCH --partition=devel
#SBATCH --clusters=all
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mark.pritchard@ndm.ox.ac.uk

export JULIA_NUM_THREADS=16
module load Julia/1.8.5-linux-x86_64

n_rounds=4

for n in {1..2}
do
	julia scripts/rsvanalysis.jl "$n" "$n_rounds" &
done

wait


