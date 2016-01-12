#!/bin/bash

# Submit jobs. Assume 16 cores per proc. Make sure line numbers are the same.
for (( i = 1; i<= 16 ; i++ ))
do
    sed -i "/--nodes=/c #SBATCH --nodes=$i                      # number of nodes requested (n)" slurm_omp.sh
    sed -i "/--ntasks=/c #SBATCH --ntasks=$i                    # required number of CPUs (n)" slurm_omp.sh
    sbatch slurm_omp.sh
done

