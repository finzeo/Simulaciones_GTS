#! /bin/bash

#SBATCH --mail-type=ALL --mail-user=fainzeo@gmail.com
#SBATCH -J GTStest --no-requeue

## Time format is: DD-HH or DD-HH:MM:SS

## TRUE HPC RUN (D-HH)
#SBATCH -t 1-00 --nodes=4 --ntasks-per-node=20

## HPC RUN USING MIX NODES
#||SBATCH -t 2-00:00:00 --ntasks=10

## TEST RUN
#||SBATCH -t 0-00:05:00 --nodes=1 --partition=test --ntasks=4

echo -n "Starting: " ; date
cd $SLURM_SUBMIT_DIR 
mpiexec -genv OMP_NUM_THREADS 1 ./cs_solver
echo -n "Ending: " ; date
