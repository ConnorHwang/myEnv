#!/bin/bash

#SBATCH -J jet2
#SBATCH -o jet2.out
#SBATCH -N 12
#SBATCH -n 384
#SBATCH -t 900:00:00           
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hanul@stanford.edu

WORKDIR='/fastscratch/hanul/jfm2020_response/n6epsilon1/'
export WORKDIR

### ---------------------------------------
### BEGINNING OF EXECUTION
### ---------------------------------------

echo The master node of this job is `hostname`
echo This job runs on the following nodes:
echo `scontrol show hostname $SLURM_JOB_NODELIST`
echo "Starting at `date`"
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NPROCS processors."
echo "Current working directory is `echo $WORKDIR`"
echo
echo Output from code
echo ----------------
echo "module purge & load"
module purge
module load intel/18.0.2.199
#module load apps/matlab/R2020a
module load openmpi/1.10.7
echo "Run my commends"
# Run your commends
qcc -source -D_MPI=1 uniform.c
mpicc -Wall -std=c99 -O2 _uniform.c -o uniform -lm
#./compile
mpiexec -n 384 uniform
