#!/bin/bash

#SBATCH -J e1
#SBATCH -o e1.out
#SBATCH -p pdebug
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH -t 01:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hanul@stanford.edu

WORKDIR='/p/lustre2/hwang8/helmholtz/setupcase'
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

# module
module load ensight

# Run your commends
ensight -batch -p ./ensight_run1.py
#srun -n 32 /g/g92/hwang8/cascade/helmholtz_vof/nextgen/bin/charles_vof.exe -i ./charles_vof.in
