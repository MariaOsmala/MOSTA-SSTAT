#!/bin/bash -l
#SBATCH --job-name=SSTAT_restTrue_artificial
#SBATCH --output=outs/SSTAT_restTrue_artificial_%A_%a.txt
#SBATCH --error=errs/SSTAT_restTrue_artificial_%A_%a.txt
#SBATCH --account=project_2007567
#SBATCH --partition=small
#SBATCH --time=3-00:00:00 #2days
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --gres=nvme:10
#SBATCH --array=0 #-1000


#time runs out

# sacct --format JobID%-20,State -j 18656531 | grep CANCELLED

srun run_SSTAT_restTrue_artificial.sh ${SLURM_ARRAY_TASK_ID}

seff $SLURM_JOBID



