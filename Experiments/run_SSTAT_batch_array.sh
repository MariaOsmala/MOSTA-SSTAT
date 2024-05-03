#!/bin/bash -l
#SBATCH --job-name=SSTAT
#SBATCH --output=outs/SSTAT_%A_%a.txt
#SBATCH --error=errs/SSTAT_%A_%a.txt
#SBATCH --account=project_2006472
#SBATCH --partition=small
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --gres=nvme:20 
#SBATCH --array=601-1000  #0-1000


#One can have only 300 scripts submitted
#21118511_[301-600]  #OK  
#21118410_[0-300] #OK
#21118692_601-1000
#sacct --format JobID%-20,State -j 20450264 | grep CANCELLED

#srun run_SSTAT.sh ${SLURM_ARRAY_TASK_ID} #fromYimeng_version1
#srun run_SSTAT_union.sh ${SLURM_ARRAY_TASK_ID} #union of version1 and version2
srun run_SSTAT_version2.2.sh ${SLURM_ARRAY_TASK_ID} #version2.2

seff $SLURM_JOBID


