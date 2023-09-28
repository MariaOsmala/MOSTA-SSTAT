#!/bin/bash -l
#SBATCH --job-name=SSTAT_artificial
#SBATCH --output=outs/SSTAT_artificial_%A_%a.txt
#SBATCH --error=errs/SSTAT_artificial_%A_%a.txt
#SBATCH --account=project_2007567
#SBATCH --partition=small
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --gres=nvme:20 
#SBATCH --array=0 #-1000

#sacct -j 17389178 -o state%20,jobid%20 | grep FAILED 418
#sacct -j 17389381 -o state%20,jobid%20 | grep FAILED 770
#sacct -j 17389583 -o state%20,jobid%20 | grep FAILED 

#One can have only 300 scripts submitted
# run the analysis command

#There should be 8000 motif-pairs analysed in a single array job
srun run_SSTAT_artificial.sh ${SLURM_ARRAY_TASK_ID}

seff $SLURM_JOBID


#
