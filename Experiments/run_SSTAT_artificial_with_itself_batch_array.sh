#!/bin/bash -l
#SBATCH --job-name=SSTAT_itself_artificial
#SBATCH --output=outs/SSTAT_itself_artificial_%A_%a.txt
#SBATCH --error=errs/SSTAT_itself_artificial_%A_%a.txt
#SBATCH --account=project_2006203
#SBATCH --partition=small
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --gres=nvme:20 
#SBATCH --array=0

#One can have only 300 scripts submitted
# run the analysis command

#There should be 8000 motif-pairs analysed in a single array job
srun run_SSTAT_artificial_with_itself.sh ${SLURM_ARRAY_TASK_ID}

seff $SLURM_JOBID

