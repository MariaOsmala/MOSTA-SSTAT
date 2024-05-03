#!/bin/bash -l
#SBATCH --job-name=SSTAT_true_vs_artificial
#SBATCH --output=outs/SSTAT_true_vs_artificial_%A_%a.txt
#SBATCH --error=errs/SSTAT_true_vs_artificial_%A_%a.txt
#SBATCH --account=project_2006203
#SBATCH --partition=small
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --gres=nvme:1
#SBATCH --array=201-393 # 0-393

srun run_SSTAT_true_vs_corresponding_artificial.sh ${SLURM_ARRAY_TASK_ID}

seff $SLURM_JOBID

#21354174_[0]  
#21354374_[1-200]
#21354380_[201-393]