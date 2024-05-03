#!/bin/bash -l
#SBATCH --job-name=SSTAT_itself
#SBATCH --output=outs/SSTAT_itself.txt
#SBATCH --error=errs/SSTAT_itself.txt
#SBATCH --account=project_2006203
#SBATCH --partition=small
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G

#21091913


#srun run_SSTAT_with_itself.sh ${SLURM_ARRAY_TASK_ID}
#srun run_SSTAT_with_itself_fromYimeng_version2.sh 
#srun run_SSTAT_with_itself_union.sh #DONE
srun run_SSTAT_with_itself_version2.2.sh #DONE

seff $SLURM_JOBID


