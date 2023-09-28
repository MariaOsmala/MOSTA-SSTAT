#!/bin/bash -l
#SBATCH --job-name=SSTAT_true_vs_artificial
#SBATCH --output=outs/SSTAT_true_vs_artificial_%A_%a.txt
#SBATCH --error=errs/SSTAT_true_vs_artificial_%A_%a.txt
#SBATCH --account=project_2007567
#SBATCH --partition=small
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --gres=nvme:1
#SBATCH --array=25,26 #0-103

#18687052

srun run_SSTAT_true_vs_corresponding_artificial.sh ${SLURM_ARRAY_TASK_ID}

seff $SLURM_JOBID



#18681373_[0]     small SSTAT_tr mariaosm PD       0:00      1 (Priority)
#18681305_0     small true_vs_ mariaosm  R       2:37      1 r18c05

# 
# 18679396_0     small SSTAT_it mariaosm  R    1:00:50      1 r18c48
#         
# 18674431_0     small SSTAT_re mariaosm  R    3:55:50      1 r18c48
# 18676855_0     small restTrue mariaosm  R    2:54:29      1 r18c35 FINISHED
#         
# 18675181_0     small SSTAT_ar mariaosm  R    3:28:42      1 r18c46
# 18656993_0     small artifici mariaosm  R   23:00:47      1 r18c38 finished?