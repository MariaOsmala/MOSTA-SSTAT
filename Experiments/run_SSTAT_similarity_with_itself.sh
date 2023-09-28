#!/bin/bash -l
#SBATCH --job-name=SSTAT_itself
#SBATCH --output=outs/SSTAT_itself_%A_%a.txt
#SBATCH --error=errs/SSTAT_itself_%A_%a.txt
#SBATCH --account=project_2006203
#SBATCH --partition=small
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --array=0
#TIMEOUT
#32,39,41,43-45,48,52,64,82-83,93-96,100-103,107-108,126,130-131,155-156,158-159,196-198,201,219-220,222,228,248,249


#One can have only 300 scripts submitted
# run the analysis command

#There should be 8000 motif-pairs analysed in a single array job
srun run_SSTAT_with_itself.sh ${SLURM_ARRAY_TASK_ID}

seff $SLURM_JOBID


#25,26,30-32,37-39,41-45,48,50,51,64-66,75,77,79,82,93-97,99-112,130,135,151,153,155,156,158-162,174,175,180-182,191-193,195-198,201,210,215,219,220,222,227,228,231,234,235,245,248,249,273,286-288,293,296,305,308,315-316,335-337,339-340,359,361-363,365,370,376,390-391,402-405,411,412,423,439,441,442,444-455,460,462-467,475,487,506,583

#746,766,767,782,795,797,798,803,806,850,936,506
