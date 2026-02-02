#!/bin/sh
#SBATCH --chdir=./
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=15G
#SBATCH --partition=nocona
#SBATCH --time=48:00:00
#SBATCH --job-name=sim_ratio
#SBATCH --output=_ratio_sim_out
#SBATCH --error=_ratio_sim_error
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=egyllenh@ttu.edu
#SBATCH --array=1-675

# environment with SLiM5
source activate slim5-env

# picks a migrant/rep combo based on job array ID
param_set=$( head -n${SLURM_ARRAY_TASK_ID} $dir/params_v1_v2.txt | tail -n1 )
arr_param=(${param_set// / })

# run SLiM with relevant variables from array
slim -d coeff=${arr_param[0]} -d freq=${arr_param[1]} -d asym=${arr_param[2]} \
     -d mig=${arr_param[3]} -d rep=${arr_param[4]} \
     fst_div_ratio.slim
