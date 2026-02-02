#!/bin/sh
#SBATCH --chdir=./
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --partition=condo
#SBATCH --time=48:00:00
#SBATCH --job-name=eems_noZ_bevi
#SBATCH --output=_eems_noZ_out
#SBATCH --error=_eems_noZ_error
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=egyllenhaal@unm.edu
#SBATCH --array=1-20

# submitted from /lustre/scratch/egyllenh/greenbul/eems/
dir=$SLURM_SUBMIT_DIR
cd $dir

module load miniconda3
source activate eems-env

# set up a different seen, one off, for each chain
seed=$((708-${SLURM_ARRAY_TASK_ID}))

# run eems for each of four chains in a job array
~/eems/runeems_snps/src/runeems_snps --params params/params_noZ_chain${SLURM_ARRAY_TASK_ID} --seed ${seed}
