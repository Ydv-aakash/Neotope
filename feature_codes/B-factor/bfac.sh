#!/bin/sh
#SBATCH -N 1
#SBATCH --ntasks-per-node=17
#SBATCH --time=1:00:00
#SBATCH --job-name=bfac_t
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out
#SBATCH --partition=standard
cd $SLURM_SUBMIT_DIR
###ENVIRONMENT
module load python/conda-python/3.7_new
###running
python /home/19bt60r19/model_running/bfactor_new_test.py
