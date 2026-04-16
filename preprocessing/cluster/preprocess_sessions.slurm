#!/bin/bash
#SBATCH --job-name=dyadic_preproc
#SBATCH --account=???
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:20:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=???
#SBATCH --output=/logs/preproc_%A_%a.out
#SBATCH --error=/logs/preproc_%A_%a.err
#SBATCH --array=1-247

module purge
module load matlab/R2024b

MATLAB_SCRIPT_DIR=/projects/m/mkislin/DyadicInteractions/preprocessing
INPUT_DIR=/projects/m/mkislin/DyadicInteractions/files2run
OUTPUT_DIR=/projects/m/mkislin/DyadicInteractions/preproc_outputs

#mkdir -p "$OUTPUT_DIR"

matlab -singleCompThread -nodisplay -nosplash -r "addpath(genpath('$MATLAB_SCRIPT_DIR')); run_preprocessing_array('$INPUT_DIR','$OUTPUT_DIR',$SLURM_ARRAY_TASK_ID); exit"

