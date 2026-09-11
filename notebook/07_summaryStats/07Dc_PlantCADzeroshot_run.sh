#!/bin/bash

#SBATCH --time=00:20:00  # walltime limit (HH:MM:SS)
#SBATCH --account=buckler_lab_panand
#SBATCH --nodes=1  # number of nodes
#SBATCH --exclusive
#SBATCH --partition=atlas  # standard node(s)
#SBATCH --job-name="PlantCAD zeroshotscore"
#SBATCH --mail-user=sh2246@cornell.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

echo "SLURM job id: ${SLURM_JOB_ID}"

source ~/miniconda3/bin/activate cbsuenv

# single job over the whole combined FASTA/npz - mirrors 07Db's logit job, no array splitting
COMBINED_FASTA=$1
Rscript src/6_PlantCAD_logit2zeroShot.R output/PlantCAD_logits/PlantCAD_embedding.npz "${COMBINED_FASTA}" output/PlantCAD_zeroShot/PlantCAD_zeroShotScores.txt
