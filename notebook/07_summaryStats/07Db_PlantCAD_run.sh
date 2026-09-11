#!/bin/bash

#SBATCH --time=24:00:00  # walltime limit (HH:MM:SS)
#SBATCH --account=buckler_lab_panand
#SBATCH --nodes=1  # number of nodes
#SBATCH --partition=gpu-a100  # gpu node
#SBATCH --gres=gpu:a100:1 # request 1 GPU
#SBATCH --job-name="PlantCAD logit"
#SBATCH --mail-user=sh2246@cornell.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

echo "SLURM job id: ${SLURM_JOB_ID}"

source ~/miniconda3/bin/activate /project/buckler_lab_panand/ana.berthel/conda/pytorch

# single job over the whole combined FASTA - no per-OG array splitting (unlike ESM2's
# per-OG protein MSAs), since PlantCAD's input is one combined nucleotide FASTA and
# src/5_PlantCAD_logits.py already loops over every record in it.
COMBINED_FASTA=$1
python src/5_PlantCAD_logits.py -input "${COMBINED_FASTA}" -output output/PlantCAD_logits/PlantCAD_embedding.npz -model 'kuleshov-group/PlantCaduceus_l32'
