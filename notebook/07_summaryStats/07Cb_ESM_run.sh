#!/bin/bash

#SBATCH --time=02:00:00  # walltime limit (HH:MM:SS)
#SBATCH --account=buckler_lab_panand
#SBATCH --nodes=1  # number of nodes
#SBATCH --partition=gpu-a100  # gpu node
#SBATCH --gres=gpu:a100:1 # request 1 GPU
#SBATCH --job-name="ESM logit"
#SBATCH --mail-user=sh2246@cornell.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
 
# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

echo "All jobs in this array have:"
echo "- SLURM array job id: ${SLURM_ARRAY_JOB_ID}"
echo "- SLURM array task count: ${SLURM_ARRAY_TASK_COUNT}"
echo "- SLURM array starting task: ${SLURM_ARRAY_TASK_MIN}"
echo "- SLURM array ending task: ${SLURM_ARRAY_TASK_MAX}"
echo "This job in the array has:"
echo "- SLURM job id: ${SLURM_JOB_ID}"
echo "- SLURM array task id: ${SLURM_ARRAY_TASK_ID}"


source ~/miniconda3/bin/activate /project/buckler_lab_panand/ana.berthel/conda/pytorch
ml parallel

id=`printf "%04d\n" ${SLURM_ARRAY_TASK_ID}`
cat cmd/subcmd/esm${id} |parallel -j 1 "{}"


# id=`ls /90daydata/buckler_lab_panand/aimee.schulz/panand/output/orthofinderProteinMSAs/ |head -n ${SLURM_ARRAY_TASK_ID}|tail -n1|sed 's/.fa//g'`

# python src/4_ESM_logits.py -input /project/90daydata/buckler_lab_panand/aimee.schulz/panand/output/orthofinderProteinMSAs/${id}.fa -output output/ESM_logits/${id}_ESM_embedding.npz 

# python src/4_ESM_logits.py -input /project/90daydata/buckler_lab_panand/aimee.schulz/panand/output/orthofinderProteinMSAs/OG0000508.fa -output output/ESM_logits/OG0000508_ESM_logits.npz 
