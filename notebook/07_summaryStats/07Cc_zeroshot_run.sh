#!/bin/bash

#SBATCH --time=00:20:00  # walltime limit (HH:MM:SS)
#SBATCH --account=buckler_lab_panand
#SBATCH --nodes=1  # number of nodes
#SBATCH --exclusive
#SBATCH --partition=atlas  # standard node(s)
#SBATCH --job-name="zeroshotscore"
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

source ~/miniconda3/bin/activate cbsuenv
ml parallel

id=`printf "%04d\n" ${SLURM_ARRAY_TASK_ID}`
cat cmd/subcmd/zeroshot${id} |parallel -j 40 "{}"
