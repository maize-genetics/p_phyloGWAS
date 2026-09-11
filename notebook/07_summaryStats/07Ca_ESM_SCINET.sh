#on SCINET atlas
# ESM2 zero-shot scoring is GPU-bound and stays SCINET-only (see WORKFLOW.md); this just
# documents the collaborator project path as one variable instead of repeating it, and
# clarifies which zero-shot conversion script is authoritative (src/10_logit2zeroShot.R -
# NOT src/4_ESM_logits_to_zero_shot.py, which has a real bug: see
# archived/4_ESM_logits_to_zero_shot.py's header note).
SCINET_PROJECT_DIR=/project/90daydata/buckler_lab_panand/aimee.schulz/panand

find ${SCINET_PROJECT_DIR}/output/orthofinderProteinMSAs/ -type f |parallel --dryrun -j 40 "python src/4_ESM_logits.py -input /{} -output output/ESM_logits/{/.}_ESM_embedding.npz -model 'facebook/esm2_t33_650M_UR50D'" > cmd/ESM_logits.cmd

mkdir cmd/subcmd
split -d -a 4 -l 40 cmd/ESM_logits.cmd cmd/subcmd/esm

mkdir output/ESM_logits

sbatch -A buckler_lab_panand --array=202-499 slurm/ESM_run.sh -e log/


find ${SCINET_PROJECT_DIR}/output/orthofinderProteinMSAs -type f |parallel --dryrun -j 40 "Rscript src/10_logit2zeroShot.R output/ESM_logits/{/.}_ESM_embedding.npz {} output/ESM_zeroShot/{/.}_ESM_zeroShotScores.txt" > cmd/ESM_zeroShot.cmd

split -d -a 4 -l 40 cmd/ESM_zeroShot.cmd cmd/subcmd/zeroshot
mkdir output/ESM_zeroShot
sbatch -A buckler_lab_panand --array=0-499 slurm/zeroshot_run.sh -e log/


#additionalOG
find ${SCINET_PROJECT_DIR}/output/orthofinderProteinMSAs_additionalOGs/ -type f |parallel --dryrun -j 40 "python src/4_ESM_logits.py -input /{} -output output/ESM_logits_additionalOGs/{/.}_ESM_embedding.npz -model 'facebook/esm2_t33_650M_UR50D'" > cmd/ESM_logits_additionalOGs.cmd
split -d -a 4 -l 40 cmd/ESM_logits_additionalOGs.cmd cmd/subcmd/esm_additionalOGs

mkdir output/ESM_logits_additionalOGs
sbatch -A buckler_lab_panand --array=0-55 slurm/ESM_run.sh -e log/


find ${SCINET_PROJECT_DIR}/output/orthofinderProteinMSAs_additionalOGs -type f |parallel --dryrun -j 40 "Rscript src/10_logit2zeroShot.R output/ESM_logits_additionalOGs/{/.}_ESM_embedding.npz {} output/ESM_zeroShot_additionalOGs/{/.}_ESM_zeroShotScores.txt" > cmd/ESM_zeroShot_additionalOGs.cmd

split -d -a 4 -l 40 cmd/ESM_zeroShot_additionalOGs.cmd cmd/subcmd/zeroshot_additionalOGs
mkdir output/ESM_zeroShot_additionalOGs
sbatch -A buckler_lab_panand --array=0-55 slurm/zeroshot_run.sh -e log/
awk FNR!=1 output/ESM_zeroShot_additionalOGs/*.txt > output/combined_ESM_zeroShotScores_additionalOGs.txt