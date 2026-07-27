#on SCINET atlas
# PlantCAD zero-shot scoring, mirroring the ESM2 pipeline (07Ca/07Cb/07Cc) with one structural
# difference: PlantCAD's input is a single combined nucleotide FASTA (patch-scaffolds
# build-msa --single-fasta, no --protein-mode - see notebook/04_msaGeneration/README.md), not
# one file per OG like ESM2's protein MSAs, so this runs as 2 single SLURM jobs (logit
# extraction, then zero-shot conversion) rather than array jobs over per-OG cmd files.
SCINET_PROJECT_DIR=/project/90daydata/buckler_lab_panand/aimee.schulz/panand
# TODO: fill in the real path once patch-scaffolds build-msa --single-fasta has been run for
# the CDS/nucleotide case (see notebook/04_msaGeneration/README.md's step 0)
COMBINED_FASTA=${SCINET_PROJECT_DIR}/output/orthofinderCDSMSAs_singleFasta/allOGs_combined.fa

mkdir -p output/PlantCAD_logits output/PlantCAD_zeroShot

logits_job=$(sbatch --parsable -A buckler_lab_panand notebook/07_summaryStats/07Db_PlantCAD_run.sh "${COMBINED_FASTA}")
sbatch -A buckler_lab_panand --dependency=afterok:${logits_job} notebook/07_summaryStats/07Dc_PlantCADzeroshot_run.sh "${COMBINED_FASTA}"
