# Archived: original SCINET/SLURM-array version of premature-stop scoring, run against a
# collaborator's own SCINET project directory across two independent OG batches (main +
# additionalOGs). Superseded by ../07Ba_PMS_run.sh, a local GNU-parallel runner, now that
# output/orthofinderProteinMSAs_fullset_20250710/ has the full, current OG set on disk locally.
# The SLURM array template this script's sbatch calls referenced (slurm/PMS_run.sh) is a
# byte-identical duplicate of notebook/slurm/PMS_run.sh and was removed rather than archived.
find /project/90daydata/buckler_lab_panand/aimee.schulz/panand/output/orthofinderProteinMSAs/ -type f |parallel --dryrun -j 40 "python src/11_find_premature_stops.py {} output/prematureStop/{/.}_prematureStops.txt" > cmd/PMS.cmd

mkdir cmd/subcmd
split -d -a 4 -l 40 cmd/PMS.cmd cmd/subcmd/pms

sbatch -A buckler_lab_panand --array=0-499 slurm/PMS_run.sh -e log/
awk FNR!=1 output/prematureStop/*.txt > output/combined_PMS.txt


find /90daydata/buckler_lab_panand/aimee.schulz/panand/output/orthofinderProteinMSAs_additionalOGs -type f |parallel --dryrun -j 40 "python src/11_find_premature_stops.py {} output/prematureStop_additionalOGs/{/.}_prematureStops.txt" > cmd/PMS_additionalOGs.cmd
split -d -a 4 -l 40 cmd/PMS_additionalOGs.cmd cmd/subcmd/pms_additionalOGs

mkdir output/prematureStop_additionalOGs
sbatch -A buckler_lab_panand --array=0-55 slurm/PMS_run.sh -e log/
awk FNR!=1 output/prematureStop_additionalOGs/*.txt > output/combined_PMS_additionalOGs.txt
