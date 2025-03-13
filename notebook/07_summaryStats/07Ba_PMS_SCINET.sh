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
