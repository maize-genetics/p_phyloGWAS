# run p_phyloGWAS/notebook/05 seqFilesRearrange


#on SCINET atlas
find output/CDSPerOG_noATGFiltered -type f |parallel --dryrun -j 4 "mafft --ep 0 --genafpair --maxiterate 1000 {} > output/MSAPerOG_noATGFiltered/{/.}.fasta" > cmd/mafft.cmd

mkdir cmd/subcmd
split -d -a 4 -l 40 cmd/mafft.cmd cmd/subcmd/mafft

mkdir output/MSAPerOG_noATGFiltered

sbatch -A buckler_lab_panand --array=0-1146 slurm/mafft_run.sh -e log/

