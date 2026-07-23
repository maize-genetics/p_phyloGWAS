
# updated 2024.05.21
# Charlie's 32 representative assemblies (34 in an earlier, superseded version)
rsync -av cbsublfs1:/data1/users/coh22/poaceae_tfbs/og_consensus ./

# annotation: scp -r cbsublfs1:/data4/users/zrm22/HelixerRuns/annotations/ data/
# genomes: /data1/users/mcs368/panand_genomes/genomes

mkdir output/aa_poaceaeRepAssemblies

find data/og_consensus/rep_annotations/ -name '*helixer.gff' -type f |cut -d "/" -f 4| sed 's/_helixer.gff//g' |parallel -j 25 'gffread -g data/og_consensus/rep_assemblies/{}.fa -y output/aa_poaceaeRepAssemblies/{}.aa.fa data/og_consensus/rep_annotations/{}_helixer.gff'

orthofinder -S diamond -I 1.5 -t 30 -a 12 -M msa -f output/aa_poaceaeRepAssemblies/ >output/orthoFinder_poaceaeRun_stdout.log 2> output/orthoFinder_poaceaeRun_stderr.log

orthofinder -t 20 -fg output/Results_Mar10 -M msa > output/orthoFinder_stdout.log 2> output/orthoFinder_stderr.log

# continues in 03B_OGFilter.ipynb (OG filtering)
