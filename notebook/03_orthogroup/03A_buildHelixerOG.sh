
# updated 2024.05.21
# Charlie's 32 representative assemblies (34 in an earlier, superseded version)
rsync -av cbsublfs1:/data1/users/coh22/poaceae_tfbs/og_consensus ./

mkdir output/aa_poaceaeRepAssemblies

find data/og_consensus/rep_annotations/ -name '*helixer.gff' -type f |cut -d "/" -f 4| sed 's/_helixer.gff//g' |parallel -j 25 'gffread -g data/og_consensus/rep_assemblies/{}.fa -y output/aa_poaceaeRepAssemblies/{}.aa.fa data/og_consensus/rep_annotations/{}_helixer.gff'

orthofinder -S diamond -I 1.5 -t 30 -a 12 -M msa -f output/aa_poaceaeRepAssemblies/ >output/orthoFinder_poaceaeRun_stdout.log 2> output/orthoFinder_poaceaeRun_stderr.log

# continues in 03B_OGFilter.ipynb (OG filtering)
