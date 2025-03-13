
# updated 2024.05.21
# Charlie's 34 representive assemblies
rsync -av cbsublfs1:/data1/users/coh22/poaceae_tfbs/og_consensus ./

# annotation: scp -r cbsublfs1:/data4/users/zrm22/HelixerRuns/annotations/ data/
# genomes: /data1/users/mcs368/panand_genomes/genomes

mkdir output/aa_poaceaeRepAssemblies

find data/og_consensus/rep_annotations/ -name '*helixer.gff' -type f |cut -d "/" -f 4| sed 's/_helixer.gff//g' |parallel -j 25 'gffread -g data/og_consensus/rep_assemblies/{}.fa -y output/aa_poaceaeRepAssemblies/{}.aa.fa data/og_consensus/rep_annotations/{}_helixer.gff' 

orthofinder -S diamond -I 1.5 -t 30 -a 12 -M msa -f output/aa_poaceaeRepAssemblies/ >output/orthoFinder_poaceaeRun_stdout.log 2> output/orthoFinder_poaceaeRun_stderr.log

orthofinder -t 20 -fg output/Results_Mar10 -M msa > output/orthoFinder_stdout.log 2> output/orthoFinder_stderr.log


# run 09_OG_comparison.ipynb to do OG filtering
# reconstruct ancestral sequence
mkdir output/poaceaeHelixOGMSA_plusAnc_toAdd
cat output/poaceaeHelixerOG_filtered_v2_toRecover.txt | parallel -j 35 'Rscript src/S02_phangornAncestralSeqReconstruct_AA.R --input output/OrthoFinder/Results_Jun06/MultipleSequenceAlignments/{}.fa --output output/poaceaeHelixOGMSA_plusAnc_toAdd/{}'

# extract the ancestral sequence
mkdir output/poaceaeHelixOG_ancSeq_toAdd/
find output/poaceaeHelixOGMSA_plusAnc_toAdd/ -type f| cut -d / -f 3| sort |parallel -j 40 'tail -n 2 output/poaceaeHelixOGMSA_plusAnc_toAdd/{} > output/poaceaeHelixOG_ancSeq_toAdd/{}'
# rename
mkdir output/poaceaeHelixOG_ancSeq_renamed_toAdd
find output/poaceaeHelixOG_ancSeq_toAdd/ -type f| cut -d / -f 3|sort|parallel -j 40 "cat output/poaceaeHelixOG_ancSeq_toAdd/{} | sed 's/MRCA/{= s/.fa// =}/g' > output/poaceaeHelixOG_ancSeq_renamed_toAdd/{}"
# add back
# cat output/poaceaeHelixOG_ancSeq_renamed_toAdd/* > output/poaceaeHelixerOG_ancSeq_filtered_v2_toAdd.fa
find output/poaceaeHelixOG_ancSeq_renamed_toAdd/ -type f| sort | xargs cat > output/poaceaeHelixerOG_ancSeq_filtered_toAdd.fa

# remove gap
cat output/poaceaeHelixerOG_ancSeq_filtered_toAdd.fa |sed 's/-//g' > output/poaceaeHelixerOG_ancSeq_filtered_toAdd_gapRemoved.fa

/programs/seqkit-0.15.0/seqkit seq -m 1 poaceaeHelixerOG_ancSeq_filtered_toAdd_gapRemoved.fa > poaceaeHelixerOG_ancSeq_filtered_toAdd_gapRemoved_noempty.fa