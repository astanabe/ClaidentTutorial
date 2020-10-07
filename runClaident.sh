if test -z $THREADS; then
export THREADS=8 || exit $?
fi
# Demultiplex Type A
clsplitseq \
--runname=ClaidentTutorial \
--forwardprimerfile=forwardprimer.fasta \
--reverseprimerfile=reverseprimer.fasta \
--truncateN=enable \
--index1file=index1.fasta \
--index2file=index2.fasta \
--minqualtag=30 \
--compress=xz \
--numthreads=$THREADS \
--seqnamestyle=other \
01_RawSequences/Undemultiplexed_R1.fastq.xz \
01_RawSequences/Undemultiplexed_I1.fastq.xz \
01_RawSequences/Undemultiplexed_I2.fastq.xz \
01_RawSequences/Undemultiplexed_R2.fastq.xz \
02a_DemultiplexedSequences
# Demultiplex Type B
for s in `ls 01_RawSequences/Blank??_R1.fastq.xz 01_RawSequences/Sample??_R1.fastq.xz | grep -o -P '[A-Z][a-z]+\d\d'`
do clsplitseq \
--runname=ClaidentTutorial \
--indexname=$s \
--forwardprimerfile=forwardprimer.fasta \
--reverseprimerfile=reverseprimer.fasta \
--truncateN=enable \
--compress=xz \
--numthreads=$THREADS \
--seqnamestyle=other \
--append \
01_RawSequences/$s\_R1.fastq.xz \
01_RawSequences/$s\_R2.fastq.xz \
02b_DemultiplexedSequences
done
# Compare Type A and B
rm -f TypeA.txt TypeB.txt
for f in 02a_DemultiplexedSequences/*.fastq.xz
do xz -dc $f | grep -c -P '^\+\r?\n?$' >> TypeA.txt
done
for f in 02b_DemultiplexedSequences/*.fastq.xz
do xz -dc $f | grep -c -P '^\+\r?\n?$' >> TypeB.txt
done
diff -u TypeA.txt TypeB.txt
# Concatenate pairs
clconcatpair \
--compress=xz \
--numthreads=$THREADS \
02a_DemultiplexedSequences \
03_ConcatenatedSequences
# Filfer out low quality sequences
clfilterseqv \
--minlen=100 \
--maxlen=250 \
--maxnee=1.0 \
--maxnN=0 \
--compress=xz \
--numthreads=$THREADS \
03_ConcatenatedSequences \
04_FilteredSequences
# Denoise using DADA2
cldenoiseseqd \
--pool=pseudo \
--numthreads=$THREADS \
04_FilteredSequences \
05_DenoisedSequences
# Remove chimeras using UCHIME3
clremovechimera \
--uchimedenovo=3 \
--referencedb=cdu12s \
--numthreads=$THREADS \
05_DenoisedSequences \
06_NonchimericSequences
# Subtract contamination
clsubtractcontam \
--blank='Blank\d\d' \
06_NonchimericSequences \
07_DecontaminatedSequences
# Cluster remaining sequences
clclassseqv \
--minident=0.99 \
--strand=plus \
--numthreads=$THREADS \
07_DecontaminatedSequences \
08_ClusteredSequences
# Make final output folder
mkdir -p 09_FinalResults
# Assign taxonomy based on QCauto method
clmakecachedb \
--bdb=animals_mt_species \
--numthreads=8 \
08_ClusteredSequences/clustered.fasta \
09_FinalResults/cachedb
clidentseq \
--method=QC \
--bdb=09_FinalResults/cachedb \
--numthreads=8 \
08_ClusteredSequences/clustered.fasta \
09_FinalResults/neighborhoods_qc.txt
classigntax \
--taxdb=animals_mt_species \
09_FinalResults/neighborhoods_qc.txt \
09_FinalResults/taxonomy_qc.txt
# Assign taxonomy based on 1-NN method
clidentseq \
--method=1,95% \
--bdb=09_FinalResults/cachedb \
--numthreads=8 \
08_ClusteredSequences/clustered.fasta \
09_FinalResults/neighborhoods_1nn.txt
classigntax \
--taxdb=animals_mt_species \
 --minnsupporter=1 \
09_FinalResults/neighborhoods_1nn.txt \
09_FinalResults/taxonomy_1nn.txt
# Merge 2 taxonomic assignment results
clmergeassign \
--preferlower \
--priority=descend \
09_FinalResults/taxonomy_qc.txt \
09_FinalResults/taxonomy_1nn.txt \
09_FinalResults/taxonomy_merged.txt
# Fill blank cells of taxonomic assignment
clfillassign \
09_FinalResults/taxonomy_merged.txt \
09_FinalResults/taxonomy_merged_filled.txt
# Make OTU-based community data maxrix
clsumclass \
--output=matrix \
08_ClusteredSequences/clustered.otu.gz \
09_FinalResults/sample_otu_matrix.txt
# Filter out non-Actinopterygii OTUs
clfiltersum \
--taxfile=09_FinalResults/taxonomy_merged_filled.txt \
--taxfilter=superclass:Actinopterygii \
09_FinalResults/sample_otu_matrix.txt \
09_FinalResults/sample_otu_matrix_bonyfishes.txt
# Make species-based community data matrix
clsumtaxa \
--output=matrix \
--taxfile=09_FinalResults/taxonomy_merged_filled.txt \
--taxrank=species \
--numbering=disable \
09_FinalResults/sample_otu_matrix_bonyfishes.txt \
09_FinalResults/sample_species_matrix_bonyfishes.txt
# Make top-10 species community data matrix for barplot
clsumtaxa \
--output=column \
--taxfile=09_FinalResults/taxonomy_merged_filled.txt \
--taxrank=species \
--topN=10 \
--numbering=enable \
09_FinalResults/sample_otu_matrix_bonyfishes.txt \
09_FinalResults/sample_10species_matrix_bonyfishes.txt
# Make top-10 families community data matrix for barplot
clsumtaxa \
--output=column \
--taxfile=09_FinalResults/taxonomy_merged_filled.txt \
--taxrank=family \
--topN=10 \
--numbering=enable \
09_FinalResults/sample_otu_matrix_bonyfishes.txt \
09_FinalResults/sample_10family_matrix_bonyfishes.txt
