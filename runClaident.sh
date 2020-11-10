export THREADS=32

# Move previous analysis results
mkdir -p previous
mv 02a_DemultiplexedSequences 02b_DemultiplexedSequences 03_ConcatenatedSequences 04_FilteredSequences 05_DenoisedSequences 06_NonchimericSequences 07_NonhoppedSequences 08_DecontaminatedSequences 09_ClusteredSequences 10_ClaidentResults 11_RAnalysisResults previous/

# Demultiplex Type A (If you have undemultiplexed FASTQ files)
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
01_RawSequences/Undemultiplexed_R1_001.fastq.xz \
01_RawSequences/Undemultiplexed_I1_001.fastq.xz \
01_RawSequences/Undemultiplexed_I2_001.fastq.xz \
01_RawSequences/Undemultiplexed_R2_001.fastq.xz \
02a_DemultiplexedSequences

# Demultiplex Type B (If FASTQ files have been already demultiplexed)
for s in `ls 01_RawSequences/Blank??_R1_001.fastq.xz 01_RawSequences/Sample??_R1_001.fastq.xz | grep -o -P '[A-Z][a-z]+\d\d'`
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
01_RawSequences/$s\_R1_001.fastq.xz \
01_RawSequences/$s\_R2_001.fastq.xz \
02b_DemultiplexedSequences
done

# Compare Type A and B
rm -f TypeA.txt TypeB.txt

cd 02a_DemultiplexedSequences

for f in *.fastq.xz
do echo $f >> ../TypeA.txt; xz -dc $f | grep -c -P '^\+\r?\n?$' >> ../TypeA.txt
done

cd ../02b_DemultiplexedSequences

for f in *.fastq.xz
do echo $f >> ../TypeB.txt; xz -dc $f | grep -c -P '^\+\r?\n?$' >> ../TypeB.txt
done

cd ..

diff -u TypeA.txt TypeB.txt

# Concatenate pairs
clconcatpairv \
--compress=xz \
--numthreads=$THREADS \
02a_DemultiplexedSequences \
03_ConcatenatedSequences

# Calculate FASTQ statistics
clcalcfastqstatv \
--mode=2 \
03_ConcatenatedSequences \
03_ConcatenatedSequences/fastq_eestats2.txt

# Filfer out low quality sequences
clfilterseqv \
--maxqual=41 \
--minlen=100 \
--maxlen=250 \
--maxnee=2.0 \
--maxnNs=0 \
--compress=xz \
--numthreads=$THREADS \
03_ConcatenatedSequences \
04_FilteredSequences

# Denoise using DADA2
cldenoiseseqd \
--pool=enable \
--numthreads=$THREADS \
04_FilteredSequences \
05_DenoisedSequences

# Remove chimeras using UCHIME3
clremovechimev \
--mode=both \
--uchimedenovo=3 \
--referencedb=cdu12s \
--numthreads=$THREADS \
05_DenoisedSequences \
06_NonchimericSequences

# Eliminate index-hopping
clremovecontam \
--index1file=index1.fasta \
--index2file=index2.fasta \
--mode=eliminate \
06_NonchimericSequences \
07_NonhoppedSequences

# Eliminate contamination
clremovecontam \
--blanklist=blanklist.txt \
--mode=eliminate \
07_NonhoppedSequences \
08_DecontaminatedSequences

# Cluster remaining sequences
clclassseqv \
--minident=0.99 \
--strand=plus \
--numthreads=$THREADS \
08_DecontaminatedSequences \
09_ClusteredSequences

# Make final output folder
mkdir -p 10_ClaidentResults

# Assign taxonomy based on QCauto method
clmakecachedb \
--blastdb=animals_mt_species \
--numthreads=$THREADS \
09_ClusteredSequences/clustered.fasta \
10_ClaidentResults/cachedb

clidentseq \
--method=QC \
--blastdb=10_ClaidentResults/cachedb \
--numthreads=$THREADS \
09_ClusteredSequences/clustered.fasta \
10_ClaidentResults/neighborhoods_qc.txt

classigntax \
--taxdb=animals_mt_species \
10_ClaidentResults/neighborhoods_qc.txt \
10_ClaidentResults/taxonomy_qc.tsv

# Assign taxonomy based on 1-NN method
clidentseq \
--method=1,95% \
--blastdb=10_ClaidentResults/cachedb \
--numthreads=$THREADS \
09_ClusteredSequences/clustered.fasta \
10_ClaidentResults/neighborhoods_1nn.txt

classigntax \
--taxdb=animals_mt_species \
 --minnsupporter=1 \
10_ClaidentResults/neighborhoods_1nn.txt \
10_ClaidentResults/taxonomy_1nn.tsv

# Merge 2 taxonomic assignment results
clmergeassign \
--preferlower \
--priority=descend \
10_ClaidentResults/taxonomy_qc.tsv \
10_ClaidentResults/taxonomy_1nn.tsv \
10_ClaidentResults/taxonomy_merged.tsv

# Fill blank cells of taxonomic assignment
clfillassign \
10_ClaidentResults/taxonomy_merged.tsv \
10_ClaidentResults/taxonomy_merged_filled.tsv

# Filter out non-Actinopterygii/Sarcopterygii OTUs
clfiltersum \
--taxfile=10_ClaidentResults/taxonomy_merged_filled.tsv \
--includetaxa=superclass,Actinopterygii,superclass,Sarcopterygii \
09_ClusteredSequences/clustered.tsv \
10_ClaidentResults/sample_otu_matrix_fishes.tsv

# Make species-based community data matrix
clsumtaxa \
--tableformat=matrix \
--taxfile=10_ClaidentResults/taxonomy_merged_filled.tsv \
--targetrank=species \
--numbering=disable \
10_ClaidentResults/sample_otu_matrix_fishes.tsv \
10_ClaidentResults/sample_species_matrix_fishes.tsv

# Make top-50 species community data matrix for barplot
clsumtaxa \
--tableformat=column \
--taxfile=10_ClaidentResults/taxonomy_merged_filled.tsv \
--targetrank=species \
--topN=50 \
--numbering=enable \
10_ClaidentResults/sample_otu_matrix_fishes.tsv \
10_ClaidentResults/sample_top50species_matrix_fishes.tsv

# Make top-50 families community data matrix for barplot
clsumtaxa \
--tableformat=column \
--taxfile=10_ClaidentResults/taxonomy_merged_filled.tsv \
--targetrank=family \
--topN=50 \
--numbering=enable \
10_ClaidentResults/sample_otu_matrix_fishes.tsv \
10_ClaidentResults/sample_top50family_matrix_fishes.tsv
