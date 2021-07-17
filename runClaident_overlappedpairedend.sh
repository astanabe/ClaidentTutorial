# Set number of processor cores used for computation
export THREADS=32

# Move previous analysis results
mkdir -p previous

if test -e OverlappedPairedEnd_03_ConcatenatedSequences; then
mv \
OverlappedPairedEnd_03_ConcatenatedSequences \
previous/
fi

if test -e OverlappedPairedEnd_04_FilteredSequences; then
mv \
OverlappedPairedEnd_04_FilteredSequences \
previous/
fi

if test -e OverlappedPairedEnd_05_DenoisedSequences; then
mv \
OverlappedPairedEnd_05_DenoisedSequences \
previous/
fi

if test -e OverlappedPairedEnd_06_NonchimericSequences; then
mv \
OverlappedPairedEnd_06_NonchimericSequences \
previous/
fi

if test -e OverlappedPairedEnd_07_NonhoppedSequences; then
mv \
OverlappedPairedEnd_07_NonhoppedSequences \
previous/
fi

if test -e OverlappedPairedEnd_08_DecontaminatedSequences; then
mv \
OverlappedPairedEnd_08_DecontaminatedSequences \
previous/
fi

if test -e OverlappedPairedEnd_09_ClusteredSequences; then
mv \
OverlappedPairedEnd_09_ClusteredSequences \
previous/
fi

if test -e OverlappedPairedEnd_10_ClaidentResults; then
mv \
OverlappedPairedEnd_10_ClaidentResults \
previous/
fi

if test -e OverlappedPairedEnd_11_RAnalysisResults; then
mv OverlappedPairedEnd_11_RAnalysisResults \
previous/
fi

# Demultiplex Type A (If you have undemultiplexed FASTQ files)
# --seqnamestyle=illumina should be used for real Illumina outputs.
if ! test -e PairedEnd_02a_DemultiplexedSequences; then
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
PairedEnd_02a_DemultiplexedSequences
fi

# Demultiplex Type B (If FASTQ files have been already demultiplexed)
# --seqnamestyle=illumina should be used for real Illumina outputs.
if ! test -e PairedEnd_02b_DemultiplexedSequences; then
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
PairedEnd_02b_DemultiplexedSequences
done
fi

# Compare Type A and B
rm -f PairedEnd_TypeA.txt PairedEnd_TypeB.txt

cd PairedEnd_02a_DemultiplexedSequences

for f in *.fastq.xz
do echo $f >> ../PairedEnd_TypeA.txt; xz -dc $f | grep -c -P '^\+\r?\n?$' >> ../PairedEnd_TypeA.txt
done

cd ../PairedEnd_02b_DemultiplexedSequences

for f in *.fastq.xz
do echo $f >> ../PairedEnd_TypeB.txt; xz -dc $f | grep -c -P '^\+\r?\n?$' >> ../PairedEnd_TypeB.txt
done

cd ..

diff -u PairedEnd_TypeA.txt PairedEnd_TypeB.txt

# Concatenate pairs
clconcatpairv \
--mode=ovl \
--compress=xz \
--numthreads=$THREADS \
PairedEnd_02a_DemultiplexedSequences \
OverlappedPairedEnd_03_ConcatenatedSequences

# Calculate FASTQ statistics
clcalcfastqstatv \
--mode=2 \
OverlappedPairedEnd_03_ConcatenatedSequences \
OverlappedPairedEnd_03_ConcatenatedSequences/fastq_eestats2.txt

# Apply filtering out low quality sequences
clfilterseqv \
--maxqual=41 \
--minlen=100 \
--maxlen=250 \
--maxnee=2.0 \
--maxnNs=0 \
--compress=xz \
--numthreads=$THREADS \
OverlappedPairedEnd_03_ConcatenatedSequences \
OverlappedPairedEnd_04_FilteredSequences

# Denoise using DADA2
cldenoiseseqd \
--pool=enable \
--numthreads=$THREADS \
OverlappedPairedEnd_04_FilteredSequences \
OverlappedPairedEnd_05_DenoisedSequences

# Remove chimeras using UCHIME3
clremovechimev \
--mode=both \
--uchimedenovo=3 \
--referencedb=cdu12s \
--numthreads=$THREADS \
OverlappedPairedEnd_05_DenoisedSequences \
OverlappedPairedEnd_06_NonchimericSequences

# Eliminate index-hopping
# This step cannot apply to TypeB demultiplexed sequences
clremovecontam \
--index1file=index1.fasta \
--index2file=index2.fasta \
--mode=eliminate \
OverlappedPairedEnd_06_NonchimericSequences \
OverlappedPairedEnd_07_NonhoppedSequences

# Eliminate contamination
# Note that this process is incompatible with normalization of concentration/sequencing depth.
# Do not apply this process in such cases.
clremovecontam \
--blanklist=blanklist.txt \
--mode=eliminate \
OverlappedPairedEnd_07_NonhoppedSequences \
OverlappedPairedEnd_08_DecontaminatedSequences

# Cluster remaining sequences
# Note that this step is meaningless on this data because additional clustering has no effect.
clclassseqv \
--minident=0.99 \
--strand=plus \
--numthreads=$THREADS \
OverlappedPairedEnd_08_DecontaminatedSequences \
OverlappedPairedEnd_09_ClusteredSequences

# Make final output folder
mkdir -p OverlappedPairedEnd_10_ClaidentResults

# Assign taxonomy based on QCauto method using animals_mt_species
clmakecachedb \
--blastdb=animals_mt_species \
--numthreads=$THREADS \
OverlappedPairedEnd_09_ClusteredSequences/clustered.fasta \
OverlappedPairedEnd_10_ClaidentResults/cachedb_species

clidentseq \
--method=QC \
--blastdb=OverlappedPairedEnd_10_ClaidentResults/cachedb_species \
--numthreads=$THREADS \
OverlappedPairedEnd_09_ClusteredSequences/clustered.fasta \
OverlappedPairedEnd_10_ClaidentResults/neighborhoods_qc_species.txt

classigntax \
--taxdb=animals_mt_species \
OverlappedPairedEnd_10_ClaidentResults/neighborhoods_qc_species.txt \
OverlappedPairedEnd_10_ClaidentResults/taxonomy_qc_species.tsv

# Assign taxonomy based on (95%-)1-NN method using animals_mt_species
clidentseq \
--method=1,95% \
--blastdb=OverlappedPairedEnd_10_ClaidentResults/cachedb_species \
--numthreads=$THREADS \
OverlappedPairedEnd_09_ClusteredSequences/clustered.fasta \
OverlappedPairedEnd_10_ClaidentResults/neighborhoods_1nn_species.txt

classigntax \
--taxdb=animals_mt_species \
--minnsupporter=1 \
OverlappedPairedEnd_10_ClaidentResults/neighborhoods_1nn_species.txt \
OverlappedPairedEnd_10_ClaidentResults/taxonomy_1nn_species.tsv

# Assign taxonomy based on QCauto method using animals_mt_species_wsp
clmakecachedb \
--blastdb=animals_mt_species_wsp \
--numthreads=$THREADS \
OverlappedPairedEnd_09_ClusteredSequences/clustered.fasta \
OverlappedPairedEnd_10_ClaidentResults/cachedb_species_wsp

clidentseq \
--method=QC \
--blastdb=OverlappedPairedEnd_10_ClaidentResults/cachedb_species_wsp \
--numthreads=$THREADS \
OverlappedPairedEnd_09_ClusteredSequences/clustered.fasta \
OverlappedPairedEnd_10_ClaidentResults/neighborhoods_qc_species_wsp.txt

classigntax \
--taxdb=animals_mt_species_wsp \
OverlappedPairedEnd_10_ClaidentResults/neighborhoods_qc_species_wsp.txt \
OverlappedPairedEnd_10_ClaidentResults/taxonomy_qc_species_wsp.tsv

# Assign taxonomy based on (95%-)1-NN method using animals_mt_species_wsp
clidentseq \
--method=1,95% \
--blastdb=OverlappedPairedEnd_10_ClaidentResults/cachedb_species_wsp \
--numthreads=$THREADS \
OverlappedPairedEnd_09_ClusteredSequences/clustered.fasta \
OverlappedPairedEnd_10_ClaidentResults/neighborhoods_1nn_species_wsp.txt

classigntax \
--taxdb=animals_mt_species_wsp \
--minnsupporter=1 \
OverlappedPairedEnd_10_ClaidentResults/neighborhoods_1nn_species_wsp.txt \
OverlappedPairedEnd_10_ClaidentResults/taxonomy_1nn_species_wsp.tsv

# Assign taxonomy based on QCauto method using animals_mt_species_wosp
clmakecachedb \
--blastdb=animals_mt_species_wosp \
--numthreads=$THREADS \
OverlappedPairedEnd_09_ClusteredSequences/clustered.fasta \
OverlappedPairedEnd_10_ClaidentResults/cachedb_species_wosp

clidentseq \
--method=QC \
--blastdb=OverlappedPairedEnd_10_ClaidentResults/cachedb_species_wosp \
--numthreads=$THREADS \
OverlappedPairedEnd_09_ClusteredSequences/clustered.fasta \
OverlappedPairedEnd_10_ClaidentResults/neighborhoods_qc_species_wosp.txt

classigntax \
--taxdb=animals_mt_species_wosp \
OverlappedPairedEnd_10_ClaidentResults/neighborhoods_qc_species_wosp.txt \
OverlappedPairedEnd_10_ClaidentResults/taxonomy_qc_species_wosp.tsv

# Assign taxonomy based on (95%-)1-NN method using animals_mt_species_wosp
clidentseq \
--method=1,95% \
--blastdb=OverlappedPairedEnd_10_ClaidentResults/cachedb_species_wosp \
--numthreads=$THREADS \
OverlappedPairedEnd_09_ClusteredSequences/clustered.fasta \
OverlappedPairedEnd_10_ClaidentResults/neighborhoods_1nn_species_wosp.txt

classigntax \
--taxdb=animals_mt_species_wosp \
--minnsupporter=1 \
OverlappedPairedEnd_10_ClaidentResults/neighborhoods_1nn_species_wosp.txt \
OverlappedPairedEnd_10_ClaidentResults/taxonomy_1nn_species_wosp.tsv

# Merge 6 taxonomic assignment results
# Note that merge of QCauto results and (95%-)1-NN results has no effects in many cases because (95%-)1-NN results are always consistent to QCauto results excluding the case when there is no 95% or more similar reference sequences to the query.
# However, merge of results using different reference database is often useful.
clmergeassign \
--preferlower \
--priority=descend \
OverlappedPairedEnd_10_ClaidentResults/taxonomy_qc_species_wosp.tsv \
OverlappedPairedEnd_10_ClaidentResults/taxonomy_qc_species.tsv \
OverlappedPairedEnd_10_ClaidentResults/taxonomy_qc_species_wsp.tsv \
OverlappedPairedEnd_10_ClaidentResults/taxonomy_1nn_species_wosp.tsv \
OverlappedPairedEnd_10_ClaidentResults/taxonomy_1nn_species.tsv \
OverlappedPairedEnd_10_ClaidentResults/taxonomy_1nn_species_wsp.tsv \
OverlappedPairedEnd_10_ClaidentResults/taxonomy_merged.tsv

# Fill blank cells of taxonomic assignment
clfillassign \
OverlappedPairedEnd_10_ClaidentResults/taxonomy_merged.tsv \
OverlappedPairedEnd_10_ClaidentResults/taxonomy_merged_filled.tsv

# Filter out non-Actinopterygii/Sarcopterygii OTUs
clfiltersum \
--taxfile=OverlappedPairedEnd_10_ClaidentResults/taxonomy_merged_filled.tsv \
--includetaxa=superclass,Actinopterygii,superclass,Sarcopterygii \
OverlappedPairedEnd_09_ClusteredSequences/clustered.tsv \
OverlappedPairedEnd_10_ClaidentResults/sample_otu_matrix_fishes.tsv

# Make top-50 species community data matrix for barplot
clsumtaxa \
--tableformat=column \
--taxfile=OverlappedPairedEnd_10_ClaidentResults/taxonomy_merged_filled.tsv \
--targetrank=species \
--topN=50 \
--numbering=enable \
OverlappedPairedEnd_10_ClaidentResults/sample_otu_matrix_fishes.tsv \
OverlappedPairedEnd_10_ClaidentResults/sample_top50species_nreads_fishes.tsv

# Make top-50 families community data matrix for barplot
clsumtaxa \
--tableformat=column \
--taxfile=OverlappedPairedEnd_10_ClaidentResults/taxonomy_merged_filled.tsv \
--targetrank=family \
--topN=50 \
--numbering=enable \
OverlappedPairedEnd_10_ClaidentResults/sample_otu_matrix_fishes.tsv \
OverlappedPairedEnd_10_ClaidentResults/sample_top50family_nreads_fishes.tsv

# Make species-based community data matrix for heatmap
clsumtaxa \
--tableformat=column \
--taxfile=OverlappedPairedEnd_10_ClaidentResults/taxonomy_merged_filled.tsv \
--targetrank=species \
--numbering=enable \
OverlappedPairedEnd_10_ClaidentResults/sample_otu_matrix_fishes.tsv \
OverlappedPairedEnd_10_ClaidentResults/sample_species_nreads_fishes.tsv

# Make family-based community data matrix for heatmap
clsumtaxa \
--tableformat=column \
--taxfile=OverlappedPairedEnd_10_ClaidentResults/taxonomy_merged_filled.tsv \
--targetrank=family \
--numbering=enable \
OverlappedPairedEnd_10_ClaidentResults/sample_otu_matrix_fishes.tsv \
OverlappedPairedEnd_10_ClaidentResults/sample_family_nreads_fishes.tsv

# Run R
Rscript runR_overlappedpairedend.R

# Remove cachedb
rm -r OverlappedPairedEnd_10_ClaidentResults/cachedb_species*
