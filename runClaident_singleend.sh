# Set number of processor cores used for computation
export THREADS=32

# Move previous analysis results
mkdir -p previous

if test -e SingleEnd_02a_DemultiplexedSequences; then
mv \
SingleEnd_02a_DemultiplexedSequences \
previous/
fi

if test -e SingleEnd_02b_DemultiplexedSequences; then
mv \
SingleEnd_02b_DemultiplexedSequences \
previous/
fi

if test -e SingleEnd_03_FilteredSequences; then
mv \
SingleEnd_03_FilteredSequences \
previous/
fi

if test -e SingleEnd_04_DenoisedSequences; then
mv \
SingleEnd_04_DenoisedSequences \
previous/
fi

if test -e SingleEnd_05_NonchimericSequences; then
mv \
SingleEnd_05_NonchimericSequences \
previous/
fi

if test -e SingleEnd_06_NonhoppedSequences; then
mv \
SingleEnd_06_NonhoppedSequences \
previous/
fi

if test -e SingleEnd_07_DecontaminatedSequences; then
mv \
SingleEnd_07_DecontaminatedSequences \
previous/
fi

if test -e SingleEnd_08_ClusteredSequences; then
mv \
SingleEnd_08_ClusteredSequences \
previous/
fi

if test -e SingleEnd_09_ClaidentResults; then
mv \
SingleEnd_09_ClaidentResults \
previous/
fi

if test -e SingleEnd_10_RAnalysisResults; then
mv \
SingleEnd_10_RAnalysisResults \
previous/
fi

# Demultiplex Type A (If you have undemultiplexed FASTQ files)
# --seqnamestyle=illumina should be used for real Illumina outputs.
clsplitseq \
--runname=ClaidentTutorial \
--forwardprimerfile=forwardprimer.fasta \
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
SingleEnd_02a_DemultiplexedSequences

# Demultiplex Type B (If FASTQ files have been already demultiplexed)
# --seqnamestyle=illumina should be used for real Illumina outputs.
for s in `ls 01_RawSequences/Blank??_R1_001.fastq.xz 01_RawSequences/Sample??_R1_001.fastq.xz | grep -o -P '[A-Z][a-z]+\d\d'`
do clsplitseq \
--runname=ClaidentTutorial \
--indexname=$s \
--forwardprimerfile=forwardprimer.fasta \
--truncateN=enable \
--compress=xz \
--numthreads=$THREADS \
--seqnamestyle=other \
--append \
01_RawSequences/$s\_R1_001.fastq.xz \
SingleEnd_02b_DemultiplexedSequences
done

# Compare Type A and B
rm -f SingleEnd_TypeA.txt SingleEnd_TypeB.txt

cd SingleEnd_02a_DemultiplexedSequences

for f in *.fastq.xz
do echo $f >> ../SingleEnd_TypeA.txt; xz -dc $f | grep -c -P '^\+\r?\n?$' >> ../SingleEnd_TypeA.txt
done

cd ../SingleEnd_02b_DemultiplexedSequences

for f in *.fastq.xz
do echo $f >> ../SingleEnd_TypeB.txt; xz -dc $f | grep -c -P '^\+\r?\n?$' >> ../SingleEnd_TypeB.txt
done

cd ..

diff -u SingleEnd_TypeA.txt SingleEnd_TypeB.txt

# Calculate FASTQ statistics
clcalcfastqstatv \
--mode=2 \
SingleEnd_02a_DemultiplexedSequences \
SingleEnd_02a_DemultiplexedSequences/fastq_eestats2.txt

# Apply sequence trimming and filtering out low quality sequences
# Trimmed length determined based on fastq_eestats2.txt
clfilterseqv \
--maxqual=41 \
--minlen=120 \
--maxlen=120 \
--maxnee=2.0 \
--maxnNs=0 \
--compress=xz \
--numthreads=$THREADS \
SingleEnd_02a_DemultiplexedSequences \
SingleEnd_03_FilteredSequences

# Denoise using DADA2
cldenoiseseqd \
--pool=enable \
--numthreads=$THREADS \
SingleEnd_03_FilteredSequences \
SingleEnd_04_DenoisedSequences

# Remove chimeras using UCHIME3
clremovechimev \
--mode=both \
--uchimedenovo=3 \
--referencedb=cdu12s \
--numthreads=$THREADS \
SingleEnd_04_DenoisedSequences \
SingleEnd_05_NonchimericSequences

# Eliminate index-hopping
# This step cannot apply to TypeB demultiplexed sequences and/or single index sequences
clremovecontam \
--index1file=index1.fasta \
--index2file=index2.fasta \
--mode=eliminate \
SingleEnd_05_NonchimericSequences \
SingleEnd_06_NonhoppedSequences

# Eliminate contamination
# Note that this process is incompatible with normalization of concentration/sequencing depth.
# Do not apply this process in such cases.
clremovecontam \
--blanklist=blanklist.txt \
--mode=eliminate \
SingleEnd_06_NonhoppedSequences \
SingleEnd_07_DecontaminatedSequences

# Cluster remaining sequences
# Note that this step is meaningless on this data because additional clustering has no effect.
clclassseqv \
--minident=0.99 \
--strand=plus \
--numthreads=$THREADS \
SingleEnd_07_DecontaminatedSequences \
SingleEnd_08_ClusteredSequences

# Make final output folder
mkdir -p SingleEnd_09_ClaidentResults

# Assign taxonomy based on QCauto method using animals_mt_species
clmakecachedb \
--blastdb=animals_mt_species \
--numthreads=$THREADS \
SingleEnd_08_ClusteredSequences/clustered.fasta \
SingleEnd_09_ClaidentResults/cachedb_species

clidentseq \
--method=QC \
--blastdb=SingleEnd_09_ClaidentResults/cachedb_species \
--numthreads=$THREADS \
SingleEnd_08_ClusteredSequences/clustered.fasta \
SingleEnd_09_ClaidentResults/neighborhoods_qc_species.txt

classigntax \
--taxdb=animals_mt_species \
SingleEnd_09_ClaidentResults/neighborhoods_qc_species.txt \
SingleEnd_09_ClaidentResults/taxonomy_qc_species.tsv

# Assign taxonomy based on (95%-)1-NN method using animals_mt_species
clidentseq \
--method=1,95% \
--blastdb=SingleEnd_09_ClaidentResults/cachedb_species \
--numthreads=$THREADS \
SingleEnd_08_ClusteredSequences/clustered.fasta \
SingleEnd_09_ClaidentResults/neighborhoods_1nn_species.txt

classigntax \
--taxdb=animals_mt_species \
 --minnsupporter=1 \
SingleEnd_09_ClaidentResults/neighborhoods_1nn_species.txt \
SingleEnd_09_ClaidentResults/taxonomy_1nn_species.tsv

# Assign taxonomy based on QCauto method using animals_mt_species_wsp
clmakecachedb \
--blastdb=animals_mt_species_wsp \
--numthreads=$THREADS \
SingleEnd_08_ClusteredSequences/clustered.fasta \
SingleEnd_09_ClaidentResults/cachedb_species_wsp

clidentseq \
--method=QC \
--blastdb=SingleEnd_09_ClaidentResults/cachedb_species_wsp \
--numthreads=$THREADS \
SingleEnd_08_ClusteredSequences/clustered.fasta \
SingleEnd_09_ClaidentResults/neighborhoods_qc_species_wsp.txt

classigntax \
--taxdb=animals_mt_species_wsp \
SingleEnd_09_ClaidentResults/neighborhoods_qc_species_wsp.txt \
SingleEnd_09_ClaidentResults/taxonomy_qc_species_wsp.tsv

# Assign taxonomy based on (95%-)1-NN method using animals_mt_species_wsp
clidentseq \
--method=1,95% \
--blastdb=SingleEnd_09_ClaidentResults/cachedb_species_wsp \
--numthreads=$THREADS \
SingleEnd_08_ClusteredSequences/clustered.fasta \
SingleEnd_09_ClaidentResults/neighborhoods_1nn_species_wsp.txt

classigntax \
--taxdb=animals_mt_species_wsp \
 --minnsupporter=1 \
SingleEnd_09_ClaidentResults/neighborhoods_1nn_species_wsp.txt \
SingleEnd_09_ClaidentResults/taxonomy_1nn_species_wsp.tsv

# Assign taxonomy based on QCauto method using animals_mt_species_wosp
clmakecachedb \
--blastdb=animals_mt_species_wosp \
--numthreads=$THREADS \
SingleEnd_08_ClusteredSequences/clustered.fasta \
SingleEnd_09_ClaidentResults/cachedb_species_wosp

clidentseq \
--method=QC \
--blastdb=SingleEnd_09_ClaidentResults/cachedb_species_wosp \
--numthreads=$THREADS \
SingleEnd_08_ClusteredSequences/clustered.fasta \
SingleEnd_09_ClaidentResults/neighborhoods_qc_species_wosp.txt

classigntax \
--taxdb=animals_mt_species_wosp \
SingleEnd_09_ClaidentResults/neighborhoods_qc_species_wosp.txt \
SingleEnd_09_ClaidentResults/taxonomy_qc_species_wosp.tsv

# Assign taxonomy based on (95%-)1-NN method using animals_mt_species_wosp
clidentseq \
--method=1,95% \
--blastdb=SingleEnd_09_ClaidentResults/cachedb_species_wosp \
--numthreads=$THREADS \
SingleEnd_08_ClusteredSequences/clustered.fasta \
SingleEnd_09_ClaidentResults/neighborhoods_1nn_species_wosp.txt

classigntax \
--taxdb=animals_mt_species_wosp \
 --minnsupporter=1 \
SingleEnd_09_ClaidentResults/neighborhoods_1nn_species_wosp.txt \
SingleEnd_09_ClaidentResults/taxonomy_1nn_species_wosp.tsv

# Merge 6 taxonomic assignment results
# Note that merge of QCauto results and (95%-)1-NN results is meaningless in many cases because (95%-)1-NN results are always consistent to QCauto results excluding the case when there is no 95% or more similar reference sequences.
# However, merge of results using different reference database is often useful.
clmergeassign \
--preferlower \
--priority=descend \
SingleEnd_09_ClaidentResults/taxonomy_qc_species_wosp.tsv \
SingleEnd_09_ClaidentResults/taxonomy_qc_species.tsv \
SingleEnd_09_ClaidentResults/taxonomy_qc_species_wsp.tsv \
SingleEnd_09_ClaidentResults/taxonomy_1nn_species_wosp.tsv \
SingleEnd_09_ClaidentResults/taxonomy_1nn_species.tsv \
SingleEnd_09_ClaidentResults/taxonomy_1nn_species_wsp.tsv \
SingleEnd_09_ClaidentResults/taxonomy_merged.tsv

# Fill blank cells of taxonomic assignment
clfillassign \
SingleEnd_09_ClaidentResults/taxonomy_merged.tsv \
SingleEnd_09_ClaidentResults/taxonomy_merged_filled.tsv

# Filter out non-Actinopterygii/Sarcopterygii OTUs
clfiltersum \
--taxfile=SingleEnd_09_ClaidentResults/taxonomy_merged_filled.tsv \
--includetaxa=superclass,Actinopterygii,superclass,Sarcopterygii \
SingleEnd_08_ClusteredSequences/clustered.tsv \
SingleEnd_09_ClaidentResults/sample_otu_matrix_fishes.tsv

# Make species-based community data matrix
clsumtaxa \
--tableformat=matrix \
--taxfile=SingleEnd_09_ClaidentResults/taxonomy_merged_filled.tsv \
--targetrank=species \
--numbering=disable \
SingleEnd_09_ClaidentResults/sample_otu_matrix_fishes.tsv \
SingleEnd_09_ClaidentResults/sample_species_matrix_fishes.tsv

# Make top-50 species community data matrix for barplot
clsumtaxa \
--tableformat=column \
--taxfile=SingleEnd_09_ClaidentResults/taxonomy_merged_filled.tsv \
--targetrank=species \
--topN=50 \
--numbering=enable \
SingleEnd_09_ClaidentResults/sample_otu_matrix_fishes.tsv \
SingleEnd_09_ClaidentResults/sample_top50species_nreads_fishes.tsv

# Make top-50 families community data matrix for barplot
clsumtaxa \
--tableformat=column \
--taxfile=SingleEnd_09_ClaidentResults/taxonomy_merged_filled.tsv \
--targetrank=family \
--topN=50 \
--numbering=enable \
SingleEnd_09_ClaidentResults/sample_otu_matrix_fishes.tsv \
SingleEnd_09_ClaidentResults/sample_top50family_nreads_fishes.tsv

# Make species-based community data matrix for heatmap
clsumtaxa \
--tableformat=column \
--taxfile=SingleEnd_09_ClaidentResults/taxonomy_merged_filled.tsv \
--targetrank=species \
--numbering=enable \
SingleEnd_09_ClaidentResults/sample_otu_matrix_fishes.tsv \
SingleEnd_09_ClaidentResults/sample_species_nreads_fishes.tsv

# Make family-based community data matrix for heatmap
clsumtaxa \
--tableformat=column \
--taxfile=SingleEnd_09_ClaidentResults/taxonomy_merged_filled.tsv \
--targetrank=family \
--numbering=enable \
SingleEnd_09_ClaidentResults/sample_otu_matrix_fishes.tsv \
SingleEnd_09_ClaidentResults/sample_family_nreads_fishes.tsv
