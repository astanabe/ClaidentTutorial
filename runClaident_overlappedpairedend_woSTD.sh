# Set number of processor cores used for computation
export THREADS=32

# Move previous analysis results
mkdir -p previous

if test -e OverlappedPairedEnd_woSTD_03_ConcatenatedSequences; then
mv \
OverlappedPairedEnd_woSTD_03_ConcatenatedSequences \
previous/
fi

if test -e OverlappedPairedEnd_woSTD_04_FilteredSequences; then
mv \
OverlappedPairedEnd_woSTD_04_FilteredSequences \
previous/
fi

if test -e OverlappedPairedEnd_woSTD_05_DenoisedSequences; then
mv \
OverlappedPairedEnd_woSTD_05_DenoisedSequences \
previous/
fi

if test -e OverlappedPairedEnd_woSTD_06_NonchimericSequences; then
mv \
OverlappedPairedEnd_woSTD_06_NonchimericSequences \
previous/
fi

if test -e OverlappedPairedEnd_woSTD_07_NonhoppedSequences; then
mv \
OverlappedPairedEnd_woSTD_07_NonhoppedSequences \
previous/
fi

if test -e OverlappedPairedEnd_woSTD_08_DecontaminatedSequences; then
mv \
OverlappedPairedEnd_woSTD_08_DecontaminatedSequences \
previous/
fi

if test -e OverlappedPairedEnd_woSTD_09_ClusteredSequences; then
mv \
OverlappedPairedEnd_woSTD_09_ClusteredSequences \
previous/
fi

if test -e OverlappedPairedEnd_woSTD_10_ClaidentResults; then
mv \
OverlappedPairedEnd_woSTD_10_ClaidentResults \
previous/
fi

if test -e OverlappedPairedEnd_woSTD_11_RAnalysisResults; then
mv OverlappedPairedEnd_woSTD_11_RAnalysisResults \
previous/
fi

# Demultiplex Type A (If you have undemultiplexed FASTQ files)
# --seqnamestyle=illumina should be used for real Illumina outputs.
if ! test -e PairedEnd_woSTD_02a_DemultiplexedSequences; then
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
--seqnamestyle=illumina \
01a_RawSequences_woSTD/Undetermined_S0_L001_R1_001.fastq.gz \
01a_RawSequences_woSTD/Undetermined_S0_L001_I1_001.fastq.gz \
01a_RawSequences_woSTD/Undetermined_S0_L001_I2_001.fastq.gz \
01a_RawSequences_woSTD/Undetermined_S0_L001_R2_001.fastq.gz \
PairedEnd_woSTD_02a_DemultiplexedSequences
fi

# Demultiplex Type B (If FASTQ files have been already demultiplexed)
# --seqnamestyle=illumina should be used for real Illumina outputs.
if ! test -e PairedEnd_woSTD_02b_DemultiplexedSequences; then
cltruncprimer \
--runname=ClaidentTutorial \
--forwardprimerfile=forwardprimer.fasta \
--reverseprimerfile=reverseprimer.fasta \
--truncateN=enable \
--index1file=index1.fasta \
--index2file=index2.fasta \
--compress=xz \
--numthreads=$THREADS \
--seqnamestyle=illumina \
01a_RawSequences_woSTD/Sample??_R?_001.fastq.gz \
01a_RawSequences_woSTD/Blank??_R?_001.fastq.gz \
PairedEnd_woSTD_02b_DemultiplexedSequences
fi

# Compare Type A and B
rm -f PairedEnd_TypeA.txt PairedEnd_TypeB.txt

cd PairedEnd_woSTD_02a_DemultiplexedSequences

for f in *.fastq.xz
do echo $f >> ../PairedEnd_TypeA.txt; xz -dc $f | grep -c -P '^\+\r?\n?$' >> ../PairedEnd_TypeA.txt
done

cd ../PairedEnd_woSTD_02b_DemultiplexedSequences

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
PairedEnd_woSTD_02a_DemultiplexedSequences \
OverlappedPairedEnd_woSTD_03_ConcatenatedSequences

# Calculate FASTQ statistics
clcalcfastqstatv \
--mode=2 \
OverlappedPairedEnd_woSTD_03_ConcatenatedSequences \
OverlappedPairedEnd_woSTD_03_ConcatenatedSequences/fastq_eestats2.txt

# Apply filtering out low quality sequences
clfilterseqv \
--maxqual=41 \
--minlen=100 \
--maxlen=250 \
--maxnee=2.0 \
--maxnNs=0 \
--compress=xz \
--numthreads=$THREADS \
OverlappedPairedEnd_woSTD_03_ConcatenatedSequences \
OverlappedPairedEnd_woSTD_04_FilteredSequences

# Denoise using DADA2
cldenoiseseqd \
--pool=pseudo \
--numthreads=$THREADS \
OverlappedPairedEnd_woSTD_04_FilteredSequences \
OverlappedPairedEnd_woSTD_05_DenoisedSequences

# Remove chimeras using UCHIME3
clremovechimev \
--mode=both \
--uchimedenovo=3 \
--referencedb=cdu12s \
--numthreads=$THREADS \
OverlappedPairedEnd_woSTD_05_DenoisedSequences \
OverlappedPairedEnd_woSTD_06_NonchimericSequences

# Eliminate index-hopping
# This step cannot apply to TypeB demultiplexed sequences
clremovecontam \
--test=thompson \
--ignoresamplelist=blanklist.txt \
--index1file=index1.fasta \
--index2file=index2.fasta \
--numthreads=$THREADS \
OverlappedPairedEnd_woSTD_06_NonchimericSequences \
OverlappedPairedEnd_woSTD_07_NonhoppedSequences

# Eliminate contamination
# Note that this process is incompatible with normalization of concentration/sequencing depth.
# Do not apply this process in such cases.
clremovecontam \
--test=thompson \
--blanklist=blanklist.txt \
--numthreads=$THREADS \
OverlappedPairedEnd_woSTD_07_NonhoppedSequences \
OverlappedPairedEnd_woSTD_08_DecontaminatedSequences

# Cluster remaining sequences
# Note that this step is meaningless on this data because additional clustering has no effect.
clclassseqv \
--minident=1.0 \
--strand=plus \
--numthreads=$THREADS \
OverlappedPairedEnd_woSTD_08_DecontaminatedSequences \
OverlappedPairedEnd_woSTD_09_ClusteredSequences

# Make final output folder
mkdir -p OverlappedPairedEnd_woSTD_10_ClaidentResults

# Assign taxonomy based on QCauto method using animals_mt_species
clmakecachedb \
--blastdb=animals_mt_species \
--numthreads=$THREADS \
OverlappedPairedEnd_woSTD_09_ClusteredSequences/clustered.fasta \
OverlappedPairedEnd_woSTD_10_ClaidentResults/cachedb_species

clidentseq \
--method=QC \
--blastdb=OverlappedPairedEnd_woSTD_10_ClaidentResults/cachedb_species \
--numthreads=$THREADS \
OverlappedPairedEnd_woSTD_09_ClusteredSequences/clustered.fasta \
OverlappedPairedEnd_woSTD_10_ClaidentResults/neighborhoods_qc_species.txt

classigntax \
--taxdb=animals_mt_species \
OverlappedPairedEnd_woSTD_10_ClaidentResults/neighborhoods_qc_species.txt \
OverlappedPairedEnd_woSTD_10_ClaidentResults/taxonomy_qc_species.tsv

# Assign taxonomy based on (95%-)3-NN method using animals_mt_species
clidentseq \
--method=3,95% \
--blastdb=OverlappedPairedEnd_woSTD_10_ClaidentResults/cachedb_species \
--numthreads=$THREADS \
OverlappedPairedEnd_woSTD_09_ClusteredSequences/clustered.fasta \
OverlappedPairedEnd_woSTD_10_ClaidentResults/neighborhoods_3nn_species.txt

classigntax \
--taxdb=animals_mt_species \
--minnsupporter=3 \
OverlappedPairedEnd_woSTD_10_ClaidentResults/neighborhoods_3nn_species.txt \
OverlappedPairedEnd_woSTD_10_ClaidentResults/taxonomy_3nn_species.tsv

# Assign taxonomy based on QCauto method using animals_mt_species_wsp
clmakecachedb \
--blastdb=animals_mt_species_wsp \
--numthreads=$THREADS \
OverlappedPairedEnd_woSTD_09_ClusteredSequences/clustered.fasta \
OverlappedPairedEnd_woSTD_10_ClaidentResults/cachedb_species_wsp

clidentseq \
--method=QC \
--blastdb=OverlappedPairedEnd_woSTD_10_ClaidentResults/cachedb_species_wsp \
--numthreads=$THREADS \
OverlappedPairedEnd_woSTD_09_ClusteredSequences/clustered.fasta \
OverlappedPairedEnd_woSTD_10_ClaidentResults/neighborhoods_qc_species_wsp.txt

classigntax \
--taxdb=animals_mt_species_wsp \
OverlappedPairedEnd_woSTD_10_ClaidentResults/neighborhoods_qc_species_wsp.txt \
OverlappedPairedEnd_woSTD_10_ClaidentResults/taxonomy_qc_species_wsp.tsv

# Assign taxonomy based on (95%-)3-NN method using animals_mt_species_wsp
clidentseq \
--method=3,95% \
--blastdb=OverlappedPairedEnd_woSTD_10_ClaidentResults/cachedb_species_wsp \
--numthreads=$THREADS \
OverlappedPairedEnd_woSTD_09_ClusteredSequences/clustered.fasta \
OverlappedPairedEnd_woSTD_10_ClaidentResults/neighborhoods_3nn_species_wsp.txt

classigntax \
--taxdb=animals_mt_species_wsp \
--minnsupporter=3 \
OverlappedPairedEnd_woSTD_10_ClaidentResults/neighborhoods_3nn_species_wsp.txt \
OverlappedPairedEnd_woSTD_10_ClaidentResults/taxonomy_3nn_species_wsp.tsv

# Assign taxonomy based on QCauto method using animals_mt_species_wosp
clmakecachedb \
--blastdb=animals_mt_species_wosp \
--numthreads=$THREADS \
OverlappedPairedEnd_woSTD_09_ClusteredSequences/clustered.fasta \
OverlappedPairedEnd_woSTD_10_ClaidentResults/cachedb_species_wosp

clidentseq \
--method=QC \
--blastdb=OverlappedPairedEnd_woSTD_10_ClaidentResults/cachedb_species_wosp \
--numthreads=$THREADS \
OverlappedPairedEnd_woSTD_09_ClusteredSequences/clustered.fasta \
OverlappedPairedEnd_woSTD_10_ClaidentResults/neighborhoods_qc_species_wosp.txt

classigntax \
--taxdb=animals_mt_species_wosp \
OverlappedPairedEnd_woSTD_10_ClaidentResults/neighborhoods_qc_species_wosp.txt \
OverlappedPairedEnd_woSTD_10_ClaidentResults/taxonomy_qc_species_wosp.tsv

# Assign taxonomy based on (95%-)3-NN method using animals_mt_species_wosp
clidentseq \
--method=3,95% \
--blastdb=OverlappedPairedEnd_woSTD_10_ClaidentResults/cachedb_species_wosp \
--numthreads=$THREADS \
OverlappedPairedEnd_woSTD_09_ClusteredSequences/clustered.fasta \
OverlappedPairedEnd_woSTD_10_ClaidentResults/neighborhoods_3nn_species_wosp.txt

classigntax \
--taxdb=animals_mt_species_wosp \
--minnsupporter=3 \
OverlappedPairedEnd_woSTD_10_ClaidentResults/neighborhoods_3nn_species_wosp.txt \
OverlappedPairedEnd_woSTD_10_ClaidentResults/taxonomy_3nn_species_wosp.tsv

# Merge 6 taxonomic assignment results
# Note that merge of QCauto results and (95%-)3-NN results has no effects in many cases because (95%-)3-NN results are always consistent to QCauto results excluding the case when there is no 95% or more similar reference sequences to the query.
# However, merge of results using different reference database is often useful.
clmergeassign \
--preferlower \
--priority=descend \
OverlappedPairedEnd_woSTD_10_ClaidentResults/taxonomy_qc_species_wosp.tsv \
OverlappedPairedEnd_woSTD_10_ClaidentResults/taxonomy_qc_species.tsv \
OverlappedPairedEnd_woSTD_10_ClaidentResults/taxonomy_qc_species_wsp.tsv \
OverlappedPairedEnd_woSTD_10_ClaidentResults/taxonomy_3nn_species_wosp.tsv \
OverlappedPairedEnd_woSTD_10_ClaidentResults/taxonomy_3nn_species.tsv \
OverlappedPairedEnd_woSTD_10_ClaidentResults/taxonomy_3nn_species_wsp.tsv \
OverlappedPairedEnd_woSTD_10_ClaidentResults/taxonomy_merged.tsv

# Fill blank cells of taxonomic assignment
clfillassign \
--fullfill=enable \
OverlappedPairedEnd_woSTD_10_ClaidentResults/taxonomy_merged.tsv \
OverlappedPairedEnd_woSTD_10_ClaidentResults/taxonomy_merged_filled.tsv

# Filter out non-fish OTUs
clfiltersum \
--taxfile=OverlappedPairedEnd_woSTD_10_ClaidentResults/taxonomy_merged_filled.tsv \
--includetaxa=class,Hyperoartia,class,Myxini,class,Chondrichthyes \
--includetaxa=superclass,Actinopterygii,order,Coelacanthiformes \
--includetaxa=subclass,Dipnomorpha \
OverlappedPairedEnd_woSTD_09_ClusteredSequences/clustered.tsv \
OverlappedPairedEnd_woSTD_10_ClaidentResults/sample_otu_matrix_fishes.tsv

# Extracting non-fish OTUs
clfiltersum \
--taxfile=OverlappedPairedEnd_woSTD_10_ClaidentResults/taxonomy_merged_filled.tsv \
--excludetaxa=class,Hyperoartia,class,Myxini,class,Chondrichthyes \
--excludetaxa=superclass,Actinopterygii,order,Coelacanthiformes \
--excludetaxa=subclass,Dipnomorpha \
OverlappedPairedEnd_woSTD_09_ClusteredSequences/clustered.tsv \
OverlappedPairedEnd_woSTD_10_ClaidentResults/sample_otu_matrix_nonfishes.tsv

# Plot word cloud
clplotwordcloud \
--taxfile=OverlappedPairedEnd_woSTD_10_ClaidentResults/taxonomy_merged_filled.tsv \
--targetrank=family,species \
--numthreads=$THREADS \
OverlappedPairedEnd_woSTD_10_ClaidentResults/sample_otu_matrix_fishes.tsv \
OverlappedPairedEnd_woSTD_10_ClaidentResults/wordcloud

# Make top-50 species community data matrix for barplot
clsumtaxa \
--tableformat=column \
--taxfile=OverlappedPairedEnd_woSTD_10_ClaidentResults/taxonomy_merged_filled.tsv \
--targetrank=species \
--topN=50 \
--numbering=enable \
OverlappedPairedEnd_woSTD_10_ClaidentResults/sample_otu_matrix_fishes.tsv \
OverlappedPairedEnd_woSTD_10_ClaidentResults/sample_top50species_nreads_fishes.tsv

# Make top-50 families community data matrix for barplot
clsumtaxa \
--tableformat=column \
--taxfile=OverlappedPairedEnd_woSTD_10_ClaidentResults/taxonomy_merged_filled.tsv \
--targetrank=family \
--topN=50 \
--numbering=enable \
OverlappedPairedEnd_woSTD_10_ClaidentResults/sample_otu_matrix_fishes.tsv \
OverlappedPairedEnd_woSTD_10_ClaidentResults/sample_top50family_nreads_fishes.tsv

# Make species-based community data matrix for heatmap
clsumtaxa \
--tableformat=column \
--taxfile=OverlappedPairedEnd_woSTD_10_ClaidentResults/taxonomy_merged_filled.tsv \
--targetrank=species \
--numbering=enable \
OverlappedPairedEnd_woSTD_10_ClaidentResults/sample_otu_matrix_fishes.tsv \
OverlappedPairedEnd_woSTD_10_ClaidentResults/sample_species_nreads_fishes.tsv

# Make family-based community data matrix for heatmap
clsumtaxa \
--tableformat=column \
--taxfile=OverlappedPairedEnd_woSTD_10_ClaidentResults/taxonomy_merged_filled.tsv \
--targetrank=family \
--numbering=enable \
OverlappedPairedEnd_woSTD_10_ClaidentResults/sample_otu_matrix_fishes.tsv \
OverlappedPairedEnd_woSTD_10_ClaidentResults/sample_family_nreads_fishes.tsv

# Coverage-based rarefaction
clrarefysum \
--minpcov=0.99 \
--minntotalseqsample=500 \
--nreplicate=4 \
--numthreads=$THREADS \
OverlappedPairedEnd_woSTD_10_ClaidentResults/sample_otu_matrix_fishes.tsv \
OverlappedPairedEnd_woSTD_10_ClaidentResults/sample_otu_matrix_fishes_rarefied

# Run R
Rscript runR_overlappedpairedend_woSTD.R

# Remove cachedb
rm -r OverlappedPairedEnd_woSTD_10_ClaidentResults/cachedb_species*
