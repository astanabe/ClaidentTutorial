# Set number of processor cores used for computation
if test "$(uname)" = 'Darwin'; then
export THREADS=`sysctl -n hw.logicalcpu_max`
elif test "$(expr substr $(uname -s) 1 5)" = 'Linux'; then
export THREADS=`grep -c processor /proc/cpuinfo`
else
echo "Your platform ($(uname -a)) is not supported."
exit 1
fi

# Backup previous analysis results
mkdir -p previous

if test -e NonoverlappedPairedEnd_woSTD_03_TrimmedSequences; then
mv \
NonoverlappedPairedEnd_woSTD_03_TrimmedSequences \
previous/
fi

if test -e NonoverlappedPairedEnd_woSTD_04_ConcatenatedSequences; then
mv \
NonoverlappedPairedEnd_woSTD_04_ConcatenatedSequences \
previous/
fi

if test -e NonoverlappedPairedEnd_woSTD_05_FilteredSequences; then
mv \
NonoverlappedPairedEnd_woSTD_05_FilteredSequences \
previous/
fi

if test -e NonoverlappedPairedEnd_woSTD_06_DenoisedSequences; then
mv \
NonoverlappedPairedEnd_woSTD_06_DenoisedSequences \
previous/
fi

if test -e NonoverlappedPairedEnd_woSTD_07_NonchimericSequences; then
mv \
NonoverlappedPairedEnd_woSTD_07_NonchimericSequences \
previous/
fi

if test -e NonoverlappedPairedEnd_woSTD_08_NonhoppedSequences; then
mv \
NonoverlappedPairedEnd_woSTD_08_NonhoppedSequences \
previous/
fi

if test -e NonoverlappedPairedEnd_woSTD_09_DecontaminatedSequences; then
mv \
NonoverlappedPairedEnd_woSTD_09_DecontaminatedSequences \
previous/
fi

if test -e NonoverlappedPairedEnd_woSTD_10_ClusteredSequences; then
mv \
NonoverlappedPairedEnd_woSTD_10_ClusteredSequences \
previous/
fi

if test -e NonoverlappedPairedEnd_woSTD_11_ClaidentResults; then
mv \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults \
previous/
fi

if test -e NonoverlappedPairedEnd_woSTD_12_RAnalysisResults; then
mv \
NonoverlappedPairedEnd_woSTD_12_RAnalysisResults \
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
rm -f PairedEnd_woSTD_TypeA.txt PairedEnd_woSTD_TypeB.txt

cd PairedEnd_woSTD_02a_DemultiplexedSequences

for f in *.fastq.xz
do echo $f >> ../PairedEnd_woSTD_TypeA.txt; xz -dc $f | grep -c -P '^\+\r?\n?$' >> ../PairedEnd_woSTD_TypeA.txt
done

cd ../PairedEnd_woSTD_02b_DemultiplexedSequences

for f in *.fastq.xz
do echo $f >> ../PairedEnd_woSTD_TypeB.txt; xz -dc $f | grep -c -P '^\+\r?\n?$' >> ../PairedEnd_woSTD_TypeB.txt
done

cd ..

diff -u PairedEnd_woSTD_TypeA.txt PairedEnd_woSTD_TypeB.txt

# Calculate FASTQ statistics
clcalcfastqstatv \
--mode=2 \
"PairedEnd_woSTD_02a_DemultiplexedSequences/*.forward.fastq.xz" \
PairedEnd_woSTD_02a_DemultiplexedSequences/forward_fastq_eestats2.txt

clcalcfastqstatv \
--mode=2 \
"PairedEnd_woSTD_02a_DemultiplexedSequences/*.reverse.fastq.xz" \
PairedEnd_woSTD_02a_DemultiplexedSequences/reverse_fastq_eestats2.txt

# Apply sequence trimming
# Trimmed length determined based on forward_fastq_eestats2.txt and reverse_fastq_eestats2.txt
clfilterseqv \
--minlen=120 \
--maxlen=120 \
--maxnNs=0 \
--compress=xz \
--numthreads=$THREADS \
"PairedEnd_woSTD_02a_DemultiplexedSequences/*.forward.fastq.xz" \
NonoverlappedPairedEnd_woSTD_03_TrimmedSequences

clfilterseqv \
--minlen=115 \
--maxlen=115 \
--maxnNs=0 \
--compress=xz \
--numthreads=$THREADS \
--append \
"PairedEnd_woSTD_02a_DemultiplexedSequences/*.reverse.fastq.xz" \
NonoverlappedPairedEnd_woSTD_03_TrimmedSequences

# Concatenate pairs as nonoverlapped paired-end
# Because this data is not nonoverlapped paired-end but overlapped paired-end, this procedure is inappropreate.
# This is just demonstration.
clconcatpairv \
--mode=non \
--compress=xz \
--numthreads=$THREADS \
NonoverlappedPairedEnd_woSTD_03_TrimmedSequences \
NonoverlappedPairedEnd_woSTD_04_ConcatenatedSequences

# Filfer out low quality sequences
clfilterseqv \
--maxqual=41 \
--maxnee=2.0 \
--maxnNs=0 \
--compress=xz \
--numthreads=$THREADS \
NonoverlappedPairedEnd_woSTD_04_ConcatenatedSequences \
NonoverlappedPairedEnd_woSTD_05_FilteredSequences

# Denoise using DADA2
cldenoiseseqd \
--pool=pseudo \
--numthreads=$THREADS \
NonoverlappedPairedEnd_woSTD_05_FilteredSequences \
NonoverlappedPairedEnd_woSTD_06_DenoisedSequences

# Remove chimeras using UCHIME3
# Do not apply reference-based chimera removal (Do not use "both" or "ref" for --mode).
clremovechimev \
--mode=denovo \
--uchimedenovo=3 \
--numthreads=$THREADS \
NonoverlappedPairedEnd_woSTD_06_DenoisedSequences \
NonoverlappedPairedEnd_woSTD_07_NonchimericSequences

# Eliminate index-hopping
# This step cannot apply to TypeB demultiplexed sequences
clremovecontam \
--test=thompson \
--ignoresamplelist=blanklist.txt \
--index1file=index1.fasta \
--index2file=index2.fasta \
--numthreads=$THREADS \
NonoverlappedPairedEnd_woSTD_07_NonchimericSequences \
NonoverlappedPairedEnd_woSTD_08_NonhoppedSequences

# Eliminate contamination
# Note that this process is incompatible with normalization of concentration/sequencing depth.
# Do not apply this process in such cases.
clremovecontam \
--test=thompson \
--blanklist=blanklist.txt \
--numthreads=$THREADS \
NonoverlappedPairedEnd_woSTD_08_NonhoppedSequences \
NonoverlappedPairedEnd_woSTD_09_DecontaminatedSequences

# Cluster remaining sequences
# Note that this step is meaningless on this data because additional clustering has no effect.
clclassseqv \
--minident=1.0 \
--strand=plus \
--paddinglen=16 \
--numthreads=$THREADS \
NonoverlappedPairedEnd_woSTD_09_DecontaminatedSequences \
NonoverlappedPairedEnd_woSTD_10_ClusteredSequences

# Divide forward and reverse sequences
cldivseq \
--query=ACGTACGTACGTACGT \
--border=both \
NonoverlappedPairedEnd_woSTD_10_ClusteredSequences/clustered.fasta \
NonoverlappedPairedEnd_woSTD_10_ClusteredSequences/reverse_revcomp.fasta \
NonoverlappedPairedEnd_woSTD_10_ClusteredSequences/forward.fasta

# Make final output folder
mkdir -p NonoverlappedPairedEnd_woSTD_11_ClaidentResults

# Assign taxonomy based on QCauto method using animals_12S_species and forward region
clmakecachedb \
--blastdb=animals_12S_species \
--numthreads=$THREADS \
NonoverlappedPairedEnd_woSTD_10_ClusteredSequences/forward.fasta \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/forwardcachedb_species

clidentseq \
--method=QC \
--blastdb=NonoverlappedPairedEnd_woSTD_11_ClaidentResults/forwardcachedb_species \
--numthreads=$THREADS \
NonoverlappedPairedEnd_woSTD_10_ClusteredSequences/forward.fasta \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/forwardneighborhoods_qc_species.txt

classigntax \
--taxdb=animals_12S_species \
--maxpopposer=0.05 \
--minsoratio=19 \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/forwardneighborhoods_qc_species.txt \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/forwardtaxonomy_qc_species.tsv

# Assign taxonomy based on (95%-)3-NN method using animals_12S_species and forward region
clidentseq \
--method=3,95% \
--blastdb=NonoverlappedPairedEnd_woSTD_11_ClaidentResults/forwardcachedb_species \
--numthreads=$THREADS \
NonoverlappedPairedEnd_woSTD_10_ClusteredSequences/forward.fasta \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/forwardneighborhoods_3nn_species.txt

classigntax \
--taxdb=animals_12S_species \
--minnsupporter=3 \
--maxpopposer=0.05 \
--minsoratio=19 \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/forwardneighborhoods_3nn_species.txt \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/forwardtaxonomy_3nn_species.tsv

# Assign taxonomy based on QCauto method using animals_12S_species_wsp and forward region
clmakecachedb \
--blastdb=animals_12S_species_wsp \
--numthreads=$THREADS \
NonoverlappedPairedEnd_woSTD_10_ClusteredSequences/forward.fasta \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/forwardcachedb_species_wsp

clidentseq \
--method=QC \
--blastdb=NonoverlappedPairedEnd_woSTD_11_ClaidentResults/forwardcachedb_species_wsp \
--numthreads=$THREADS \
NonoverlappedPairedEnd_woSTD_10_ClusteredSequences/forward.fasta \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/forwardneighborhoods_qc_species_wsp.txt

classigntax \
--taxdb=animals_12S_species_wsp \
--maxpopposer=0.05 \
--minsoratio=19 \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/forwardneighborhoods_qc_species_wsp.txt \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/forwardtaxonomy_qc_species_wsp.tsv

# Assign taxonomy based on (95%-)3-NN method using animals_12S_species_wsp and forward region
clidentseq \
--method=3,95% \
--blastdb=NonoverlappedPairedEnd_woSTD_11_ClaidentResults/forwardcachedb_species_wsp \
--numthreads=$THREADS \
NonoverlappedPairedEnd_woSTD_10_ClusteredSequences/forward.fasta \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/forwardneighborhoods_3nn_species_wsp.txt

classigntax \
--taxdb=animals_12S_species_wsp \
--minnsupporter=3 \
--maxpopposer=0.05 \
--minsoratio=19 \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/forwardneighborhoods_3nn_species_wsp.txt \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/forwardtaxonomy_3nn_species_wsp.tsv

# Assign taxonomy based on QCauto method using animals_12S_species_wosp and forward region
clmakecachedb \
--blastdb=animals_12S_species_wosp \
--numthreads=$THREADS \
NonoverlappedPairedEnd_woSTD_10_ClusteredSequences/forward.fasta \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/forwardcachedb_species_wosp

clidentseq \
--method=QC \
--blastdb=NonoverlappedPairedEnd_woSTD_11_ClaidentResults/forwardcachedb_species_wosp \
--numthreads=$THREADS \
NonoverlappedPairedEnd_woSTD_10_ClusteredSequences/forward.fasta \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/forwardneighborhoods_qc_species_wosp.txt

classigntax \
--taxdb=animals_12S_species_wosp \
--maxpopposer=0.05 \
--minsoratio=19 \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/forwardneighborhoods_qc_species_wosp.txt \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/forwardtaxonomy_qc_species_wosp.tsv

# Assign taxonomy based on (95%-)3-NN method using animals_12S_species_wosp and forward region
clidentseq \
--method=3,95% \
--blastdb=NonoverlappedPairedEnd_woSTD_11_ClaidentResults/forwardcachedb_species_wosp \
--numthreads=$THREADS \
NonoverlappedPairedEnd_woSTD_10_ClusteredSequences/forward.fasta \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/forwardneighborhoods_3nn_species_wosp.txt

classigntax \
--taxdb=animals_12S_species_wosp \
--minnsupporter=3 \
--maxpopposer=0.05 \
--minsoratio=19 \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/forwardneighborhoods_3nn_species_wosp.txt \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/forwardtaxonomy_3nn_species_wosp.tsv

# Assign taxonomy based on QCauto method using animals_12S_species and reverse region
clmakecachedb \
--blastdb=animals_12S_species \
--numthreads=$THREADS \
NonoverlappedPairedEnd_woSTD_10_ClusteredSequences/reverse_revcomp.fasta \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/reversecachedb_species

clidentseq \
--method=QC \
--blastdb=NonoverlappedPairedEnd_woSTD_11_ClaidentResults/reversecachedb_species \
--numthreads=$THREADS \
NonoverlappedPairedEnd_woSTD_10_ClusteredSequences/reverse_revcomp.fasta \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/reverseneighborhoods_qc_species.txt

classigntax \
--taxdb=animals_12S_species \
--maxpopposer=0.05 \
--minsoratio=19 \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/reverseneighborhoods_qc_species.txt \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/reversetaxonomy_qc_species.tsv

# Assign taxonomy based on (95%-)3-NN method using animals_12S_species and reverse region
clidentseq \
--method=3,95% \
--blastdb=NonoverlappedPairedEnd_woSTD_11_ClaidentResults/reversecachedb_species \
--numthreads=$THREADS \
NonoverlappedPairedEnd_woSTD_10_ClusteredSequences/reverse_revcomp.fasta \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/reverseneighborhoods_3nn_species.txt

classigntax \
--taxdb=animals_12S_species \
--minnsupporter=3 \
--maxpopposer=0.05 \
--minsoratio=19 \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/reverseneighborhoods_3nn_species.txt \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/reversetaxonomy_3nn_species.tsv

# Assign taxonomy based on QCauto method using animals_12S_species_wsp and reverse region
clmakecachedb \
--blastdb=animals_12S_species_wsp \
--numthreads=$THREADS \
NonoverlappedPairedEnd_woSTD_10_ClusteredSequences/reverse_revcomp.fasta \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/reversecachedb_species_wsp

clidentseq \
--method=QC \
--blastdb=NonoverlappedPairedEnd_woSTD_11_ClaidentResults/reversecachedb_species_wsp \
--numthreads=$THREADS \
NonoverlappedPairedEnd_woSTD_10_ClusteredSequences/reverse_revcomp.fasta \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/reverseneighborhoods_qc_species_wsp.txt

classigntax \
--taxdb=animals_12S_species_wsp \
--maxpopposer=0.05 \
--minsoratio=19 \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/reverseneighborhoods_qc_species_wsp.txt \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/reversetaxonomy_qc_species_wsp.tsv

# Assign taxonomy based on (95%-)3-NN method using animals_12S_species_wsp and reverse region
clidentseq \
--method=3,95% \
--blastdb=NonoverlappedPairedEnd_woSTD_11_ClaidentResults/reversecachedb_species_wsp \
--numthreads=$THREADS \
NonoverlappedPairedEnd_woSTD_10_ClusteredSequences/reverse_revcomp.fasta \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/reverseneighborhoods_3nn_species_wsp.txt

classigntax \
--taxdb=animals_12S_species_wsp \
--minnsupporter=3 \
--maxpopposer=0.05 \
--minsoratio=19 \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/reverseneighborhoods_3nn_species_wsp.txt \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/reversetaxonomy_3nn_species_wsp.tsv

# Assign taxonomy based on QCauto method using animals_12S_species_wosp and reverse region
clmakecachedb \
--blastdb=animals_12S_species_wosp \
--numthreads=$THREADS \
NonoverlappedPairedEnd_woSTD_10_ClusteredSequences/reverse_revcomp.fasta \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/reversecachedb_species_wosp

clidentseq \
--method=QC \
--blastdb=NonoverlappedPairedEnd_woSTD_11_ClaidentResults/reversecachedb_species_wosp \
--numthreads=$THREADS \
NonoverlappedPairedEnd_woSTD_10_ClusteredSequences/reverse_revcomp.fasta \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/reverseneighborhoods_qc_species_wosp.txt

classigntax \
--taxdb=animals_12S_species_wosp \
--maxpopposer=0.05 \
--minsoratio=19 \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/reverseneighborhoods_qc_species_wosp.txt \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/reversetaxonomy_qc_species_wosp.tsv

# Assign taxonomy based on (95%-)3-NN method using animals_12S_species_wosp and reverse region
clidentseq \
--method=3,95% \
--blastdb=NonoverlappedPairedEnd_woSTD_11_ClaidentResults/reversecachedb_species_wosp \
--numthreads=$THREADS \
NonoverlappedPairedEnd_woSTD_10_ClusteredSequences/reverse_revcomp.fasta \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/reverseneighborhoods_3nn_species_wosp.txt

classigntax \
--taxdb=animals_12S_species_wosp \
--minnsupporter=3 \
--maxpopposer=0.05 \
--minsoratio=19 \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/reverseneighborhoods_3nn_species_wosp.txt \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/reversetaxonomy_3nn_species_wosp.tsv

# Merge 6 taxonomic assignment results
# Note that merge of QCauto results and (95%-)3-NN results has no effects in many cases because (95%-)3-NN results are always consistent to QCauto results excluding the case when there is no 95% or more similar reference sequences to the query.
# However, merge of results using different reference database is often useful.
clmergeassign \
--preferlower \
--priority=descend \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/forwardtaxonomy_qc_species_wosp.tsv \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/forwardtaxonomy_qc_species.tsv \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/forwardtaxonomy_qc_species_wsp.tsv \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/forwardtaxonomy_3nn_species_wosp.tsv \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/forwardtaxonomy_3nn_species.tsv \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/forwardtaxonomy_3nn_species_wsp.tsv \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/forwardtaxonomy_merged.tsv

clmergeassign \
--preferlower \
--priority=descend \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/reversetaxonomy_qc_species_wosp.tsv \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/reversetaxonomy_qc_species.tsv \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/reversetaxonomy_qc_species_wsp.tsv \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/reversetaxonomy_3nn_species_wosp.tsv \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/reversetaxonomy_3nn_species.tsv \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/reversetaxonomy_3nn_species_wsp.tsv \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/reversetaxonomy_merged.tsv

clmergeassign \
--preferlower \
--priority=equal \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/forwardtaxonomy_merged.tsv \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/reversetaxonomy_merged.tsv \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/taxonomy_merged.tsv \

# Fill blank cells of taxonomic assignment
clfillassign \
--fullfill=enable \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/taxonomy_merged.tsv \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/taxonomy_merged_filled.tsv

# Filter out non-fish OTUs
clfiltersum \
--taxfile=NonoverlappedPairedEnd_woSTD_11_ClaidentResults/taxonomy_merged_filled.tsv \
--includetaxa=class,Hyperoartia,class,Myxini,class,Chondrichthyes \
--includetaxa=superclass,Actinopterygii,order,Coelacanthiformes \
--includetaxa=subclass,Dipnomorpha \
NonoverlappedPairedEnd_woSTD_10_ClusteredSequences/clustered.tsv \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/sample_otu_matrix_fishes.tsv

# Extract non-fish OTUs
clfiltersum \
--taxfile=NonoverlappedPairedEnd_woSTD_11_ClaidentResults/taxonomy_merged_filled.tsv \
--excludetaxa=class,Hyperoartia,class,Myxini,class,Chondrichthyes \
--excludetaxa=superclass,Actinopterygii,order,Coelacanthiformes \
--excludetaxa=subclass,Dipnomorpha \
NonoverlappedPairedEnd_woSTD_10_ClusteredSequences/clustered.tsv \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/sample_otu_matrix_nonfishes.tsv

# Plot word cloud
# Note that this command requires Google Chrome or Chromium browser
clplotwordcloud \
--taxfile=NonoverlappedPairedEnd_woSTD_11_ClaidentResults/taxonomy_merged_filled.tsv \
--targetrank=family,species \
--numthreads=$THREADS \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/sample_otu_matrix_fishes.tsv \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/wordcloud

# Make top-50 species community data matrix for barplot
clsumtaxa \
--tableformat=column \
--taxfile=NonoverlappedPairedEnd_woSTD_11_ClaidentResults/taxonomy_merged_filled.tsv \
--targetrank=species \
--topN=50 \
--numbering=enable \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/sample_otu_matrix_fishes.tsv \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/sample_top50species_nreads_fishes.tsv

# Make top-50 families community data matrix for barplot
clsumtaxa \
--tableformat=column \
--taxfile=NonoverlappedPairedEnd_woSTD_11_ClaidentResults/taxonomy_merged_filled.tsv \
--targetrank=family \
--topN=50 \
--numbering=enable \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/sample_otu_matrix_fishes.tsv \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/sample_top50family_nreads_fishes.tsv

# Make species-based community data matrix for heatmap
clsumtaxa \
--tableformat=column \
--taxfile=NonoverlappedPairedEnd_woSTD_11_ClaidentResults/taxonomy_merged_filled.tsv \
--targetrank=species \
--numbering=enable \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/sample_otu_matrix_fishes.tsv \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/sample_species_nreads_fishes.tsv

# Make family-based community data matrix for heatmap
clsumtaxa \
--tableformat=column \
--taxfile=NonoverlappedPairedEnd_woSTD_11_ClaidentResults/taxonomy_merged_filled.tsv \
--targetrank=family \
--numbering=enable \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/sample_otu_matrix_fishes.tsv \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/sample_family_nreads_fishes.tsv

# Coverage-based rarefaction
clrarefysum \
--minpcov=0.99 \
--minntotalseqsample=500 \
--nreplicate=4 \
--numthreads=$THREADS \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/sample_otu_matrix_fishes.tsv \
NonoverlappedPairedEnd_woSTD_11_ClaidentResults/sample_otu_matrix_fishes_rarefied

# Remove cachedb
rm -r NonoverlappedPairedEnd_woSTD_11_ClaidentResults/*cachedb_species*
