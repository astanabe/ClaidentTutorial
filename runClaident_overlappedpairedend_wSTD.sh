# Set number of processor cores used for computation
if test "$(uname)" = 'Darwin'; then
export THREADS=`sysctl -n hw.logicalcpu_max`
elif test "$(expr substr $(uname -s) 1 5)" = 'Linux'; then
export THREADS=`grep -c processor /proc/cpuinfo`
else
echo "Your platform ($(uname -a)) is not supported."
exit 1
fi

# Move previous analysis results
mkdir -p previous

if test -e OverlappedPairedEnd_wSTD_03_ConcatenatedSequences; then
mv \
OverlappedPairedEnd_wSTD_03_ConcatenatedSequences \
previous/
fi

if test -e OverlappedPairedEnd_wSTD_04_FilteredSequences; then
mv \
OverlappedPairedEnd_wSTD_04_FilteredSequences \
previous/
fi

if test -e OverlappedPairedEnd_wSTD_05_DenoisedSequences; then
mv \
OverlappedPairedEnd_wSTD_05_DenoisedSequences \
previous/
fi

if test -e OverlappedPairedEnd_wSTD_06_NonchimericSequences1; then
mv \
OverlappedPairedEnd_wSTD_06_NonchimericSequences1 \
previous/
fi

if test -e OverlappedPairedEnd_wSTD_07_STDClusteredSequences; then
mv \
OverlappedPairedEnd_wSTD_07_STDClusteredSequences \
previous/
fi

if test -e OverlappedPairedEnd_wSTD_08_NonchimericSequences2; then
mv \
OverlappedPairedEnd_wSTD_08_NonchimericSequences2 \
previous/
fi

if test -e OverlappedPairedEnd_wSTD_09_NonhoppedSequences; then
mv \
OverlappedPairedEnd_wSTD_09_NonhoppedSequences \
previous/
fi

if test -e OverlappedPairedEnd_wSTD_10_DecontaminatedSequences; then
mv \
OverlappedPairedEnd_wSTD_10_DecontaminatedSequences \
previous/
fi

if test -e OverlappedPairedEnd_wSTD_11_ClusteredSequences; then
mv \
OverlappedPairedEnd_wSTD_11_ClusteredSequences \
previous/
fi

if test -e OverlappedPairedEnd_wSTD_12_ClaidentResults; then
mv \
OverlappedPairedEnd_wSTD_12_ClaidentResults \
previous/
fi

if test -e OverlappedPairedEnd_wSTD_13_RAnalysisResults; then
mv OverlappedPairedEnd_wSTD_13_RAnalysisResults \
previous/
fi

# Demultiplex Type A (If you have undemultiplexed FASTQ files)
# --seqnamestyle=illumina should be used for real Illumina outputs.
if ! test -e PairedEnd_wSTD_02a_DemultiplexedSequences; then
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
01b_RawSequences_wSTD/Undetermined_S0_L001_R1_001.fastq.gz \
01b_RawSequences_wSTD/Undetermined_S0_L001_I1_001.fastq.gz \
01b_RawSequences_wSTD/Undetermined_S0_L001_I2_001.fastq.gz \
01b_RawSequences_wSTD/Undetermined_S0_L001_R2_001.fastq.gz \
PairedEnd_wSTD_02a_DemultiplexedSequences
fi

# Demultiplex Type B (If FASTQ files have been already demultiplexed)
# --seqnamestyle=illumina should be used for real Illumina outputs.
if ! test -e PairedEnd_wSTD_02b_DemultiplexedSequences; then
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
01b_RawSequences_wSTD/Sample??_R?_001.fastq.gz \
01b_RawSequences_wSTD/Blank??_R?_001.fastq.gz \
PairedEnd_wSTD_02b_DemultiplexedSequences
fi

# Compare Type A and B
rm -f PairedEnd_wSTD_TypeA.txt PairedEnd_wSTD_TypeB.txt

cd PairedEnd_wSTD_02a_DemultiplexedSequences

for f in *.fastq.xz
do echo $f >> ../PairedEnd_wSTD_TypeA.txt; xz -dc $f | grep -c -P '^\+\r?\n?$' >> ../PairedEnd_wSTD_TypeA.txt
done

cd ../PairedEnd_wSTD_02b_DemultiplexedSequences

for f in *.fastq.xz
do echo $f >> ../PairedEnd_wSTD_TypeB.txt; xz -dc $f | grep -c -P '^\+\r?\n?$' >> ../PairedEnd_wSTD_TypeB.txt
done

cd ..

diff -u PairedEnd_wSTD_TypeA.txt PairedEnd_wSTD_TypeB.txt

# Concatenate pairs
clconcatpairv \
--mode=ovl \
--compress=xz \
--numthreads=$THREADS \
PairedEnd_wSTD_02a_DemultiplexedSequences \
OverlappedPairedEnd_wSTD_03_ConcatenatedSequences

# Calculate FASTQ statistics
clcalcfastqstatv \
--mode=2 \
OverlappedPairedEnd_wSTD_03_ConcatenatedSequences \
OverlappedPairedEnd_wSTD_03_ConcatenatedSequences/fastq_eestats2.txt

# Apply filtering out low quality sequences
clfilterseqv \
--maxqual=41 \
--minlen=100 \
--maxlen=250 \
--maxnee=2.0 \
--maxnNs=0 \
--compress=xz \
--numthreads=$THREADS \
OverlappedPairedEnd_wSTD_03_ConcatenatedSequences \
OverlappedPairedEnd_wSTD_04_FilteredSequences

# Denoise using DADA2
cldenoiseseqd \
--pool=pseudo \
--numthreads=$THREADS \
OverlappedPairedEnd_wSTD_04_FilteredSequences \
OverlappedPairedEnd_wSTD_05_DenoisedSequences

# Remove chimeras using UCHIME3
clremovechimev \
--mode=denovo \
--uchimedenovo=3 \
OverlappedPairedEnd_wSTD_05_DenoisedSequences \
OverlappedPairedEnd_wSTD_06_NonchimericSequences1

# Cluster internal standard sequences to otus
clclusterstdv \
--standardseq=standard.fasta \
--minident=0.9 \
--numthreads=$THREADS \
OverlappedPairedEnd_wSTD_06_NonchimericSequences1 \
OverlappedPairedEnd_wSTD_07_STDClusteredSequences

# Remove chimeras using UCHIME3
clremovechimev \
--mode=ref \
--referencedb=cdu12s \
--addtoref=OverlappedPairedEnd_wSTD_07_STDClusteredSequences/stdvariations.fasta \
--numthreads=$THREADS \
OverlappedPairedEnd_wSTD_07_STDClusteredSequences \
OverlappedPairedEnd_wSTD_08_NonchimericSequences2

# Eliminate index-hopping
# This step cannot apply to TypeB demultiplexed sequences
clremovecontam \
--test=thompson \
--ignoresamplelist=blanklist.txt \
--index1file=index1.fasta \
--index2file=index2.fasta \
--numthreads=$THREADS \
OverlappedPairedEnd_wSTD_08_NonchimericSequences2 \
OverlappedPairedEnd_wSTD_09_NonhoppedSequences

# Eliminate contamination
# Note that this process can be applied to read-depth-normalized library data because of internal standard.
clremovecontam \
--test=thompson \
--blanklist=blanklist.txt \
--ignoreotuseq=standard.fasta \
--stdconctable=stdconctable.tsv \
--solutionvoltable=solutionvoltable.tsv \
--watervoltable=watervoltable.tsv \
--numthreads=$THREADS \
OverlappedPairedEnd_wSTD_09_NonhoppedSequences \
OverlappedPairedEnd_wSTD_10_DecontaminatedSequences

# Cluster remaining sequences
# Note that this step is meaningless on this data because additional clustering has no effect.
clclassseqv \
--minident=1.0 \
--strand=plus \
--numthreads=$THREADS \
--ignoreotuseq=standard.fasta \
OverlappedPairedEnd_wSTD_10_DecontaminatedSequences \
OverlappedPairedEnd_wSTD_11_ClusteredSequences

# Make final output folder
mkdir -p OverlappedPairedEnd_wSTD_12_ClaidentResults

# Assign taxonomy based on QCauto method using animals_mt_species
clmakecachedb \
--blastdb=animals_mt_species \
--ignoreotuseq=standard.fasta \
--numthreads=$THREADS \
OverlappedPairedEnd_wSTD_11_ClusteredSequences/clustered.fasta \
OverlappedPairedEnd_wSTD_12_ClaidentResults/cachedb_species

clidentseq \
--method=QC \
--blastdb=OverlappedPairedEnd_wSTD_12_ClaidentResults/cachedb_species \
--ignoreotuseq=standard.fasta \
--numthreads=$THREADS \
OverlappedPairedEnd_wSTD_11_ClusteredSequences/clustered.fasta \
OverlappedPairedEnd_wSTD_12_ClaidentResults/neighborhoods_qc_species.txt

classigntax \
--taxdb=animals_mt_species \
OverlappedPairedEnd_wSTD_12_ClaidentResults/neighborhoods_qc_species.txt \
OverlappedPairedEnd_wSTD_12_ClaidentResults/taxonomy_qc_species.tsv

# Assign taxonomy based on (95%-)3-NN method using animals_mt_species
clidentseq \
--method=3,95% \
--blastdb=OverlappedPairedEnd_wSTD_12_ClaidentResults/cachedb_species \
--ignoreotuseq=standard.fasta \
--numthreads=$THREADS \
OverlappedPairedEnd_wSTD_11_ClusteredSequences/clustered.fasta \
OverlappedPairedEnd_wSTD_12_ClaidentResults/neighborhoods_3nn_species.txt

classigntax \
--taxdb=animals_mt_species \
--minnsupporter=3 \
OverlappedPairedEnd_wSTD_12_ClaidentResults/neighborhoods_3nn_species.txt \
OverlappedPairedEnd_wSTD_12_ClaidentResults/taxonomy_3nn_species.tsv

# Assign taxonomy based on QCauto method using animals_mt_species_wsp
clmakecachedb \
--blastdb=animals_mt_species_wsp \
--ignoreotuseq=standard.fasta \
--numthreads=$THREADS \
OverlappedPairedEnd_wSTD_11_ClusteredSequences/clustered.fasta \
OverlappedPairedEnd_wSTD_12_ClaidentResults/cachedb_species_wsp

clidentseq \
--method=QC \
--blastdb=OverlappedPairedEnd_wSTD_12_ClaidentResults/cachedb_species_wsp \
--ignoreotuseq=standard.fasta \
--numthreads=$THREADS \
OverlappedPairedEnd_wSTD_11_ClusteredSequences/clustered.fasta \
OverlappedPairedEnd_wSTD_12_ClaidentResults/neighborhoods_qc_species_wsp.txt

classigntax \
--taxdb=animals_mt_species_wsp \
OverlappedPairedEnd_wSTD_12_ClaidentResults/neighborhoods_qc_species_wsp.txt \
OverlappedPairedEnd_wSTD_12_ClaidentResults/taxonomy_qc_species_wsp.tsv

# Assign taxonomy based on (95%-)3-NN method using animals_mt_species_wsp
clidentseq \
--method=3,95% \
--blastdb=OverlappedPairedEnd_wSTD_12_ClaidentResults/cachedb_species_wsp \
--ignoreotuseq=standard.fasta \
--numthreads=$THREADS \
OverlappedPairedEnd_wSTD_11_ClusteredSequences/clustered.fasta \
OverlappedPairedEnd_wSTD_12_ClaidentResults/neighborhoods_3nn_species_wsp.txt

classigntax \
--taxdb=animals_mt_species_wsp \
--minnsupporter=3 \
OverlappedPairedEnd_wSTD_12_ClaidentResults/neighborhoods_3nn_species_wsp.txt \
OverlappedPairedEnd_wSTD_12_ClaidentResults/taxonomy_3nn_species_wsp.tsv

# Assign taxonomy based on QCauto method using animals_mt_species_wosp
clmakecachedb \
--blastdb=animals_mt_species_wosp \
--ignoreotuseq=standard.fasta \
--numthreads=$THREADS \
OverlappedPairedEnd_wSTD_11_ClusteredSequences/clustered.fasta \
OverlappedPairedEnd_wSTD_12_ClaidentResults/cachedb_species_wosp

clidentseq \
--method=QC \
--blastdb=OverlappedPairedEnd_wSTD_12_ClaidentResults/cachedb_species_wosp \
--ignoreotuseq=standard.fasta \
--numthreads=$THREADS \
OverlappedPairedEnd_wSTD_11_ClusteredSequences/clustered.fasta \
OverlappedPairedEnd_wSTD_12_ClaidentResults/neighborhoods_qc_species_wosp.txt

classigntax \
--taxdb=animals_mt_species_wosp \
OverlappedPairedEnd_wSTD_12_ClaidentResults/neighborhoods_qc_species_wosp.txt \
OverlappedPairedEnd_wSTD_12_ClaidentResults/taxonomy_qc_species_wosp.tsv

# Assign taxonomy based on (95%-)3-NN method using animals_mt_species_wosp
clidentseq \
--method=3,95% \
--blastdb=OverlappedPairedEnd_wSTD_12_ClaidentResults/cachedb_species_wosp \
--ignoreotuseq=standard.fasta \
--numthreads=$THREADS \
OverlappedPairedEnd_wSTD_11_ClusteredSequences/clustered.fasta \
OverlappedPairedEnd_wSTD_12_ClaidentResults/neighborhoods_3nn_species_wosp.txt

classigntax \
--taxdb=animals_mt_species_wosp \
--minnsupporter=3 \
OverlappedPairedEnd_wSTD_12_ClaidentResults/neighborhoods_3nn_species_wosp.txt \
OverlappedPairedEnd_wSTD_12_ClaidentResults/taxonomy_3nn_species_wosp.tsv

# Merge 6 taxonomic assignment results
# Note that merge of QCauto results and (95%-)3-NN results has no effects in many cases because (95%-)3-NN results are always consistent to QCauto results excluding the case when there is no 95% or more similar reference sequences to the query.
# However, merge of results using different reference database is often useful.
clmergeassign \
--preferlower \
--priority=descend \
OverlappedPairedEnd_wSTD_12_ClaidentResults/taxonomy_qc_species_wosp.tsv \
OverlappedPairedEnd_wSTD_12_ClaidentResults/taxonomy_qc_species.tsv \
OverlappedPairedEnd_wSTD_12_ClaidentResults/taxonomy_qc_species_wsp.tsv \
OverlappedPairedEnd_wSTD_12_ClaidentResults/taxonomy_3nn_species_wosp.tsv \
OverlappedPairedEnd_wSTD_12_ClaidentResults/taxonomy_3nn_species.tsv \
OverlappedPairedEnd_wSTD_12_ClaidentResults/taxonomy_3nn_species_wsp.tsv \
OverlappedPairedEnd_wSTD_12_ClaidentResults/taxonomy_merged.tsv

# Fill blank cells of taxonomic assignment
clfillassign \
--fullfill=enable \
OverlappedPairedEnd_wSTD_12_ClaidentResults/taxonomy_merged.tsv \
OverlappedPairedEnd_wSTD_12_ClaidentResults/taxonomy_merged_filled.tsv

# Extract standard OTUs
clfiltersum \
--otuseq=standard.fasta \
OverlappedPairedEnd_wSTD_11_ClusteredSequences/clustered.tsv \
OverlappedPairedEnd_wSTD_12_ClaidentResults/sample_otu_matrix_standard.tsv

# Filter out non-fish OTUs
clfiltersum \
--negativeotuseq=standard.fasta \
--taxfile=OverlappedPairedEnd_wSTD_12_ClaidentResults/taxonomy_merged_filled.tsv \
--includetaxa=class,Hyperoartia,class,Myxini,class,Chondrichthyes \
--includetaxa=superclass,Actinopterygii,order,Coelacanthiformes \
--includetaxa=subclass,Dipnomorpha \
OverlappedPairedEnd_wSTD_11_ClusteredSequences/clustered.tsv \
OverlappedPairedEnd_wSTD_12_ClaidentResults/sample_otu_matrix_fishes.tsv

# Extract non-fish OTUs
clfiltersum \
--negativeotuseq=standard.fasta \
--taxfile=OverlappedPairedEnd_wSTD_12_ClaidentResults/taxonomy_merged_filled.tsv \
--excludetaxa=class,Hyperoartia,class,Myxini,class,Chondrichthyes \
--excludetaxa=superclass,Actinopterygii,order,Coelacanthiformes \
--excludetaxa=subclass,Dipnomorpha \
OverlappedPairedEnd_wSTD_11_ClusteredSequences/clustered.tsv \
OverlappedPairedEnd_wSTD_12_ClaidentResults/sample_otu_matrix_nonfishes.tsv

# Plot word cloud
clplotwordcloud \
--taxfile=OverlappedPairedEnd_wSTD_12_ClaidentResults/taxonomy_merged_filled.tsv \
--targetrank=family,species \
--numthreads=$THREADS \
OverlappedPairedEnd_wSTD_12_ClaidentResults/sample_otu_matrix_fishes.tsv \
OverlappedPairedEnd_wSTD_12_ClaidentResults/wordcloud

# Estimating DNA concentrations
clestimateconc \
--stdtable=OverlappedPairedEnd_wSTD_12_ClaidentResults/sample_otu_matrix_standard.tsv \
--stdconctable=stdconctable.tsv \
--solutionvoltable=solutionvoltable.tsv \
--watervoltable=watervoltable.tsv \
--numthreads=$THREADS \
OverlappedPairedEnd_wSTD_12_ClaidentResults/sample_otu_matrix_fishes.tsv \
OverlappedPairedEnd_wSTD_12_ClaidentResults/sample_otu_matrix_fishes_concentration.tsv

# Make top-50 species community data matrix for barplot
clsumtaxa \
--tableformat=column \
--taxfile=OverlappedPairedEnd_wSTD_12_ClaidentResults/taxonomy_merged_filled.tsv \
--targetrank=species \
--topN=50 \
--numbering=enable \
OverlappedPairedEnd_wSTD_12_ClaidentResults/sample_otu_matrix_fishes_concentration.tsv \
OverlappedPairedEnd_wSTD_12_ClaidentResults/sample_top50species_nreads_fishes_concentration.tsv

# Make top-50 families community data matrix for barplot
clsumtaxa \
--tableformat=column \
--taxfile=OverlappedPairedEnd_wSTD_12_ClaidentResults/taxonomy_merged_filled.tsv \
--targetrank=family \
--topN=50 \
--numbering=enable \
OverlappedPairedEnd_wSTD_12_ClaidentResults/sample_otu_matrix_fishes_concentration.tsv \
OverlappedPairedEnd_wSTD_12_ClaidentResults/sample_top50family_nreads_fishes_concentration.tsv

# Make species-based community data matrix for heatmap
clsumtaxa \
--tableformat=column \
--taxfile=OverlappedPairedEnd_wSTD_12_ClaidentResults/taxonomy_merged_filled.tsv \
--targetrank=species \
--numbering=enable \
OverlappedPairedEnd_wSTD_12_ClaidentResults/sample_otu_matrix_fishes_concentration.tsv \
OverlappedPairedEnd_wSTD_12_ClaidentResults/sample_species_nreads_fishes_concentration.tsv

# Make family-based community data matrix for heatmap
clsumtaxa \
--tableformat=column \
--taxfile=OverlappedPairedEnd_wSTD_12_ClaidentResults/taxonomy_merged_filled.tsv \
--targetrank=family \
--numbering=enable \
OverlappedPairedEnd_wSTD_12_ClaidentResults/sample_otu_matrix_fishes_concentration.tsv \
OverlappedPairedEnd_wSTD_12_ClaidentResults/sample_family_nreads_fishes_concentration.tsv

# Coverage-based rarefaction
clrarefysum \
--minpcov=0.99 \
--minntotalseqsample=500 \
--nreplicate=4 \
--numthreads=$THREADS \
OverlappedPairedEnd_wSTD_11_ClusteredSequences/clustered.tsv \
OverlappedPairedEnd_wSTD_12_ClaidentResults/sample_otu_matrix_all_rarefied

# Extract standard OTUs from rarefied tables
for n in `seq -w 1 4`
do clfiltersum \
--otuseq=standard.fasta \
OverlappedPairedEnd_wSTD_12_ClaidentResults/sample_otu_matrix_all_rarefied-r$n.tsv \
OverlappedPairedEnd_wSTD_12_ClaidentResults/sample_otu_matrix_standard_rarefied-r$n.tsv
done

# Extract fish OTUs
for n in `seq -w 1 4`
do clfiltersum \
--negativeotuseq=standard.fasta \
--taxfile=OverlappedPairedEnd_wSTD_12_ClaidentResults/taxonomy_merged_filled.tsv \
--includetaxa=class,Hyperoartia,class,Myxini,class,Chondrichthyes \
--includetaxa=superclass,Actinopterygii,order,Coelacanthiformes \
--includetaxa=subclass,Dipnomorpha \
OverlappedPairedEnd_wSTD_12_ClaidentResults/sample_otu_matrix_all_rarefied-r$n.tsv \
OverlappedPairedEnd_wSTD_12_ClaidentResults/sample_otu_matrix_fishes_rarefied-r$n.tsv
done

# Estimating DNA concentrations of rarefied tables
for n in `seq -w 1 4`
do clestimateconc \
--stdtable=OverlappedPairedEnd_wSTD_12_ClaidentResults/sample_otu_matrix_standard_rarefied-r$n.tsv \
--stdconctable=stdconctable.tsv \
--solutionvoltable=solutionvoltable.tsv \
--watervoltable=watervoltable.tsv \
--numthreads=$THREADS \
OverlappedPairedEnd_wSTD_12_ClaidentResults/sample_otu_matrix_fishes_rarefied-r$n.tsv \
OverlappedPairedEnd_wSTD_12_ClaidentResults/sample_otu_matrix_fishes_rarefied-r$n\_concentration.tsv
done

# Run R
Rscript runR_overlappedpairedend_wSTD.R

# Remove cachedb
rm -r OverlappedPairedEnd_wSTD_12_ClaidentResults/cachedb_species*
