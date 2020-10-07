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
--numthreads=8 \
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
--indedxname=$s \
--forwardprimerfile=forwardprimer.fasta \
--reverseprimerfile=reverseprimer.fasta \
--truncateN=enable \
--compress=xz \
--numthreads=8 \
--seqnamestyle=other \
01_RawSequences/$s\_R1.fastq.xz \
01_RawSequences/$s\_R2.fastq.xz \
02b_DemultiplexedSequences
done
# Denoising using DADA2

