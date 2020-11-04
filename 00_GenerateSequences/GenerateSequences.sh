# Download taxonomy dump file and extract
mkdir -p taxonomy
cd taxonomy
wget -c ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzf taxdump.tar.gz
chmod 644 taxdump.tar.gz
rm taxdump.tar.gz
cd ..
# Install ecoPCR
tar -xzf ecoPCR-0.8.0.tar.gz
cd ecoPCR-0.8.0/src
make -j4
make ecoisundertaxon
cp ecofind ecogrep ecoisundertaxon ecoPCR ../..
cd ../tools
cp *.py ../..
cd ../..
# Add taxonomic information to 12S reference sequences
perl adddummytaxid2fasta.pl < 12Sreferences.fasta > 12Sreferences.ecoPCR
# Make ecoPCR database
./ecoPCRFormat.py -f -t ./taxonomy -n 12Sreferences 12Sreferences.ecoPCR
# Run ecoPCR and convert to amplicon FASTA file
./ecoPCR -d 12Sreferences -e 5 -l 150 -L 250 GTCGGTAAAACTCGTGCCAGC CATAGTGGGGTATCTAATCCCAGTTTG | perl convertecopcr2fasta.pl -fp=GTCGGTAAAACTCGTGCCAGC -rp=CATAGTGGGGTATCTAATCCCAGTTTG | perl -npe 's/^>.+; />/' > 12Sbarcodes.fasta
# Select representative sequences
vsearch --fasta_width 0 --notrunclabels --threads 32 --cluster_fast 12Sbarcodes.fasta --qmask none --id 0.95 --strand plus --centroids 12Sbarcodes_representatives.fasta
# Make fish community references for each sample
perl makecommref.pl 20 4 50 10 12Sbarcodes_representatives.fasta
# Download ART simulator
wget -c https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier2016.06.05linux64.tgz
# Extract ART simulator
tar -xzf artbinmountrainier2016.06.05linux64.tgz
# Generate simulated sequences
ls Sample??.fasta | grep -o -P '^[^\.]+' | xargs -L 1 -P 32 -I {} sh -c './art_bin_MountRainier/art_illumina --amplicon --seqSys MSv1 --in {}.fasta --len 144 --paired --noALN --fcov 500 --out {}_'
ls Blank??.fasta | grep -o -P '^[^\.]+' | xargs -L 1 -P 32 -I {} sh -c './art_bin_MountRainier/art_illumina --amplicon --seqSys MSv1 --in {}.fasta --len 144 --paired --noALN --fcov 50 --out {}_'
# Make demultiplexed FASTQ
ls Sample??.fasta Blank??.fasta | grep -o -P '^[^\.]+' | xargs -L 1 -P 32 -I {} sh -c 'perl addNNNNNN.pl {}_1.fq > {}_R1_001.fastq; perl addNNNNNN.pl {}_2.fq > {}_R2_001.fastq'
# Make output folder
mkdir -p ../01_RawSequences
# Make 2 index sequence files
perl makeindexfastq.pl ../index1.fasta Sample??_1.fq Blank??_1.fq | xz -c9e > ../01_RawSequences/Undemultiplexed_I1_001.fastq.xz &
perl makeindexfastq.pl ../index2.fasta Sample??_2.fq Blank??_2.fq | xz -c9e > ../01_RawSequences/Undemultiplexed_I2_001.fastq.xz &
# Make Undemultiplexed files
sh -c 'cat Sample??_R1_001.fastq Blank??_R1_001.fastq | xz -c9e > ../01_RawSequences/Undemultiplexed_R1_001.fastq.xz' &
sh -c 'cat Sample??_R2_001.fastq Blank??_R2_001.fastq | xz -c9e > ../01_RawSequences/Undemultiplexed_R2_001.fastq.xz' &
wait
# Make demultiplexed files
clsplitseq \
--runname=TEMP \
--primername=TEMP \
--index1file=../index1.fasta \
--index2file=../index2.fasta \
--minqualtag=0 \
--compress=disable \
--numthreads=32 \
--seqnamestyle=nochange \
../01_RawSequences/Undemultiplexed_R1_001.fastq.xz \
../01_RawSequences/Undemultiplexed_I1_001.fastq.xz \
../01_RawSequences/Undemultiplexed_I2_001.fastq.xz \
../01_RawSequences/Undemultiplexed_R2_001.fastq.xz \
TEMP
# Move demultiplexed files
cd TEMP
ls *.fastq | perl -nle '$fn=$_;m/TEMP__(.+)__TEMP\.(forward|reverse)\.fastq/;$sn=$1;$fr=$2;if($sn=~/^[ACGT]+\+[ACGT]+$/){unlink($fn)}else{if($fr eq "forward"){rename($fn,"$sn\_R1_001.fastq")}else{rename($fn,"$sn\_R2_001.fastq")}}'
ls *fastq | xargs -L 1 -P 32 xz -9e
mv *.fastq.xz ../../01_RawSequences/
cd ..
rm -rf TEMP
