# Download taxonomy dump file and extract
mkdir -p taxonomy
cd taxonomy
wget -c https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
wget -c https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz.md5
md5sum -c taxdump.tar.gz.md5 || exit $?
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
python2.7 ./ecoPCRFormat.py -f -t ./taxonomy -n 12Sreferences 12Sreferences.ecoPCR
# Run ecoPCR and convert to amplicon FASTA file
./ecoPCR -d 12Sreferences -e 5 -l 150 -L 250 GTCGGTAAAACTCGTGCCAGC CATAGTGGGGTATCTAATCCCAGTTTG | perl convertecopcr2fasta.pl -fp=GTCGGTAAAACTCGTGCCAGC -rp=CATAGTGGGGTATCTAATCCCAGTTTG | perl -npe 's/^>.+; />/' > 12Sbarcodes.fasta
# Select representative sequences
vsearch --fasta_width 0 --notrunclabels --threads 32 --cluster_fast 12Sbarcodes.fasta --qmask none --id 0.9 --strand plus --centroids 12Sbarcodes_representatives.fasta
# Make fish community references for each sample
perl makecommref.pl 20 4 50 10 12Sbarcodes_representatives.fasta
# Download ART simulator
wget -c https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier2016.06.05linux64.tgz
# Extract ART simulator
tar -xzf artbinmountrainier2016.06.05linux64.tgz
# Generate simulated sequences
ls Sample??.fasta | grep -o -P '^[^\.]+' | xargs -P 32 -I {} sh -c './art_bin_MountRainier/art_illumina --amplicon --seqSys MSv1 --in {}.fasta --len 144 --paired --noALN --fcov 250 --out {}_'
ls Blank??.fasta | grep -o -P '^[^\.]+' | xargs -P 32 -I {} sh -c './art_bin_MountRainier/art_illumina --amplicon --seqSys MSv1 --in {}.fasta --len 144 --paired --noALN --fcov 25 --out {}_'
# Make demultiplexed FASTQ
ls Sample??.fasta Blank??.fasta | grep -o -P '^[^\.]+' | xargs -P 32 -I {} sh -c 'perl addNNNNNN.pl {}_1.fq > {}_R1_001.fastq; perl addNNNNNN.pl {}_2.fq > {}_R2_001.fastq'
# Make output folder
mkdir -p ../01a_RawSequences_woSTD
# Make 2 index sequence files
perl makeindexfastq.pl ../index1.fasta Sample??_1.fq Blank??_1.fq | perl illuminaseqnamestyle.pl | gzip -c9 > ../01_RawSequences_woSTD/Undetermined_S0_L001_I1_001.fastq.gz &
perl makeindexfastq.pl ../index2.fasta Sample??_2.fq Blank??_2.fq | perl illuminaseqnamestyle.pl | gzip -c9 > ../01_RawSequences_woSTD/Undetermined_S0_L001_I2_001.fastq.gz &
# Make undemultiplexed files
sh -c 'cat Sample??_R1_001.fastq Blank??_R1_001.fastq | perl illuminaseqnamestyle.pl | gzip -c9 > ../01_RawSequences_woSTD/Undetermined_S0_L001_R1_001.fastq.gz' &
sh -c 'cat Sample??_R2_001.fastq Blank??_R2_001.fastq | perl illuminaseqnamestyle.pl | gzip -c9 > ../01_RawSequences_woSTD/Undetermined_S0_L001_R2_001.fastq.gz' &
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
--seqnamestyle=illumina \
../01_RawSequences_woSTD/Undetermined_S0_L001_R1_001.fastq.gz \
../01_RawSequences_woSTD/Undetermined_S0_L001_I1_001.fastq.gz \
../01_RawSequences_woSTD/Undetermined_S0_L001_I2_001.fastq.gz \
../01_RawSequences_woSTD/Undetermined_S0_L001_R2_001.fastq.gz \
TEMP
# Move demultiplexed files
cd TEMP
ls *.fastq | perl -nle '$fn=$_;m/TEMP__(.+)__TEMP\.(forward|reverse)\.fastq/;$sn=$1;$fr=$2;if($sn=~/^[ACGT]+\+[ACGT]+$/){unlink($fn)}else{if($fr eq "forward"){rename($fn,"$sn\_R1_001.fastq")}else{rename($fn,"$sn\_R2_001.fastq")}}'
ls *.fastq | xargs -L 1 -P 32 gzip -9
mv *.fastq.gz ../../01_RawSequences_woSTD/
cd ..
rm -rf TEMP

# Generate simulated standard sequences
./art_bin_MountRainier/art_illumina --amplicon --seqSys MSv1 --in MiFish_STD_01.fasta --len 144 --paired --noALN --fcov 100 --out MiFish_STD_01_
./art_bin_MountRainier/art_illumina --amplicon --seqSys MSv1 --in MiFish_STD_02.fasta --len 144 --paired --noALN --fcov 200 --out MiFish_STD_02_
./art_bin_MountRainier/art_illumina --amplicon --seqSys MSv1 --in MiFish_STD_03.fasta --len 144 --paired --noALN --fcov 400 --out MiFish_STD_03_
./art_bin_MountRainier/art_illumina --amplicon --seqSys MSv1 --in MiFish_STD_04-2.fasta --len 144 --paired --noALN --fcov 800 --out MiFish_STD_04-2_
ls Sample??.fasta Blank??.fasta | grep -o -P '^[^\.]+' | xargs -P 32 -I {} sh -c 'cat MiFish_STD_*_1.fq >> {}_1.fq; cat MiFish_STD_*_2.fq >> {}_2.fq'
# Make demultiplexed FASTQ
ls Sample??.fasta Blank??.fasta | grep -o -P '^[^\.]+' | xargs -P 32 -I {} sh -c 'perl addNNNNNN.pl {}_1.fq > {}_R1_001.fastq; perl addNNNNNN.pl {}_2.fq > {}_R2_001.fastq'
# Make output folder
mkdir -p ../01b_RawSequences_wSTD
# Make 2 index sequence files
perl makeindexfastq.pl ../index1.fasta Sample??_1.fq Blank??_1.fq | perl illuminaseqnamestyle.pl | gzip -c9 > ../01_RawSequences_wSTD/Undetermined_S0_L001_I1_001.fastq.gz &
perl makeindexfastq.pl ../index2.fasta Sample??_2.fq Blank??_2.fq | perl illuminaseqnamestyle.pl | gzip -c9 > ../01_RawSequences_wSTD/Undetermined_S0_L001_I2_001.fastq.gz &
# Make undemultiplexed files
sh -c 'cat Sample??_R1_001.fastq Blank??_R1_001.fastq | perl illuminaseqnamestyle.pl | gzip -c9 > ../01_RawSequences_wSTD/Undetermined_S0_L001_R1_001.fastq.gz' &
sh -c 'cat Sample??_R2_001.fastq Blank??_R2_001.fastq | perl illuminaseqnamestyle.pl | gzip -c9 > ../01_RawSequences_wSTD/Undetermined_S0_L001_R2_001.fastq.gz' &
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
--seqnamestyle=illumina \
../01_RawSequences_wSTD/Undetermined_S0_L001_R1_001.fastq.gz \
../01_RawSequences_wSTD/Undetermined_S0_L001_I1_001.fastq.gz \
../01_RawSequences_wSTD/Undetermined_S0_L001_I2_001.fastq.gz \
../01_RawSequences_wSTD/Undetermined_S0_L001_R2_001.fastq.gz \
TEMP
# Move demultiplexed files
cd TEMP
ls *.fastq | perl -nle '$fn=$_;m/TEMP__(.+)__TEMP\.(forward|reverse)\.fastq/;$sn=$1;$fr=$2;if($sn=~/^[ACGT]+\+[ACGT]+$/){unlink($fn)}else{if($fr eq "forward"){rename($fn,"$sn\_R1_001.fastq")}else{rename($fn,"$sn\_R2_001.fastq")}}'
ls *.fastq | xargs -L 1 -P 32 gzip -9
mv *.fastq.gz ../../01_RawSequences_wSTD/
cd ..
rm -rf TEMP

# Delete temporary files
rm *.fastq *.fq Sample??.fasta Blank??.fasta
