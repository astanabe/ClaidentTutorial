# Install several packages
sudo apt install libgsl-dev libmath-random-mt-perl python2.7
# Install ecoPCR
tar -xzf ecoPCR-0.8.0.tar.gz
cd ecoPCR-0.8.0/src
make -j4
make ecoisundertaxon
cp ecofind ecogrep ecoisundertaxon ecoPCR ../..
cd ../tools
cp *.py ../..
cd ../..
# Download and compile Simera
wget -c -O Simera-master.tar.gz https://github.com/bnichols1979/Simera/archive/refs/heads/master.tar.gz
tar -xzf Simera-master.tar.gz
patch Simera-master/Simera.h < Simera.diff
cd Simera-master
make clean
make
cd ..
# Download ART simulator
wget -c https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier2016.06.05linux64.tgz
# Extract ART simulator
tar -xzf artbinmountrainier2016.06.05linux64.tgz

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

# Perform chimera formation simulation
ls Sample??.fasta | grep -o -P '^[^\.]+' | xargs -P 32 -I {} sh -c './Simera-master/Simera -i {}.fasta -o {}_simeraout -n 30 -s 10000 -l 0.00005 -c 10000 -f GTCGGTAAAACTCGTGCCAGC -r CATAGTGGGGTATCTAATCCCAGTTTG'
ls Blank??.fasta | grep -o -P '^[^\.]+' | xargs -P 32 -I {} sh -c './Simera-master/Simera -i {}.fasta -o {}_simeraout -n 30 -s 100 -l 0.00005 -c 10000 -f GTCGGTAAAACTCGTGCCAGC -r CATAGTGGGGTATCTAATCCCAGTTTG'
# Replace abundance format
ls Sample??.fasta Blank??.fasta | grep -o -P '^[^\.]+' | xargs -P 32 -I {} sh -c 'perl -i -npe "if(/^>/){s/_(\\d+)\\n/\\;size=\$1\\n/;}" {}_simeraout/samp_all_seqs.fa'
# Rereplicate sequences
ls Sample??.fasta Blank??.fasta | grep -o -P '^[^\.]+' | xargs -P 32 -I {} sh -c 'vsearch --fasta_width 0 --notrunclabels --rereplicate {}_simeraout/samp_all_seqs.fa --sizein --xsize --output {}_simeraout/samp_rerep.fasta'
# Generate simulated sequences
ls Sample??.fasta Blank??.fasta | grep -o -P '^[^\.]+' | xargs -P 32 -I {} sh -c './art_bin_MountRainier/art_illumina --amplicon --seqSys MSv1 --in {}_simeraout/samp_rerep.fasta --len 144 --paired --noALN --fcov 1 --out {}_'
# Make demultiplexed FASTQ
ls Sample??.fasta Blank??.fasta | grep -o -P '^[^\.]+' | xargs -P 32 -I {} sh -c 'perl addNNNNNN.pl {}_1.fq > {}_R1_001.fastq; perl addNNNNNN.pl {}_2.fq > {}_R2_001.fastq'
# Make output folder
mkdir -p ../01a_RawSequences_woSTD
# Make 2 index sequence files
perl makeindexfastq.pl ../index1.fasta Sample??_1.fq Blank??_1.fq | perl illuminaseqnamestyle.pl | gzip -c9 > ../01a_RawSequences_woSTD/Undetermined_S0_L001_I1_001.fastq.gz &
perl makeindexfastq.pl ../index2.fasta Sample??_2.fq Blank??_2.fq | perl illuminaseqnamestyle.pl | gzip -c9 > ../01a_RawSequences_woSTD/Undetermined_S0_L001_I2_001.fastq.gz &
# Make undemultiplexed files
sh -c 'cat Sample??_R1_001.fastq Blank??_R1_001.fastq | perl illuminaseqnamestyle.pl | gzip -c9 > ../01a_RawSequences_woSTD/Undetermined_S0_L001_R1_001.fastq.gz' &
sh -c 'cat Sample??_R2_001.fastq Blank??_R2_001.fastq | perl illuminaseqnamestyle.pl | gzip -c9 > ../01a_RawSequences_woSTD/Undetermined_S0_L001_R2_001.fastq.gz' &
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
../01a_RawSequences_woSTD/Undetermined_S0_L001_R1_001.fastq.gz \
../01a_RawSequences_woSTD/Undetermined_S0_L001_I1_001.fastq.gz \
../01a_RawSequences_woSTD/Undetermined_S0_L001_I2_001.fastq.gz \
../01a_RawSequences_woSTD/Undetermined_S0_L001_R2_001.fastq.gz \
TEMP
# Move demultiplexed files
cd TEMP
ls *.fastq | perl -nle '$fn=$_;m/TEMP__(.+)__TEMP\.(forward|reverse)\.fastq/;$sn=$1;$fr=$2;if($sn=~/^[ACGT]+\+[ACGT]+$/){unlink($fn)}else{if($fr eq "forward"){rename($fn,"$sn\_R1_001.fastq")}else{rename($fn,"$sn\_R2_001.fastq")}}'
ls *.fastq | xargs -L 1 -P 32 gzip -9
mv *.fastq.gz ../../01a_RawSequences_woSTD/
cd ..
rm -rf TEMP

# Make "wSTD" sequences
ls Sample??.fasta | grep -o -P '^[^\.]+' | xargs -P 32 -I {} sh -c 'cat {}.fasta MiFish_STD_*.fasta > {}_wSTD.fasta'
ls Blank??.fasta | grep -o -P '^[^\.]+' | xargs -P 32 -I {} sh -c 'cat {}.fasta > {}_wSTD.fasta; cat MiFish_STD_*.fasta | perl -npe "if(/^>/){s/\$/0/}" >> {}_wSTD.fasta'
# Perform chimera formation simulation
ls Sample??_wSTD.fasta | grep -o -P '^[^\.]+' | xargs -P 32 -I {} sh -c './Simera-master/Simera -i {}.fasta -o {}_simeraout -n 30 -s 20000 -l 0.00005 -c 10000 -f GTCGGTAAAACTCGTGCCAGC -r CATAGTGGGGTATCTAATCCCAGTTTG'
ls Blank??_wSTD.fasta | grep -o -P '^[^\.]+' | xargs -P 32 -I {} sh -c './Simera-master/Simera -i {}.fasta -o {}_simeraout -n 30 -s 10000 -l 0.00005 -c 10000 -f GTCGGTAAAACTCGTGCCAGC -r CATAGTGGGGTATCTAATCCCAGTTTG'
# Replace abundance format
ls Sample??_wSTD.fasta Blank??_wSTD.fasta | grep -o -P '^[^\.]+' | xargs -P 32 -I {} sh -c 'perl -i -npe "if(/^>/){s/_(\\d+)\\n/\\;size=\$1\\n/;}" {}_simeraout/samp_all_seqs.fa'
# Rereplicate sequences
ls Sample??_wSTD.fasta Blank??_wSTD.fasta | grep -o -P '^[^\.]+' | xargs -P 32 -I {} sh -c 'vsearch --fasta_width 0 --notrunclabels --rereplicate {}_simeraout/samp_all_seqs.fa --sizein --xsize --output {}_simeraout/samp_rerep.fasta'
# Generate simulated sequences
ls Sample??_wSTD.fasta Blank??_wSTD.fasta | grep -o -P '^[^\.]+' | xargs -P 32 -I {} sh -c './art_bin_MountRainier/art_illumina --amplicon --seqSys MSv1 --in {}_simeraout/samp_rerep.fasta --len 144 --paired --noALN --fcov 1 --out {}_'
# Make demultiplexed FASTQ
ls Sample??_wSTD.fasta Blank??_wSTD.fasta | grep -o -P '^[^\.]+' | xargs -P 32 -I {} sh -c 'perl addNNNNNN.pl {}_1.fq > {}_R1_001.fastq; perl addNNNNNN.pl {}_2.fq > {}_R2_001.fastq'
# Make output folder
mkdir -p ../01b_RawSequences_wSTD
# Make temporary index files
perl -npe 'if(/^>/){s/\n$/_wSTD\n/}' ../index1.fasta > index1_wSTD.fasta &
perl -npe 'if(/^>/){s/\n$/_wSTD\n/}' ../index2.fasta > index2_wSTD.fasta &
wait
# Make 2 index sequence files
perl makeindexfastq.pl index1_wSTD.fasta Sample??_wSTD_1.fq Blank??_wSTD_1.fq | perl illuminaseqnamestyle.pl | gzip -c9 > ../01b_RawSequences_wSTD/Undetermined_S0_L001_I1_001.fastq.gz &
perl makeindexfastq.pl index2_wSTD.fasta Sample??_wSTD_2.fq Blank??_wSTD_2.fq | perl illuminaseqnamestyle.pl | gzip -c9 > ../01b_RawSequences_wSTD/Undetermined_S0_L001_I2_001.fastq.gz &
# Make undemultiplexed files
sh -c 'cat Sample??_wSTD_R1_001.fastq Blank??_wSTD_R1_001.fastq | perl illuminaseqnamestyle.pl | gzip -c9 > ../01b_RawSequences_wSTD/Undetermined_S0_L001_R1_001.fastq.gz' &
sh -c 'cat Sample??_wSTD_R2_001.fastq Blank??_wSTD_R2_001.fastq | perl illuminaseqnamestyle.pl | gzip -c9 > ../01b_RawSequences_wSTD/Undetermined_S0_L001_R2_001.fastq.gz' &
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
../01b_RawSequences_wSTD/Undetermined_S0_L001_R1_001.fastq.gz \
../01b_RawSequences_wSTD/Undetermined_S0_L001_I1_001.fastq.gz \
../01b_RawSequences_wSTD/Undetermined_S0_L001_I2_001.fastq.gz \
../01b_RawSequences_wSTD/Undetermined_S0_L001_R2_001.fastq.gz \
TEMP
# Move demultiplexed files
cd TEMP
ls *.fastq | perl -nle '$fn=$_;m/TEMP__(.+)__TEMP\.(forward|reverse)\.fastq/;$sn=$1;$fr=$2;if($sn=~/^[ACGT]+\+[ACGT]+$/){unlink($fn)}else{if($fr eq "forward"){rename($fn,"$sn\_R1_001.fastq")}else{rename($fn,"$sn\_R2_001.fastq")}}'
ls *.fastq | xargs -L 1 -P 32 gzip -9
mv *.fastq.gz ../../01b_RawSequences_wSTD/
cd ..
rm -rf TEMP

# Delete temporary files
rm *.fastq *.fq Sample??.fasta Blank??.fasta Sample??_wSTD.fasta Blank??_wSTD.fasta
rm -r *_simeraout
