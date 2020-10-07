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
./ecoPCR -d 12Sreferences -e 5 -l 150 -L 250 GTCGGTAAAACTCGTGCCAGC CATAGTGGGGTATCTAATCCCAGTTTG | perl convertecopcr2fasta.pl | perl -npe 's/^>.+; />/' > 12Sbarcodes.fasta
# Make fish community references for each sample
perl makecommref.pl 20 4 50 10 12Sbarcodes.fasta
# Download ART simulator
wget -c https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier2016.06.05linux64.tgz
# Extract ART simulator
tar -xzf artbinmountrainier2016.06.05linux64.tgz
# Generate simulated sequences
for f in `ls Sample??.fasta | grep -o -P '^[^\.]+'`
do ./art_bin_MountRainier/art_illumina --amplicon --seqSys MSv3 --in $f.fasta --len 144 --paired --noALN --fcov 1000 --out $f\_
done
for f in `ls Blank??.fasta | grep -o -P '^[^\.]+'`
do ./art_bin_MountRainier/art_illumina --amplicon --seqSys MSv3 --in $f.fasta --len 144 --paired --noALN --fcov 10 --out $f\_
done
# Make demultiplexed FASTQ
mkdir -p ../01_RawSequences
for f in `ls Sample??.fasta | grep -o -P '^[^\.]+'`
do perl addNNNNNN.pl $f\_1.fq | xz -c9 > ../01_RawSequences/$f\_R1.fastq.xz
perl addNNNNNN.pl $f\_2.fq | xz -c9 > ../01_RawSequences/$f\_R2.fastq.xz
done
for f in `ls Blank??.fasta | grep -o -P '^[^\.]+'`
do perl addNNNNNN.pl $f\_1.fq | xz -c9 > ../01_RawSequences/$f\_R1.fastq.xz
perl addNNNNNN.pl $f\_2.fq | xz -c9 > ../01_RawSequences/$f\_R2.fastq.xz
done
# Join FASTQ
xz -dc ../01_RawSequences/Sample??_R1.fastq.xz ../01_RawSequences/Blank??_R1.fastq.xz | xz -c9 > ../01_RawSequences/Undemultiplexed_R1.fastq.xz &
xz -dc ../01_RawSequences/Sample??_R2.fastq.xz ../01_RawSequences/Blank??_R2.fastq.xz | xz -c9 > ../01_RawSequences/Undemultiplexed_R2.fastq.xz &
# Make 2 index sequence files
perl makeindexfastq.pl ../index1.fasta Sample??_1.fq Blank??_1.fq | xz -c9 > ../01_RawSequences/Undemultiplexed_I1.fastq.xz &
perl makeindexfastq.pl ../index2.fasta Sample??_2.fq Blank??_2.fq | xz -c9 > ../01_RawSequences/Undemultiplexed_I2.fastq.xz &
wait
