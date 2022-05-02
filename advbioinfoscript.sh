#!/bin/bash #

##INSTALLING DEPENDENCIES

cd ~

wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x ./ Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh 
source ~/.bashrc 
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install samtools
conda install bwa
conda install freebayes
conda install picard
conda install bedtools
conda install trimmomatic
conda install fastqc
conda install vcflib

##ORGANISING DATA

mkdir ngs_pipeline
mkdir ngs_pipeline/dnaseq

cd ngs_pipeline/dnaseq
mkdir data results

ls -lF

cd ~/ngs_pipeline/dnaseq/data
mkdir untrimmed_fastq
mkdir trimmed_fastq

##DOWNLOADING INPUT FILES

wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R1.fastq.qz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R2.fastq.qz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/annotation.bed

mv NGS0001.R1.fastq.qz ~/ngs_pipeline/dnaseq/data/untrimmed_fastq
mv NGS0001.R2.fastq.qz ~/ngs_pipeline/dnaseq/data/untrimmed_fastq
mv annotation.bed ~/ngs_pipeline/dnaseq/data

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
mv hg19.fa.gz ~/ngs_pipeline/dnaseq/data/

##QUALITY ASSSESSMENT 

cd ~/ngs_pipeline/dnaseq/data/untrimmed_fastq
zcat NGS0001.R1.fastq.qz > NGS0001.R1.fastq
zcat NGS0001.R2.fastq.qz > NGS0001.R2.fastq

fastqc -t 4 *.fastq

mkdir ~/ngs_pipeline/dnaseq/results/fastqc_untrimmed_reads
mv *fastqc* ~/ngs_pipeline/dnaseq/results/fastqc_untrimmed_reads/
ls -lh ~/ngs_pipeline/dnaseq/results/fastqc_untrimmed_reads/

##moved .html files onto personal computer through FileZilla to view them and check for quality
cd ~/ngs_pipeline/dnaseq/results
rm -r fastqc_untrimmed_reads

##TRIMMING DATA

cd ~/ngs_pipeline/dnaseq/data/untrimmed_fastq

trimmomatic PE  \
-threads 4 \
-phred33 \
/home/ubuntu/ngs_pipeline/dnaseq/data/untrimmed_fastq/NGS0001.R1.fastq /home/ubuntu/ngs_pipeline/dnaseq/data/untrimmed_fastq/NGS0001.R2.fastq \
-baseout /home/ubuntu/ngs_pipeline/dnaseq/data/trimmed_fastq/NGS0001_trimmed_R \
ILLUMINACLIP:/home/ubuntu/miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:30:10 \
TRAILING:25 MINLEN:50

rm NGS0001.R1.fastq 
rm NGS0001.R2.fastq

##QUALITY ASSESSMENT OF TRIMMED DATA

cd ~/ngs_pipeline/dnaseq/data/trimmed_fastq
fastqc -t 4 NGS0001_trimmed_R_1P
fastqc -t 4 NGS0001_trimmed_R_2P

mkdir ~/ngs_pipeline/dnaseq/results/fastqc_trimmed_reads
mv *fastqc* ~/ngs_pipeline/dnaseq/results/fastqc_trimmed_reads/

##moved .html files onto personal computer through FileZilla to view them and check for quality
cd ~/ngs_pipeline/dnaseq/results
rm -r fastqc_trimmed_reads

##INDEXING REFERENCE FILE AND ALIGNING TRIMMED READS TO REFERENCE 

mkdir -p ~/ngs_pipeline/dnaseq/data/reference
mv ~/ngs_pipeline/dnaseq/data/hg19.fa.gz ~/ngs_pipeline/dnaseq/data/reference/
bwa index ~/ngs_pipeline/dnaseq/data/reference/hg19.fa.gz

mkdir ~/ngs_pipeline/dnaseq/data/aligned_data
bwa mem -t 4 -v 1 -R '@RG\tID:HWI-D0011.50.H7AP8ADXX.1.NGS0001\tSM:WES01\tPL:ILLUMINA\tLB:nextera-ngs0001-blood\tDT:2017-02-23\tPU:HWI-D00119' -I 250,50  ~/ngs_pipeline/dnaseq/data/reference/hg19.fa.gz ~/ngs_pipeline/dnaseq/data/trimmed_fastq/NGS0001_trimmed_R_1P ~/ngs_pipeline/dnaseq/data/trimmed_fastq/NGS0001_trimmed_R_2P > ~/ngs_pipeline/dnaseq/data/aligned_data/NGS0001.sam

rm -r ~/ngs_pipeline/dnaseq/data/untrimmed_fastq

cd ~/ngs_pipeline/dnaseq/data/aligned_data
samtools view -h -b NGS0001.sam > NGS0001.bam
samtools sort NGS0001.bam > NGS0001_sorted.bam
samtools index NGS0001_sorted.bam

rm -r ~/ngs_pipeline/dnaseq/data/trimmed_fastq
rm NGS0001.sam
rm NGS0001.bam

##DUPLICATE MARKING AND FILTERING

picard MarkDuplicates I=NGS0001_sorted.bam O=NGS0001_sorted_marked.bam M=marked_dup_metrics.txt
samtools index NGS0001_sorted_marked.bam

samtools view -F 1796  -q 20 -o NGS0001_sorted_filtered.bam NGS0001_sorted_marked.bam
samtools index NGS0001_sorted_filtered.bam
rm NGS0001_sorted.bam
rm NGS0001_sorted_marked.bam
rm NGS0001_sorted.bam.bai
rm NGS0001_sorted_marked.bam.bai

##ALIGNMENT STATS

samtools flagstat NGS0001_sorted_filtered.bam
samtools idxstats NGS0001_sorted_filtered.bam
picard CollectInsertSizeMetrics I=NGS0001_sorted_filtered.bam O=insertsizemetrics.txt H=insertsizehistogram.pdf

samtools depth -a NGS0001_sorted_filtered.bam > depthofcovergae.txt

rm depthofcovergae.txt
rm insertsizehistogram.pdf
rm insertsizemetrics.txt
rm marked_dup_metrics.txt


##VARIANT CALLING

zcat ~/ngs_pipeline/dnaseq/data/reference/hg19.fa.gz > ~/ngs_pipeline/dnaseq/data/reference/hg19.fa 

samtools faidx ~/ngs_pipeline/dnaseq/data/reference/hg19.fa

freebayes --bam ~/ngs_pipeline/dnaseq/data/aligned_data/NGS0001_sorted_filtered.bam --fasta-reference ~/ngs_pipeline/dnaseq/data/reference/hg19.fa --vcf ~/ngs_pipeline/dnaseq/results/NGS0001.vcf

rm -r ~/ngs_pipeline/dnaseq/data/reference
rm NGS0001_sorted_filtered.bam
rm NGS0001_sorted_filtered.bam.bai

bgzip ~/ngs_pipeline/dnaseq/results/NGS0001.vcf

tabix -p vcf ~/ngs_pipeline/dnaseq/results/NGS0001.vcf.gz

##QUALITY FILTERING CALLED VARIANTS
vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" ~/ngs_pipeline/dnaseq/results/NGS0001.vcf.gz > ~/ngs_pipeline/dnaseq/results/NGS0001_filtered.vcf

bedtools intersect -header -wa -a ~/ngs_pipeline/dnaseq/results/NGS0001_filtered.vcf -b ~/ngs_pipeline/dnaseq/data/annotation.bed > ~/ngs_pipeline/dnaseq/results/NGS0001_filtered_R.vcf

bgzip ~/ngs_pipeline/dnaseq/results/NGS0001_filtered_R.vcf

tabix -p vcf ~/ngs_pipeline/dnaseq/results/NGS0001_filtered_R.vcf.gz

##ANNOTATING VARIANTS THROUGH ANNOVAR
cd ~
##assuming tar file has been imported into the home directory through filezilla
tar -zxvf annovar.latest.tar.gz
cd annovar
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar knownGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp31a_interpro humandb/

./convert2annovar.pl -format vcf4 ~/ngs_pipeline/dnaseq/results/NGS0001_filtered_R.vcf.gz > ~/ngs_pipeline/dnaseq/results/NGS0001_filtered_annotated.avinput


./table_annovar.pl ~/ngs_pipeline/dnaseq/results/NGS0001_filtered_annotated.avinput humandb/ -buildver hg19 -out ~/ngs_pipeline/dnaseq/results/NGS0001_filtered_annotated -remove -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro -operation g,g,f,f,f -otherinfo -nastring . -csvout

##ANNOTATING VARIANTS THROUGH SNPEFF
#Download latest version
#wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip

#Unzip file
#unzip snpEff_latest_core.zip

#sudo add-apt-repository ppa:openjdk-r/ppa
#sudo apt-get update
#sudo apt install openjdk-11-jdk

#java -jar snpEff.jar download GRCh38.76
#java -Xmx8g -jar snpEff.jar GRCh38.76 ~/ngs_pipeline/dnaseq/results/NGS0001_filtered_R.vcf.gz ~/ngs_pipeline/dnaseq/results/NGS0001_filtered_annotated.ann.vcf 

##VARIANT PRIORITISATION
#cd ~/ngs_pipeline/dnaseq/results

