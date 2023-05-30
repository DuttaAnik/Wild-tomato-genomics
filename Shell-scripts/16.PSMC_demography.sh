#!/bin/bash

# Anik Dutta

# Demographic analysis

# Use the pairwise Markovian sequential coalescent (PSMC) to estimate effective population size through time. Source: https://github.com/lh3/psmc

# The steps below follow [Al Ludington's PSMC workflow](https://github.com/a-lud/snakes-demographic-history) fairly closely.

# Choose samples and prepare inputs

# PSMC takes diploid sequences from an individual as input. These can be generated using mappings (i.e., bam files) and individual variant calls using samtools/bcftools.

# Individual inputs will be generated for a single male from each population: C069, C074, C005, C012, C025, C021, C036, C046, C052, C060, C079, C081, C093, C101, C104, C109, C122, C126

# Generating consensus variants per individual randomly chosen from each population for PSMC analysis.

bcftools mpileup -C 50 -q 30 -Q 25 -Ou -f S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C069.bam | bcftools call -c | bcftools filter --SnpGap 5 -i "DP>=5 & DP<=50" | bcftools view --exclude-types indels | bcftools sort --temp-dir /gxfs_work1/cau/suaph275/Demography/C069/temp_C069 -Oz -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C069.vcf.gz

bcftools mpileup -C 50 -q 10 -Ou -f S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C069.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=10 & DP<=60" | bcftools view --exclude-types indels | bcftools sort --temp-dir /gxfs_work1/cau/suaph275/Demography/C069/temp_C069 -Oz -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C069.vcf.gz
bcftools mpileup -C 50 -q 10 -Ou -f S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C074.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=10 & DP<=60" | bcftools view --exclude-types indels | bcftools sort --temp-dir /gxfs_work1/cau/suaph275/Demography/C069/temp_C069 -Oz -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C074.vcf.gz

bcftools mpileup -C 50 -q 10 -Ou -f S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C005.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=10 & DP<=60" | bcftools view --exclude-types indels | bcftools sort --temp-dir /gxfs_work1/cau/suaph275/Demography/C005/temp_C005 -Oz -o /gxfs_work1/cau/suaph275/Demography/C005/C_2931_C005.vcf.gz
bcftools mpileup -C 50 -q 10 -Ou -f S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C012.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=10 & DP<=60" | bcftools view --exclude-types indels | bcftools sort --temp-dir /gxfs_work1/cau/suaph275/Demography/C005/temp_C005 -Oz -o /gxfs_work1/cau/suaph275/Demography/C005/C_2931_C012.vcf.gz

bcftools mpileup -C 50 -q 10 -Ou -f S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C025.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=10 & DP<=60" | bcftools view --exclude-types indels | bcftools sort --temp-dir /gxfs_work1/cau/suaph275/Demography/C025/temp_C025 -Oz -o /gxfs_work1/cau/suaph275/Demography/C025/C_2747_C025.vcf.gz
bcftools mpileup -C 50 -q 10 -Ou -f S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C021.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=10 & DP<=60" | bcftools view --exclude-types indels | bcftools sort --temp-dir /gxfs_work1/cau/suaph275/Demography/C025/temp_C025 -Oz -o /gxfs_work1/cau/suaph275/Demography/C025/C_2747_C021.vcf.gz

bcftools mpileup -C 50 -q 10 -Ou -f S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C036.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=10 & DP<=60" | bcftools view --exclude-types indels | bcftools sort --temp-dir /gxfs_work1/cau/suaph275/Demography/C036/temp_C036 -Oz -o /gxfs_work1/cau/suaph275/Demography/C036/SC_2932_C036.vcf.gz
bcftools mpileup -C 50 -q 10 -Ou -f S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C046.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=10 & DP<=60" | bcftools view --exclude-types indels | bcftools sort --temp-dir /gxfs_work1/cau/suaph275/Demography/C036/temp_C036 -Oz -o /gxfs_work1/cau/suaph275/Demography/C036/SC_2932_C046.vcf.gz

bcftools mpileup -C 50 -q 10 -Ou -f S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C052.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=10 & DP<=60" | bcftools view --exclude-types indels | bcftools sort --temp-dir /gxfs_work1/cau/suaph275/Demography/C052/temp_C052 -Oz -o /gxfs_work1/cau/suaph275/Demography/C052/C_1963_C052.vcf.gz
bcftools mpileup -C 50 -q 10 -Ou -f S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C060.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=10 & DP<=60" | bcftools view --exclude-types indels | bcftools sort --temp-dir /gxfs_work1/cau/suaph275/Demography/C052/temp_C052 -Oz -o /gxfs_work1/cau/suaph275/Demography/C052/C_1963_C060.vcf.gz

bcftools mpileup -C 50 -q 10 -Ou -f S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C079.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=10 & DP<=60" | bcftools view --exclude-types indels | bcftools sort --temp-dir /gxfs_work1/cau/suaph275/Demography/C079/temp_C079 -Oz -o /gxfs_work1/cau/suaph275/Demography/C079/SH_4330_C079.vcf.gz
bcftools mpileup -C 50 -q 10 -Ou -f S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C081.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=10 & DP<=60" | bcftools view --exclude-types indels | bcftools sort --temp-dir /gxfs_work1/cau/suaph275/Demography/C079/temp_C079 -Oz -o /gxfs_work1/cau/suaph275/Demography/C079/SH_4330_C081.vcf.gz

bcftools mpileup -C 50 -q 10 -Ou -f S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C093.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=10 & DP<=60" | bcftools view --exclude-types indels | bcftools sort --temp-dir /gxfs_work1/cau/suaph275/Demography/C093/temp_C093 -Oz -o /gxfs_work1/cau/suaph275/Demography/C093/C_3111_C093.vcf.gz
bcftools mpileup -C 50 -q 10 -Ou -f S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C101.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=10 & DP<=60" | bcftools view --exclude-types indels | bcftools sort --temp-dir /gxfs_work1/cau/suaph275/Demography/C093/temp_C093 -Oz -o /gxfs_work1/cau/suaph275/Demography/C093/C_3111_C101.vcf.gz

bcftools mpileup -C 50 -q 10 -Ou -f S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C109.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=10 & DP<=60" | bcftools view --exclude-types indels | bcftools sort --temp-dir /gxfs_work1/cau/suaph275/Demography/C109/temp_C109 -Oz -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C109.vcf.gz
bcftools mpileup -C 50 -q 10 -Ou -f S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C104.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=10 & DP<=60" | bcftools view --exclude-types indels | bcftools sort --temp-dir /gxfs_work1/cau/suaph275/Demography/C109/temp_C109 -Oz -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C104.vcf.gz

bcftools mpileup -C 50 -q 10 -Ou -f S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C122.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=10 & DP<=60" | bcftools view --exclude-types indels | bcftools sort --temp-dir /gxfs_work1/cau/suaph275/Demography/C122/temp_C122 -Oz -o /gxfs_work1/cau/suaph275/Demography/C122/SH_4117_C122.vcf.gz
bcftools mpileup -C 50 -q 10 -Ou -f S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C126.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=10 & DP<=60" | bcftools view --exclude-types indels | bcftools sort --temp-dir /gxfs_work1/cau/suaph275/Demography/C122/temp_C122 -Oz -o /gxfs_work1/cau/suaph275/Demography/C122/SH_4117_C126.vcf.gz

# Index VCF files

tabix -p vcf /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C069.vcf.gz
tabix -p vcf /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C074.vcf.gz
tabix -p vcf /gxfs_work1/cau/suaph275/Demography/C005/C_2931_C005.vcf.gz
tabix -p vcf /gxfs_work1/cau/suaph275/Demography/C005/C_2931_C012.vcf.gz
tabix -p vcf /gxfs_work1/cau/suaph275/Demography/C025/C_2747_C025.vcf.gz
tabix -p vcf /gxfs_work1/cau/suaph275/Demography/C025/C_2747_C021.vcf.gz
tabix -p vcf /gxfs_work1/cau/suaph275/Demography/C036/SC_2932_C036.vcf.gz
tabix -p vcf /gxfs_work1/cau/suaph275/Demography/C036/SC_2932_C046.vcf.gz
tabix -p vcf /gxfs_work1/cau/suaph275/Demography/C052/C_1963_C052.vcf.gz
tabix -p vcf /gxfs_work1/cau/suaph275/Demography/C052/C_1963_C060.vcf.gz
tabix -p vcf /gxfs_work1/cau/suaph275/Demography/C079/SH_4330_C079.vcf.gz
tabix -p vcf /gxfs_work1/cau/suaph275/Demography/C079/SH_4330_C081.vcf.gz
tabix -p vcf /gxfs_work1/cau/suaph275/Demography/C093/C_3111_C093.vcf.gz
tabix -p vcf /gxfs_work1/cau/suaph275/Demography/C093/C_3111_C101.vcf.gz
tabix -p vcf /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C109.vcf.gz
tabix -p vcf /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C104.vcf.gz
tabix -p vcf /gxfs_work1/cau/suaph275/Demography/C122/SH_4117_C122.vcf.gz
tabix -p vcf /gxfs_work1/cau/suaph275/Demography/C122/SH_4117_C126.vcf.gz

# Get diploid sequence

module load bcftools vcftools

bcftools view /gxfs_work1/cau/suaph275/Demography/C005/C_2931_C005.vcf.gz | vcfutils.pl vcf2fq | gzip > /gxfs_work1/cau/suaph275/Demography/C005/C_2931_C005.fastq.gz
bcftools view /gxfs_work1/cau/suaph275/Demography/C005/C_2931_C012.vcf.gz | vcfutils.pl vcf2fq | gzip > /gxfs_work1/cau/suaph275/Demography/C005/C_2931_C012.fastq.gz

bcftools view /gxfs_work1/cau/suaph275/Demography/C025/C_2747_C025.vcf.gz | vcfutils.pl vcf2fq | gzip > /gxfs_work1/cau/suaph275/Demography/C025/C_2747_C025.fastq.gz
bcftools view /gxfs_work1/cau/suaph275/Demography/C025/C_2747_C021.vcf.gz | vcfutils.pl vcf2fq | gzip > /gxfs_work1/cau/suaph275/Demography/C025/C_2747_C021.fastq.gz

bcftools view /gxfs_work1/cau/suaph275/Demography/C036/SC_2932_C036.vcf.gz | vcfutils.pl vcf2fq | gzip > /gxfs_work1/cau/suaph275/Demography/C036/SC_2932_C036.fastq.gz
bcftools view /gxfs_work1/cau/suaph275/Demography/C036/SC_2932_C046.vcf.gz | vcfutils.pl vcf2fq | gzip > /gxfs_work1/cau/suaph275/Demography/C036/SC_2932_C046.fastq.gz

bcftools view /gxfs_work1/cau/suaph275/Demography/C052/C_1963_C052.vcf.gz | vcfutils.pl vcf2fq | gzip > /gxfs_work1/cau/suaph275/Demography/C052/C_1963_C052.fastq.gz
bcftools view /gxfs_work1/cau/suaph275/Demography/C052/C_1963_C060.vcf.gz | vcfutils.pl vcf2fq | gzip > /gxfs_work1/cau/suaph275/Demography/C052/C_1963_C060.fastq.gz

bcftools view /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C069.vcf.gz | vcfutils.pl vcf2fq | gzip > /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C069.fastq.gz
bcftools view /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C074.vcf.gz | vcfutils.pl vcf2fq | gzip > /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C074.fastq.gz

bcftools view /gxfs_work1/cau/suaph275/Demography/C079/SH_4330_C079.vcf.gz | vcfutils.pl vcf2fq | gzip > /gxfs_work1/cau/suaph275/Demography/C079/SH_4330_C079.fastq.gz
bcftools view /gxfs_work1/cau/suaph275/Demography/C079/SH_4330_C081.vcf.gz | vcfutils.pl vcf2fq | gzip > /gxfs_work1/cau/suaph275/Demography/C079/SH_4330_C081.fastq.gz

bcftools view /gxfs_work1/cau/suaph275/Demography/C093/C_3111_C093.vcf.gz | vcfutils.pl vcf2fq | gzip > /gxfs_work1/cau/suaph275/Demography/C093/C_3111_C093.fastq.gz
bcftools view /gxfs_work1/cau/suaph275/Demography/C093/C_3111_C101.vcf.gz | vcfutils.pl vcf2fq | gzip > /gxfs_work1/cau/suaph275/Demography/C093/C_3111_C101.fastq.gz

bcftools view /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C109.vcf.gz | vcfutils.pl vcf2fq | gzip > /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C109.fastq.gz
bcftools view /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C104.vcf.gz | vcfutils.pl vcf2fq | gzip > /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C104.fastq.gz

bcftools view /gxfs_work1/cau/suaph275/Demography/C122/SH_4117_C122.vcf.gz | vcfutils.pl vcf2fq | gzip > /gxfs_work1/cau/suaph275/Demography/C122/SH_4117_C122.fastq.gz
bcftools view /gxfs_work1/cau/suaph275/Demography/C122/SH_4117_C126.vcf.gz | vcfutils.pl vcf2fq | gzip > /gxfs_work1/cau/suaph275/Demography/C122/SH_4117_C126.fastq.gz

# Convert into PSMCFA format

fq2psmcfa -q 20 /gxfs_work1/cau/suaph275/Demography/C005/C_2931_C005.fastq.gz > /gxfs_work1/cau/suaph275/Demography/C005/C_2931_C005.psmcfa
fq2psmcfa -q 20 /gxfs_work1/cau/suaph275/Demography/C005/C_2931_C012.fastq.gz > /gxfs_work1/cau/suaph275/Demography/C005/C_2931_C012.psmcfa

fq2psmcfa -q 20 /gxfs_work1/cau/suaph275/Demography/C025/C_2747_C025.fastq.gz > /gxfs_work1/cau/suaph275/Demography/C025/C_2747_C025.psmcfa
fq2psmcfa -q 20 /gxfs_work1/cau/suaph275/Demography/C025/C_2747_C021.fastq.gz > /gxfs_work1/cau/suaph275/Demography/C025/C_2747_C021.psmcfa

fq2psmcfa -q 20 /gxfs_work1/cau/suaph275/Demography/C036/SC_2932_C036.fastq.gz > /gxfs_work1/cau/suaph275/Demography/C036/SC_2932_C036.psmcfa
fq2psmcfa -q 20 /gxfs_work1/cau/suaph275/Demography/C036/SC_2932_C046.fastq.gz > /gxfs_work1/cau/suaph275/Demography/C036/SC_2932_C046.psmcfa

fq2psmcfa -q 20 /gxfs_work1/cau/suaph275/Demography/C052/C_1963_C052.fastq.gz > /gxfs_work1/cau/suaph275/Demography/C052/C_1963_C052.psmcfa
fq2psmcfa -q 20 /gxfs_work1/cau/suaph275/Demography/C052/C_1963_C060.fastq.gz > /gxfs_work1/cau/suaph275/Demography/C052/C_1963_C060.psmcfa

fq2psmcfa -q 20 /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C069.fastq.gz > /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C069.psmcfa
fq2psmcfa -q 20 /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C074.fastq.gz > /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C074.psmcfa

fq2psmcfa -q 20 /gxfs_work1/cau/suaph275/Demography/C079/SH_4330_C079.fastq.gz > /gxfs_work1/cau/suaph275/Demography/C079/SH_4330_C079.psmcfa
fq2psmcfa -q 20 /gxfs_work1/cau/suaph275/Demography/C079/SH_4330_C081.fastq.gz > /gxfs_work1/cau/suaph275/Demography/C079/SH_4330_C081.psmcfa

fq2psmcfa -q 20 /gxfs_work1/cau/suaph275/Demography/C093/C_3111_C093.fastq.gz > /gxfs_work1/cau/suaph275/Demography/C093/C_3111_C093.psmcfa
fq2psmcfa -q 20 /gxfs_work1/cau/suaph275/Demography/C093/C_3111_C101.fastq.gz > /gxfs_work1/cau/suaph275/Demography/C093/C_3111_C101.psmcfa

fq2psmcfa -q 20 /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C109.fastq.gz > /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C109.psmcfa
fq2psmcfa -q 20 /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C104.fastq.gz > /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C104.psmcfa

fq2psmcfa -q 20 /gxfs_work1/cau/suaph275/Demography/C122/SH_4117_C122.fastq.gz > /gxfs_work1/cau/suaph275/Demography/C122/SH_4117_C122.psmcfa
fq2psmcfa -q 20 /gxfs_work1/cau/suaph275/Demography/C122/SH_4117_C126.fastq.gz > /gxfs_work1/cau/suaph275/Demography/C122/SH_4117_C126.psmcfa

# Run PSMC using the time segment pattern "4+25*2+4+6" that is specified by the '-p' flag in PSMC analysis

#psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C005/C_2931_C005_1.psmc /gxfs_work1/cau/suaph275/Demography/C005/C_2931_C005.psmcfa This was the initial setting.
#psmc -N25 -t15 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C005/C_2931_C005_2.psmc /gxfs_work1/cau/suaph275/Demography/C005/C_2931_C005.psmcfa This was the initial setting.

psmc -N25 -t5 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C005/C_2931_C005_3.psmc /gxfs_work1/cau/suaph275/Demography/C005/C_2931_C005.psmcfa
psmc -N25 -t5 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C005/C_2931_C012_3.psmc /gxfs_work1/cau/suaph275/Demography/C005/C_2931_C012.psmcfa

psmc -N25 -t5 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C025/C_2747_C025_3.psmc /gxfs_work1/cau/suaph275/Demography/C025/C_2747_C025.psmcfa
psmc -N25 -t5 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C025/C_2747_C021_3.psmc /gxfs_work1/cau/suaph275/Demography/C025/C_2747_C021.psmcfa

psmc -N25 -t5 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C036/SC_2932_C036_3.psmc /gxfs_work1/cau/suaph275/Demography/C036/SC_2932_C036.psmcfa
psmc -N25 -t5 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C036/SC_2932_C046_3.psmc /gxfs_work1/cau/suaph275/Demography/C036/SC_2932_C046.psmcfa

psmc -N25 -t5 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C052/C_1963_C052_3.psmc /gxfs_work1/cau/suaph275/Demography/C052/C_1963_C052.psmcfa
psmc -N25 -t5 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C052/C_1963_C060_3.psmc /gxfs_work1/cau/suaph275/Demography/C052/C_1963_C060.psmcfa

psmc -N25 -t5 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C069_3.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C069.psmcfa
psmc -N25 -t5 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C074_3.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C074.psmcfa

psmc -N25 -t5 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C079/SH_4330_C079_3.psmc /gxfs_work1/cau/suaph275/Demography/C079/SH_4330_C079.psmcfa
psmc -N25 -t5 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C079/SH_4330_C081_3.psmc /gxfs_work1/cau/suaph275/Demography/C079/SH_4330_C081.psmcfa

psmc -N25 -t5 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C093/C_3111_C093_3.psmc /gxfs_work1/cau/suaph275/Demography/C093/C_3111_C093.psmcfa
psmc -N25 -t5 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C093/C_3111_C101_3.psmc /gxfs_work1/cau/suaph275/Demography/C093/C_3111_C101.psmcfa

psmc -N25 -t5 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C109_3.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C109.psmcfa
psmc -N25 -t5 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C104_3.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C104.psmcfa

psmc -N25 -t5 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C122/SH_4117_C122_3.psmc /gxfs_work1/cau/suaph275/Demography/C122/SH_4117_C122.psmcfa
psmc -N25 -t5 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C122/SH_4117_C126_3.psmc /gxfs_work1/cau/suaph275/Demography/C122/SH_4117_C126.psmcfa

# Re-Run PSMC using the model "4+30*2+4+6+10"

psmc -N25 -t5 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C005/C_2931_C005_4.psmc /gxfs_work1/cau/suaph275/Demography/C005/C_2931_C005.psmcfa
psmc -N25 -t5 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C005/C_2931_C012_4.psmc /gxfs_work1/cau/suaph275/Demography/C005/C_2931_C012.psmcfa

psmc -N25 -t5 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C025/C_2747_C025_4.psmc /gxfs_work1/cau/suaph275/Demography/C025/C_2747_C025.psmcfa
psmc -N25 -t5 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C025/C_2747_C021_4.psmc /gxfs_work1/cau/suaph275/Demography/C025/C_2747_C021.psmcfa

psmc -N25 -t5 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C036/SC_2932_C036_4.psmc /gxfs_work1/cau/suaph275/Demography/C036/SC_2932_C036.psmcfa
psmc -N25 -t5 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C036/SC_2932_C046_4.psmc /gxfs_work1/cau/suaph275/Demography/C036/SC_2932_C046.psmcfa

psmc -N25 -t5 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C052/C_1963_C052_4.psmc /gxfs_work1/cau/suaph275/Demography/C052/C_1963_C052.psmcfa
psmc -N25 -t5 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C052/C_1963_C060_4.psmc /gxfs_work1/cau/suaph275/Demography/C052/C_1963_C060.psmcfa

psmc -N25 -t5 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C069_4.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C069.psmcfa
psmc -N25 -t5 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C074_4.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C074.psmcfa

psmc -N25 -t5 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C079/SH_4330_C079_4.psmc /gxfs_work1/cau/suaph275/Demography/C079/SH_4330_C079.psmcfa
psmc -N25 -t5 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C079/SH_4330_C081_4.psmc /gxfs_work1/cau/suaph275/Demography/C079/SH_4330_C081.psmcfa

psmc -N25 -t5 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C093/C_3111_C093_4.psmc /gxfs_work1/cau/suaph275/Demography/C093/C_3111_C093.psmcfa
psmc -N25 -t5 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C093/C_3111_C101_4.psmc /gxfs_work1/cau/suaph275/Demography/C093/C_3111_C101.psmcfa

psmc -N25 -t5 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C109_4.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C109.psmcfa
psmc -N25 -t5 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C104_4.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C104.psmcfa

psmc -N25 -t5 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C122/SH_4117_C122_4.psmc /gxfs_work1/cau/suaph275/Demography/C122/SH_4117_C122.psmcfa
psmc -N25 -t5 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C122/SH_4117_C126_4.psmc /gxfs_work1/cau/suaph275/Demography/C122/SH_4117_C126.psmcfa

# Generate 'split' psmcfa inputs for bootstrapping

/gxfs_home/cau/suaph275/psmc/utils/splitfa /gxfs_work1/cau/suaph275/Demography/C005/C_2931_C005.psmc > /gxfs_work1/cau/suaph275/Demography/C005/C_2931_C005-split.psmcfa
/gxfs_home/cau/suaph275/psmc/utils/splitfa /gxfs_work1/cau/suaph275/Demography/C025/C_2747_C025.psmc > /gxfs_work1/cau/suaph275/Demography/C025/C_2747_C025-split.psmcfa
/gxfs_home/cau/suaph275/psmc/utils/splitfa /gxfs_work1/cau/suaph275/Demography/C036/SC_2932_C036.psmc > /gxfs_work1/cau/suaph275/Demography/C036/SC_2932_C036-split.psmcfa
/gxfs_home/cau/suaph275/psmc/utils/splitfa /gxfs_work1/cau/suaph275/Demography/C052/C_1963_C052.psmc > /gxfs_work1/cau/suaph275/Demography/C052/C_1963_C052-split.psmcfa
/gxfs_home/cau/suaph275/psmc/utils/splitfa /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C069.psmc > /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C069-split.psmcfa
/gxfs_home/cau/suaph275/psmc/utils/splitfa /gxfs_work1/cau/suaph275/Demography/C079/SH_4330_C079.psmc > /gxfs_work1/cau/suaph275/Demography/C079/SH_4330_C079-split.psmcfa
/gxfs_home/cau/suaph275/psmc/utils/splitfa /gxfs_work1/cau/suaph275/Demography/C093/C_3111_C093.psmc > /gxfs_work1/cau/suaph275/Demography/C093/C_3111_C093-split.psmcfa
/gxfs_home/cau/suaph275/psmc/utils/splitfa /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C109.psmc > /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C109-split.psmcfa
/gxfs_home/cau/suaph275/psmc/utils/splitfa /gxfs_work1/cau/suaph275/Demography/C122/SH_4117_C122.psmc > /gxfs_work1/cau/suaph275/Demography/C122/SH_4117_C122-split.psmcfa

# Perform 100 bootstrapping analysis

seq 100 | xargs -i echo /gxfs_home/cau/suaph275/psmc/psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C005/bootstrap/C_2931_C005.round-{}.psmc /gxfs_work1/cau/suaph275/Demography/C005/C_2931_C005-split.psmcfa | sh
seq 100 | xargs -i echo /gxfs_home/cau/suaph275/psmc/psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C025/bootstrap/C_2747_C025.round-{}.psmc /gxfs_work1/cau/suaph275/Demography/C025/C_2747_C025-split.psmcfa | sh
seq 100 | xargs -i echo /gxfs_home/cau/suaph275/psmc/psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C036/bootstrap/SC_2932_C036.round-{}.psmc /gxfs_work1/cau/suaph275/Demography/C036/SC_2932_C036-split.psmcfa | sh
seq 100 | xargs -i echo /gxfs_home/cau/suaph275/psmc/psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C052/bootstrap/C_1963_C052.round-{}.psmc /gxfs_work1/cau/suaph275/Demography/C052/C_1963_C052-split.psmcfa | sh
seq 100 | xargs -i echo /gxfs_home/cau/suaph275/psmc/psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C069/bootstrap/C_1958_C069.round-{}.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C069-split.psmcfa | sh
seq 100 | xargs -i echo /gxfs_home/cau/suaph275/psmc/psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C079/bootstrap/SH_4330_C079.round-{}.psmc /gxfs_work1/cau/suaph275/Demography/C079/SH_4330_C079-split.psmcfa | sh
seq 100 | xargs -i echo /gxfs_home/cau/suaph275/psmc/psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C093/bootstrap/C_3111_C093.round-{}.psmc /gxfs_work1/cau/suaph275/Demography/C093/C_3111_C093-split.psmcfa | sh
seq 100 | xargs -i echo /gxfs_home/cau/suaph275/psmc/psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C109/bootstrap/SC_4107_C109.round-{}.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C109-split.psmcfa | sh
seq 100 | xargs -i echo /gxfs_home/cau/suaph275/psmc/psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C122/bootstrap/SH_4117_C122.round-{}.psmc /gxfs_work1/cau/suaph275/Demography/C122/SH_4117_C122-split.psmcfa | sh

# Plot PSMC results
# Concatenate all the results into one file to plot all samples on one figure

cat C005/C_2931_C005_1.psmc C025/C_2747_C021_1.psmc C036/SC_2932_C036_1.psmc C052/C_1963_C052_1.psmc C069/C_1958_C069_1.psmc C079/SH_4330_C079_1.psmc C093/C_3111_C093_1.psmc C109/SC_4107_C109_1.psmc C122/SH_4117_C122_1.psmc > combined.psmc

psmc_plot.pl -u 1e-08 -g 5 -M 'C2931,C2747,SC2932,C1963,C1958,SH4330,C3111,SC4107,SH4117' Cpmbined combined.psmc 

cat C005/C_2931_C005_2.psmc C025/C_2747_C021_2.psmc C036/SC_2932_C036_2.psmc C052/C_1963_C052_2.psmc C069/C_1958_C069_2.psmc C079/SH_4330_C079_2.psmc C093/C_3111_C093_2.psmc C109/SC_4107_C109_2.psmc C122/SH_4117_C122_2.psmc > combined_2.psmc
psmc_plot.pl -u 1e-08 -g 5 -M 'C2931,C2747,SC2932,C1963,C1958,SH4330,C3111,SC4107,SH4117' Cpmbined_2 combined_2.psmc 

cat C005/C_2931_C012_2.psmc C025/C_2747_C025_2.psmc C036/SC_2932_C046_2.psmc C052/C_1963_C060_2.psmc C069/C_1958_C074_2.psmc C079/SH_4330_C081_2.psmc C093/C_3111_C101_2.psmc C109/SC_4107_C104_2.psmc C122/SH_4117_C126_2.psmc > combined_1_1.psmc
psmc_plot.pl -u 1e-08 -g 5 -M 'C2931,C2747,SC2932,C1963,C1958,SH4330,C3111,SC4107,SH4117' Cpmbined_1_1 combined_1_1.psmc 

cat C005/C_2931_C012_1.psmc C025/C_2747_C025_1.psmc C036/SC_2932_C046_1.psmc C052/C_1963_C060_1.psmc C069/C_1958_C074_1.psmc C079/SH_4330_C081_1.psmc C093/C_3111_C101_1.psmc C109/SC_4107_C104_1.psmc C122/SH_4117_C126_1.psmc > combined_2_2.psmc
psmc_plot.pl -u 1e-08 -g 5 -M 'C2931,C2747,SC2932,C1963,C1958,SH4330,C3111,SC4107,SH4117' Cpmbined_2_2 combined_2_2.psmc 

cat C005/C_2931_C005_1.psmc C025/C_2747_C025_1.psmc C036/SC_2932_C036_1.psmc C052/C_1963_C052_1.psmc C069/C_1958_C069_1.psmc C079/SH_4330_C079_1.psmc C093/C_3111_C093_1.psmc C109/SC_4107_C109_1.psmc C122/SH_4117_C122_1.psmc > combined_1_1.psmc
psmc_plot.pl -u 1e-08 -g 5 -M 'C2931,C2747,SC2932,C1963,C1958,SH4330,C3111,SC4107,SH4117' Combined_1_1 combined_1_1.psmc

cat C005/C_2931_C005_2.psmc C025/C_2747_C025_2.psmc C036/SC_2932_C036_2.psmc C052/C_1963_C052_2.psmc C069/C_1958_C069_2.psmc C079/SH_4330_C079_2.psmc C093/C_3111_C093_2.psmc C109/SC_4107_C109_2.psmc C122/SH_4117_C122_2.psmc > combined_1_2.psmc
psmc_plot.pl -u 1e-08 -g 5 -M 'C2931,C2747,SC2932,C1963,C1958,SH4330,C3111,SC4107,SH4117' Combined_1_2 combined_1_2.psmc

cat C005/C_2931_C012_2.psmc C025/C_2747_C021_2.psmc C036/SC_2932_C046_2.psmc C052/C_1963_C060_2.psmc C069/C_1958_C074_2.psmc C079/SH_4330_C081_2.psmc C093/C_3111_C101_2.psmc C109/SC_4107_C104_2.psmc C122/SH_4117_C126_2.psmc > combined_2_2.psmc
psmc_plot.pl -u 1e-08 -g 5 -M 'C2931,C2747,SC2932,C1963,C1958,SH4330,C3111,SC4107,SH4117' Combined_2_2 combined_2_2.psmc

cat C005/C_2931_C012_1.psmc C025/C_2747_C021_1.psmc C036/SC_2932_C046_1.psmc C052/C_1963_C060_1.psmc C069/C_1958_C074_1.psmc C079/SH_4330_C081_1.psmc C093/C_3111_C101_1.psmc C109/SC_4107_C104_1.psmc C122/SH_4117_C126_1.psmc > combined_2_1.psmc
psmc_plot.pl -u 1e-08 -g 5 -M 'C2931,C2747,SC2932,C1963,C1958,SH4330,C3111,SC4107,SH4117' Combined_2_1 combined_2_1.psmc

# Convert EPS output to PDF (an example)

epstopdf.pl Combined_2_1.eps

## I have re-run the whole script for different individuals from each population to see if the demographic pattern holds true or not. You can follow this process from top to bottom for any sample.

## Additional samples that I used for PSMC analysis

bcftools mpileup -C 50 -q 10 -Ou -f S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C105.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=10 & DP<=60" | bcftools view --exclude-types indels | bcftools sort --temp-dir /gxfs_work1/cau/suaph275/Demography/C109/temp_C105 -Oz -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C105.vcf.gz
bcftools mpileup -C 50 -q 10 -Ou -f S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C106.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=10 & DP<=60" | bcftools view --exclude-types indels | bcftools sort --temp-dir /gxfs_work1/cau/suaph275/Demography/C109/temp_C106 -Oz -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C106.vcf.gz

bcftools mpileup -C 50 -q 10 -Ou -f S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C107.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=10 & DP<=60" | bcftools view --exclude-types indels | bcftools sort --temp-dir /gxfs_work1/cau/suaph275/Demography/C109/temp_C107 -Oz -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C107.vcf.gz
bcftools mpileup -C 50 -q 10 -Ou -f S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C108.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=10 & DP<=60" | bcftools view --exclude-types indels | bcftools sort --temp-dir /gxfs_work1/cau/suaph275/Demography/C109/temp_C108 -Oz -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C108.vcf.gz

bcftools mpileup -C 50 -q 10 -Ou -f S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C110.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=10 & DP<=60" | bcftools view --exclude-types indels | bcftools sort --temp-dir /gxfs_work1/cau/suaph275/Demography/C109/temp_C110 -Oz -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C110.vcf.gz
bcftools mpileup -C 50 -q 10 -Ou -f S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C111.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=10 & DP<=60" | bcftools view --exclude-types indels | bcftools sort --temp-dir /gxfs_work1/cau/suaph275/Demography/C109/temp_C111 -Oz -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C111.vcf.gz

bcftools mpileup -C 50 -q 10 -Ou -f S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C112.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=10 & DP<=60" | bcftools view --exclude-types indels | bcftools sort --temp-dir /gxfs_work1/cau/suaph275/Demography/C109/temp_C112 -Oz -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C112.vcf.gz
bcftools mpileup -C 50 -q 10 -Ou -f S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C113.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=10 & DP<=60" | bcftools view --exclude-types indels | bcftools sort --temp-dir /gxfs_work1/cau/suaph275/Demography/C109/temp_C113 -Oz -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C113.vcf.gz

bcftools mpileup -C 50 -q 10 -Ou -f S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C114.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=10 & DP<=60" | bcftools view --exclude-types indels | bcftools sort --temp-dir /gxfs_work1/cau/suaph275/Demography/C109/temp_C114 -Oz -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C114.vcf.gz

tabix -p vcf /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C105.vcf.gz
tabix -p vcf /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C106.vcf.gz
tabix -p vcf /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C107.vcf.gz
tabix -p vcf /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C108.vcf.gz
tabix -p vcf /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C110.vcf.gz
tabix -p vcf /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C111.vcf.gz
tabix -p vcf /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C112.vcf.gz
tabix -p vcf /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C113.vcf.gz
tabix -p vcf /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C114.vcf.gz

bcftools view /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C105.vcf.gz | vcfutils.pl vcf2fq | gzip > /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C105.fastq.gz
bcftools view /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C106.vcf.gz | vcfutils.pl vcf2fq | gzip > /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C106.fastq.gz

bcftools view /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C107.vcf.gz | vcfutils.pl vcf2fq | gzip > /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C107.fastq.gz
bcftools view /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C108.vcf.gz | vcfutils.pl vcf2fq | gzip > /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C108.fastq.gz

bcftools view /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C110.vcf.gz | vcfutils.pl vcf2fq | gzip > /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C110.fastq.gz
bcftools view /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C111.vcf.gz | vcfutils.pl vcf2fq | gzip > /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C111.fastq.gz

bcftools view /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C112.vcf.gz | vcfutils.pl vcf2fq | gzip > /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C112.fastq.gz
bcftools view /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C113.vcf.gz | vcfutils.pl vcf2fq | gzip > /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C113.fastq.gz

bcftools view /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C114.vcf.gz | vcfutils.pl vcf2fq | gzip > /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C114.fastq.gz

fq2psmcfa -q 20 /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C105.fastq.gz > /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C105.psmcfa
fq2psmcfa -q 20 /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C106.fastq.gz > /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C106.psmcfa

fq2psmcfa -q 20 /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C107.fastq.gz > /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C107.psmcfa
fq2psmcfa -q 20 /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C108.fastq.gz > /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C108.psmcfa

fq2psmcfa -q 20 /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C110.fastq.gz > /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C110.psmcfa
fq2psmcfa -q 20 /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C111.fastq.gz > /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C111.psmcfa

fq2psmcfa -q 20 /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C112.fastq.gz > /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C112.psmcfa
fq2psmcfa -q 20 /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C113.fastq.gz > /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C113.psmcfa

fq2psmcfa -q 20 /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C114.fastq.gz > /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C114.psmcfa

# -t5 -r5 -p "4+25*2+4+6" is 1
# -t15 -r5 -p "4+25*2+4+6" is 2
# -t5 -r5 -p "4+30*2+4+6+10" is 3
# -t15 -r5 -p "4+30*2+4+6+10" is 4


psmc -N25 -t5 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C105_1.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C105.psmcfa
psmc -N25 -t5 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C106_1.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C106.psmcfa

psmc -N25 -t5 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C107_1.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C107.psmcfa
psmc -N25 -t5 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C108_1.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C108.psmcfa

psmc -N25 -t5 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C110_1.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C110.psmcfa
psmc -N25 -t5 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C111_1.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C111.psmcfa

psmc -N25 -t5 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C112_1.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C112.psmcfa
psmc -N25 -t5 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C113_1.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C113.psmcfa

psmc -N25 -t5 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C114_1.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C114.psmcfa

psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C105_2.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C105.psmcfa
psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C106_2.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C106.psmcfa

psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C107_2.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C107.psmcfa
psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C108_2.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C108.psmcfa

psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C110_2.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C110.psmcfa
psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C111_2.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C111.psmcfa

psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C112_2.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C112.psmcfa
psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C113_2.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C113.psmcfa

psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C114_2.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C114.psmcfa

psmc -N25 -t5 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C105_3.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C105.psmcfa
psmc -N25 -t5 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C106_3.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C106.psmcfa

psmc -N25 -t5 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C107_3.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C107.psmcfa
psmc -N25 -t5 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C108_3.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C108.psmcfa

psmc -N25 -t5 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C110_3.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C110.psmcfa
psmc -N25 -t5 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C111_3.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C111.psmcfa

psmc -N25 -t5 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C112_3.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C112.psmcfa
psmc -N25 -t5 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C113_3.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C113.psmcfa

psmc -N25 -t5 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C114_3.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C114.psmcfa

psmc -N25 -t15 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C105_4.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C105.psmcfa
psmc -N25 -t15 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C106_4.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C106.psmcfa

psmc -N25 -t15 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C107_4.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C107.psmcfa
psmc -N25 -t15 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C108_4.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C108.psmcfa

psmc -N25 -t15 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C110_4.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C110.psmcfa
psmc -N25 -t15 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C111_4.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C111.psmcfa

psmc -N25 -t15 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C112_4.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C112.psmcfa
psmc -N25 -t15 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C113_4.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C113.psmcfa

psmc -N25 -t15 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C114_4.psmc /gxfs_work1/cau/suaph275/Demography/C109/SC_4107_C114.psmcfa



## Additional samples that I used for PSMC (re-ran everything using different time segment patterns)

bcftools mpileup -C 50 -q 10 -Ou -f S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C061.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=10 & DP<=60" | bcftools view --exclude-types indels | bcftools sort --temp-dir /gxfs_work1/cau/suaph275/Demography/C069/temp_C061 -Oz -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C061.vcf.gz
bcftools mpileup -C 50 -q 10 -Ou -f S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C062.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=10 & DP<=60" | bcftools view --exclude-types indels | bcftools sort --temp-dir /gxfs_work1/cau/suaph275/Demography/C069/temp_C062 -Oz -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C062.vcf.gz

bcftools mpileup -C 50 -q 10 -Ou -f S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C063.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=10 & DP<=60" | bcftools view --exclude-types indels | bcftools sort --temp-dir /gxfs_work1/cau/suaph275/Demography/C069/temp_C063 -Oz -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C063.vcf.gz
bcftools mpileup -C 50 -q 10 -Ou -f S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C064.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=10 & DP<=60" | bcftools view --exclude-types indels | bcftools sort --temp-dir /gxfs_work1/cau/suaph275/Demography/C069/temp_C064 -Oz -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C064.vcf.gz

bcftools mpileup -C 50 -q 10 -Ou -f S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C067.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=10 & DP<=60" | bcftools view --exclude-types indels | bcftools sort --temp-dir /gxfs_work1/cau/suaph275/Demography/C069/temp_C067 -Oz -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C067.vcf.gz
bcftools mpileup -C 50 -q 10 -Ou -f S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C070.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=10 & DP<=60" | bcftools view --exclude-types indels | bcftools sort --temp-dir /gxfs_work1/cau/suaph275/Demography/C069/temp_C070 -Oz -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C070.vcf.gz

bcftools mpileup -C 50 -q 10 -Ou -f S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C071.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=10 & DP<=60" | bcftools view --exclude-types indels | bcftools sort --temp-dir /gxfs_work1/cau/suaph275/Demography/C069/temp_C071 -Oz -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C071.vcf.gz
bcftools mpileup -C 50 -q 10 -Ou -f S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C072.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=10 & DP<=60" | bcftools view --exclude-types indels | bcftools sort --temp-dir /gxfs_work1/cau/suaph275/Demography/C069/temp_C072 -Oz -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C072.vcf.gz

bcftools mpileup -C 50 -q 10 -Ou -f S_chilense_Hirise.fasta /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/C073.bam | bcftools call -c | bcftools filter --SnpGap 10 -i "DP>=10 & DP<=60" | bcftools view --exclude-types indels | bcftools sort --temp-dir /gxfs_work1/cau/suaph275/Demography/C069/temp_C073 -Oz -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C073.vcf.gz

tabix -p vcf /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C061.vcf.gz
tabix -p vcf /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C062.vcf.gz
tabix -p vcf /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C063.vcf.gz
tabix -p vcf /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C064.vcf.gz
tabix -p vcf /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C067.vcf.gz
tabix -p vcf /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C070.vcf.gz
tabix -p vcf /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C071.vcf.gz
tabix -p vcf /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C072.vcf.gz
tabix -p vcf /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C073.vcf.gz

bcftools view /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C061.vcf.gz | vcfutils.pl vcf2fq | gzip > /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C061.fastq.gz
bcftools view /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C062.vcf.gz | vcfutils.pl vcf2fq | gzip > /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C062.fastq.gz

bcftools view /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C063.vcf.gz | vcfutils.pl vcf2fq | gzip > /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C063.fastq.gz
bcftools view /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C064.vcf.gz | vcfutils.pl vcf2fq | gzip > /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C064.fastq.gz

bcftools view /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C067.vcf.gz | vcfutils.pl vcf2fq | gzip > /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C067.fastq.gz
bcftools view /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C070.vcf.gz | vcfutils.pl vcf2fq | gzip > /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C070.fastq.gz

bcftools view /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C071.vcf.gz | vcfutils.pl vcf2fq | gzip > /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C071.fastq.gz
bcftools view /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C072.vcf.gz | vcfutils.pl vcf2fq | gzip > /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C072.fastq.gz

bcftools view /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C073.vcf.gz | vcfutils.pl vcf2fq | gzip > /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C073.fastq.gz

fq2psmcfa -q 20 /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C061.fastq.gz > /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C061.psmcfa
fq2psmcfa -q 20 /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C062.fastq.gz > /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C062.psmcfa

fq2psmcfa -q 20 /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C063.fastq.gz > /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C063.psmcfa
fq2psmcfa -q 20 /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C064.fastq.gz > /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C064.psmcfa

fq2psmcfa -q 20 /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C067.fastq.gz > /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C067.psmcfa
fq2psmcfa -q 20 /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C070.fastq.gz > /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C070.psmcfa

fq2psmcfa -q 20 /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C071.fastq.gz > /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C071.psmcfa
fq2psmcfa -q 20 /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C072.fastq.gz > /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C072.psmcfa

fq2psmcfa -q 20 /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C073.fastq.gz > /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C073.psmcfa

# -t5 -r5 -p "4+25*2+4+6" is 1
# -t15 -r5 -p "4+25*2+4+6" is 2
# -t5 -r5 -p "4+30*2+4+6+10" is 3
# -t15 -r5 -p "4+30*2+4+6+10" is 4


psmc -N25 -t5 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C061_1.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C061.psmcfa
psmc -N25 -t5 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C062_1.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C062.psmcfa

psmc -N25 -t5 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C063_1.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C063.psmcfa
psmc -N25 -t5 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C064_1.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C064.psmcfa

psmc -N25 -t5 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C067_1.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C067.psmcfa
psmc -N25 -t5 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C070_1.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C070.psmcfa

psmc -N25 -t5 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C071_1.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C071.psmcfa
psmc -N25 -t5 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C072_1.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C072.psmcfa

psmc -N25 -t5 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C073_1.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C073.psmcfa

psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C061_2.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C061.psmcfa
psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C062_2.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C062.psmcfa

psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C063_2.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C063.psmcfa
psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C064_2.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C064.psmcfa

psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C067_2.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C067.psmcfa
psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C070_2.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C070.psmcfa

psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C071_2.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C071.psmcfa
psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C072_2.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C072.psmcfa

psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C073_2.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C073.psmcfa

psmc -N25 -t5 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C061_3.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C061.psmcfa
psmc -N25 -t5 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C062_3.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C062.psmcfa

psmc -N25 -t5 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C063_3.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C063.psmcfa
psmc -N25 -t5 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C064_3.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C064.psmcfa

psmc -N25 -t5 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C067_3.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C067.psmcfa
psmc -N25 -t5 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C070_3.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C070.psmcfa

psmc -N25 -t5 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C071_3.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C071.psmcfa
psmc -N25 -t5 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C072_3.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C072.psmcfa

psmc -N25 -t5 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C073_3.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C073.psmcfa

psmc -N25 -t15 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C061_4.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C061.psmcfa
psmc -N25 -t15 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C062_4.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C062.psmcfa

psmc -N25 -t15 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C063_4.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C063.psmcfa
psmc -N25 -t15 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C064_4.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C064.psmcfa

psmc -N25 -t15 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C067_4.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C067.psmcfa
psmc -N25 -t15 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C070_4.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C070.psmcfa

psmc -N25 -t15 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C071_4.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C071.psmcfa
psmc -N25 -t15 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C072_4.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C072.psmcfa

psmc -N25 -t15 -r5 -p "4+30*2+4+6+10" -o /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C073_4.psmc /gxfs_work1/cau/suaph275/Demography/C069/C_1958_C073.psmcfa






