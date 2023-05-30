#!/bin/bash

# GWAS script for wild tomato (Anik Dutta)

# I ran the GWAS using only 12 largest scaffolds. So, first I separated the combined VCF file into individual scaffold files and then combined the 12 scaffolds into one VCF

# Get the names of the scaffolds and the VCF header

grep -v "^#" input.vcf | cut -f1 | sort | uniq > scaff_names​
grep "^#" input.vcf > vcf_header

# Make VCF per scaffold and add the vcf header (dont forget to do that)

while read line
do
grep -v "^#" input.vcf  | awk -v scaff=$line '{ if ( $1==scaff ) print }' | cat vcf_header - > $line.vcf
done < scaff_names

# MAking a combined VCF with 12 scaffolds. First, use bgzip to make .gz vcf files. I used the vcf of scaffolds produced by the above codes.

bgzip Scaffold_1__2_contigs__length_78070536.vcf

# Then use Tabix to index each .gz file

tabix Scaffold_2__1_contigs__length_61731525.vcf.gz

# Then concat using bcftools
bcftools concat -o merged.vcf.gz Scaffold_1__2_contigs__length_78070536.vcf.gz Scaffold_2__1_contigs__length_61731525.vcf.gz Scaffold_3__3_contigs__length_114446166.vcf.gz Scaffold_4__2_contigs__length_60228994.vcf.gz Scaffold_5__2_contigs__length_88504823.vcf.gz Scaffold_6__2_contigs__length_83591330.vcf.gz Scaffold_7__2_contigs__length_80977755.vcf.gz Scaffold_8__2_contigs__length_76899493.vcf.gz Scaffold_9__1_contigs__length_70652527.vcf.gz Scaffold_10__2_contigs__length_69708562.vcf.gz Scaffold_11__3_contigs__length_75158268.vcf.gz Scaffold_12__3_contigs__length_71445296.vcf.gz

# Install the package vcf2gwas from https://github.com/frankvogt/vcf2gwas

# Activate the conda environment to run it: conda activate GWAS

# Without PCA
vcf2gwas -v merged.vcf -pf Flower_gwas_data.csv -ap -lmm 2 -nl -fs 12 -T 16

# With PCA (first 3 components -c 3)

vcf2gwas -v merged.vcf -pf FT_gwas_November.csv -ap -cf PCA -c 3 -lmm 1 -nq -nl -fs 12 -T 12

# 01/03/2023 using pca and mean FT from all the replication that Severin gave.

vcf2gwas -v merged.vcf -pf Mean_FT_GWAS.csv -ap -cf PCA -c 3 -lmm 2 -nq -nl -fs 12 -T 12

# Without PCA and with maf (-q o.05)
vcf2gwas -v merged.vcf -pf Mean_FT_GWAS.csv -ap -q 0.05 -lmm 2 -nq -nl -fs 12 -T 12

# With PCA and with maf
vcf2gwas -v merged.vcf -pf Mean_FT_GWAS.csv -ap -cf PCA -c 3 -q 0.05 -lmm 2 -nq -nl -fs 12 -T 12

