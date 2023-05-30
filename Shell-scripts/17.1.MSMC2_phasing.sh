#!/bin/bash

# Anik Dutta

# Phasing the vcf files generated from Bcftools to generate input files for msmc2 analysis

## first create index files for all the vcf files in each individual population directory

cd /gxfs_work1/cau/suaph275/MSMC_analysis/1958_pop

for file in *.vcf.gz 
do 
tabix -p vcf $file 
done

## Create a list of vcf files for each scaffold separately. It contains all the individuals for that specific scaffold

cat /gxfs_work1/cau/suaph275/MSMC_analysis/data_chr.txt | while read Scaffold; do ls ${Scaffold}*.vcf.gz > ${Scaffold}.list; done

## Merge all the vcf from all individuals into one vcf file per scaffold

cat /gxfs_work1/cau/suaph275/MSMC_analysis/data_chr.txt | while read Scaffold; do bcftools merge -l ${Scaffold}.list -Oz -o ${Scaffold}.unphase.vcf.gz; done

## Keep biallelic SNPs

cat /gxfs_work1/cau/suaph275/MSMC_analysis/data_chr.txt | while read Scaffold; do bcftools view -M 2 -O z ${Scaffold}.unphase.vcf.gz > ${Scaffold}.unphase.bial.vcf.gz; done


## Run beagle to phase the vcf per scaffold. Beagle did not run for one individual and one scaffold. Same with shapeit 2. So, I had to merge all individuals together into one file.

cat /gxfs_work1/cau/suaph275/MSMC_analysis/data_chr.txt | while read Scaffold; do java -jar /gxfs_home/cau/suaph275/beagle.19Mar22.0da.jar gt=${Scaffold}.unphase.bial.vcf.gz out=${Scaffold}.phased window=5 overlap=2 iterations=30 burnin=10 gp=true err=0.001 ; done


## Index all the phased vcf files 
cat /gxfs_work1/cau/suaph275/MSMC_analysis/data_chr.txt | while read Scaffold ; do  tabix -p vcf ${Scaffold}.phased.vcf.gz; done

## Create vcf files for each individual and each scaffold

cat data_chr.txt | while read Scaffold ; do  bcftools view -s C104 ${Scaffold}.phased.vcf.gz -Oz -o /gxfs_work1/cau/suaph275/MSMC_analysis/4107_pop/${Scaffold}.C104.phased.vcf.gz; done

cat data_chr.txt | while read Scaffold ; do  bcftools view -s C106 ${Scaffold}.phased.vcf.gz -Oz -o /gxfs_work1/cau/suaph275/MSMC_analysis/4107_pop/${Scaffold}.C106.phased.vcf.gz; done

cat data_chr.txt | while read Scaffold ; do  bcftools view -s C108 ${Scaffold}.phased.vcf.gz -Oz -o /gxfs_work1/cau/suaph275/MSMC_analysis/4107_pop/${Scaffold}.C108.phased.vcf.gz; done


cat data_chr.txt | while read Scaffold ; do  bcftools view -s C118 ${Scaffold}.phased.vcf.gz -Oz -o /gxfs_work1/cau/suaph275/MSMC_analysis/4117_pop/${Scaffold}.C118.phased.vcf.gz; done

cat data_chr.txt | while read Scaffold ; do  bcftools view -s C121 ${Scaffold}.phased.vcf.gz -Oz -o /gxfs_work1/cau/suaph275/MSMC_analysis/4117_pop/${Scaffold}.C121.phased.vcf.gz; done

cat data_chr.txt | while read Scaffold ; do  bcftools view -s C123 ${Scaffold}.phased.vcf.gz -Oz -o /gxfs_work1/cau/suaph275/MSMC_analysis/4117_pop/${Scaffold}.C123.phased.vcf.gz; done








