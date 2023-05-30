#!/bin/bash

# Anik Dutta

# Creating SNPEff database and running it for variant annotation (https://genomics.sschmeier.com/ngs-voi/index.html and helminth PDF)

# First configure the snpEff config file. Inside the config file, write the following after the # Databases & Genomes section

# Solunum_chilense

Solunum_chilense.genome : Solunum_chilense

# Then, create 'data' folder inside the snpEff folder. Copy the reference genome and gff file and give the names sequences.fa and genes.gff

mkdir /gxfs_home/cau/suaph275/snpEff/data/Solunum_chilense
cp $WORK/Chilense_Hiseq/S_chilense_Hirise.fasta /gxfs_home/cau/suaph275/snpEff/data/Solunum_chilense/sequences.fa
cp $WORK/Chilense_Hiseq/S_chilense_Hirise.gff3 /gxfs_home/cau/suaph275/snpEff/data/Solunum_chilense/genes.gff

# Build the database

java -jar snpEff.jar build -c snpEff.config -gff3 -v Solunum_chilense

# Run the annotation

java -Xmx50G -jar /gxfs_home/cau/suaph275/snpEff/snpEff.jar -v Solunum_chilense /gxfs_work1/cau/suaph275/Bowtie2_results/genotyped_vcf/dustmasked_header_maf0.01_min_meanDP2_maxmiss_0.5.recode.vcf > SNP_eff.vcf

# Create the file with list of variants of interest

cat SNP_eff.vcf | /gxfs_home/cau/suaph275/snpEff/scripts/vcfEffOnePerLine.pl | java -Xmx50G -jar /gxfs_home/cau/suaph275/snpEff/SnpSift.jar extractFields - CHROM POS REF ALT AF "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" > SNP_filtered.eff.txt