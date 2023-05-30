#!/bin/bash

# Anik Dutta

# This bash script is used to extract specific information from a variant call format (VCF) file. 
# Each command in the script extracts a particular quality metric and writes the output to a separate file.

# egrep -v "^#" : The egrep command with -v option is used to print all lines not matching the given pattern. Here it's used to ignore lines starting with '#'. These lines are usually headers in a VCF file.

# cut -f 8 : The cut command is used to extract columns from a file. Here it extracts the 8th column from the VCF file, which is the INFO field that contains various quality metrics.

# sed 's/^.*;DP=\([0-9]*\);.*$/\1/' : The sed command is used here for substitution. It finds the DP (Read Depth) value in the INFO field and outputs it.

# > /gxfs_work1/cau/suaph275/Chilense_Hiseq/genotyped_vcf/GATK_quality_check/depth_66mio.txt : The output is redirected to a text file.

# This command is repeated for various quality metrics like QD (Quality by Depth), FS (Fisher Strand Bias), SOR (Symmetric Odds Ratio), MQ (Mapping Quality), and MQRankSum (Mapping Quality Rank Sum).

egrep -v "^#" Chilense_99.genotyped.SNP.Clean.vcf.recode.vcf | cut -f 8 | sed 's/^.*;DP=\([0-9]*\);.*$/\1/' > /gxfs_work1/cau/suaph275/Chilense_Hiseq/genotyped_vcf/GATK_quality_check/depth_66mio.txt

egrep -v "^#" Chilense_99.genotyped.SNP.Clean.vcf.recode.vcf | cut -f 8 | sed 's/^.*;QD=\([0-9]*.[0-9]*\);.*$/\1/' > /gxfs_work1/cau/suaph275/Chilense_Hiseq/genotyped_vcf/GATK_quality_check/QD_66mio.txt

egrep -v "^#" Chilense_99.genotyped.SNP.Clean.vcf.recode.vcf | cut -f 8 | sed 's/^.*;FS=\([0-9]*.[0-9]*\);.*$/\1/' > /gxfs_work1/cau/suaph275/Chilense_Hiseq/genotyped_vcf/GATK_quality_check/FS_66mio.txt

egrep -v "^#" Chilense_99.genotyped.SNP.Clean.vcf.recode.vcf | cut -f 8 | sed 's/^.*;SOR=\([0-9]*.[0-9]*\);.*$/\1/' > /gxfs_work1/cau/suaph275/Chilense_Hiseq/genotyped_vcf/GATK_quality_check/SOR_66mio.txt

egrep -v "^#" Chilense_99.genotyped.SNP.Clean.vcf.recode.vcf | cut -f 8 | sed 's/^.*;MQ=\([0-9].[0-9]*\);.*$/\1/' > /gxfs_work1/cau/suaph275/Chilense_Hiseq/genotyped_vcf/GATK_quality_check/MQ_66mio.txt

egrep -v "^#" Chilense_99.genotyped.SNP.Clean.vcf.recode.vcf | cut -f 8 | sed 's/^.*;MQRankSum=\(\-\{0,1\}[0-9]\{1,\}.[0-9]*\);.*$/\1/' > /gxfs_work1/cau/suaph275/Chilense_Hiseq/genotyped_vcf/GATK_quality_check/MQ_rank_66mio.txt
