#!/bin/bash

# Anik Dutta

# Estimate chromosome length from a vcf file

# Store the input VCF file name in a variable
# input_file=$1 if you want to use this command, then you dont need to put the vcf name in the below command. you can use "$input_file". but then execute the code like ""./estimate_chromosome_length.sh input.vcf > output.txt". Here specify the input vcf.

chromosomes=$(grep -v "^#" "Chilense_99.genotyped.SNP.Clean.MAF.MAXMiss.MinDP.recode.dustmasked.headers.vcf" | awk '{print $1}' | sort | uniq)

for chromosome in $chromosomes; do
  max_pos=$(grep -v "^#" "Chilense_99.genotyped.SNP.Clean.MAF.MAXMiss.MinDP.recode.dustmasked.headers.vcf" | awk -v c=$chromosome '$1 == c {print $2}' | sort -nr | head -1)
  echo "Chromosome $chromosome has length $max_pos"
done

#execute
chmod +x estimated_chr_length.sh

sh estimated_chr_length.sh > chr_length.txt