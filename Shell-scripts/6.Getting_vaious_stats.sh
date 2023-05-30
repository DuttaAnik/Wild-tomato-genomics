#!/bin/bash

# Anik Dutta

### Getting statistics and coverage from raw reads and bam files (nifa-crown-rust/trimmed_reads) ###

## Collecting stats with samtools. Change the working directory to the folder where the bam files are.

cd /gxfs_work1/cau/suaph275/Bowtie2_results/Bam_files

(for F in `ls *.bam`; do echo $F; samtools flagstat -@24 $F; echo -en "\\n"; done;) > mapping_stats.txt

# Optional: Mosdepth to calculate coverage depth for BAM files. Here it's running in fast mode, and computing coverage in 100,000 base pair intervals.

for F in `ls *.bam`; do echo $F; mosdepth -n --fast-mode --d4 --by 100000 /gxfs_work1/cau/suaph275/Chilense_Hiseq/Bam_RG_files/MosDepth/"$F".wgs $F; done

## This script is used to generate coverage statistics for a set of BAM files (main code used to generate the figure)

# For each sequence in the FASTA file, print the sequence name and its length and then sort by the sequence name.

cat S_chilense_Hirise.fasta | bioawk -c fastx '{ print $name, length($seq) }' |  sort -k1,1 > S_chilense_Hirise.scaffoldlen

# Use AWK to process the scaffold lengths file. For each line in the file, 
# print the sequence name, 0, and the sequence length, separated by tabs.
# This creates a BED (Browser Extensible Data) file with the sequence name 
# and length for each sequence in the FASTA file.

awk '{print $1"\t" 0 "\t" $2}' S_chilense_Hirise.scaffoldlen > S_chilense_Hirise.scaffoldlen.bed

# Create an array 'input' containing the paths to all BAM files in a directory.
input=(/gxfs_work1/cau/suaph275/Bowtie2_results/Bam_files/*.bam)

## Use 'bedtools' to convert the BAM file to BED format, then sort the BED file by the first and second fields.
# Use 'bedtools coverage' to calculate the mean coverage for each entry in the BED file, 
# using the scaffold lengths file for the genome file (-g) and the scaffold lengths BED file for the BED file (-a).
# The output is written to a file with the same name as the BAM file, but with a .cov extension.
# For each BAM file in the 'input' array:

for ((i=0;i<=${#input[@]};i++))
do
  bedtools bamtobed -i "${input[i]}" | sort -k1,1 -k2,2n -S 75% | bedtools coverage -mean -sorted -g S_chilense_Hirise.scaffoldlen -a S_chilense_Hirise.scaffoldlen.bed -b - > $(basename "${input[i]%.rm.sorted.bam}.cov")
done

### Breadth of coverage

# count the length of the reference genome
bowtie2-inspect -s S_chilense_Hirise | awk '{ FS = "\t" } ; BEGIN{L=0}; {L=L+$3}; END{print L}'

seqkit stats S_chilense_Hirise.fasta

file                     format  type  num_seqs        sum_len  min_len    avg_len      max_len

S_chilense_Hirise.fasta  FASTA   DNA      3,563  1,100,934,140    4,999  308,990.8  114,446,168

## Getting number of singletons and putting it together with the sample names

grep 'singletons' mapping_stats.txt | cut -d ' ' -f1 > number_singleton.txt # -d is the delimitor and -f1 is the column number
paste bam_list.txt number_singleton.txt > final.txt

# Remove SNPs within 5 bases of an INDEL

# B=1025_230_255.recode.vcf

# bcftools filter -g 5 $B >filter.1.vcf

# Using VCF tools to extract different statistics from a hard filtered VCF file

vcftools --vcf input.vcf --het --out het
vcftools --vcf input.vcf --missing-site --out miss_site
vcftools --vcf input.vcf --missing-indv --out miss_indv
vcftools --vcf input.vcf --site-quality --out site_qual
vcftools --vcf input.vcf --depth --out depth_ind
vcftools --vcf input.vcf --freq2 --max-alleles 2 --out allele_freq
vcftools --vcf input.vcf --site-mean-depth --out mean_depth_site

# Singleton stats from the VCf file. Use the raw Clean vcf before applying maf filtering for singleton and non-ref variants counting. If you apply MAF then you cannot count the singletons. Here, I accidentally used the filtered VCF file. So, did not use this result to produce the figures. Rather I used the last command with bcftools stats

for i in 1958 1963 2747 2931 2932 3111 4107 4117 4330 
do 
full_pop=""pop."$i" 
module load vcftools 
vcftools --vcf Chilense_99.genotyped.SNP.Clean.MAF.MAXMiss.MinDP.recode.dustmasked.headers.vcf --keep "$i".txt --singletons --out "$full_pop".singleton 
done

# Remove the first line of the file (sed 1d), extract columns 3 and 5 (cut -f 3,5),
# remove lines starting with 'D' (grep -v '^D'), sort the remaining lines, 
# and count the unique entries (uniq -c). The output is redirected to a new file 4117.singletone.count. This contains the number of singletons for each individual in that population.

sed 1d pop.4117.singletons | cut -f 3,5 | grep -v '^D' | sort | uniq -c > 4117.singletone.count

# Use bcftools to get different statistics per population. Here, the VCF is not hard filtered. Hence, this was used for singletons count. This is the real result. 

for i in 1958_pop 1963_pop 2747_pop 2931_pop 2932_pop 3111_pop 4107_pop 4117_pop 4330_pop
do
bcftools stats -S "$i".txt /gxfs_work1/cau/suaph275/Clean_VCF/ABBA_pop_stats/Chilense_99.genotyped.SNP.Clean.vcf.recode.vcf > "$i".stats
done

