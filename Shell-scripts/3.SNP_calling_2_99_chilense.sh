####################################################################
# This is a continuation of the previous script for SNP calling. script use to combine gvcfs, genotype vcfs, SNP calling, filtering.
# Anik Dutta (02/04/2022)     
####################################################################


# CombineGVCFs from the previous script. Assign the input directory to a variable and set the base name for finding the input files

REF="/gxfs_work1/cau/suaph275/Chilense_Hiseq/S_chilense_Hirise"
BASE="/gxfs_work1/cau/suaph275/Chilense_Hiseq"

# The first step is to make a list of all the g.vcf files generated with the HaplotypeCaller.
find $BASE/gvcf_files -name "*.g.vcf" > Chilense_99.g.vcf.list

# Then, it is possible to join all the files with the tool CombineGVCFs

module load gatk

echo -e '#!/bin/bash' > SNP_calling_99_bwa.sh
echo "gatk --java-options \"-Xmx120G\" CombineGVCFs -R "$REF".fasta \
-V $BASE/Chilense_99.g.vcf.list  \
-O $BASE/combined_gvcf/Chilense_99.g.vcf" >> SNP_calling_99_bwa.sh

####################################################################
# script use to map to reference target species #     
####################################################################
# run GATK4 GenotypeGVCFs on each file individually

## run GenotypeGVCF

echo "gatk --java-options \"-Xmx120G\" GenotypeGVCFs -R "$REF".fasta \
	   --verbosity ERROR \
	   -V $BASE/combined_gvcf/Chilense_99.g.vcf \
	   --max-alternate-alleles 3 \
	   -O $BASE/genotyped_vcf/Chilense_99.genotypedGVC.vcf" >> SNP_calling_99_bwa.sh


# SelectVariants. INDELs and SNPs
echo "gatk SelectVariants \
	    -R "$REF".fasta \
	    -V $BASE/genotyped_vcf/Chilense_99.genotypedGVC.vcf \
	    --select-type-to-include INDEL \
		--select-type-to-include SNP \
	    -O $BASE/genotyped_vcf/Chilense_99.genotyped.SNP.and.INDEL.vcf" >> SNP_calling_99_bwa.sh

# Keep only Biallelic SNPs
echo "gatk SelectVariants -R /gxfs_work1/cau/suaph275/Chilense_Hiseq/S_chilense_Hirise.fasta -V /gxfs_work1/cau/suaph275/Chilense_Hiseq/genotyped_vcf/Chilense_99.genotyped.SNP.and.INDEL.vcf --restrict-alleles-to BIALLELIC -O /gxfs_work1/cau/suaph275/Chilense_Hiseq/genotyped_vcf/Chilense_99.genotyped.BIALLELIC.SNP.and.INDEL.vcf" >> SNP_calling_99_bwa.sh

# Remove SNPs within 5 bp of an indel

echo "bcftools filter -g 5 Chilense_99.genotyped.BIALLELIC.SNP.and.INDEL.vcf > Chilense_99.genotyped.BIALLELIC.SNP.and.INDEL.5bp.vcf" >> SNP_calling_99_bwa.sh

### IMPORTANT set the minimum genotyping rate
# Indel plus SNPS

"gatk --java-options "-Xmx120G" VariantFiltration -R /gxfs_work1/cau/suaph275/Chilense_Hiseq/S_chilense_Hirise.fasta    -V /gxfs_work1/cau/suaph275/Chilense_Hiseq/genotyped_vcf/Chilense_99.genotyped.BIALLELIC.SNP.and.INDEL.5bp.vcf    --filter-name "nSamples" --filter-expression "AN < 5"  --filter-name "QDFilter" --filter-expression "QD < 2.0"  --filter-name "MQ" --filter-expression "MQ < 40.0"  --filter-name "QUAL30" --filter-expression "QUAL < 30.0" --filter-name "SOR"  --filter-expression "SOR > 3.0"  --filter-name "FS"  --filter-expression "FS > 60.0"     --filter-name "Low_depth3" --filter-expression "DP < 3" --filter-name "ReadPosRankSum_filter" --filter-expression "ReadPosRankSum < -8.0"   --filter-name "ReadPosRankSum_filter"  --filter-expression "ReadPosRankSum > 8.0"  --filter-name "MQRankSum_filter"   --filter-expression "MQRankSum < -12.5"   --filter-name "MQRankSum_filter"   --filter-expression "MQRankSum > 12.5"  -O /gxfs_work1/cau/suaph275/Chilense_Hiseq/genotyped_vcf/Chilense_99.genotyped.SNP.and.INDEL.5bp.Filter.Applied.vcf" >> SNP_calling_99.sh

# Remove filtered variants

echo "gatk SelectVariants -variant /gxfs_work1/cau/suaph275/Chilense_Hiseq/genotyped_vcf/Chilense_99.genotyped.SNP.and.INDEL.5bp.Filter.Applied.vcf --reference /gxfs_work1/cau/suaph275/Chilense_Hiseq/S_chilense_Hirise.fasta --output /gxfs_work1/cau/suaph275/Chilense_Hiseq/genotyped_vcf/Chilense_99.genotyped.SNP.and.INDEL.5bp.Filter.Removed.vcf --exclude-filtered true" >> SNP_calling_99.sh

sbatch --job-name=SNP_calling --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=32 --qos=long --time=10-00:00:00 --mem=120G --error=job.%J.err --output=job.%J.out --mail-type=END,FAIL --mail-user=anik.dutta@phytomed.uni-kiel.de --partition=cluster SNP_calling_99.sh

# hard filtering using vcftools

vcftools --vcf /gxfs_work1/cau/suaph275/Chilense_Hiseq/genotyped_vcf/Chilense_99.genotyped.SNP.and.INDEL.5bp.Filter.Removed.vcf --remove-indels --recode --recode-INFO-all --out /gxfs_work1/cau/suaph275/Chilense_Hiseq/genotyped_vcf/Chilense_99.genotyped.SNP.Clean.vcf

# Number of SNPS=66984130

vcftools --vcf Chilense_99.genotyped.SNP.Clean.vcf.recode.vcf --recode --recode-INFO-all --maf 0.03 --max-missing 0.5 --min-meanDP 2 --out Chilense_99.genotyped.SNP.Clean.MAF.MAXMiss.MinDP

# Number of SNPS=32074624

## dustmask snps in the low complexity regions
dustmasker -in scaffolds.reduced.fa -out scaffolds.reduced.dustmasked.list -outfmt acclist
awk 'split($1,a,">") {print a[2]"\t"$2-1"\t"$3}' scaffolds.reduced.dustmasked.list > dustmasked.bed #also need to remove -1 values from the bed file

bedtools subtract -a genotyped_vcf/Chilense_99.genotyped.SNP.Clean.MAF.MAXMiss.MinDP.recode.vcf -b dustmasked.bed > genotyped_vcf/Chilense_99.genotyped.SNP.Clean.MAF.MAXMiss.MinDP.recode.dustmasked.vcf
grep "#" Chilense_99.genotyped.SNP.Clean.MAF.MAXMiss.MinDP.recode.vcf > header
cat header Chilense_99.genotyped.SNP.Clean.MAF.MAXMiss.MinDP.recode.dustmasked.vcf > Chilense_99.genotyped.SNP.Clean.MAF.MAXMiss.MinDP.recode.dustmasked.headers.vcf

# Number of SNPS=30834951

This is the final VCF file that I used for all downstream analysis.

## INDEL select and filter

gatk SelectVariants --variant /gxfs_work1/cau/suaph275/Chilense_Hiseq/genotyped_vcf/Chilense_99.genotyped.BIALLELIC.SNP.and.INDEL.vcf --output /gxfs_work1/cau/suaph275/Chilense_Hiseq/genotyped_vcf/INDEL.vcf --reference /gxfs_work1/cau/suaph275/Bowtie2_results/S_chilense_Hirise.fasta --select-type-to-include INDEL

# Number of INDELs=11775300

gatk VariantFiltration -V /gxfs_work1/cau/suaph275/Chilense_Hiseq/genotyped_vcf/INDEL.vcf -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "FS > 200.0" --filter-name "FS200" -filter "SOR > 10.0" --filter-name "SOR10" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8.0" -O /gxfs_work1/cau/suaph275/Chilense_Hiseq/genotyped_vcf/INDEL.vcf.Filters.Added.vcf 

gatk SelectVariants -variant /gxfs_work1/cau/suaph275/Chilense_Hiseq/genotyped_vcf/INDEL.vcf.Filters.Added.vcf --reference /gxfs_work1/cau/suaph275/Bowtie2_results/S_chilense_Hirise.fasta --output /gxfs_work1/cau/suaph275/Chilense_Hiseq/genotyped_vcf/INDEL.vcf.Filters.Removed.vcf --exclude-filtered true

# Number of INDELs=11061744


## To include invariant sites (monomorphic variants for estimating pi (nucleotide diversity) and dxy (absolute measures of divergence). You have to use the combined g.vcf file that was created with one of the above codes. I deleted this file since it takes a lot of space. Many do not follow this procedure but I followed according to some other population genomics paper in plants.
 
gatk --java-options "-Xmx120G" GenotypeGVCFs -R /gxfs_work1/cau/suaph275/Chilense_Hiseq/S_chilense_Hirise.fasta    --verbosity ERROR    -V /gxfs_work1/cau/suaph275/Chilense_Hiseq/combined_gvcf/Chilense_99.g.vcf    -all-sites    -O /gxfs_work1/cau/suaph275/Chilense_Hiseq/genotyped_vcf/Chilense_99.Invariant.genotypedGVC.vcf


vcftools --vcf Chilense_99.Invariant.genotypedGVC.vcf --max-maf 0 --remove-indels --max-missing 0.8 --recode --stdout | bgzip -c > test_invariant.vcf.gz

vcftools --vcf Chilense_99.Invariant.genotypedGVC.vcf --mac 1 --max-missing 0.8 --minQ 30 --recode --stdout | bgzip -c > test_variant.vcf.gz

tabix -p vcf test_invariant.vcf.gz
tabix -p vcf test_variant.vcf.gz
bcftools concat --allow-overlaps test_variant.vcf.gz test_invariant.vcf.gz -O z -o test_filtered.vcf.gz
bcftools query -f '%POS\n' test_filtered.vcf.gz | wc -l

# Number of total variants=709114385 (monomorphic+polymorphic)









